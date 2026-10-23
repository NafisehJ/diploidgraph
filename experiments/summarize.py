"""Print the README results table from the JSON files written by reproduce.sh.

    python experiments/summarize.py experiments/results
"""

from __future__ import annotations

import json
import pathlib
import sys


def load(path: pathlib.Path):
    return json.loads(path.read_text()) if path.exists() else None


def fmt(x) -> str:
    if x is None:
        return "-"
    if isinstance(x, bool):
        return "passed" if x else "FAILED"
    if isinstance(x, float):
        return f"{x:,.2f}"
    return f"{x:,}"


def main() -> None:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "experiments/results")
    loci = sorted(p.name[: -len("_stats.json")] for p in out.glob("*_stats.json"))
    if not loci:
        raise SystemExit(f"no *_stats.json files in {out}")
    stats = {l: load(out / f"{l}_stats.json") for l in loci}
    sup = {l: load(out / f"{l}_unsupported.json") for l in loci}
    pairs = {l: load(out / f"{l}_pairs.json") for l in loci}
    order = {l: load(out / f"{l}_order.json") for l in loci}

    def row(label, f):
        cells = []
        for l in loci:
            try:
                cells.append(f(stats[l], sup[l]))
            except (KeyError, TypeError, ZeroDivisionError):
                cells.append("-")
        return f"| {label} | " + " | ".join(cells) + " |"

    def pct(n, d):
        return f"{100 * n / d:.2f}%"

    rows = [
        "| | " + " | ".join(loci) + " |",
        "|---|" + "---|" * len(loci),
        row("Input vertices / edges", lambda s, u: f"{fmt(s['n_vertices'])} / {fmt(s['n_edges'])}"),
        row("Haplotype walks (distinct; source-to-sink)", lambda s, u: f"{s['n_haplotype_walks']} "
            f"({s['n_distinct_haplotype_walks']}; {s['n_full_haplotype_walks']})"),
        row("P-nodes (max branches)", lambda s, u: f"{fmt(s['n_pnodes'])} ({s['max_pnode_branches']})"),
        row("Non-series-parallel core: vertices / edges", lambda s, u: f"{fmt(s['core_vertices'])} / {fmt(s['core_edges'])}"),
        row("Diploid graph vertices / edges", lambda s, u: f"{fmt(s['diploid_vertices'])} / {fmt(s['diploid_edges'])}"),
        row("Diploid vertices per input vertex", lambda s, u: f"{s['diploid_vertices'] / s['n_vertices']:.2f}"),
        row("Naive pairing of observed walks, paired vertices", lambda s, u: fmt(s["naive_ordered_observed_nodes"])),
        row("Fewer vertices than naive (x)", lambda s, u: f"{s['naive_ordered_observed_nodes'] / s['diploid_vertices']:,.0f}"),
        row("log10 of ordered pairs of all source-to-sink paths", lambda s, u: f"{s['log10_naive_ordered_all_paths']:.0f}"),
        row("Added vertices used by pairs of complete observed haplotypes",
            lambda s, u: pct(u["all_new"]["supported_full"], u["all_new"]["total"])),
        row("Added vertices used by any observed haplotype (incl. partial walks)",
            lambda s, u: pct(u["all_new"]["total"] - u["all_new"]["still_unsupported"], u["all_new"]["total"])),
        row("Decompose + build time (s)", lambda s, u: f"{s['time_decompose_s'] + s['time_build_s']:.2f}"),
        row("Peak memory (MB)", lambda s, u: f"{s['peak_rss_mb']:,.0f}"),
        row("Exactness check (weighted count = N^2)", lambda s, u: f"{fmt(s['exactness_passed'])} ({s['exactness_check']})"),
    ]
    extra = []
    for l in loci:
        pr, od = pairs[l], order[l]
        if pr:
            extra.append(f"{l}: sampled pairs decoded exactly once: observed {pr['observed']['decoded_exactly_once']}/"
                         f"{pr['observed']['pairs']}, random {pr['random']['decoded_exactly_once']}/{pr['random']['pairs']}")
        if od:
            lo, med, hi = od["vertices_min_median_max"]
            extra.append(f"{l}: diploid vertices over {od['orders']} vertex orders: GFA order {od['gfa_order_vertices']:,}, "
                         f"min {lo:,}, median {med:,.0f}, max {hi:,}")
    print("\n".join(rows))
    if extra:
        print()
        print("\n".join(extra))


if __name__ == "__main__":
    main()
