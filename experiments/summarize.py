"""Build the README results table and the figure data from the JSON files written
by reproduce.sh.

    python experiments/summarize.py experiments/results

Prints the table (the same rows as README.md) and writes
experiments/results/results_table.md and experiments/results/real_loci.json
(the input of plot_real.py).
"""

from __future__ import annotations

import json
import math
import pathlib
import sys

ORDER = ["BRCA1", "MHC"]


def load(path: pathlib.Path):
    return json.loads(path.read_text()) if path.exists() else None


def n(x) -> str:
    return f"{x:,.0f}"


def main() -> None:
    out = pathlib.Path(sys.argv[1] if len(sys.argv) > 1 else "experiments/results")
    found = [p.name[: -len("_stats.json")] for p in out.glob("*_stats.json")]
    if not found:
        raise SystemExit(f"no *_stats.json files in {out}")
    loci = [l for l in ORDER if l in found] + sorted(l for l in found if l not in ORDER)
    J = {l: {k: load(out / f"{l}_{k}.json") for k in ("stats", "product", "pairs", "order", "unsupported")}
         for l in loci}

    def row(label, f):
        cells = []
        for l in loci:
            try:
                cells.append(f(J[l]))
            except (KeyError, TypeError, ZeroDivisionError):
                cells.append("not run")
        return f"| {label} | " + " | ".join(cells) + " |"

    def core_support(j):
        u = j["unsupported"]
        sup = sum(u[k]["supported_full"] for k in ("fork", "pair", "step"))
        tot = sum(u[k]["total"] for k in ("fork", "pair", "step"))
        return f"{100 * sup / tot:.2f}%" if tot else "no core"

    def all_support(j):
        a = j["unsupported"]["all_new"]
        return f"{100 * a['supported_full'] / a['total']:.2f}%"

    rows = [
        "| | " + " | ".join(loci) + " |",
        "|---|" + "---|" * len(loci),
        row("Input vertices / edges", lambda j: f"{n(j['stats']['n_vertices'])} / {n(j['stats']['n_edges'])}"),
        row("Haplotype walks: all / distinct source-to-sink",
            lambda j: f"{j['stats']['n_haplotype_walks']} / {j['stats']['n_distinct_full_haplotype_walks']}"),
        row("Dropped links / walks with mixed orientation",
            lambda j: f"{j['stats']['dropped_links']} / {j['stats']['unoriented_haplotypes']}"),
        row("P-nodes (max branches)", lambda j: f"{n(j['stats']['n_pnodes'])} ({j['stats']['max_pnode_branches']})"),
        row("Non-series-parallel core: vertices / edges",
            lambda j: f"{n(j['stats']['core_vertices'])} / {n(j['stats']['core_edges'])}"),
        row("**Diploid graph vertices / edges** (GFA vertex order)",
            lambda j: f"**{n(j['stats']['diploid_vertices'])} / {n(j['stats']['diploid_edges'])}**"),
        row("Diploid vertices per input vertex",
            lambda j: f"{j['stats']['diploid_vertices'] / j['stats']['n_vertices']:.2f}"),
        row(f"Diploid vertices over {J[loci[0]]['order']['orders'] if J[loci[0]]['order'] else 'several'} vertex orders: min / median / max",
            lambda j: " / ".join(n(x) for x in j["order"]["vertices_min_median_max"])),
        row("Product-graph ablation vertices over the same orders: min / median / max",
            lambda j: " / ".join(n(x) for x in j["order"]["product_vertices_min_median_max"])
            + (" (lower bound)" if j["order"]["product_any_capped"] else "")),
        row("Product ablation / diploid graph, median over orders",
            lambda j: f"{j['order']['median_ratio_product_over_ours']:.1f}x"),
        row("Naive pairing of distinct complete haplotypes, paired positions",
            lambda j: n(j["stats"]["naive_ordered_observed_nodes"])),
        row("Ordered pairs of source-to-sink paths represented",
            lambda j: f"about 10^{j['stats']['log10_st_paths'] * 2:.0f}"),
        row("Sampled pairs decoded exactly once: observed haplotypes; random-successor paths",
            lambda j: f"{j['pairs']['observed']['decoded_exactly_once']}/{j['pairs']['observed']['pairs']}; "
                      f"{j['pairs']['random']['decoded_exactly_once']}/{j['pairs']['random']['pairs']}"),
        row("Whole-graph check (weighted walk count = N^2)",
            lambda j: ("passed" if j["stats"]["exactness_passed"] else "FAILED")
            + (" (exact)" if j["stats"]["exactness_check"] == "exact" else " (modulo two primes)")),
        row("Added vertices used by pairs of complete observed haplotypes: all / core only",
            lambda j: f"{all_support(j)} / {core_support(j)}"),
        row("Decompose + build time", lambda j: f"{j['stats']['time_decompose_s'] + j['stats']['time_build_s']:.2f} s"),
        row("Peak memory (whole process)", lambda j: f"{j['stats']['peak_rss_mb']:,.0f} MB"),
    ]
    table = "\n".join(rows)
    print(table)
    (out / "results_table.md").write_text(table + "\n")

    fig = {"source": "written by experiments/summarize.py from reproduce.sh results"}
    for l in loci:
        j = J[l]
        try:
            fig[l] = {
                "input_nodes": j["stats"]["n_vertices"],
                "diploid_nodes": j["stats"]["diploid_vertices"],
                "product_nodes": j["product"]["product_vertices"],
                "naive_observed_nodes": j["stats"]["naive_ordered_observed_nodes"],
                "haplotypes": j["stats"]["n_distinct_full_haplotype_walks"],
                "build_s": j["stats"]["time_decompose_s"] + j["stats"]["time_build_s"],
                "ratio_median": j["order"]["median_ratio_product_over_ours"],
            }
        except (KeyError, TypeError):
            pass
    (out / "real_loci.json").write_text(json.dumps(fig, indent=2) + "\n")
    print(f"\nwrote {out / 'results_table.md'} and {out / 'real_loci.json'}")


if __name__ == "__main__":
    main()
