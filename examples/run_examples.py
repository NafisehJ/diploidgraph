"""Small worked examples: input graph, decomposition, diploid graph, exactness check.

    python examples/run_examples.py            # writes examples/out/
    python examples/run_examples.py --png      # also renders PNGs (needs Graphviz `dot`)
"""

from __future__ import annotations

import argparse
import pathlib
import shutil
import subprocess

from diploidgraph import SequenceGraph, build_diploid, decompose, naive_pair_count
from diploidgraph.export import diploid_graph_dot, sequence_graph_dot
from diploidgraph.generators import random_sp_dag
from diploidgraph.validate import brute_force_check

OUT = pathlib.Path(__file__).parent / "out"


def examples():
    yield "01_bubble", "One bubble (a SNP): s -> {a, b} -> t.", SequenceGraph.from_edges(
        [("s", "a"), ("a", "t"), ("s", "b"), ("b", "t")])
    yield "02_three_alleles", "One bubble with three alleles.", SequenceGraph.from_edges(
        [("s", "a"), ("a", "t"), ("s", "b"), ("b", "t"), ("s", "c"), ("c", "t")])
    yield "03_two_bubbles", "Two bubbles in series: pairs differ in either, both, or none.", SequenceGraph.from_edges(
        [("s", "a"), ("s", "b"), ("a", "m"), ("b", "m"), ("m", "c"), ("m", "d"), ("c", "t"), ("d", "t")])
    yield "04_nested", "A bubble nested inside one allele of an outer bubble (an insertion carrying a SNP).", SequenceGraph.from_edges(
        [("s", "a"), ("a", "b"), ("a", "c"), ("b", "d"), ("c", "d"), ("d", "t"), ("s", "t")])
    yield "05_random_sp", "Random nested series-parallel graph (seed 7).", random_sp_dag(n_ops=10, seed=7)
    yield "06_w_motif", "Smallest non-series-parallel graph (W-motif) in series with a bubble; core handled by pair states.", SequenceGraph.from_edges(
        [("s", "a"), ("s", "b"), ("a", "m"), ("b", "m"),                  # bubble
         ("m", "c"), ("m", "d"), ("c", "d"), ("c", "t"), ("d", "t")])     # W-motif: diagonal c -> d


def render(dot_path: pathlib.Path) -> None:
    subprocess.run(["dot", "-Tpng", str(dot_path), "-o", str(dot_path.with_suffix(".png"))], check=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--png", action="store_true")
    args = ap.parse_args()
    if args.png and not shutil.which("dot"):
        raise SystemExit("Graphviz `dot` not found; install graphviz or drop --png")
    OUT.mkdir(exist_ok=True)

    # Het vertices are ordered pairs, so the like-for-like naive count is N(N-1) (ordered).
    rows = ["| Example | Description | \\|V\\| | \\|E\\| | s-t paths N | P-nodes | \\|V(G_D)\\| | new vertices | naive N(N-1) | exact |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    for name, desc, g in examples():
        s, t = g.terminals()
        n = g.count_paths(s, t)
        dec = decompose(g)
        (OUT / f"{name}_input.dot").write_text(sequence_graph_dot(g))
        gd, _ = build_diploid(g, dec)
        exact = brute_force_check(g, gd)["exact"]
        (OUT / f"{name}_diploid.dot").write_text(diploid_graph_dot(gd))
        row = [gd.num_vertices(), gd.num_het_vertices(), naive_pair_count(n, ordered=True), "yes" if exact else "NO"]
        rows.append(f"| {name} | {desc} | {g.num_vertices()} | {g.num_edges()} | {n} | {len(dec.nodes('P'))} | " + " | ".join(map(str, row)) + " |")
        if args.png:
            for f in OUT.glob(f"{name}_*.dot"):
                render(f)

    table = "\n".join(rows)
    (OUT / "summary.md").write_text(table + "\n")
    print(table)


if __name__ == "__main__":
    main()
