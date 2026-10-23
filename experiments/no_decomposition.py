"""Compare the diploid graph with the same construction without decomposition (ablation).

Both represent exactly the same set of haplotype pairs (every ordered pair of
source-to-sink paths, each once) and both are counted in vertices of a graph.
The baseline applies the pair-state construction to every vertex of the input
(no series-parallel decomposition); the diploid graph applies it only inside the
non-series-parallel core and handles bubbles with P-node vertices. The difference
is the saving due to the decomposition.

    python experiments/no_decomposition.py BRCA1.gfa --json brca1_nodecomp.json
"""

from __future__ import annotations

import argparse
import json
import time

from diploidgraph import decompose, diploid_size
from diploidgraph.decompose import trivial_decomposition
from diploidgraph.gfa import read_gfa


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--max-pair-states", type=int, default=30_000_000,
                    help="stop the baseline after this many pair states (its size is then a lower bound)")
    ap.add_argument("--json")
    args = ap.parse_args()

    g = read_gfa(args.gfa).graph
    t0 = time.perf_counter()
    ours = diploid_size(g, decompose(g))
    t1 = time.perf_counter()
    base = diploid_size(g, trivial_decomposition(g), max_pair_states=args.max_pair_states)
    t2 = time.perf_counter()
    out = {
        "input": args.gfa,
        "input_vertices": g.num_vertices(),
        "diploid_vertices": ours["vertices"],
        "diploid_edges": ours["edges"],
        "no_decomposition_vertices": base["vertices"],
        "no_decomposition_edges": base["edges"],
        "no_decomposition_is_lower_bound": base["capped"],
        "no_decomposition_over_diploid": base["vertices"] / ours["vertices"],
        "time_diploid_s": t1 - t0,
        "time_no_decomposition_s": t2 - t1,
    }
    text = json.dumps(out, indent=2)
    print(text)
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
