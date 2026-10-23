"""How much does the diploid graph size depend on the vertex order?

In the non-series-parallel core the haplotype that is earlier in a fixed
topological order moves first; ties between topological orders are broken by
the order of vertices in the input. The P-node part does not depend on order.
This script reports the size for the GFA order and for random vertex orders.

    python experiments/order_sensitivity.py BRCA1.gfa --orders 20 --json brca1_order.json
"""

from __future__ import annotations

import argparse
import json
import statistics
import time

from diploidgraph import decompose, diploid_size
from diploidgraph.generators import shuffled
from diploidgraph.gfa import read_gfa

KINDS = ("het", "fork", "pair", "step")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--orders", type=int, default=20)
    ap.add_argument("--json")
    args = ap.parse_args()

    t0 = time.perf_counter()
    g = read_gfa(args.gfa).graph
    runs = []
    for k in range(args.orders + 1):
        h = g if k == 0 else shuffled(g, seed=k)
        size = diploid_size(h, decompose(h))
        runs.append({"order": "gfa" if k == 0 else f"random-{k}", "vertices": size["vertices"],
                     "edges": size["edges"], **{kk: size[kk] for kk in KINDS}})
    v = [r["vertices"] for r in runs]
    out = {
        "input": args.gfa,
        "orders": len(runs),
        "gfa_order_vertices": runs[0]["vertices"],
        "vertices_min_median_max": [min(v), statistics.median(v), max(v)],
        "relative_spread": (max(v) - min(v)) / min(v),
        "het_constant": len({r["het"] for r in runs}) == 1,
        "runs": runs,
        "time_s": time.perf_counter() - t0,
    }
    text = json.dumps(out, indent=2)
    print(json.dumps({k: out[k] for k in out if k != "runs"}, indent=2))
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
