"""Per-pair exactness check on a real graph.

Builds the diploid graph and, for sampled ordered pairs of haplotypes, counts the
walks that decode to each pair; an exact diploid graph gives exactly 1 for every
pair. Pairs are drawn from (a) the observed source-to-sink haplotypes and
(b) random source-to-sink paths (uniform random successor at every vertex),
which include recombinant pairs that no sample carries.

    python experiments/pair_check.py BRCA1.gfa --pairs 500 --json brca1_pairs.json
"""

from __future__ import annotations

import argparse
import json
import random
import time

from diploidgraph import build_diploid
from diploidgraph.gfa import iter_gfa_walks, read_gfa
from diploidgraph.validate import PairDecoder


def random_path(g, s, t, rng):
    path = [s]
    while path[-1] != t:
        path.append(rng.choice(sorted(g.succ[path[-1]], key=str)))
    return path


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--pairs", type=int, default=500, help="pairs of each kind (observed, random)")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--json")
    args = ap.parse_args()
    rng = random.Random(args.seed)

    t0 = time.perf_counter()
    g = read_gfa(args.gfa).graph
    s, t = g.terminals()
    gd, dec = build_diploid(g)
    decoder = PairDecoder(gd, dec)
    haps = [w for w in iter_gfa_walks(args.gfa) if w[0] == s and w[-1] == t]

    results = {}
    for kind in ("observed", "random"):
        counts = {}
        if kind == "observed" and not haps:
            results[kind] = {"pairs": 0, "decoded_exactly_once": 0, "other_counts": {}}
            continue
        for _ in range(args.pairs):
            if kind == "observed":
                P, Q = rng.choice(haps), rng.choice(haps)
            else:
                P, Q = random_path(g, s, t, rng), random_path(g, s, t, rng)
            c = decoder.count(P, Q)
            counts[c] = counts.get(c, 0) + 1
        results[kind] = {"pairs": args.pairs, "decoded_exactly_once": counts.get(1, 0),
                         "other_counts": {str(k): v for k, v in counts.items() if k != 1}}
    out = {"input": args.gfa, "observed_haplotypes": len(haps), **results,
           "passed": all(r["decoded_exactly_once"] == r["pairs"] for r in results.values()),
           "note": None if haps else "no source-to-sink haplotype walks in the GFA; only random pairs checked",
           "time_s": time.perf_counter() - t0}
    text = json.dumps(out, indent=2)
    print(text)
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
