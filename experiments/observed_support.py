"""How many new vertices of the diploid graph are used by pairs of observed haplotypes?

Answers: "Of the diploid nodes added, how many are contained in a sample and how
many aren't?" A new vertex is *supported* if the diploid-graph walk of some
ordered pair of distinct observed haplotypes (full source-to-sink walks in the
GFA) passes through it; the others represent recombinant pairs only.

    python experiments/observed_support.py BRCA1.gfa --json brca1_support.json
"""

from __future__ import annotations

import argparse
import json
import time

from diploidgraph import decompose, diploid_size
from diploidgraph.gfa import iter_gfa_walks, read_gfa
from diploidgraph.stats import peak_rss_mb
from diploidgraph.support import KINDS, observed_support


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--json")
    args = ap.parse_args()

    t0 = time.perf_counter()
    g = read_gfa(args.gfa).graph
    dec = decompose(g)
    total = diploid_size(g, dec)
    t1 = time.perf_counter()
    supported, counts = observed_support(g, iter_gfa_walks(args.gfa), dec)
    t2 = time.perf_counter()

    out = {"input": args.gfa, **counts}
    sup_all = tot_all = 0
    for k in KINDS:
        n_sup, n_tot = len(supported[k]), total[k]
        out[f"{k}_supported"] = n_sup
        out[f"{k}_total"] = n_tot
        sup_all += n_sup
        tot_all += n_tot
    out["new_supported"] = sup_all
    out["new_total"] = tot_all
    out["new_supported_fraction"] = sup_all / tot_all if tot_all else None
    out["time_size_s"] = t1 - t0
    out["time_support_s"] = t2 - t1
    out["peak_rss_mb"] = peak_rss_mb()
    text = json.dumps(out, indent=2)
    print(text)
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
