"""State-space size of the PanGenie model on the same GFA, next to ours.

    python experiments/external_state_space.py BRCA1.gfa MHC.gfa --out experiments/results/external.csv

Size comparison only (no reads, no genotyping). For each method the table
reports what it covers next to how many states it uses:

* ours: every ordered pair of source-to-sink paths, exactly once
  (|V(G_D)| comes from real_graphs.py / summarize.py);
* PanGenie: N^2 ordered pairs of panel haplotypes per bubble, bubbles closer
  than k = 31 bp combined; covers ordered pairs of panel mosaics that switch
  only between bubbles;

Definitions and proxies: see the docstring of diploidgraph/external.py.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import time

from diploidgraph.external import external_sizes
from diploidgraph.stats import peak_rss_mb


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa", nargs="+")
    ap.add_argument("-k", type=int, default=31, help="PanGenie k-mer size used to combine nearby bubbles (default 31)")
    ap.add_argument("--reference-samples", default="GRCh38,CHM13",
                    help="samples left out of the PanGenie panel (comma-separated; empty string keeps all)")
    ap.add_argument("--out", default="experiments/results/external.csv")
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in args.gfa:
        t0 = time.perf_counter()
        try:
            refs = tuple(x for x in args.reference_samples.split(",") if x)
            r = external_sizes(path, k=args.k, reference_samples=refs)
        except ValueError as err:  # e.g. a GFA without complete haplotype walks
            r = {"input": path, "error": str(err)}
        r["time_s"] = time.perf_counter() - t0
        r["peak_rss_mb"] = peak_rss_mb()
        print(json.dumps(r, indent=2))
        rows.append(r)
    fields = list(dict.fromkeys(k for r in rows for k in r))
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
