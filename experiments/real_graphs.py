"""Run the diploid construction on real pangenome subgraphs (GFA) and collect a table.

    python experiments/real_graphs.py BRCA1.gfa MHC.gfa --out experiments/results/real.csv

Input graphs (HPRC v2.0 Minigraph-Cactus, GRCh38, ec1M) and how to get GFA from
.gbz are described in experiments/README.md.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib

from diploidgraph.stats import gfa_stats


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa", nargs="+")
    ap.add_argument("--out", default="experiments/results/real.csv")
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for path in args.gfa:
        r = gfa_stats(path)
        print(json.dumps(r, indent=2))
        rows.append(r)
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
