"""Are the unsupported new vertices explained by partial haplotype walks?

Runs the observed-support computation twice in one pass: with full source-to-sink
walks only, and with partial walks added (their P-node branches and complete
block traversals). Vertices supported only in the second run are seen only
together with partial walks.

    python experiments/explain_unsupported.py MHC.gfa --json mhc_unsupported.json
"""

from __future__ import annotations

import argparse
import json
import time

from diploidgraph import decompose, diploid_size
from diploidgraph.gfa import iter_gfa_walks, read_gfa
from diploidgraph.stats import peak_rss_mb
from diploidgraph.support import KINDS, observed_support_with_partial


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--json")
    args = ap.parse_args()

    t0 = time.perf_counter()
    g = read_gfa(args.gfa).graph
    dec = decompose(g)
    total = diploid_size(g, dec)
    full, with_partial, counts = observed_support_with_partial(g, iter_gfa_walks(args.gfa), dec)

    out = {"input": args.gfa, **counts}
    sums = {"total": 0, "supported_full": 0, "explained_by_partial": 0, "still_unsupported": 0}
    for k in KINDS:
        row = {
            "total": total[k],
            "supported_full": len(full[k]),
            "explained_by_partial": len(with_partial[k] - full[k]),
        }
        row["still_unsupported"] = row["total"] - len(with_partial[k])
        out[k] = row
        for key in sums:
            sums[key] += row[key]
    out["all_new"] = sums
    out["time_s"] = time.perf_counter() - t0
    out["peak_rss_mb"] = peak_rss_mb()
    text = json.dumps(out, indent=2)
    print(text)
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")


if __name__ == "__main__":
    main()
