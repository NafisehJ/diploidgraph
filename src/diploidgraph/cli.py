"""Command line: build a diploid graph from a GFA file and print size statistics.

    diploidgraph graph.gfa [--count-only] [--json out.json]
"""

from __future__ import annotations

import argparse
import json

from .stats import gfa_stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--count-only", action="store_true",
                    help="compute the size of the diploid graph without building it (large graphs)")
    ap.add_argument("--json", help="write statistics to this file")
    args = ap.parse_args(argv)

    text = json.dumps(gfa_stats(args.gfa, count_only=args.count_only), indent=2)
    print(text)
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
