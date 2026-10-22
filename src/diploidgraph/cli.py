"""Command line: build a diploid graph from a GFA file and print size statistics.

    diploidgraph graph.gfa [--allow-partial] [--json out.json]
"""

from __future__ import annotations

import argparse
import json

from .stats import gfa_stats


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--allow-partial", action="store_true",
                    help="for non-series-parallel inputs, copy the non-SP region homozygously only")
    ap.add_argument("--json", help="write statistics to this file")
    args = ap.parse_args(argv)

    text = json.dumps(gfa_stats(args.gfa, allow_partial=args.allow_partial), indent=2)
    print(text)
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
