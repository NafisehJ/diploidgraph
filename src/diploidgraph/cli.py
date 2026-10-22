"""Command line: build a diploid graph from a GFA file and print size statistics.

    python -m diploidgraph.cli graph.gfa [--allow-partial] [--json out.json]
"""

from __future__ import annotations

import argparse
import json
import time
import tracemalloc

from .decompose import decompose
from .diploid import build_diploid, naive_pair_count
from .gfa import read_gfa
from .validate import coverage


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--allow-partial", action="store_true", help="copy R-node regions homozygously only")
    ap.add_argument("--json", help="write statistics to this file")
    args = ap.parse_args(argv)

    stats = {"input": args.gfa}
    tracemalloc.start()

    t0 = time.perf_counter()
    data = read_gfa(args.gfa)
    g = data.graph
    stats["time_read_s"] = time.perf_counter() - t0

    s, t = g.terminals()
    stats.update(
        n_vertices=g.num_vertices(),
        n_edges=g.num_edges(),
        dropped_links=data.dropped_links,
        n_haplotype_walks=len(data.haplotypes),
        n_distinct_haplotype_walks=len({w for _, w in data.haplotypes}),
        unoriented_haplotypes=data.unoriented_haplotypes,
    )

    t0 = time.perf_counter()
    dec = decompose(g)
    stats["time_decompose_s"] = time.perf_counter() - t0
    pnodes = dec.nodes("P")
    stats.update(
        is_series_parallel=dec.is_series_parallel,
        n_pnodes=len(pnodes),
        max_pnode_branches=max((len(p.children) for p in pnodes), default=0),
        n_rnodes=len(dec.nodes("R")),
    )

    t0 = time.perf_counter()
    gd, _ = build_diploid(g, dec, allow_partial=args.allow_partial)
    stats["time_build_s"] = time.perf_counter() - t0
    stats["peak_python_memory_bytes"] = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    n_paths = g.count_paths(s, t)
    represented, total = coverage(g, gd)
    stats.update(
        diploid_vertices=gd.num_vertices(),
        diploid_edges=gd.num_edges(),
        diploid_het_vertices=gd.num_het_vertices(),
        n_st_paths=str(n_paths),  # can exceed 64-bit; keep as string in JSON
        naive_vertices_all_paths=str(naive_pair_count(n_paths)),
        naive_vertices_observed=naive_pair_count(stats["n_distinct_haplotype_walks"]),
        pair_coverage=represented / total if total else None,
    )

    text = json.dumps(stats, indent=2)
    print(text)
    if args.json:
        with open(args.json, "w") as fh:
            fh.write(text + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
