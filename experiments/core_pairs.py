"""Diagnostic for non-series-parallel graphs: how many heterozygous vertices would an
exact construction need in the irreducible core?

After all series and parallel reductions, the core is a DAG whose edges each stand
for a series-parallel subgraph. For an ordered haplotype pair (P, Q), every place
where P and Q diverge inside the core is a segment from a core vertex x to a core
vertex y along two internally disjoint core paths. This script counts, over all
(x, y), the ordered pairs of internally disjoint core paths from x to y.

    python experiments/core_pairs.py BRCA1.gfa [--cap 2000000]
"""

from __future__ import annotations

import argparse
from collections import defaultdict
from typing import Dict, List, Set, Tuple

from diploidgraph import decompose
from diploidgraph.gfa import read_gfa


def core_graph(dec) -> Tuple[Dict, List]:
    (r,) = dec.nodes("R")
    succ: Dict = defaultdict(list)
    for c in r.children:
        succ[c.source].append(c.sink)
        succ.setdefault(c.sink, [])
    return succ, r.children


def topo(succ) -> List:
    indeg = defaultdict(int)
    for u in succ:
        for v in succ[u]:
            indeg[v] += 1
    order, stack = [], [u for u in succ if indeg[u] == 0]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in succ[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                stack.append(v)
    return order


def count_avoiding(succ, order, x, y, banned: Set) -> int:
    ways = {x: 1}
    for u in order:
        w = ways.get(u, 0)
        if not w or u == y:
            continue
        for v in succ[u]:
            if v not in banned:
                ways[v] = ways.get(v, 0) + w
    return ways.get(y, 0)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--cap", type=int, default=2_000_000, help="stop after enumerating this many core paths")
    args = ap.parse_args()

    dec = decompose(read_gfa(args.gfa).graph)
    if dec.is_series_parallel:
        print("graph is series-parallel: no core")
        return
    succ, edges = core_graph(dec)
    order = topo(succ)
    print(f"core: {len(succ)} vertices, {len(edges)} edges")

    total_pairs = 0
    positions = 0
    enumerated = 0
    for x in order:
        per_y = defaultdict(int)
        stack = [(x, (x,))]
        while stack:
            u, path = stack.pop()
            for v in succ[u]:
                p = path + (v,)
                enumerated += 1
                if enumerated > args.cap:
                    print(f"STOPPED: more than {args.cap} core paths; pairs so far {total_pairs}")
                    return
                inner = set(p[1:-1])
                n_q = count_avoiding(succ, order, x, v, inner)
                if len(p) == 2:
                    n_q -= 1  # q may not equal the single edge p itself
                per_y[v] += n_q
                stack.append((v, p))
        for y, k in per_y.items():
            if k:
                positions += 1
                total_pairs += k
    print(f"hetero-positions in core (x, y with a disjoint pair): {positions}")
    print(f"ordered internally disjoint core path pairs (= het vertices needed): {total_pairs}")
    print(f"core paths enumerated: {enumerated}")


if __name__ == "__main__":
    main()
