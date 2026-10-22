"""Correctness checks for a diploid graph G_D built from G.

Decoding semantics. A walk W in G_D from ("hom", s) to ("hom", t) decodes to
the set of ordered haplotype pairs (P, Q) obtained by reading W left to right:
a homozygous step ("hom", u) -> ("hom", v) appends v to both P and Q; a
heterozygous vertex ("het", p, i, j) between ("hom", x) and ("hom", y)
appends to P any x-y path inside branch i and to Q any x-y path inside
branch j.

G_D is *exact* if every ordered pair of s-t paths of G is decoded by exactly
one walk and every decoded pair is a pair of s-t paths of G.

Two checks:
* brute_force_check: enumerates everything; small graphs only.
* weighted_walk_count: a linear-time necessary condition for exactness that
  scales to real graphs. It sums, over all walks, the number of pairs each
  walk decodes to. For an exact G_D this equals N^2 (N = number of s-t paths).
"""

from __future__ import annotations

from collections import Counter, deque
from typing import Dict, List, Tuple

from .diploid import DiploidGraph
from .seqgraph import SequenceGraph


def _topo(gd: DiploidGraph) -> List:
    indeg: Dict = {v: 0 for v in gd.succ}
    for outs in gd.succ.values():
        for w in outs:
            indeg[w] += 1
    q = deque(v for v, d in indeg.items() if d == 0)
    order = []
    while q:
        u = q.popleft()
        order.append(u)
        for w in gd.succ[u]:
            indeg[w] -= 1
            if indeg[w] == 0:
                q.append(w)
    if len(order) != len(gd.succ):
        raise ValueError("diploid graph has a cycle")
    return order


def weighted_walk_count(gd: DiploidGraph) -> int:
    ways = {v: 0 for v in gd.succ}
    ways[gd.source] = 1
    for u in _topo(gd):
        w_u = ways[u]
        if not w_u:
            continue
        if u in gd.het:
            info = gd.het[u]
            w_u *= info.branch1.n_paths * info.branch2.n_paths
        for v in gd.succ[u]:
            ways[v] += w_u
    return ways[gd.sink]


def coverage(g: SequenceGraph, gd: DiploidGraph) -> Tuple[int, int]:
    """(pairs represented counted with multiplicity, N^2)."""
    s, t = g.terminals()
    n = g.count_paths(s, t)
    return weighted_walk_count(gd), n * n


def brute_force_check(g: SequenceGraph, gd: DiploidGraph) -> Dict[str, object]:
    s, t = g.terminals()
    haps = g.all_paths(s, t)
    expected = Counter((p, q) for p in haps for q in haps)

    decoded: Counter = Counter()
    stack: List[Tuple[object, List[Tuple[Tuple, Tuple]]]] = [(gd.source, [((s,), (s,))])]
    while stack:
        u, partial = stack.pop()
        if u == gd.sink:
            decoded.update(partial)
            continue
        for v in gd.succ[u]:
            if v in gd.het:
                info = gd.het[v]
                x, y = info.pnode.source, info.pnode.sink
                p_ext = g.all_paths(x, y, info.branch1.edges())
                q_ext = g.all_paths(x, y, info.branch2.edges())
                (nxt,) = gd.succ[v]
                new = [(a + pe[1:], b + qe[1:]) for a, b in partial for pe in p_ext for qe in q_ext]
                stack.append((nxt, new))
            else:
                y = v[1]
                stack.append((v, [(a + (y,), b + (y,)) for a, b in partial]))

    missing = expected - decoded
    spurious = decoded - expected
    duplicated = {k: c for k, c in decoded.items() if c > 1}
    return {
        "n_haplotypes": len(haps),
        "expected_pairs": sum(expected.values()),
        "decoded_pairs": sum(decoded.values()),
        "missing": sum(missing.values()),
        "spurious": sum(spurious.values()),
        "duplicated": len(duplicated),
        "exact": not missing and not spurious and not duplicated,
    }
