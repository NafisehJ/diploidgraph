"""Correctness checks for a diploid graph G_D built from G.

Decoding of walks is defined in diploid.py. G_D is *exact* if every ordered
pair of s-t paths of G is decoded by exactly one walk and every decoded pair is
a pair of s-t paths of G.

* brute_force_check enumerates everything; small graphs only.
* weighted_walk_count is a linear-time necessary condition for exactness: the
  number of pairs decoded, summed over all walks. It equals N^2 for an exact
  G_D (N = number of s-t paths).
"""

from __future__ import annotations

from collections import Counter, deque
from typing import Dict, List, Tuple

from .diploid import DiploidGraph, Segment
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


def weighted_walk_count(gd: DiploidGraph, modulus: int | None = None) -> int:
    """Sum over all walks of the number of pairs each walk decodes to.

    With `modulus`, the count is computed modulo that number. Path counts on
    locus-level graphs have tens of thousands of digits, so the exact count
    would hold a huge integer per vertex; comparing residues modulo large
    primes keeps memory linear.
    """
    ways = {v: 0 for v in gd.succ}
    ways[gd.source] = 1
    for u in _topo(gd):
        w_u = ways[u]
        if not w_u:
            continue
        if u in gd.gadget:
            w_u *= gd.gadget[u].weight
        if modulus:
            w_u %= modulus
        for v in gd.succ[u]:
            ways[v] = (ways[v] + w_u) % modulus if modulus else ways[v] + w_u
    return ways[gd.sink]


PRIMES = (2**61 - 1, 2**89 - 1)


def exactness_check(g: SequenceGraph, gd: DiploidGraph, n_paths: int | None = None,
                    exact_limit_bits: int = 8000) -> Dict[str, object]:
    """Check weighted walk count == N^2 (N = number of source-sink paths of G).

    Exact when N is small; otherwise modulo two Mersenne primes, which keeps one
    machine-size number per vertex instead of a number with tens of thousands of
    digits. Pass n_paths (e.g. the decomposition root's count) to avoid
    recounting paths.
    """
    if n_paths is None:
        s, t = g.terminals()
        n_paths = g.count_paths(s, t)
    n = n_paths
    if n.bit_length() <= exact_limit_bits:
        return {"method": "exact", "passed": weighted_walk_count(gd) == n * n}
    ok = all(weighted_walk_count(gd, p) == (n % p) * (n % p) % p for p in PRIMES)
    return {"method": "modulo 2^61-1 and 2^89-1", "passed": ok}


def coverage(g: SequenceGraph, gd: DiploidGraph) -> Tuple[int, int]:
    """(pairs represented counted with multiplicity, N^2)."""
    s, t = g.terminals()
    n = g.count_paths(s, t)
    return weighted_walk_count(gd), n * n


def brute_force_check(g: SequenceGraph, gd: DiploidGraph) -> Dict[str, object]:
    s, t = g.terminals()
    haps = g.all_paths(s, t)
    expected = Counter((p, q) for p in haps for q in haps)
    seg_cache: Dict[Segment, List[Tuple]] = {}

    def extensions(seg):
        if seg is None:
            return [()]
        if seg not in seg_cache:
            seg_cache[seg] = [p[1:] for p in g.all_paths(seg.source, seg.sink, seg.node.edges())]
        return seg_cache[seg]

    decoded: Counter = Counter()
    stack: List[Tuple[object, List[Tuple[Tuple, Tuple]]]] = [(gd.source, [((s,), (s,))])]
    while stack:
        u, partial = stack.pop()
        if u == gd.sink:
            decoded.update(partial)
            continue
        for v in gd.succ[u]:
            if v in gd.gadget:
                gad = gd.gadget[v]
                new = [(a + e1, b + e2) for a, b in partial
                       for e1 in extensions(gad.seg1) for e2 in extensions(gad.seg2)]
            elif u[0] == "hom" and v[0] == "hom":
                new = [(a + (v[1],), b + (v[1],)) for a, b in partial]
            else:
                new = partial
            stack.append((v, new))

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
