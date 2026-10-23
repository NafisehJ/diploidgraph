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


def exactness_check(g: SequenceGraph, gd: DiploidGraph, exact_limit_vertices: int = 50_000) -> Dict[str, object]:
    """Necessary condition for exactness: weighted walk count == N^2.

    N (the number of source-sink paths of G) is counted directly on G, so it does
    not depend on the decomposition. For graphs with at most exact_limit_vertices
    vertices the comparison is exact; for larger ones it is done modulo two
    Mersenne primes, which keeps one machine-size number per vertex.

    This check cannot see a missing pair compensated by a duplicated pair of the
    same weight; count_decodings checks individual pairs and catches that.
    """
    s, t = g.terminals()
    if g.num_vertices() <= exact_limit_vertices:
        n = g.count_paths(s, t)
        return {"method": "exact", "passed": weighted_walk_count(gd) == n * n}
    ok = True
    for p in PRIMES:
        n = g.count_paths_mod(s, t, p)
        ok = ok and weighted_walk_count(gd, p) == n * n % p
    return {"method": "modulo 2^61-1 and 2^89-1", "passed": ok}


class PairDecoder:
    """Counts the walks of G_D that decode to a given ordered pair (P, Q).

    For an exact G_D every ordered pair of source-sink paths of G gives exactly 1.
    Unlike the weighted count, this is checked pair by pair, so it detects a
    missing pair even if another pair is duplicated. Used on sampled pairs of
    real haplotypes where enumerating all pairs is impossible.
    """

    def __init__(self, gd: DiploidGraph, dec):
        self.gd = gd
        self.order = {v: i for i, v in enumerate(_topo(gd))}
        self.tin: Dict[int, int] = {}
        self.tout: Dict[int, int] = {}
        self.leaf: Dict[Tuple, int] = {}
        clock = 0
        stack = [(dec.root, False)]
        while stack:
            node, done = stack.pop()
            if done:
                self.tout[id(node)] = clock
                continue
            self.tin[id(node)] = clock
            clock += 1
            if node.kind == "Q":
                self.leaf[(node.source, node.sink)] = self.tin[id(node)]
            stack.append((node, True))
            for c in reversed(node.children):
                stack.append((c, False))

    def _inside(self, edge, node) -> bool:
        k = self.leaf.get(edge)
        return k is not None and self.tin[id(node)] <= k < self.tout[id(node)]

    def _advance(self, path, pos, i, seg):
        """Index k such that seg.source + path[i+1..k] is a path inside seg.node ending at seg.sink, or None."""
        if seg is None:
            return i
        k = pos.get(seg.sink)
        if k is None or k <= i:
            return None
        prev = seg.source
        for x in range(i + 1, k + 1):
            if not self._inside((prev, path[x]), seg.node):
                return None
            prev = path[x]
        return k

    def count(self, P, Q) -> int:
        import heapq
        gd = self.gd
        if P[0] != gd.source[1] or Q[0] != gd.source[1]:
            return 0
        posP = {v: i for i, v in enumerate(P)}
        posQ = {v: i for i, v in enumerate(Q)}
        start = (gd.source, 0, 0)
        ways = {start: 1}
        heap = [(self.order[gd.source], 0, 0, start)]
        seen = {start}
        total = 0
        lastP, lastQ = len(P) - 1, len(Q) - 1
        while heap:
            _, _, _, st = heapq.heappop(heap)
            u, i, j = st
            w = ways.pop(st)
            if u == gd.sink and i == lastP and j == lastQ:
                total += w
                continue
            for v in gd.succ[u]:
                if v in gd.gadget:
                    gad = gd.gadget[v]
                    i2 = self._advance(P, posP, i, gad.seg1)
                    j2 = self._advance(Q, posQ, j, gad.seg2)
                    if i2 is None or j2 is None:
                        continue
                elif u[0] == "hom" and v[0] == "hom":
                    if i >= lastP or j >= lastQ or P[i + 1] != v[1] or Q[j + 1] != v[1]:
                        continue
                    i2, j2 = i + 1, j + 1
                else:
                    i2, j2 = i, j
                nxt = (v, i2, j2)
                ways[nxt] = ways.get(nxt, 0) + w
                if nxt not in seen:
                    seen.add(nxt)
                    heapq.heappush(heap, (self.order[v], i2, j2, nxt))
        return total


def count_decodings(gd: DiploidGraph, dec, P, Q) -> int:
    return PairDecoder(gd, dec).count(P, Q)


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
