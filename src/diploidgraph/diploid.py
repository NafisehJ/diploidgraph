"""Diploid sequence graph construction.

Reading a walk of G_D builds an ordered pair of haplotypes (P, Q). A walk starts
at ("hom", s) with P = Q = (s,) and is read edge by edge:

* an edge ("hom", u) -> ("hom", v) appends v to both P and Q;
* entering a *gadget* vertex appends its segments: segment 1 to P and segment 2
  to Q. A segment (x, y, B) means "any x-y path inside the two-terminal
  subgraph B", appended without x and with y. A missing segment appends nothing;
* all other edges append nothing.

Vertices of G_D:

* ("hom", v): both haplotypes are at v (one per vertex of G).
* ("het", p, i, j): P-node p with terminals (x, y), ordered pair of distinct
  branches i != j. Gadget with segments (x, y, B_i) and (x, y, B_j), placed
  between ("hom", x) and ("hom", y).

Series-parallel inputs need nothing else. If the input is not series-parallel,
the decomposition leaves an irreducible core: a DAG whose edges ("core edges")
each stand for a series-parallel subgraph, and whose distinct core edges are
internally vertex-disjoint. Two haplotypes that leave a core vertex x along
different core edges cannot meet again except at a core vertex. The core part
of G_D tracks such a divergence with a pair state:

* ("fork", c1, c2): at core vertex x, haplotype 1 leaves along core edge c1
  and haplotype 2 along c2 != c1. Gadget with segments (x, head(c1), c1) and
  (x, head(c2), c2).
* ("pair", u, v): haplotypes are at distinct core vertices u (hap 1), v (hap 2)
  and have not met since they diverged. Appends nothing.
* ("step", h, c, w): haplotype h moves along core edge c while the other waits
  at w. Gadget with one segment. Only the haplotype at the earlier vertex in a
  fixed topological order moves, which makes the walk of every pair unique and
  guarantees that the two haplotypes cannot share a vertex before they meet.
  When both reach the same core vertex z the walk returns to ("hom", z).

Only pair states reachable from a fork are created, so the core part has at
most n_c(n_c - 1) pair states for a core with n_c vertices.
"""

from __future__ import annotations

from collections import deque
from dataclasses import dataclass, field
from typing import Dict, Hashable, List, Optional, Set, Tuple

from .decompose import DecompNode, Decomposition, decompose
from .seqgraph import SequenceGraph, Vertex

DVertex = Tuple[Hashable, ...]


@dataclass(frozen=True)
class Segment:
    source: Vertex
    sink: Vertex
    node: DecompNode  # two-terminal subgraph the path runs inside


@dataclass
class Gadget:
    seg1: Optional[Segment]  # appended to haplotype 1
    seg2: Optional[Segment]  # appended to haplotype 2

    @property
    def weight(self) -> int:
        """Number of haplotype pairs this vertex stands for."""
        w = 1
        for s in (self.seg1, self.seg2):
            if s is not None:
                w *= s.node.n_paths
        return w


@dataclass
class DiploidGraph:
    source: DVertex
    sink: DVertex
    succ: Dict[DVertex, Set[DVertex]] = field(default_factory=dict)
    gadget: Dict[DVertex, Gadget] = field(default_factory=dict)

    def add_vertex(self, a: DVertex) -> None:
        self.succ.setdefault(a, set())

    def add_edge(self, a: DVertex, b: DVertex) -> None:
        self.succ.setdefault(a, set()).add(b)
        self.succ.setdefault(b, set())

    def num_vertices(self) -> int:
        return len(self.succ)

    def num_edges(self) -> int:
        return sum(len(o) for o in self.succ.values())

    def count_by_kind(self) -> Dict[str, int]:
        out: Dict[str, int] = {}
        for v in self.succ:
            out[v[0]] = out.get(v[0], 0) + 1
        return out

    def num_het_vertices(self) -> int:
        """All vertices other than the homozygous copy of G."""
        return len(self.succ) - self.count_by_kind().get("hom", 0)


def build_diploid(g: SequenceGraph, dec: Decomposition | None = None) -> Tuple[DiploidGraph, Decomposition]:
    dec = dec or decompose(g)
    gd = DiploidGraph(source=("hom", dec.source), sink=("hom", dec.sink))

    # Homozygous copy of G.
    for v in g.vertices:
        gd.add_vertex(("hom", v))
    for u, v in g.edges():
        gd.add_edge(("hom", u), ("hom", v))

    # P-nodes.
    for p in dec.nodes("P"):
        x, y = p.source, p.sink
        for i, b1 in enumerate(p.children):
            for j, b2 in enumerate(p.children):
                if i != j:
                    h = ("het", p.uid, i, j)
                    gd.add_edge(("hom", x), h)
                    gd.add_edge(h, ("hom", y))
                    gd.gadget[h] = Gadget(Segment(x, y, b1), Segment(x, y, b2))

    # Irreducible core (non-series-parallel inputs only).
    for r in dec.nodes("R"):
        _add_core(gd, r, {v: i for i, v in enumerate(g.topological_order())})
    return gd, dec


def _add_core(gd: DiploidGraph, r: DecompNode, rank: Dict[Vertex, int]) -> None:
    out_core: Dict[Vertex, List[DecompNode]] = {}
    for c in r.children:
        out_core.setdefault(c.source, []).append(c)
        out_core.setdefault(c.sink, [])
    for u, cs in out_core.items():
        heads = [c.sink for c in cs]
        assert len(heads) == len(set(heads)), "core has parallel edges"

    def state(a: Vertex, b: Vertex) -> DVertex:
        return ("hom", a) if a == b else ("pair", a, b)

    queue: deque = deque()

    def link(frm: DVertex, to: DVertex) -> None:
        if to[0] == "pair" and to not in gd.succ:
            queue.append(to)
        gd.add_edge(frm, to)

    for x, cs in out_core.items():
        for c1 in cs:
            for c2 in cs:
                if c1 is c2:
                    continue
                f = ("fork", c1.uid, c2.uid)
                gd.gadget[f] = Gadget(Segment(x, c1.sink, c1), Segment(x, c2.sink, c2))
                gd.add_edge(("hom", x), f)
                link(f, state(c1.sink, c2.sink))

    while queue:
        st = queue.popleft()
        _, u, v = st
        if rank[u] < rank[v]:
            h, lag, other = 1, u, v
        else:
            h, lag, other = 2, v, u
        for c in out_core[lag]:
            step = ("step", h, c.uid, other)
            seg = Segment(lag, c.sink, c)
            gd.gadget[step] = Gadget(seg, None) if h == 1 else Gadget(None, seg)
            gd.add_edge(st, step)
            link(step, state(c.sink, other) if h == 1 else state(other, c.sink))


def diploid_size(g: SequenceGraph, dec: Decomposition | None = None) -> Dict[str, int]:
    """Vertex and edge counts of G_D, by kind, without building it.

    Gives the same numbers as build_diploid (checked in the tests) but stores
    only integer-coded pair states, so it runs on graphs whose core is too large
    to materialize in Python.
    """
    dec = dec or decompose(g)
    het = sum(len(p.children) * (len(p.children) - 1) for p in dec.nodes("P"))
    counts = {"hom": g.num_vertices(), "het": het, "fork": 0, "pair": 0, "step": 0}
    n_edges = g.num_edges() + 2 * het
    rank = {v: i for i, v in enumerate(g.topological_order())}
    for r in dec.nodes("R"):
        order = sorted({c.source for c in r.children} | {c.sink for c in r.children}, key=rank.__getitem__)
        pos = {v: i for i, v in enumerate(order)}
        n = len(order)
        out: List[List[int]] = [[] for _ in range(n)]
        for c in r.children:
            out[pos[c.source]].append(pos[c.sink])
        seen: Set[int] = set()
        queue: deque = deque()
        for heads in out:
            counts["fork"] += len(heads) * (len(heads) - 1)
            for a in heads:
                for b in heads:
                    if a != b and a * n + b not in seen:
                        seen.add(a * n + b)
                        queue.append(a * n + b)
        while queue:
            u, v = divmod(queue.popleft(), n)
            lag, other, first = (u, v, True) if u < v else (v, u, False)
            for w in out[lag]:
                counts["step"] += 1
                if w != other:
                    key = (w * n + other) if first else (other * n + w)
                    if key not in seen:
                        seen.add(key)
                        queue.append(key)
        counts["pair"] += len(seen)
    n_edges += 2 * counts["fork"] + 2 * counts["step"]
    counts["vertices"] = sum(counts[k] for k in ("hom", "het", "fork", "pair", "step"))
    counts["edges"] = n_edges
    return counts


def naive_pair_count(n_haplotypes: int, ordered: bool = False) -> int:
    """Number of vertices a naive construction creates: one per haplotype pair."""
    n = n_haplotypes
    return n * (n - 1) if ordered else n * (n - 1) // 2
