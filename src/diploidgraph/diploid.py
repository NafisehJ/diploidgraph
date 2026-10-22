"""Diploid sequence graph construction (series-parallel case).

Vertices of the diploid graph G_D:

* homozygous vertex ("hom", v) for every v in V(G), meaning both haplotypes
  pass through v;
* heterozygous vertex ("het", p, i, j) for every P-node p with terminals
  (x, y) and every ordered pair of distinct branches i != j of p, meaning
  haplotype 1 goes from x to y through branch i and haplotype 2 through
  branch j (any path inside each branch).

Edges: ("hom", u) -> ("hom", v) for every (u, v) in E(G), and
("hom", x) -> ("het", p, i, j) -> ("hom", y) for every heterozygous vertex.

A P-node with k branches therefore adds k(k-1) vertices and 2k(k-1) edges, so
|V(G_D)| = |V| + sum_p k_p(k_p - 1) and |E(G_D)| = |E| + 2 sum_p k_p(k_p - 1).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, Hashable, List, Set, Tuple

from .decompose import DecompNode, Decomposition, decompose
from .seqgraph import SequenceGraph, Vertex

DVertex = Tuple[Hashable, ...]


class NonSeriesParallelError(NotImplementedError):
    """The input has an R-node; the W-motif extension is not implemented yet."""


@dataclass
class HetInfo:
    pnode: DecompNode
    branch1: DecompNode
    branch2: DecompNode


@dataclass
class DiploidGraph:
    source: DVertex
    sink: DVertex
    succ: Dict[DVertex, Set[DVertex]] = field(default_factory=dict)
    het: Dict[DVertex, HetInfo] = field(default_factory=dict)

    def add_edge(self, a: DVertex, b: DVertex) -> None:
        self.succ.setdefault(a, set()).add(b)
        self.succ.setdefault(b, set())

    def num_vertices(self) -> int:
        return len(self.succ)

    def num_edges(self) -> int:
        return sum(len(o) for o in self.succ.values())

    def num_het_vertices(self) -> int:
        return len(self.het)


def build_diploid(g: SequenceGraph, dec: Decomposition | None = None, allow_partial: bool = False) -> Tuple[DiploidGraph, Decomposition]:
    """Build G_D from a two-terminal DAG.

    If the graph is not series-parallel, raises NonSeriesParallelError unless
    allow_partial=True, in which case R-node regions are copied homozygously
    only (pairs that diverge inside them are NOT represented; use
    validate.coverage() to measure what fraction of pairs is represented).
    """
    dec = dec or decompose(g)
    if not dec.is_series_parallel and not allow_partial:
        raise NonSeriesParallelError("non-series-parallel input: W-motif extension not implemented yet")

    gd = DiploidGraph(source=("hom", dec.source), sink=("hom", dec.sink))
    for v in g.vertices:
        gd.succ.setdefault(("hom", v), set())
    for u, v in g.edges():
        gd.add_edge(("hom", u), ("hom", v))

    for p in dec.nodes("P"):
        x, y = ("hom", p.source), ("hom", p.sink)
        for i, b1 in enumerate(p.children):
            for j, b2 in enumerate(p.children):
                if i == j:
                    continue
                h = ("het", p.uid, i, j)
                gd.add_edge(x, h)
                gd.add_edge(h, y)
                gd.het[h] = HetInfo(p, b1, b2)
    return gd, dec


def naive_pair_count(n_haplotypes: int, ordered: bool = False) -> int:
    """Number of vertices a naive construction creates: one per haplotype pair."""
    n = n_haplotypes
    return n * (n - 1) if ordered else n * (n - 1) // 2
