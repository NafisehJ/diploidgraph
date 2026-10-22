"""Two-terminal DAG sequence graphs.

A *sequence graph* here is a directed acyclic graph G = (V, E) in which every
vertex carries a nucleotide string. It is *two-terminal* if it has exactly one
source s (in-degree 0) and exactly one sink t (out-degree 0); in a DAG this
implies that every vertex lies on at least one s-t path.
"""

from __future__ import annotations

import heapq
from collections import deque
from dataclasses import dataclass, field
from typing import Dict, Hashable, Iterable, Iterator, List, Optional, Set, Tuple

Vertex = Hashable
Edge = Tuple[Vertex, Vertex]


class NotTwoTerminalDAG(ValueError):
    """Raised when an input graph violates the two-terminal DAG assumption."""


@dataclass
class SequenceGraph:
    """A simple directed graph with optional vertex sequences.

    Parallel edges are not allowed (a GFA link set cannot contain them).
    """

    succ: Dict[Vertex, Set[Vertex]] = field(default_factory=dict)
    pred: Dict[Vertex, Set[Vertex]] = field(default_factory=dict)
    seq: Dict[Vertex, str] = field(default_factory=dict)

    # ------------------------------------------------------------------ build
    @classmethod
    def from_edges(cls, edges: Iterable[Edge], seq: Optional[Dict[Vertex, str]] = None) -> "SequenceGraph":
        g = cls()
        for u, v in edges:
            g.add_edge(u, v)
        for v, s in (seq or {}).items():
            g.add_vertex(v, s)
        return g

    def add_vertex(self, v: Vertex, sequence: Optional[str] = None) -> None:
        self.succ.setdefault(v, set())
        self.pred.setdefault(v, set())
        if sequence is not None:
            self.seq[v] = sequence

    def add_edge(self, u: Vertex, v: Vertex) -> None:
        if u == v:
            raise NotTwoTerminalDAG(f"self-loop on {u!r}")
        self.add_vertex(u)
        self.add_vertex(v)
        self.succ[u].add(v)
        self.pred[v].add(u)

    # ------------------------------------------------------------------ query
    @property
    def vertices(self) -> List[Vertex]:
        return list(self.succ)

    def edges(self) -> Iterator[Edge]:
        for u, outs in self.succ.items():
            for v in outs:
                yield (u, v)

    def num_vertices(self) -> int:
        return len(self.succ)

    def num_edges(self) -> int:
        return sum(len(o) for o in self.succ.values())

    def sources(self) -> List[Vertex]:
        return [v for v in self.succ if not self.pred[v]]

    def sinks(self) -> List[Vertex]:
        return [v for v in self.succ if not self.succ[v]]

    def topological_order(self) -> List[Vertex]:
        """Kahn's algorithm with ties broken by vertex insertion order, so the
        result is the same on every run (set iteration order is not).
        Raises NotTwoTerminalDAG if a cycle exists."""
        index = {v: i for i, v in enumerate(self.succ)}
        indeg = {v: len(p) for v, p in self.pred.items()}
        heap = [index[v] for v, d in indeg.items() if d == 0]
        heapq.heapify(heap)
        verts = list(self.succ)
        order: List[Vertex] = []
        while heap:
            u = verts[heapq.heappop(heap)]
            order.append(u)
            for w in self.succ[u]:
                indeg[w] -= 1
                if indeg[w] == 0:
                    heapq.heappush(heap, index[w])
        if len(order) != len(self.succ):
            raise NotTwoTerminalDAG("graph contains a directed cycle")
        return order

    def terminals(self) -> Tuple[Vertex, Vertex]:
        """Return (s, t) after checking the two-terminal DAG assumption."""
        self.topological_order()
        srcs, snks = self.sources(), self.sinks()
        if len(srcs) != 1 or len(snks) != 1:
            raise NotTwoTerminalDAG(f"expected 1 source and 1 sink, found {len(srcs)} and {len(snks)}")
        return srcs[0], snks[0]

    def count_paths(self, s: Vertex, t: Vertex) -> int:
        """Number of directed s-t paths (exact integer, linear time)."""
        count: Dict[Vertex, int] = {v: 0 for v in self.succ}
        count[s] = 1
        for u in self.topological_order():
            if count[u]:
                for w in self.succ[u]:
                    count[w] += count[u]
        return count[t]

    def all_paths(self, s: Vertex, t: Vertex, allowed_edges: Optional[Set[Edge]] = None) -> List[Tuple[Vertex, ...]]:
        """Enumerate all s-t paths. Exponential; for small validation graphs only."""
        out: List[Tuple[Vertex, ...]] = []
        stack: List[Tuple[Vertex, Tuple[Vertex, ...]]] = [(s, (s,))]
        while stack:
            u, path = stack.pop()
            if u == t:
                out.append(path)
                continue
            for w in self.succ[u]:
                if allowed_edges is None or (u, w) in allowed_edges:
                    stack.append((w, path + (w,)))
        return out
