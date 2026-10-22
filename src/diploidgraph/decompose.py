"""Series-parallel decomposition of a two-terminal DAG.

We build the decomposition tree by exhaustive series and parallel reductions
on a directed multigraph (Valdes, Tarjan & Lawler 1982). Each multigraph edge
carries a tree node describing the subgraph it stands for:

* Q-node: one edge of G.
* S-node: a series composition (ordered children, first child starts at the
  node's source, last child ends at its sink).
* P-node: a parallel composition of k >= 2 *branches* between the same pair of
  terminals. Branches are S- or Q-nodes (nested P-nodes are flattened), so the
  terminals of every P-node are distinct.
* R-node: what remains when no reduction applies. For an SP graph there is no
  R-node. For a non-SP graph we currently return a single R-node whose children
  are the irreducible core edges; splitting it into the proper R-components of
  the SPQR-tree is the next step [TODO: non-SP case].

For a two-terminal DAG, a P-node produced here is a maximal set of >= 2
internally disjoint two-terminal subgraphs sharing terminals (x, y); these are
the P-node terminal pairs ("hetero-positions") used by the diploid
construction.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from itertools import count
from typing import Dict, List, Optional, Set, Tuple

from .seqgraph import Edge, SequenceGraph, Vertex


@dataclass
class DecompNode:
    kind: str  # "Q", "S", "P" or "R"
    source: Vertex
    sink: Vertex
    children: List["DecompNode"] = field(default_factory=list)
    uid: int = -1

    # Filled by annotate():
    n_paths: int = 0  # number of source-sink paths inside this subgraph

    def edges(self) -> Set[Edge]:
        """Edges of G inside this subgraph (computed on demand, O(size))."""
        return {(q.source, q.sink) for q in self.walk() if q.kind == "Q"}

    def walk(self):
        """Post-order traversal (children before parent)."""
        stack: List[Tuple["DecompNode", bool]] = [(self, False)]
        while stack:
            node, done = stack.pop()
            if done:
                yield node
            else:
                stack.append((node, True))
                for c in reversed(node.children):
                    stack.append((c, False))

    def __repr__(self) -> str:  # compact, avoids recursion into children
        return f"{self.kind}#{self.uid}({self.source!r}->{self.sink!r}, {len(self.children)} children)"


@dataclass
class Decomposition:
    root: DecompNode
    source: Vertex
    sink: Vertex
    is_series_parallel: bool

    def nodes(self, kind: Optional[str] = None) -> List[DecompNode]:
        return [n for n in self.root.walk() if kind is None or n.kind == kind]


def decompose(g: SequenceGraph) -> Decomposition:
    s, t = g.terminals()
    ids = count()

    # Multigraph state: edge id -> (tail, head, node)
    edge_tail: Dict[int, Vertex] = {}
    edge_head: Dict[int, Vertex] = {}
    edge_node: Dict[int, DecompNode] = {}
    out_e: Dict[Vertex, Set[int]] = {v: set() for v in g.vertices}
    in_e: Dict[Vertex, Set[int]] = {v: set() for v in g.vertices}

    def new_edge(u: Vertex, v: Vertex, node: DecompNode) -> int:
        e = next(ids)
        edge_tail[e], edge_head[e], edge_node[e] = u, v, node
        out_e[u].add(e)
        in_e[v].add(e)
        return e

    def drop_edge(e: int) -> DecompNode:
        out_e[edge_tail[e]].discard(e)
        in_e[edge_head[e]].discard(e)
        del edge_tail[e], edge_head[e]
        return edge_node.pop(e)

    for u, v in g.edges():
        new_edge(u, v, DecompNode("Q", u, v))

    def parallel_at(u: Vertex) -> Set[Vertex]:
        """Merge parallel out-edges of u. Returns vertices to revisit."""
        by_head: Dict[Vertex, List[int]] = {}
        for e in out_e[u]:
            by_head.setdefault(edge_head[e], []).append(e)
        touched: Set[Vertex] = set()
        for v, es in by_head.items():
            if len(es) < 2:
                continue
            branches: List[DecompNode] = []
            for e in es:
                n = drop_edge(e)
                branches.extend(n.children if n.kind == "P" else [n])
            new_edge(u, v, DecompNode("P", u, v, branches))
            touched.update((u, v))
        return touched

    def series_at(w: Vertex) -> Set[Vertex]:
        if w in (s, t) or len(in_e[w]) != 1 or len(out_e[w]) != 1:
            return set()
        (e1,), (e2,) = in_e[w], out_e[w]
        u, v = edge_tail[e1], edge_head[e2]
        a, b = drop_edge(e1), drop_edge(e2)
        parts: List[DecompNode] = []
        for n in (a, b):
            parts.extend(n.children if n.kind == "S" else [n])
        new_edge(u, v, DecompNode("S", u, v, parts))
        return {u, v}

    work = list(g.vertices)
    queued = set(work)
    while work:
        w = work.pop()
        queued.discard(w)
        if w not in out_e:
            continue
        touched = parallel_at(w) | series_at(w)
        # A reduction only changes degrees / parallel classes at its endpoints.
        for x in touched:
            if x not in queued:
                queued.add(x)
                work.append(x)

    remaining = list(edge_node)
    if len(remaining) == 1 and edge_tail[remaining[0]] == s and edge_head[remaining[0]] == t:
        root, is_sp = edge_node[remaining[0]], True
    else:
        root = DecompNode("R", s, t, [edge_node[e] for e in sorted(remaining)])
        is_sp = False

    uid = count()
    for n in root.walk():
        n.uid = next(uid)
    annotate(root)
    return Decomposition(root=root, source=s, sink=t, is_series_parallel=is_sp)


def annotate(root: DecompNode) -> None:
    """Fill n_paths bottom-up (exact integers)."""
    for n in root.walk():
        if n.kind == "Q":
            n.n_paths = 1
        elif n.kind == "S":
            n.n_paths = 1
            for c in n.children:
                n.n_paths *= c.n_paths
        elif n.kind == "P":
            n.n_paths = sum(c.n_paths for c in n.children)
        else:  # R: count paths in the core DAG with edge multiplicities
            n.n_paths = _count_core_paths(n)


def _count_core_paths(r: DecompNode) -> int:
    succ: Dict[Vertex, List[DecompNode]] = {}
    indeg: Dict[Vertex, int] = {}
    for c in r.children:
        succ.setdefault(c.source, []).append(c)
        succ.setdefault(c.sink, [])
        indeg[c.sink] = indeg.get(c.sink, 0) + 1
        indeg.setdefault(c.source, 0)
    ways = {v: 0 for v in succ}
    ways[r.source] = 1
    frontier = [v for v, d in indeg.items() if d == 0]
    while frontier:
        u = frontier.pop()
        for c in succ[u]:
            ways[c.sink] += ways[u] * c.n_paths
            indeg[c.sink] -= 1
            if indeg[c.sink] == 0:
                frontier.append(c.sink)
    return ways[r.sink]
