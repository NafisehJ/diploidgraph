"""Synthetic test graphs."""

from __future__ import annotations

import random
from typing import Tuple

from .seqgraph import SequenceGraph


def random_sp_dag(n_ops: int, p_series: float = 0.3, max_branches: int = 3, seed: int | None = None,
                  p_deletion: float = 0.0) -> SequenceGraph:
    """Random two-terminal series-parallel DAG without parallel edges.

    Start from s=0 -> t=1. Each step picks an edge (u, v) and replaces it either
    by a path u -> w -> v (series) or by 2..max_branches parallel two-edge
    paths u -> w_i -> v (parallel bubble). With probability p_deletion a bubble
    keeps the direct edge u -> v as one of its branches (a deletion allele).
    """
    rng = random.Random(seed)
    edges = [(0, 1)]  # list + swap-remove gives O(1) random choice, reproducible by seed
    nxt = 2
    for _ in range(n_ops):
        i = rng.randrange(len(edges))
        u, v = edges[i]
        edges[i] = edges[-1]
        edges.pop()
        k = 1 if rng.random() < p_series else rng.randint(2, max_branches)
        if k > 1 and rng.random() < p_deletion:
            edges.append((u, v))  # direct edge as one branch
            k -= 1
        for _ in range(k):
            w, nxt = nxt, nxt + 1
            edges += [(u, w), (w, v)]
    return SequenceGraph.from_edges(edges)


def random_dag(n_ops: int, extra_edges: int, max_branches: int = 3, seed: int | None = None) -> SequenceGraph:
    """Random two-terminal DAG, usually not series-parallel.

    A random SP DAG plus `extra_edges` random forward edges (u before v in a
    topological order, u != sink, v != source), which keeps one source and one
    sink and creates W-motifs.
    """
    rng = random.Random(seed)
    g = random_sp_dag(n_ops, max_branches=max_branches, seed=rng.randrange(1 << 30))
    s, t = g.terminals()
    order = g.topological_order()
    added = 0
    for _ in range(50 * extra_edges + 50):
        if added == extra_edges:
            break
        i, j = sorted(rng.sample(range(len(order)), 2))
        u, v = order[i], order[j]
        if u == t or v == s or v in g.succ[u]:
            continue
        g.add_edge(u, v)
        added += 1
    return g


def w_graph() -> SequenceGraph:
    """The smallest non-SP two-terminal DAG (the W / Wheatstone-bridge pattern).

    s=u, t=v, diagonal vertices a=u', b=v':
    u->a, u->b, a->b, a->v, b->v.
    """
    return SequenceGraph.from_edges([("u", "a"), ("u", "b"), ("a", "b"), ("a", "v"), ("b", "v")])


def series_of(*graphs: Tuple[SequenceGraph, ...]) -> SequenceGraph:
    """Glue two-terminal graphs in series (sink of one = source of next)."""
    out = SequenceGraph()
    prev_sink = None
    for idx, g in enumerate(graphs):
        s, t = g.terminals()
        ren = {v: (idx, v) for v in g.vertices}
        if prev_sink is not None:
            ren[s] = prev_sink
        for a, b in g.edges():
            out.add_edge(ren[a], ren[b])
        prev_sink = ren[t]
    return out


def dense_random_dag(n: int, p: float, seed: int | None = None) -> SequenceGraph:
    """Random two-terminal DAG on vertices 0..n-1 (edges i -> j, i < j, with probability p).

    Every vertex other than 0 gets an in-edge and every vertex other than n-1 an
    out-edge, so 0 is the only source and n-1 the only sink. Unlike the other
    generators this is not built from series-parallel pieces.
    """
    rng = random.Random(seed)
    edges = {(i, j) for i in range(n) for j in range(i + 1, n) if rng.random() < p}
    for j in range(1, n):
        if not any((i, j) in edges for i in range(j)):
            edges.add((rng.randrange(j), j))
    for i in range(n - 1):
        if not any((i, j) in edges for j in range(i + 1, n)):
            edges.add((i, rng.randrange(i + 1, n)))
    return SequenceGraph.from_edges(sorted(edges))


def shuffled(g: SequenceGraph, seed: int | None = None) -> SequenceGraph:
    """The same graph with vertices and edges inserted in a random order.

    The construction breaks topological-order ties by insertion order, so this
    changes which haplotype moves first in the non-series-parallel core.
    """
    rng = random.Random(seed)
    verts = list(g.vertices)
    rng.shuffle(verts)
    edges = list(g.edges())
    rng.shuffle(edges)
    out = SequenceGraph()
    for v in verts:
        out.add_vertex(v, g.seq.get(v))
    for u, v in edges:
        out.add_edge(u, v)
    return out
