"""Synthetic test graphs."""

from __future__ import annotations

import random
from typing import Tuple

from .seqgraph import SequenceGraph


def random_sp_dag(n_ops: int, p_series: float = 0.3, max_branches: int = 3, seed: int | None = None) -> SequenceGraph:
    """Random two-terminal series-parallel DAG without parallel edges.

    Start from s=0 -> t=1. Each step picks an edge (u, v) and replaces it either
    by a path u -> w -> v (series) or by 2..max_branches parallel two-edge
    paths u -> w_i -> v (parallel bubble).
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
        for _ in range(k):
            w, nxt = nxt, nxt + 1
            edges += [(u, w), (w, v)]
    return SequenceGraph.from_edges(edges)


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
