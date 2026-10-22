"""Which new vertices of G_D are used by pairs of observed haplotypes?

Given observed haplotypes (source-to-sink walks of G), a new vertex of G_D is
*supported* if the walk of G_D for some ordered pair (h1, h2) of distinct
observed haplotypes passes through it. The other new vertices represent pairs
that arise only by recombination between the observed haplotypes.

The computation streams over the walks and never builds G_D:

* P-node vertex ("het", p, i, j) lies on the walk of (h1, h2) exactly when both
  haplotypes cross p, h1 through branch i and h2 through branch j. It is
  therefore supported iff i != j are both branches of p used by some observed
  haplotype.
* Core vertices: a *cut* is a core vertex on every source-sink path, so any two
  haplotypes are together at every cut. Each walk is split at the cuts into
  sub-walks; for each block between consecutive cuts we keep the distinct
  sub-walks and replay every ordered pair of them with the same rule as
  diploid.build_diploid.
"""

from __future__ import annotations

import hashlib
from collections import defaultdict
from typing import Dict, Hashable, Iterable, List, Sequence, Set, Tuple

from .decompose import DecompNode, Decomposition, decompose
from .seqgraph import SequenceGraph, Vertex

Key = Tuple[Hashable, ...]
KINDS = ("het", "fork", "pair", "step")


def core_cuts(core_edges: Sequence[DecompNode], rank: Dict[Vertex, int]) -> Set[Vertex]:
    """Core vertices that every source-sink path passes through (source and sink included)."""
    verts = sorted({c.source for c in core_edges} | {c.sink for c in core_edges}, key=rank.__getitem__)
    pos = {v: i for i, v in enumerate(verts)}
    cover = [0] * (len(verts) + 1)
    for c in core_edges:
        a, b = pos[c.source], pos[c.sink]
        if b - a > 1:
            cover[a + 1] += 1
            cover[b] -= 1
    cuts, run = set(), 0
    for i, v in enumerate(verts):
        run += cover[i]
        if run == 0:
            cuts.add(v)
    return cuts


def observed_support(g: SequenceGraph, walks: Iterable[Sequence[Vertex]], dec: Decomposition | None = None
                     ) -> Tuple[Dict[str, Set[Key]], Dict[str, int]]:
    """Return (supported new vertices by kind, walk counts).

    `walks` may be any iterable of vertex sequences, e.g. a generator over a GFA
    file. Only distinct walks from the source to the sink of G are used.
    """
    dec = dec or decompose(g)
    s, t = dec.source, dec.sink
    rank = {v: i for i, v in enumerate(g.topological_order())}

    parent: Dict[int, Tuple[DecompNode, int]] = {}
    leaf: Dict[Tuple[Vertex, Vertex], DecompNode] = {}
    for n in dec.root.walk():
        if n.kind == "Q":
            leaf[(n.source, n.sink)] = n
        for i, c in enumerate(n.children):
            parent[id(c)] = (n, i)

    core_edges = [c for r in dec.nodes("R") for c in r.children]
    edge_uid = {(c.source, c.sink): c.uid for c in core_edges}
    core = {c.source for c in core_edges} | {c.sink for c in core_edges}
    cuts = core_cuts(core_edges, rank) if core_edges else set()

    used_branches: Dict[int, Set[int]] = defaultdict(set)
    blocks: Dict[Vertex, Set[Tuple[Vertex, ...]]] = defaultdict(set)  # keyed by the opening cut
    seen: Set[bytes] = set()
    counts = {"walks": 0, "full": 0, "distinct_full": 0}

    for walk in walks:
        counts["walks"] += 1
        if not walk or walk[0] != s or walk[-1] != t:
            continue
        counts["full"] += 1
        digest = hashlib.blake2b("\t".join(map(str, walk)).encode(), digest_size=16).digest()
        if digest in seen:
            continue
        seen.add(digest)
        counts["distinct_full"] += 1

        crossed: Set[int] = set()
        for a, b in zip(walk, walk[1:]):
            node = leaf[(a, b)]
            while id(node) in parent:
                par, idx = parent[id(node)]
                if par.kind == "P":
                    if par.uid in crossed:
                        break  # this P-node and all its ancestors are already recorded
                    crossed.add(par.uid)
                    used_branches[par.uid].add(idx)
                node = par

        if core:
            cur: List[Vertex] = []
            for v in walk:
                if v not in core:
                    continue
                cur.append(v)
                if v in cuts:
                    if len(cur) > 1:
                        blocks[cur[0]].add(tuple(cur))
                    cur = [v]

    supported: Dict[str, Set[Key]] = {k: set() for k in KINDS}
    for uid, branches in used_branches.items():
        for i in branches:
            for j in branches:
                if i != j:
                    supported["het"].add(("het", uid, i, j))
    for subs in blocks.values():
        subs = sorted(subs)
        for A in subs:
            for B in subs:
                if A != B:
                    _replay(A, B, rank, edge_uid, supported)
    return supported, counts


def _replay(A: Tuple[Vertex, ...], B: Tuple[Vertex, ...], rank, edge_uid, supported) -> None:
    """Follow the G_D walk of the haplotype pair (A, B) through one block (same rule as build_diploid)."""
    i = j = 0
    last = len(A) - 1
    while i < last:
        if A[i + 1] == B[j + 1]:  # same core edge: homozygous copy
            i += 1
            j += 1
            continue
        supported["fork"].add(("fork", edge_uid[(A[i], A[i + 1])], edge_uid[(B[j], B[j + 1])]))
        i += 1
        j += 1
        while A[i] != B[j]:
            supported["pair"].add(("pair", A[i], B[j]))
            if rank[A[i]] < rank[B[j]]:
                supported["step"].add(("step", 1, edge_uid[(A[i], A[i + 1])], B[j]))
                i += 1
            else:
                supported["step"].add(("step", 2, edge_uid[(B[j], B[j + 1])], A[i]))
                j += 1
