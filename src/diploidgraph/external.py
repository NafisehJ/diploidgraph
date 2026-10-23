"""State-space sizes of two external genotyping models, computed from the same GFA.

This module does not run PanGenie or vg call. It counts, on the input graph and
its haplotype walks, the objects that each model's genotyping step works with,
as described in the methods papers:

* PanGenie (Ebler et al., Nat. Genet. 54:518-525, 2022, Methods, "Hidden states
  and transitions"): for each bubble v and N given haplotype paths, the hidden
  states are H_{v,i,j} for all i, j <= N, i.e. N^2 ordered pairs per bubble.
  Bubbles less than the k-mer size apart (k = 31) are combined into one bubble,
  whose alleles are the haplotype paths through the combined region. Between
  consecutive bubbles every state can move to every state (Li-Stephens
  transitions), so a decoded diploid genome is an ordered pair of mosaics of
  the panel haplotypes that switch only between bubbles.

* vg call (Hickey et al., Genome Biol. 21:35, 2020, Methods, "SV genotyping
  algorithm"): each top-level snarl is genotyped on its own from a set of
  candidate traversals; the two most supported traversals give the genotype.
  No state links one snarl to the next.

Proxies used here (state them in the paper):

* Top-level bubbles. In a two-terminal DAG, a vertex lies on every
  source-to-sink path iff no edge jumps over its position in a topological
  order (every vertex is assumed to lie on some source-to-sink path). These
  "chain vertices" c_0 = s, ..., c_m = t cut the graph into consecutive
  intervals; an interval with at least one interior vertex is a top-level
  bubble. They play the role of PanGenie's bubbles and vg's top-level snarls.
* Panel haplotypes. The P/W walks that run from the source to the sink, except
  the reference genomes (samples GRCh38 and CHM13 by default): PanGenie's panel
  consists of the assembly haplotypes. vg call candidates are counted from all
  complete walks, the references included.
* Alleles of a bubble. The distinct sub-walks of the panel haplotypes through it
  (what PanGenie's bubble VCF and a haplotype-based traversal finder would list).
  A bubble with a single allele carries no variation among the panel
  haplotypes and is not a PanGenie bubble.
"""

from __future__ import annotations

import hashlib
import math
from dataclasses import dataclass
from typing import Dict, List, Optional, Sequence, Tuple

from .gfa import iter_gfa_walks, read_gfa
from .seqgraph import SequenceGraph, Vertex
from .stats import log10_int


# --------------------------------------------------------------------------- chain
def chain_vertices(g: SequenceGraph) -> List[Vertex]:
    """Vertices on every source-to-sink path, in topological order.

    Assumes g is a two-terminal DAG (every vertex on some s-t path). Linear time:
    position i is a chain vertex iff no edge (a, b) has pos(a) < i < pos(b).
    """
    g.terminals()
    order = g.topological_order()
    pos = {v: i for i, v in enumerate(order)}
    diff = [0] * (len(order) + 1)
    for a, outs in g.succ.items():
        pa = pos[a]
        for b in outs:
            pb = pos[b]
            if pb - pa >= 2:
                diff[pa + 1] += 1
                diff[pb] -= 1
    out: List[Vertex] = []
    running = 0
    for i, v in enumerate(order):
        running += diff[i]
        if running == 0:
            out.append(v)
    return out


@dataclass
class Bubble:
    start: Vertex  # chain vertex where the bubble opens
    end: Vertex  # next chain vertex, where it closes
    start_idx: int  # index of `start` in the chain


def top_level_bubbles(g: SequenceGraph, chain: Sequence[Vertex]) -> List[Bubble]:
    """Intervals between consecutive chain vertices that contain an interior vertex."""
    order = g.topological_order()
    pos = {v: i for i, v in enumerate(order)}
    return [
        Bubble(chain[i], chain[i + 1], i)
        for i in range(len(chain) - 1)
        if pos[chain[i + 1]] - pos[chain[i]] > 1
    ]


# --------------------------------------------------------------------------- walks
@dataclass
class AlleleTable:
    """Allele id of every complete haplotype walk at every top-level bubble."""

    bubbles: List[Bubble]
    alleles: List[List[int]]  # one row per complete walk: allele id at each bubble
    n_alleles: List[int]  # distinct alleles per bubble
    n_distinct_walks: int
    chain_bp: Optional[List[int]]  # start coordinate (bp) of each chain vertex on the walks; None if sequences missing
    chain_len: Optional[List[int]]  # sequence length of each chain vertex


def _digest(parts: Sequence[str]) -> bytes:
    return hashlib.blake2b("\t".join(parts).encode(), digest_size=16).digest()


def allele_table(path: str, g: SequenceGraph, chain: Sequence[Vertex], bubbles: Sequence[Bubble],
                 exclude_samples: Sequence[str] = ()) -> AlleleTable:
    """Stream the GFA walks once; keep only complete (source-to-sink) walks,
    skipping walks of the samples in exclude_samples."""
    s, t = chain[0], chain[-1]
    cidx = {v: i for i, v in enumerate(chain)}
    by_start = {b.start_idx: k for k, b in enumerate(bubbles)}
    ids: List[Dict[bytes, int]] = [dict() for _ in bubbles]
    rows: List[List[int]] = []
    walk_digests = set()
    has_seq = all(v in g.seq for v in g.succ)
    chain_bp: Optional[List[int]] = None
    chain_len = [len(g.seq[v]) for v in chain] if has_seq else None

    for walk in iter_gfa_walks(path, exclude_samples):
        if walk[0] != s or walk[-1] != t:
            continue
        walk_digests.add(_digest([str(v) for v in walk]))
        at = [-1] * len(chain)  # index in the walk of each chain vertex
        for i, v in enumerate(walk):
            c = cidx.get(v)
            if c is not None:
                at[c] = i
        if min(at) < 0:
            raise ValueError("a complete walk misses a chain vertex; is every walk a path of the graph?")
        row = [0] * len(bubbles)
        for c, k in by_start.items():
            sub = walk[at[c] : at[c + 1] + 1]
            d = _digest([str(v) for v in sub])
            row[k] = ids[k].setdefault(d, len(ids[k]))
        rows.append(row)
        if has_seq and chain_bp is None:
            bp, acc = {}, 0
            for v in walk:
                if v in cidx:
                    bp[cidx[v]] = acc
                acc += len(g.seq[v])
            chain_bp = [bp[i] for i in range(len(chain))]

    if not rows:
        raise ValueError("no complete source-to-sink haplotype walks in the GFA")
    n_alleles = [len(d) for d in ids]
    return AlleleTable(list(bubbles), rows, n_alleles, len(walk_digests), chain_bp, chain_len)


# --------------------------------------------------------------------------- PanGenie
def pangenie_groups(tab: AlleleTable, k: int = 31) -> Optional[List[List[int]]]:
    """Combine variant bubbles (>= 2 alleles) that are less than k bp apart.

    The gap between consecutive variant bubbles is the invariant sequence
    between them: from the start of the closing chain vertex of the first to
    the end of the opening chain vertex of the second. Returns lists of bubble
    indices, or None when the GFA has no segment sequences.
    """
    var = [i for i, a in enumerate(tab.n_alleles) if a >= 2]
    if tab.chain_bp is None or tab.chain_len is None:
        return None
    groups: List[List[int]] = []
    for i in var:
        if groups:
            prev = tab.bubbles[groups[-1][-1]]
            cur = tab.bubbles[i]
            e, s = prev.start_idx + 1, cur.start_idx  # chain indices
            gap = tab.chain_bp[s] + tab.chain_len[s] - tab.chain_bp[e]
            if gap < k:
                groups[-1].append(i)
                continue
        groups.append([i])
    return groups


def _distinct_rows(rows: List[List[int]], cols: List[int]) -> int:
    """Number of distinct allele combinations of the walks over the given bubbles."""
    return len({tuple(r[c] for c in cols) for r in rows})


REFERENCE_SAMPLES = ("GRCh38", "CHM13")


def external_sizes(path: str, k: int = 31, reference_samples: Sequence[str] = REFERENCE_SAMPLES) -> Dict[str, object]:
    """All counts for one GFA, as a flat dict (see module docstring)."""
    data = read_gfa(path)
    g = data.graph
    s, t = g.terminals()
    chain = chain_vertices(g)
    bubbles = top_level_bubbles(g, chain)
    tab = allele_table(path, g, chain, bubbles)

    n_walks = len(tab.alleles)
    n_paths = g.count_paths(s, t)
    a = tab.n_alleles
    var = [x for x in a if x >= 2]

    out: Dict[str, object] = {
        "input": path,
        "n_vertices": g.num_vertices(),
        "n_complete_walks": n_walks,
        "n_distinct_complete_walks": tab.n_distinct_walks,
        "n_chain_vertices": len(chain),
        "n_top_level_bubbles": len(bubbles),
        "n_variant_bubbles": len(var),
        "alleles_per_variant_bubble_max": max(var, default=0),
        "alleles_per_variant_bubble_mean": (sum(var) / len(var)) if var else 0.0,
        # this paper: every ordered pair of s-t paths, once
        "log10_st_paths": log10_int(n_paths),
        "log10_ours_covered_ordered_pairs": 2 * log10_int(n_paths),
    }

    # vg call: independent per top-level snarl; candidates = haplotype traversals
    out.update(
        vgcall_sites=len(bubbles),
        vgcall_candidate_traversals=sum(a),
        vgcall_unordered_genotypes=sum(x * (x + 1) // 2 for x in a),
        log10_vgcall_unphased_genotype_vectors=sum(math.log10(x * (x + 1) // 2) for x in a),
    )

    # PanGenie: N^2 ordered states per (combined) bubble; panel = complete walks
    # of the assembly haplotypes (reference samples excluded)
    tab = allele_table(path, g, chain, bubbles, exclude_samples=reference_samples) if reference_samples else tab
    groups = pangenie_groups(tab, k)
    N, Nd = len(tab.alleles), tab.n_distinct_walks
    out.update(
        pangenie_reference_samples_excluded=",".join(reference_samples),
        pangenie_panel_haplotypes=N,
        pangenie_panel_distinct_haplotypes=Nd,
    )
    if groups is None:
        out.update(pangenie_k=k, pangenie_bubbles=None, pangenie_note="no segment sequences; bubbles not combined")
        return out
    group_alleles = [_distinct_rows(tab.alleles, gi) for gi in groups]
    M = len(groups)
    out.update(
        pangenie_k=k,
        pangenie_bubbles=M,
        pangenie_states=M * N * N,
        pangenie_states_if_duplicates_collapsed=M * Nd * Nd,
        pangenie_transitions=(M - 1) * N**4 + 2 * N * N if M else 0,
        pangenie_alleles_per_bubble_max=max(group_alleles, default=0),
        log10_pangenie_covered_ordered_pairs=2 * sum(math.log10(x) for x in group_alleles),
    )
    return out
