"""Minimal GFA 1.x reader for DAG pangenome subgraphs.

Reads S (segments), L (links), P and W (haplotype paths). Links with mixed
orientation (+/- or -/+) cannot be represented in a directed graph and are
counted and dropped; the caller decides whether that is acceptable.

Haplotype walks can be long (hundreds of thousands of steps on locus-level
graphs), so by default only a summary of each walk is kept: its length, first
and last vertex, and a digest for counting distinct walks. Pass
keep_walks=True to keep the full vertex sequences.
"""

from __future__ import annotations

import hashlib
import re
from dataclasses import dataclass, field
from typing import List, Optional, Tuple

from .seqgraph import SequenceGraph

_W_STEP = re.compile(r"([<>])([^<>]+)")


@dataclass
class WalkSummary:
    name: str
    length: int  # number of vertices
    first: str
    last: str
    digest: bytes


@dataclass
class GFAData:
    graph: SequenceGraph
    walks: List[WalkSummary] = field(default_factory=list)
    haplotypes: List[Tuple[str, Tuple[str, ...]]] = field(default_factory=list)  # only with keep_walks=True
    dropped_links: int = 0
    unoriented_haplotypes: int = 0  # walks that mix orientations; not summarized

    def n_distinct_walks(self) -> int:
        return len({w.digest for w in self.walks})

    def n_full_walks(self) -> int:
        """Walks that run from the graph's source to its sink."""
        s, t = self.graph.terminals()
        return sum(1 for w in self.walks if w.first == s and w.last == t)


def _orient(steps: List[Tuple[str, bool]]) -> Optional[List[str]]:
    """steps = [(segment, is_forward)]. Return a forward walk or None if mixed."""
    if all(f for _, f in steps):
        return [s for s, _ in steps]
    if not any(f for _, f in steps):
        return [s for s, _ in reversed(steps)]
    return None


def read_gfa(path: str, keep_walks: bool = False) -> GFAData:
    g = SequenceGraph()
    data = GFAData(graph=g)
    with open(path) as fh:
        for line in fh:
            f = line.rstrip("\n").split("\t")
            tag = f[0]
            if tag == "S":
                g.add_vertex(f[1], f[2] if len(f) > 2 and f[2] != "*" else None)
            elif tag == "L":
                a, oa, b, ob = f[1], f[2], f[3], f[4]
                if oa == "+" and ob == "+":
                    g.add_edge(a, b)
                elif oa == "-" and ob == "-":
                    g.add_edge(b, a)
                else:
                    data.dropped_links += 1
            elif tag in ("P", "W"):
                if tag == "P":
                    name = f[1]
                    steps = [(x[:-1], x[-1] == "+") for x in f[2].split(",") if x]
                else:
                    name = "#".join(f[1:4])
                    steps = [(m.group(2), m.group(1) == ">") for m in _W_STEP.finditer(f[6])]
                walk = _orient(steps)
                del steps
                if walk is None or not walk:
                    data.unoriented_haplotypes += 1
                    continue
                digest = hashlib.blake2b("\t".join(walk).encode(), digest_size=16).digest()
                data.walks.append(WalkSummary(name, len(walk), walk[0], walk[-1], digest))
                if keep_walks:
                    data.haplotypes.append((name, tuple(walk)))
    return data
