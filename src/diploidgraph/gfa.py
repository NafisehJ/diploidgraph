"""Minimal GFA 1.x reader for DAG pangenome subgraphs.

Reads S (segments), L (links), P and W (haplotype paths). Links with mixed
orientation (+/- or -/+) cannot be represented in a directed graph and are
counted and dropped; the caller decides whether that is acceptable.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from typing import List, Tuple

from .seqgraph import SequenceGraph

_W_STEP = re.compile(r"([<>])([^<>]+)")


@dataclass
class GFAData:
    graph: SequenceGraph
    haplotypes: List[Tuple[str, Tuple[str, ...]]] = field(default_factory=list)  # (name, forward vertex walk)
    dropped_links: int = 0
    unoriented_haplotypes: int = 0  # walks that mix orientations; kept out of `haplotypes`


def _orient(steps: List[Tuple[str, bool]]) -> Tuple[str, ...] | None:
    """steps = [(segment, is_forward)]. Return a forward walk or None if mixed."""
    if all(f for _, f in steps):
        return tuple(s for s, _ in steps)
    if not any(f for _, f in steps):
        return tuple(s for s, _ in reversed(steps))
    return None


def read_gfa(path: str) -> GFAData:
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
                if walk is None:
                    data.unoriented_haplotypes += 1
                else:
                    data.haplotypes.append((name, walk))
    return data
