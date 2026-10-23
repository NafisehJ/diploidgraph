"""Export sequence graphs and diploid graphs to Graphviz DOT for figures."""

from __future__ import annotations

from typing import Dict

from .decompose import DecompNode
from .diploid import DiploidGraph
from .seqgraph import SequenceGraph

HOM_COLOR = "#cfe3f7"
HET_COLOR = "#f9c784"
CORE_COLOR = "#b9e4c9"
PAIR_COLOR = "#e3e3e3"


def _q(x) -> str:
    return '"' + str(x).replace('"', r"\"") + '"'


def _expr(n: DecompNode) -> str:
    """Inner vertices of a two-terminal subgraph as a pattern, e.g. a(b|c)d."""
    if n.kind == "Q":
        return ""
    if n.kind == "S":
        parts = []
        for k, c in enumerate(n.children):
            parts.append(_expr(c))
            if k < len(n.children) - 1:
                parts.append(str(c.sink))
        return "".join(parts)
    if n.kind == "P":
        return "(" + "|".join(_expr(c) or "-" for c in n.children) + ")"
    return f"R{n.uid}"


def branch_label(b: DecompNode) -> str:
    """Label for a P-node branch: its inner vertices, with nested alternatives as (x|y); '-' for a direct edge."""
    return _expr(b) or "-"


def sequence_graph_dot(g: SequenceGraph, name: str = "G") -> str:
    lines = [f"digraph {name} {{", "  rankdir=LR;", f'  node [shape=circle, style=filled, fillcolor="{HOM_COLOR}"];']
    for v in g.vertices:
        lines.append(f"  {_q(v)};")
    for u, v in sorted(g.edges(), key=str):  # stable output across runs
        lines.append(f"  {_q(u)} -> {_q(v)};")
    lines.append("}")
    return "\n".join(lines)


def _seg_label(seg, with_sink: bool) -> str:
    if seg is None:
        return "\u00b7"  # middle dot: this haplotype does not move
    inner = _expr(seg.node)
    return (inner + str(seg.sink)) if with_sink else (inner or "-")


def diploid_graph_dot(gd: DiploidGraph, name: str = "GD", compact_steps: bool = True) -> str:
    """Blue ellipses: homozygous copy (v,v). Orange boxes: P-node pairs of branches.
    Green boxes: core forks (the two haplotypes separate). Grey diamonds: core pair
    states [hap1 vertex | hap2 vertex].

    With compact_steps (default), a step vertex (one haplotype moves along a core
    edge while the other waits) is drawn as a labelled arrow instead of a box:
    "h1 → t" means haplotype 1 moves to t.
    """
    lines = [f"digraph {name} {{", "  rankdir=LR;", "  node [style=filled];",
             '  edge [fontsize=10, fontcolor="#52514e"];']
    ids: Dict = {}
    hidden = set()
    for k, v in enumerate(gd.succ):
        ids[v] = f"n{k}"
        kind = v[0]
        if kind == "step" and compact_steps:
            hidden.add(v)
            continue
        if kind == "hom":
            label, shape, color = f"({v[1]},{v[1]})", "ellipse", HOM_COLOR
        elif kind == "pair":
            label, shape, color = f"[{v[1]} | {v[2]}]", "diamond", PAIR_COLOR
        else:
            gad = gd.gadget[v]
            ws = kind != "het"
            label = f"({_seg_label(gad.seg1, ws)}, {_seg_label(gad.seg2, ws)})"
            shape, color = "box", (HET_COLOR if kind == "het" else CORE_COLOR)
        lines.append(f'  n{k} [label={_q(label)}, shape={shape}, fillcolor="{color}"];')
    for a, outs in gd.succ.items():
        if a in hidden:
            continue
        for b in sorted(outs, key=str):  # stable output across runs
            if b in hidden:
                gad = gd.gadget[b]
                hap, seg = (1, gad.seg1) if gad.seg1 is not None else (2, gad.seg2)
                (c,) = gd.succ[b]
                lines.append(f"  {ids[a]} -> {ids[c]} [label={_q(f'h{hap} \u2192 {_seg_label(seg, True)}')}];")
            else:
                lines.append(f"  {ids[a]} -> {ids[b]};")
    lines.append("}")
    return "\n".join(lines)
