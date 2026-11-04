"""Cut a small real window out of a pangenome GFA (for example BRCA1 from HPRC v2.0).

A window is the subgraph between two cut vertices a and b of the input: every
source-to-sink path passes through a and then b, so the window is again a
two-terminal DAG and every haplotype walk crosses it from a to b. The script
lists candidate windows, prefers ones that are not series-parallel (a W-motif
or larger core), contain a large deletion and many bubbles, and writes the best
one as a GFA file with the original segment IDs, sequences and haplotype walks.

    python examples/extract_window.py ~/brca1/BRCA1.gfa --out examples/data/brca1_window.gfa
    python examples/extract_window.py ~/brca1/BRCA1.gfa --list 10        # only list candidates
    python examples/extract_window.py ~/brca1/MHC.gfa --list 10          # MHC: a few minutes
    python examples/extract_window.py ~/brca1/BRCA1.gfa --pick 3 --out examples/data/brca1_window.gfa

The script only reads the input; it does not change any result of the paper.
"""

from __future__ import annotations

import argparse
import re
from collections import deque
from typing import Dict, List, Optional, Set, Tuple

from diploidgraph import SequenceGraph, decompose

_W_STEP = re.compile(r"([<>])([^<>]+)")


def read_graph(path: str):
    """Segments and links only; haplotype walks are streamed later (they can be long)."""
    seq: Dict[str, str] = {}
    edges: Set[Tuple[str, str]] = set()
    n_walks = dropped = 0
    with open(path) as fh:
        for line in fh:
            tag = line[:2]
            if tag == "S\t":
                f = line.rstrip("\n").split("\t", 3)
                seq[f[1]] = f[2]
            elif tag == "L\t":
                f = line.split("\t", 6)
                if f[2] == "+" and f[4] == "+":
                    edges.add((f[1], f[3]))
                elif f[2] == "-" and f[4] == "-":
                    edges.add((f[3], f[1]))
                else:
                    dropped += 1
            elif tag in ("W\t", "P\t"):
                n_walks += 1
    return seq, edges, n_walks, dropped


def iter_walks(path: str):
    with open(path) as fh:
        for line in fh:
            if line.startswith("W\t"):
                f = line.rstrip("\n").split("\t")
                yield "W", f, [(m.group(2), m.group(1) == ">") for m in _W_STEP.finditer(f[6])]
            elif line.startswith("P\t"):
                f = line.rstrip("\n").split("\t")
                yield "P", f, [(x[:-1], x[-1] == "+") for x in f[2].split(",") if x]


def topo(vertices: List[str], succ: Dict[str, List[str]]) -> List[str]:
    indeg = {v: 0 for v in vertices}
    for u in vertices:
        for w in succ[u]:
            indeg[w] += 1
    queue = deque(v for v in vertices if indeg[v] == 0)
    order = []
    while queue:
        u = queue.popleft()
        order.append(u)
        for w in succ[u]:
            indeg[w] -= 1
            if indeg[w] == 0:
                queue.append(w)
    if len(order) != len(vertices):
        raise SystemExit("input graph is not a DAG")
    return order


def window_vertices(a: str, b: str, succ, pred) -> Set[str]:
    def reach(start, nbr, stop):
        seen, stack = {start}, [start]
        while stack:
            u = stack.pop()
            if u == stop:
                continue
            for w in nbr[u]:
                if w not in seen:
                    seen.add(w)
                    stack.append(w)
        return seen
    return reach(a, succ, b) & reach(b, pred, a)


def largest_deletion(vs: Set[str], pos: Dict[str, int], succ, seq) -> int:
    """Longest sequence (bp) skipped by one edge: max over edges (u, v) of the longest u-v path inside."""
    sub = sorted(vs, key=pos.get)
    best = 0
    for u in sub:
        # longest bp strictly between u and each vertex reachable from u
        dist: Dict[str, int] = {u: 0}
        for x in sub[sub.index(u):]:  # small window: quadratic is fine
            if x not in dist:
                continue
            for w in succ[x]:
                if w in vs:
                    d = dist[x] + (len(seq[x]) if x != u else 0)
                    if d > dist.get(w, -1):
                        dist[w] = d
        for w in succ[u]:
            if w in vs:
                best = max(best, dist.get(w, 0))
    return best


def slice_walk(kind, f, steps, a, b, seq) -> Optional[str]:
    fwd = all(x for _, x in steps)
    if not fwd and any(x for _, x in steps):
        return None  # mixed orientation
    ids = [s for s, _ in steps] if fwd else [s for s, _ in reversed(steps)]
    if a not in ids or b not in ids:
        return None
    i, j = ids.index(a), ids.index(b)
    if j < i:
        return None
    part = ids[i:j + 1]
    if kind == "P":
        segs = ",".join(s + "+" for s in part) if fwd else ",".join(s + "-" for s in reversed(part))
        return f"P\t{f[1]}\t{segs}\t*"
    before = sum(len(seq[s]) for s in ids[:i])
    after = sum(len(seq[s]) for s in ids[j + 1:])
    length = sum(len(seq[s]) for s in part)
    start = int(f[4]) + (before if fwd else after)
    text = "".join(">" + s for s in part) if fwd else "".join("<" + s for s in reversed(part))
    return "\t".join([f[0], f[1], f[2], f[3], str(start), str(start + length), text])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    ap.add_argument("--out", help="write the chosen window here")
    ap.add_argument("--min-vertices", type=int, default=10)
    ap.add_argument("--max-vertices", type=int, default=40)
    ap.add_argument("--max-paths", type=int, default=200, help="keep brute force checks fast")
    ap.add_argument("--min-deletion", type=int, default=20, help="bp skipped by one edge to count as a large deletion")
    ap.add_argument("--list", type=int, default=10, help="print this many best candidates")
    ap.add_argument("--pick", type=int, default=1, help="which candidate to write (1 = best)")
    args = ap.parse_args()

    seq, edges, n_walks, dropped = read_graph(args.gfa)
    g = SequenceGraph.from_edges(sorted(edges))
    vertices = g.vertices
    succ: Dict[str, List[str]] = {v: [] for v in vertices}
    pred: Dict[str, List[str]] = {v: [] for v in vertices}
    for u, v in edges:
        succ[u].append(v)
        pred[v].append(u)
    order = topo(vertices, succ)
    pos = {v: i for i, v in enumerate(order)}
    s, t = g.terminals()

    # cut vertices: every s-t path passes through them, i.e. paths(s,v) * paths(v,t) = paths(s,t).
    # Path counts are astronomically large on locus graphs (about 10^71057 for MHC), so they are
    # compared modulo two large primes, as in validate.weighted_walk_count.
    cuts = None
    for prime in (2**61 - 1, 2**89 - 1):
        fwd: Dict[str, int] = {}
        for v in order:
            fwd[v] = 1 if v == s else sum(fwd[u] for u in pred[v]) % prime
        bwd: Dict[str, int] = {}
        for v in reversed(order):
            bwd[v] = 1 if v == t else sum(bwd[w] for w in succ[v]) % prime
        total = fwd[t]
        ok = {v for v in order if fwd[v] * bwd[v] % prime == total}
        cuts = ok if cuts is None else cuts & ok
    cuts = [v for v in order if v in cuts]
    print(f"{args.gfa}: {len(vertices)} vertices, {len(edges)} edges, {n_walks} walks, "
          f"{dropped} dropped links, {len(cuts)} cut vertices")

    cands = []
    for i, a in enumerate(cuts):
        for b in cuts[i + 2:]:
            vs = window_vertices(a, b, succ, pred)
            if len(vs) > args.max_vertices:
                break
            if len(vs) < args.min_vertices:
                continue
            sub = SequenceGraph.from_edges([(u, w) for u in vs for w in succ[u] if w in vs])
            n_paths = sub.count_paths(a, b)
            if n_paths > args.max_paths:
                continue
            dec = decompose(sub)
            core = sum(len(r.children) for r in dec.nodes("R"))
            deletion = largest_deletion(vs, pos, succ, seq)
            score = (core > 0, deletion >= args.min_deletion, len(dec.nodes("P")), -len(vs))
            cands.append((score, a, b, len(vs), sub.num_edges(), n_paths, len(dec.nodes("P")), core, deletion))
    if not cands:
        raise SystemExit("no window found; try other --min-vertices / --max-vertices / --max-paths")
    cands.sort(key=lambda c: c[0], reverse=True)

    print(f"\n{'#':>3} {'from':>8} {'to':>8} {'|V|':>4} {'|E|':>4} {'paths':>6} {'P-nodes':>7} {'core edges':>10} {'largest deletion (bp)':>22}")
    for k, c in enumerate(cands[:args.list], 1):
        _, a, b, nv, ne, npaths, npn, core, dele = c
        print(f"{k:>3} {a:>8} {b:>8} {nv:>4} {ne:>4} {npaths:>6} {npn:>7} {core:>10} {dele:>22}")

    if not args.out:
        return
    _, a, b, *_ = cands[args.pick - 1]
    vs = window_vertices(a, b, succ, pred)
    lines = ["H\tVN:Z:1.1"]
    lines += [f"S\t{v}\t{seq[v]}" for v in sorted(vs, key=pos.get)]
    lines += [f"L\t{u}\t+\t{w}\t+\t0M" for u in sorted(vs, key=pos.get) for w in sorted(succ[u], key=pos.get) if w in vs]
    kept = 0
    for kind, f, steps in iter_walks(args.gfa):
        line = slice_walk(kind, f, steps, a, b, seq)
        if line:
            lines.append(line)
            kept += 1
    with open(args.out, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    print(f"\nwrote {args.out}: window {a} -> {b}, {len(vs)} segments, {kept} of {n_walks} walks")


if __name__ == "__main__":
    main()
