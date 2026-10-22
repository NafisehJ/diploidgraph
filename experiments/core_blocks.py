"""Split the non-series-parallel core into blocks and report their sizes.

A core vertex c is a *cut* if every source-sink path passes through it. Two
haplotypes that diverge before a cut must meet at that cut at the latest, so the
pair states of the diploid construction never span two blocks. The number of
pair states is therefore at most sum_b n_b (n_b - 1) over blocks b with n_b
vertices; the largest block decides whether a graph is practical.

    python experiments/core_blocks.py MHC.gfa
"""

from __future__ import annotations

import argparse
from typing import Dict, List

from diploidgraph import decompose
from diploidgraph.gfa import read_gfa


def core_blocks(dec) -> List[List]:
    edges = [(c.source, c.sink) for r in dec.nodes("R") for c in r.children]
    if not edges:
        return []
    succ: Dict = {}
    indeg: Dict = {}
    for u, v in edges:
        succ.setdefault(u, []).append(v)
        succ.setdefault(v, [])
        indeg[v] = indeg.get(v, 0) + 1
        indeg.setdefault(u, 0)
    order, stack = [], [u for u, d in indeg.items() if d == 0]
    while stack:
        u = stack.pop()
        order.append(u)
        for v in succ[u]:
            indeg[v] -= 1
            if indeg[v] == 0:
                stack.append(v)
    pos = {v: i for i, v in enumerate(order)}
    # c is a cut iff no core edge (a, b) has pos[a] < pos[c] < pos[b].
    cover = [0] * (len(order) + 1)
    for u, v in edges:
        if pos[v] - pos[u] > 1:
            cover[pos[u] + 1] += 1
            cover[pos[v]] -= 1
    blocks, cur, run = [], [order[0]], 0
    for i in range(1, len(order)):
        run += cover[i]
        cur.append(order[i])
        if run == 0:  # order[i] is a cut: close the block (cut belongs to both sides)
            blocks.append(cur)
            cur = [order[i]]
    return [b for b in blocks if len(b) > 2]  # blocks of 2 vertices are single core edges


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("gfa")
    args = ap.parse_args()
    dec = decompose(read_gfa(args.gfa).graph)
    blocks = core_blocks(dec)
    sizes = sorted((len(b) for b in blocks), reverse=True)
    print(f"non-SP blocks: {len(sizes)}")
    print(f"largest block sizes: {sizes[:10]}")
    print(f"bound on pair states, sum n_b(n_b-1): {sum(n * (n - 1) for n in sizes)}")


if __name__ == "__main__":
    main()
