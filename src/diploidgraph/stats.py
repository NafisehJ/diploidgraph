"""Measure one run: sizes, runtime per stage, peak memory, naive baselines."""

from __future__ import annotations

import math
import time
import tracemalloc
from typing import Dict, Optional

from .decompose import decompose
from .diploid import build_diploid, naive_pair_count
from .gfa import read_gfa
from .seqgraph import SequenceGraph
from .validate import coverage


def log10_int(x: int) -> float:
    """log10 of a possibly huge non-negative integer."""
    if x <= 0:
        return float("-inf")
    k = max(0, x.bit_length() - 64)
    return math.log10(x >> k) + k * math.log10(2)


def run_stats(g: SequenceGraph, allow_partial: bool = False, n_observed: Optional[int] = None) -> Dict[str, object]:
    out: Dict[str, object] = {"n_vertices": g.num_vertices(), "n_edges": g.num_edges()}
    s, t = g.terminals()

    tracemalloc.start()
    t0 = time.perf_counter()
    dec = decompose(g)
    t1 = time.perf_counter()
    gd, _ = build_diploid(g, dec, allow_partial=allow_partial)
    t2 = time.perf_counter()
    peak = tracemalloc.get_traced_memory()[1]
    tracemalloc.stop()

    pnodes = dec.nodes("P")
    n_paths = dec.root.n_paths
    represented, total = coverage(g, gd)
    out.update(
        is_series_parallel=dec.is_series_parallel,
        n_pnodes=len(pnodes),
        max_pnode_branches=max((len(p.children) for p in pnodes), default=0),
        n_rnodes=len(dec.nodes("R")),
        time_decompose_s=t1 - t0,
        time_build_s=t2 - t1,
        peak_python_memory_bytes=peak,
        diploid_vertices=gd.num_vertices(),
        diploid_edges=gd.num_edges(),
        diploid_het_vertices=gd.num_het_vertices(),
        log10_st_paths=log10_int(n_paths),
        log10_naive_ordered_all_paths=log10_int(naive_pair_count(n_paths, ordered=True)),
        pair_coverage=represented / total if total else None,
    )
    if n_observed is not None:
        out["naive_ordered_observed"] = naive_pair_count(n_observed, ordered=True)
    return out


def gfa_stats(path: str, allow_partial: bool = False) -> Dict[str, object]:
    t0 = time.perf_counter()
    data = read_gfa(path)
    t_read = time.perf_counter() - t0
    n_obs = len({w for _, w in data.haplotypes})
    out: Dict[str, object] = {
        "input": path,
        "time_read_s": t_read,
        "dropped_links": data.dropped_links,
        "n_haplotype_walks": len(data.haplotypes),
        "n_distinct_haplotype_walks": n_obs,
        "unoriented_haplotypes": data.unoriented_haplotypes,
    }
    out.update(run_stats(data.graph, allow_partial=allow_partial, n_observed=n_obs))
    return out
