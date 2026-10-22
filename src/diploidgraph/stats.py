"""Measure one run: sizes, runtime per stage, peak memory, naive baselines."""

from __future__ import annotations

import math
import time
import resource
import sys
from typing import Dict, Optional

from .decompose import decompose
from .diploid import build_diploid, diploid_size, naive_pair_count
from .gfa import read_gfa
from .seqgraph import SequenceGraph
from .validate import exactness_check


def log10_int(x: int) -> float:
    """log10 of a possibly huge non-negative integer."""
    if x <= 0:
        return float("-inf")
    k = max(0, x.bit_length() - 64)
    return math.log10(x >> k) + k * math.log10(2)


def peak_rss_mb() -> float:
    """Peak resident memory of this process so far, in MB (Linux/macOS)."""
    kb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return kb / 1024 / (1024 if sys.platform == "darwin" else 1)


def run_stats(g: SequenceGraph, n_observed: Optional[int] = None, count_only: bool = False) -> Dict[str, object]:
    """Sizes, timings and memory for one input graph.

    count_only=True computes the sizes of G_D without building it (for graphs
    whose non-series-parallel core is too large to build in Python); pair
    coverage is then not checked.
    """
    out: Dict[str, object] = {"n_vertices": g.num_vertices(), "n_edges": g.num_edges(), "count_only": count_only}
    g.terminals()

    t0 = time.perf_counter()
    dec = decompose(g)
    t1 = time.perf_counter()
    if count_only:
        size = diploid_size(g, dec)
        gd = None
    else:
        gd, _ = build_diploid(g, dec)
        kinds = gd.count_by_kind()
        size = {k: kinds.get(k, 0) for k in ("hom", "het", "fork", "pair", "step")}
        size.update(vertices=gd.num_vertices(), edges=gd.num_edges())
    t2 = time.perf_counter()

    pnodes = dec.nodes("P")
    n_paths = dec.root.n_paths
    out.update(
        is_series_parallel=dec.is_series_parallel,
        n_pnodes=len(pnodes),
        max_pnode_branches=max((len(p.children) for p in pnodes), default=0),
        core_vertices=len({c.source for r in dec.nodes("R") for c in r.children}
                          | {c.sink for r in dec.nodes("R") for c in r.children}),
        core_edges=sum(len(r.children) for r in dec.nodes("R")),
        time_decompose_s=t1 - t0,
        time_build_s=t2 - t1,
        peak_rss_mb=peak_rss_mb(),
        diploid_vertices=size["vertices"],
        diploid_edges=size["edges"],
        diploid_new_vertices=size["vertices"] - size["hom"],
        **{f"diploid_kind_{k}": size[k] for k in ("hom", "het", "fork", "pair", "step")},
        log10_st_paths=log10_int(n_paths),
        log10_naive_ordered_all_paths=log10_int(naive_pair_count(n_paths, ordered=True)),
    )
    if gd is not None:
        check = exactness_check(g, gd, n_paths=n_paths)
        out["exactness_check"] = check["method"]
        out["exactness_passed"] = check["passed"]
    else:
        out["exactness_check"] = None
        out["exactness_passed"] = None
    if n_observed is not None:
        out["naive_ordered_observed"] = naive_pair_count(n_observed, ordered=True)
    return out


def naive_observed_nodes(lengths) -> int:
    """Vertices of a naive diploid graph over observed walks: one chain per ordered
    pair (p, q) of distinct walks, with max(|p|, |q|) paired vertices."""
    ls = sorted(lengths)
    return 2 * sum(k * L for k, L in enumerate(ls))  # each L is the max for the k shorter walks


def gfa_stats(path: str, count_only: bool = False) -> Dict[str, object]:
    t0 = time.perf_counter()
    data = read_gfa(path)
    t_read = time.perf_counter() - t0
    distinct = {}
    for w in data.walks:
        distinct.setdefault(w.digest, w.length)
    lengths = list(distinct.values())
    out: Dict[str, object] = {
        "input": path,
        "time_read_s": t_read,
        "dropped_links": data.dropped_links,
        "n_haplotype_walks": len(data.walks),
        "n_distinct_haplotype_walks": len(distinct),
        "n_full_haplotype_walks": data.n_full_walks(),
        "unoriented_haplotypes": data.unoriented_haplotypes,
        "walk_length_min_mean_max": [min(lengths), sum(lengths) / len(lengths), max(lengths)] if lengths else None,
        "naive_ordered_observed_nodes": naive_observed_nodes(lengths),
    }
    out.update(run_stats(data.graph, n_observed=len(distinct), count_only=count_only))
    return out
