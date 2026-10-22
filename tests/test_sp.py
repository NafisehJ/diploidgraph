import pytest

from diploidgraph import NonSeriesParallelError, NotTwoTerminalDAG, SequenceGraph, build_diploid, decompose
from diploidgraph.generators import random_sp_dag, series_of, w_graph
from diploidgraph.validate import brute_force_check, coverage, weighted_walk_count


def bubble():
    # s -> a -> t, s -> b -> t  (one SNP-like bubble)
    return SequenceGraph.from_edges([("s", "a"), ("a", "t"), ("s", "b"), ("b", "t")])


def test_input_validation():
    with pytest.raises(NotTwoTerminalDAG):
        SequenceGraph.from_edges([(1, 2), (2, 3), (3, 1)]).terminals()
    with pytest.raises(NotTwoTerminalDAG):
        SequenceGraph.from_edges([(1, 2), (1, 3)]).terminals()  # two sinks


def test_single_bubble():
    g = bubble()
    dec = decompose(g)
    assert dec.is_series_parallel
    (p,) = dec.nodes("P")
    assert (p.source, p.sink) == ("s", "t") and len(p.children) == 2
    gd, _ = build_diploid(g, dec)
    assert gd.num_het_vertices() == 2  # (a,b) and (b,a)
    assert gd.num_vertices() == 4 + 2
    assert brute_force_check(g, gd)["exact"]


def test_nested_bubbles_size_formula():
    # bubble in series with a 3-branch bubble whose first branch has a nested bubble
    g = SequenceGraph.from_edges([
        (0, 1), (0, 2), (1, 3), (2, 3),                 # P(0,3), k=2
        (3, 4), (4, 6), (3, 5), (5, 6), (3, 7), (7, 6),  # P(3,6), k=3
        (4, 8), (8, 6),                                  # extra branch inside: P(4,6), k=2
    ])
    dec = decompose(g)
    ks = sorted(len(p.children) for p in dec.nodes("P"))
    assert ks == [2, 2, 3]
    gd, _ = build_diploid(g, dec)
    assert gd.num_het_vertices() == sum(k * (k - 1) for k in ks)
    assert gd.num_vertices() == g.num_vertices() + sum(k * (k - 1) for k in ks)
    assert gd.num_edges() == g.num_edges() + 2 * sum(k * (k - 1) for k in ks)
    assert brute_force_check(g, gd)["exact"]


@pytest.mark.parametrize("seed", range(40))
def test_random_sp_exact(seed):
    g = random_sp_dag(n_ops=12, max_branches=3, seed=seed)
    gd, dec = build_diploid(g)
    assert dec.is_series_parallel
    report = brute_force_check(g, gd)
    assert report["exact"], report
    represented, total = coverage(g, gd)
    assert represented == total


@pytest.mark.parametrize("seed", range(5))
def test_large_sp_weighted_count(seed):
    g = random_sp_dag(n_ops=2000, seed=seed)
    gd, dec = build_diploid(g)
    s, t = g.terminals()
    assert dec.root.n_paths == g.count_paths(s, t)
    assert weighted_walk_count(gd) == g.count_paths(s, t) ** 2


def test_non_sp_detected():
    g = series_of(bubble(), w_graph())
    dec = decompose(g)
    assert not dec.is_series_parallel
    with pytest.raises(NonSeriesParallelError):
        build_diploid(g, dec)
    gd, _ = build_diploid(g, dec, allow_partial=True)
    represented, total = coverage(g, gd)
    assert represented < total  # W-motif pairs are not represented yet
