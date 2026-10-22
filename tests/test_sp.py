import pytest

from diploidgraph import NotTwoTerminalDAG, SequenceGraph, build_diploid, decompose
from diploidgraph.generators import random_dag, random_sp_dag, series_of, w_graph
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


def test_w_graph_exact():
    g = w_graph()
    dec = decompose(g)
    assert not dec.is_series_parallel
    gd, _ = build_diploid(g, dec)
    report = brute_force_check(g, gd)
    assert report["exact"], report
    assert report["expected_pairs"] == 9


def test_bubble_then_w_exact():
    g = series_of(bubble(), w_graph())
    gd, dec = build_diploid(g)
    assert not dec.is_series_parallel
    assert brute_force_check(g, gd)["exact"]


@pytest.mark.parametrize("seed", range(60))
def test_random_non_sp_exact(seed):
    g = random_dag(n_ops=8, extra_edges=1 + seed % 4, seed=seed)
    gd, dec = build_diploid(g)
    report = brute_force_check(g, gd)
    assert report["exact"], (seed, report)


@pytest.mark.parametrize("seed", range(5))
def test_large_non_sp_weighted_count(seed):
    g = random_dag(n_ops=1500, extra_edges=40, seed=seed)
    gd, dec = build_diploid(g)
    assert not dec.is_series_parallel
    s, t = g.terminals()
    n = g.count_paths(s, t)
    assert dec.root.n_paths == n
    assert weighted_walk_count(gd) == n * n


def test_core_size_bound():
    g = random_dag(n_ops=300, extra_edges=15, seed=11)
    gd, dec = build_diploid(g)
    (r,) = dec.nodes("R")
    n_core = len({c.source for c in r.children} | {c.sink for c in r.children})
    assert gd.count_by_kind().get("pair", 0) <= n_core * (n_core - 1)
