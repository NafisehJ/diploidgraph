"""Exactness on inputs built differently from the main generators, and checks
that catch errors the weighted count alone would miss."""

import itertools
import random

import pytest

from diploidgraph import SequenceGraph, build_diploid, decompose, diploid_size
from diploidgraph.diploid import Gadget
from diploidgraph.generators import dense_random_dag, random_dag, random_sp_dag, shuffled
from diploidgraph.gfa import iter_gfa_walks
from diploidgraph.support import observed_support
from diploidgraph.validate import PairDecoder, brute_force_check, exactness_check, weighted_walk_count


def _all_pairs_once(g, gd, dec):
    s, t = g.terminals()
    paths = g.all_paths(s, t)
    dec_ = PairDecoder(gd, dec)
    return all(dec_.count(P, Q) == 1 for P, Q in itertools.product(paths, paths))


@pytest.mark.parametrize("seed", range(80))
def test_dense_random_dags_shuffled(seed):
    """Dense random DAGs (mostly not series-parallel), in a random vertex order."""
    rng = random.Random(seed)
    g = shuffled(dense_random_dag(rng.randint(3, 9), rng.choice([0.3, 0.5, 0.7]), seed=seed), seed=seed)
    gd, dec = build_diploid(g)
    assert brute_force_check(g, gd)["exact"], seed
    assert _all_pairs_once(g, gd, dec), seed
    size = diploid_size(g, dec)
    assert size["vertices"] == gd.num_vertices() and size["edges"] == gd.num_edges()


@pytest.mark.parametrize("seed", range(40))
def test_sp_with_deletions_and_four_alleles(seed):
    g = random_sp_dag(n_ops=9, max_branches=4, p_deletion=0.5, seed=seed)
    gd, dec = build_diploid(g)
    assert dec.is_series_parallel
    assert brute_force_check(g, gd)["exact"], seed


def test_deletion_generator_makes_direct_branches():
    found = 0
    for seed in range(40):
        g = random_sp_dag(n_ops=9, max_branches=4, p_deletion=0.5, seed=seed)
        for p in decompose(g).nodes("P"):
            found += any(c.kind == "Q" for c in p.children)
    assert found > 20


@pytest.mark.parametrize("seed", range(15))
def test_same_graph_any_order_is_exact(seed):
    base = random_dag(n_ops=8, extra_edges=3, seed=seed)
    for k in range(4):
        g = shuffled(base, seed=100 * seed + k)
        gd, dec = build_diploid(g)
        assert brute_force_check(g, gd)["exact"]


def test_pair_decoder_catches_swap_that_weighted_count_misses():
    """Replace (i, j) by a second copy of (j, i): weighted count still equals N^2."""
    g = SequenceGraph.from_edges([("s", "a"), ("a", "t"), ("s", "b"), ("b", "t"), ("s", "c"), ("c", "t")])
    gd, dec = build_diploid(g)
    assert exactness_check(g, gd)["passed"]
    h_ij = ("het", dec.nodes("P")[0].uid, 0, 1)
    gad = gd.gadget[h_ij]
    gd.gadget[h_ij] = Gadget(gad.seg2, gad.seg1)  # now a duplicate of (1, 0); (0, 1) is gone
    assert exactness_check(g, gd)["passed"]  # the necessary condition cannot see it
    assert not brute_force_check(g, gd)["exact"]
    assert not _all_pairs_once(g, gd, dec)


def test_modular_check_on_large_graph_agrees_with_exact():
    g = random_dag(n_ops=600, extra_edges=20, seed=2)
    gd, _ = build_diploid(g)
    assert exactness_check(g, gd)["method"] == "exact"
    assert exactness_check(g, gd, exact_limit_vertices=0)["passed"]


def test_single_vertex_graph():
    g = SequenceGraph()
    g.add_vertex("v")
    dec = decompose(g)
    assert dec.is_series_parallel and dec.root.n_paths == 1
    gd, _ = build_diploid(g, dec)
    assert gd.num_vertices() == 1 and weighted_walk_count(gd) == 1


def test_gfa_blank_lines_and_walks_off_graph(tmp_path):
    gfa = tmp_path / "x.gfa"
    gfa.write_text("S\t1\tA\n\nS\t2\tC\nL\t1\t+\t2\t+\t0M\n\nP\tok\t1+,2+\t*\nP\tbad\t2+,1+\t*\n")
    walks = list(iter_gfa_walks(str(gfa)))
    assert walks == [["1", "2"], ["2", "1"]]
    g = SequenceGraph.from_edges([("1", "2")])
    _, counts = observed_support(g, walks)
    assert counts["not_in_graph"] == 1 and counts["distinct_full"] == 1


@pytest.mark.parametrize("seed", range(30))
def test_product_baseline_is_exact_and_not_smaller(seed):
    """The product-graph baseline represents the same pairs (exactly once each)."""
    from diploidgraph.decompose import trivial_decomposition
    rng = random.Random(seed)
    g = random_dag(n_ops=7, extra_edges=seed % 3, seed=seed) if seed % 2 else \
        dense_random_dag(rng.randint(3, 8), 0.5, seed=seed)
    tdec = trivial_decomposition(g)
    gd_prod, _ = build_diploid(g, tdec)
    assert brute_force_check(g, gd_prod)["exact"]
    size = diploid_size(g, tdec)
    assert size["vertices"] == gd_prod.num_vertices() and not size["capped"]
    assert diploid_size(g)["vertices"] <= size["vertices"]


def test_size_cap():
    from diploidgraph.decompose import trivial_decomposition
    g = random_sp_dag(200, seed=1)
    size = diploid_size(g, trivial_decomposition(g), max_pair_states=100)
    assert size["capped"]
