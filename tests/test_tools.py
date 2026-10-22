from diploidgraph import SequenceGraph, build_diploid
from diploidgraph.export import branch_label, diploid_graph_dot, sequence_graph_dot
from diploidgraph.generators import random_sp_dag
from diploidgraph.gfa import read_gfa
from diploidgraph.stats import log10_int, run_stats


def test_branch_labels_nested():
    g = SequenceGraph.from_edges([("s", "a"), ("a", "b"), ("a", "c"), ("b", "d"), ("c", "d"), ("d", "t"), ("s", "t")])
    gd, dec = build_diploid(g)
    outer = next(p for p in dec.nodes("P") if (p.source, p.sink) == ("s", "t"))
    assert sorted(branch_label(b) for b in outer.children) == ["-", "a(b|c)d"]
    assert "digraph" in sequence_graph_dot(g) and "digraph" in diploid_graph_dot(gd)


def test_gfa_reader(tmp_path):
    gfa = tmp_path / "x.gfa"
    gfa.write_text(
        "S\t1\tA\nS\t2\tC\nS\t3\tG\nS\t4\tT\n"
        "L\t1\t+\t2\t+\t0M\nL\t1\t+\t3\t+\t0M\nL\t2\t+\t4\t+\t0M\nL\t4\t-\t3\t-\t0M\n"
        "P\tref\t1+,2+,4+\t*\n"
        "W\tHG1\t1\tchr\t0\t3\t<4<3<1\n"
    )
    d = read_gfa(str(gfa), keep_walks=True)
    assert d.graph.num_vertices() == 4 and d.graph.num_edges() == 4
    assert (3, 4) not in d.graph.edges() and ("3", "4") in set(d.graph.edges())
    walks = dict(d.haplotypes)
    assert walks["ref"] == ("1", "2", "4")
    assert walks["HG1#1#chr"] == ("1", "3", "4")  # reverse walk normalised to forward
    assert d.n_distinct_walks() == 2 and d.n_full_walks() == 2
    assert [w.length for w in d.walks] == [3, 3]


def test_stats_and_log10():
    assert abs(log10_int(10**300) - 300) < 1e-9
    r = run_stats(random_sp_dag(50, seed=3))
    assert r["pair_coverage"] == 1.0
    assert r["diploid_vertices"] == r["n_vertices"] + r["diploid_new_vertices"]


def test_naive_observed_nodes():
    from diploidgraph.stats import naive_observed_nodes
    ls = [3, 5, 4]
    brute = sum(max(a, b) for i, a in enumerate(ls) for j, b in enumerate(ls) if i != j)
    assert naive_observed_nodes(ls) == brute
