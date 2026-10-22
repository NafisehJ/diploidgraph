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
    assert r["exactness_passed"] is True
    assert r["diploid_vertices"] == r["n_vertices"] + r["diploid_new_vertices"]


def test_naive_observed_nodes():
    from diploidgraph.stats import naive_observed_nodes
    ls = [3, 5, 4]
    brute = sum(max(a, b) for i, a in enumerate(ls) for j, b in enumerate(ls) if i != j)
    assert naive_observed_nodes(ls) == brute


def test_sizes_do_not_depend_on_hash_seed(tmp_path):
    import json
    import os
    import subprocess
    import sys

    code = (
        "import json\n"
        "from diploidgraph import SequenceGraph, diploid_size\n"
        "from diploidgraph.generators import random_dag\n"
        "g0 = random_dag(n_ops=400, extra_edges=12, seed=4)\n"
        "g = SequenceGraph.from_edges([(f'v{u}', f'v{v}') for u, v in g0.edges()])\n"
        "print(json.dumps(diploid_size(g)))\n"
    )
    outs = set()
    for seed in ("1", "2", "3"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
        outs.add(r.stdout.strip())
    assert len(outs) == 1, outs


def test_modular_check_agrees_and_detects_errors():
    from diploidgraph import build_diploid
    from diploidgraph.generators import random_dag
    from diploidgraph.validate import PRIMES, exactness_check, weighted_walk_count
    g = random_dag(n_ops=600, extra_edges=20, seed=2)
    gd, _ = build_diploid(g)
    exact = weighted_walk_count(gd)
    for p in PRIMES:
        assert weighted_walk_count(gd, p) == exact % p
    assert exactness_check(g, gd, exact_limit_bits=0)["passed"]
    # remove one heterozygous vertex: both checks must fail
    h = next(v for v in gd.succ if v[0] == "het")
    for outs in gd.succ.values():
        outs.discard(h)
    del gd.succ[h]
    assert not exactness_check(g, gd)["passed"]
    assert not exactness_check(g, gd, exact_limit_bits=0)["passed"]


def test_example_figures_do_not_depend_on_hash_seed():
    import os
    import subprocess
    import sys

    code = (
        "from diploidgraph import SequenceGraph, build_diploid\n"
        "from diploidgraph.export import diploid_graph_dot\n"
        "g = SequenceGraph.from_edges([('s','a'),('s','b'),('a','m'),('b','m'),"
        "('m','c'),('m','d'),('c','d'),('c','t'),('d','t')])\n"
        "print(diploid_graph_dot(build_diploid(g)[0]))\n"
    )
    outs = set()
    for seed in ("1", "2", "3", "4"):
        env = dict(os.environ, PYTHONHASHSEED=seed)
        r = subprocess.run([sys.executable, "-c", code], env=env, capture_output=True, text=True, check=True)
        outs.add(r.stdout)
    assert len(outs) == 1
