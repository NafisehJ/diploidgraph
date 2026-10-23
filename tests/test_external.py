import itertools
import math
import random

import pytest

from diploidgraph.external import chain_vertices, external_sizes, top_level_bubbles
from diploidgraph.generators import random_dag, random_sp_dag


def brute_chain(g):
    s, t = g.terminals()
    paths = g.all_paths(s, t)
    order = g.topological_order()
    on_all = set(paths[0]).intersection(*map(set, paths[1:]))
    return [v for v in order if v in on_all]


@pytest.mark.parametrize("seed", range(40))
def test_chain_vertices_match_brute_force(seed):
    rng = random.Random(seed)
    g = random_sp_dag(rng.randint(2, 9), seed=seed) if seed % 2 else random_dag(rng.randint(2, 7), rng.randint(0, 4), seed=seed)
    assert chain_vertices(g) == brute_chain(g)


def write_gfa(tmp_path, segs, links, walks, name="g.gfa"):
    lines = ["H\tVN:Z:1.1"]
    lines += [f"S\t{s}\t{seq}" for s, seq in segs]
    lines += [f"L\t{a}\t+\t{b}\t+\t0M" for a, b in links]
    for i, w in enumerate(walks):
        lines.append(f"P\th{i}\t{','.join(x + '+' for x in w)}\t*")
    p = tmp_path / name
    p.write_text("\n".join(lines) + "\n")
    return str(p)


def two_snp_graph(tmp_path, spacer_len):
    # s -> {a1|a2} -> m -> {b1|b2} -> t ; m carries the spacer between the two SNPs
    segs = [("s", "A" * 40), ("a1", "C"), ("a2", "G"), ("m", "T" * spacer_len), ("b1", "C"), ("b2", "G"), ("t", "A" * 40)]
    links = [("s", "a1"), ("s", "a2"), ("a1", "m"), ("a2", "m"), ("m", "b1"), ("m", "b2"), ("b1", "t"), ("b2", "t")]
    walks = [["s", "a1", "m", "b1", "t"], ["s", "a2", "m", "b2", "t"], ["s", "a1", "m", "b1", "t"]]
    return write_gfa(tmp_path, segs, links, walks)


def test_far_apart_bubbles_stay_separate(tmp_path):
    r = external_sizes(two_snp_graph(tmp_path, 100), k=31)
    assert r["n_complete_walks"] == 3 and r["n_distinct_complete_walks"] == 2
    assert r["n_top_level_bubbles"] == 2 and r["pangenie_bubbles"] == 2
    assert r["pangenie_states"] == 2 * 3 * 3
    assert r["pangenie_states_if_duplicates_collapsed"] == 2 * 2 * 2
    # mosaics switch between the bubbles: 2 * 2 haplotypes, i.e. all 4 s-t paths
    assert r["log10_pangenie_covered_ordered_pairs"] == pytest.approx(2 * math.log10(4))
    assert r["log10_ours_covered_ordered_pairs"] == pytest.approx(2 * math.log10(4))


def test_close_bubbles_are_combined(tmp_path):
    # gap = shared sequence between the SNPs = len(m) = 5 < 31
    r = external_sizes(two_snp_graph(tmp_path, 5), k=31)
    assert r["pangenie_bubbles"] == 1
    # combined bubble: alleles are the observed sub-walks only (a1..b1, a2..b2)
    assert r["pangenie_alleles_per_bubble_max"] == 2
    assert r["log10_pangenie_covered_ordered_pairs"] == pytest.approx(2 * math.log10(2))
    assert r["pangenie_states"] == 1 * 3 * 3


def test_gap_threshold_is_strict(tmp_path):
    # gap = len(m) exactly; combined iff gap < k
    assert external_sizes(two_snp_graph(tmp_path, 31), k=31)["pangenie_bubbles"] == 2
    assert external_sizes(two_snp_graph(tmp_path, 30), k=31)["pangenie_bubbles"] == 1


def test_bubble_without_variation_is_not_a_pangenie_bubble(tmp_path):
    segs = [("s", "A"), ("x", "C"), ("y", "G"), ("t", "A")]
    links = [("s", "x"), ("s", "y"), ("x", "t"), ("y", "t")]
    walks = [["s", "x", "t"], ["s", "x", "t"]]
    r = external_sizes(write_gfa(tmp_path, segs, links, walks))
    assert r["n_top_level_bubbles"] == 1 and r["n_variant_bubbles"] == 0
    assert r["pangenie_bubbles"] == 0 and r["pangenie_states"] == 0


@pytest.mark.parametrize("seed", range(15))
def test_coverage_bounds_on_random_graphs(tmp_path, seed):
    """PanGenie's mosaics are s-t paths, so they never exceed ours; observed walks are among them."""
    rng = random.Random(seed)
    g = random_dag(rng.randint(3, 7), rng.randint(0, 4), seed=seed)
    s, t = g.terminals()
    paths = g.all_paths(s, t)
    walks = [list(map(str, rng.choice(paths))) for _ in range(rng.randint(1, 6))]
    segs = [(str(v), "ACGT"[i % 4] * rng.randint(1, 40)) for i, v in enumerate(g.vertices)]
    links = [(str(a), str(b)) for a, b in g.edges()]
    r = external_sizes(write_gfa(tmp_path, segs, links, walks, name=f"r{seed}.gfa"))
    n_distinct = len({tuple(w) for w in walks})
    assert r["n_distinct_complete_walks"] == n_distinct
    assert r["log10_pangenie_covered_ordered_pairs"] <= r["log10_ours_covered_ordered_pairs"] + 1e-9
    assert r["log10_pangenie_covered_ordered_pairs"] >= 2 * math.log10(n_distinct) - 1e-9


def test_reference_samples_left_out_of_pangenie_panel(tmp_path):
    # one SNP; two assembly haplotypes (a1, a2) and a reference walk (a1)
    segs = [("s", "A" * 40), ("a1", "C"), ("a2", "G"), ("t", "A" * 40)]
    links = [("s", "a1"), ("s", "a2"), ("a1", "t"), ("a2", "t")]
    lines = ["H\tVN:Z:1.1"] + [f"S\t{a}\t{b}" for a, b in segs] + [f"L\t{a}\t+\t{b}\t+\t0M" for a, b in links]
    lines += ["W\tHG1\t1\tchr\t0\t81\t>s>a1>t", "W\tHG1\t2\tchr\t0\t81\t>s>a2>t", "W\tGRCh38\t0\tchr\t0\t81\t>s>a1>t"]
    p = tmp_path / "ref.gfa"
    p.write_text("\n".join(lines) + "\n")
    r = external_sizes(str(p))
    assert r["n_complete_walks"] == 3
    assert r["pangenie_panel_haplotypes"] == 2 and r["pangenie_states"] == 4
    r_all = external_sizes(str(p), reference_samples=())
    assert r_all["pangenie_panel_haplotypes"] == 3 and r_all["pangenie_states"] == 9
