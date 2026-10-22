import random

import pytest

from diploidgraph import build_diploid
from diploidgraph.generators import random_dag, random_sp_dag
from diploidgraph.support import observed_support


def brute_support(g, gd, observed):
    """New vertices on some G_D walk that decodes to an ordered pair of distinct observed haplotypes."""
    s, t = g.terminals()
    obs = set(observed)
    ext_cache = {}

    def ext(seg):
        if seg is None:
            return [()]
        if seg not in ext_cache:
            ext_cache[seg] = [p[1:] for p in g.all_paths(seg.source, seg.sink, seg.node.edges())]
        return ext_cache[seg]

    found = set()
    stack = [(gd.source, [((s,), (s,))], ())]
    while stack:
        u, partial, on_walk = stack.pop()
        if u == gd.sink:
            if any(a != b and a in obs and b in obs for a, b in partial):
                found.update(on_walk)
            continue
        for v in gd.succ[u]:
            if v in gd.gadget:
                gad = gd.gadget[v]
                new = [(a + e1, b + e2) for a, b in partial for e1 in ext(gad.seg1) for e2 in ext(gad.seg2)]
            elif u[0] == "hom" and v[0] == "hom":
                new = [(a + (v[1],), b + (v[1],)) for a, b in partial]
            else:
                new = partial
            # keep only partial pairs that can still become an observed pair (prefix check)
            new = [(a, b) for a, b in new if any(h[:len(a)] == a for h in obs) and any(h[:len(b)] == b for h in obs)]
            if new:
                stack.append((v, new, on_walk + ((v,) if v[0] != "hom" else ())))
    return found


@pytest.mark.parametrize("seed", range(40))
def test_support_matches_brute_force(seed):
    rng = random.Random(seed)
    g = random_dag(n_ops=7, extra_edges=seed % 4, seed=seed) if seed % 3 else random_sp_dag(9, seed=seed)
    gd, dec = build_diploid(g)
    s, t = g.terminals()
    paths = g.all_paths(s, t)
    k = min(len(paths), 2 + seed % 5)
    observed = rng.sample(paths, k)
    # duplicates and a fragment must be ignored
    walks = observed + observed[:1] + [observed[0][:-1]]
    supported, counts = observed_support(g, walks, dec)
    fast = set().union(*supported.values())
    assert fast == brute_support(g, gd, observed), seed
    assert counts["distinct_full"] == k
    assert fast <= set(gd.succ)  # every supported key is a real vertex of G_D


def test_all_paths_observed_supports_everything():
    g = random_dag(n_ops=6, extra_edges=2, seed=3)
    gd, dec = build_diploid(g)
    s, t = g.terminals()
    supported, _ = observed_support(g, g.all_paths(s, t), dec)
    new_vertices = {v for v in gd.succ if v[0] != "hom"}
    assert set().union(*supported.values()) == new_vertices


@pytest.mark.parametrize("seed", range(20))
def test_partial_walks_only_add(seed):
    from diploidgraph.support import observed_support_with_partial
    rng = random.Random(seed)
    g = random_dag(n_ops=10, extra_edges=1 + seed % 3, seed=seed)
    gd, dec = build_diploid(g)
    s, t = g.terminals()
    paths = g.all_paths(s, t)
    observed = rng.sample(paths, min(len(paths), 3))
    others = [p for p in paths if p not in observed]
    partial = [p[1:] for p in others[:2]] + [p[:-1] for p in others[2:4]]
    full_only, with_partial, counts = observed_support_with_partial(g, observed + partial, dec)
    plain, _ = observed_support(g, observed, dec)
    assert full_only == plain
    for k in full_only:
        assert full_only[k] <= with_partial[k] <= {v for v in gd.succ if v[0] == k}
    # making the partial walks complete can only add more
    complete, _ = observed_support(g, observed + others[:4], dec)
    for k in full_only:
        assert with_partial[k] <= complete[k]
