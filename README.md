# diploidgraph

[![tests](https://github.com/NafisehJ/spqr-diploid-genotyper/actions/workflows/tests.yml/badge.svg)](https://github.com/NafisehJ/spqr-diploid-genotyper/actions/workflows/tests.yml)

Compact **diploid sequence graphs** from two-terminal DAG sequence graphs.

A haploid sequence graph encodes haplotypes as source-to-sink paths. A diploid
sequence graph encodes **ordered pairs of haplotypes** as walks: every pair of
source-to-sink paths, including recombinant pairs that no sample carries.
`diploidgraph` keeps homozygous regions as a single copy of the input graph and
adds vertices only where two haplotypes can differ: one vertex per ordered pair
of branches at each parallel structure (P-node) of a series-parallel
decomposition, and, for graphs that are not series-parallel, pair states inside
the irreducible core that contains their W-motifs.

On two HPRC v2.0 pangenome loci (BRCA1 and MHC):

* **Compact:** 2.3x (BRCA1) and 3.1x (MHC) fewer vertices than the product-graph
  baseline, which represents exactly the same haplotype pairs without the
  decomposition. Pairing the observed haplotypes one by one would need about
  4 x 10^4 to 5 x 10^4 times more vertices (different units; see below).
* **Exact:** every ordered pair of haplotypes is represented exactly once, checked
  by brute force on small graphs, pair by pair on sampled real haplotype pairs,
  and by a whole-graph path count.
* **Fast:** 0.03 s for BRCA1 and under 10 s for MHC on a laptop.

<img src="experiments/figures/real_loci.png" width="100%">

Paper: N. Jafarzadeh, J. M. Eizenga, B. Paten. *An Efficient Graph Algorithm for
Diploid Local Ancestry Inference.* bioRxiv 2025.
[doi:10.1101/2025.07.05.662656](https://doi.org/10.1101/2025.07.05.662656)

## Results on real pangenome loci

HPRC v2.0 Minigraph-Cactus subgraphs (GRCh38-based, `ec1M`). Both inputs are
two-terminal DAGs with no preprocessing: no link or haplotype walk was dropped.
Measured on a laptop (WSL Ubuntu 22.04, aarch64, 8 GB RAM), single run.

| | BRCA1 | MHC |
|---|---|---|
| Input vertices / edges | 2,704 / 3,828 | 349,475 / 482,220 |
| Haplotype walks: all / distinct source-to-sink | 465 / 453 | 466 / 462 |
| Dropped links / walks with mixed orientation | 0 / 0 | 0 / 0 |
| P-nodes (max branches) | 982 (3) | 114,540 (4) |
| Non-series-parallel core: vertices / edges | 203 / 340 | 28,664 / 45,545 |
| **Diploid graph vertices / edges** | **8,250 / 12,592** | **887,929 / 1,361,520** |
| Diploid vertices per input vertex | 3.05 | 2.54 |
| Product-graph baseline vertices / edges (same pairs) | 19,342 / pending | 2,709,299 / 3,352,592 |
| Naive pairing of observed walks, paired positions | 354,168,814 | 46,485,747,514 |
| Ordered pairs of source-to-sink paths represented | about 10^570 | about 10^71057 |
| Diploid vertices over 21 vertex orders (min / median / max) | 7,678 / 8,158 / 8,842 | pending |
| Sampled pairs decoded exactly once (observed; random) | 100/100; 100/100 | pending |
| Whole-graph check (weighted walk count = N^2) | passed (exact) | passed (modulo two large primes) |
| Added vertices used by pairs of complete observed haplotypes | 100% | 99.86% |
| Decompose + build time | 0.03 s | 9.6 s |
| Peak memory (whole process) | 27 MB | 1.2 GB |

How to read the comparisons:

* **Product-graph baseline** is the like-for-like comparison: it applies the same
  pair-state construction to every vertex of the input instead of only to the
  core. It represents exactly the same pairs, and its size is counted in the same
  unit. The difference is the saving due to the decomposition; it is larger where
  bubbles are nested or have long branches.
* **Naive pairing** stores one chain of paired positions for every ordered pair of
  distinct observed walks. It is not in the same unit as the diploid graph: one
  diploid vertex can stand for many haplotype pairs, while one naive position
  stands for one position of one pair. It mainly shows how the cost of listing
  pairs grows with the number of haplotypes.
* **Vertex order.** The P-node part does not depend on the vertex order. In the
  core, the haplotype that is earlier in a fixed topological order moves first,
  so the core part depends on the order; the table gives the range over the GFA
  order and 20 random orders.
* **Observed haplotypes.** On a Minigraph-Cactus graph every bubble branch is
  used by some input haplotype, so the P-node part is fully used by construction;
  the figure is informative for the core only. The 0.14% of MHC vertices not used
  by complete haplotype pairs are used by the 4 walks that cover the region only
  partly.

Reading the MHC GFA takes about 50 s in pure Python and is not included in the
build time.

## Reproduce the results

```bash
bash reproduce.sh           # BRCA1 and MHC: download, convert, build, verify, measure (about 15 min)
bash reproduce.sh --quick   # BRCA1 only (under a minute)
```

The script downloads the two HPRC graphs, converts them to GFA with
[vg](https://github.com/vgteam/vg) (only needed once), builds and verifies the
diploid graphs, measures which added vertices are used by the observed
haplotypes, and prints the table above. Sizes are reproduced exactly; times and
memory depend on the machine. If the GFA files already exist, point to them with
`DATA=/path/to/folder bash reproduce.sh` and vg is not needed.

## Install

```bash
git clone https://github.com/NafisehJ/spqr-diploid-genotyper.git
cd spqr-diploid-genotyper
python -m pip install -e ".[test]"
python -m pytest
```

Python 3.10 or newer. No dependencies for the core package; `matplotlib` for
plots and Graphviz `dot` for example figures are optional.

## Quick start

```bash
diploidgraph graph.gfa --json stats.json                 # build, verify, report sizes
diploidgraph graph.gfa --count-only --json stats.json    # sizes only, for very large graphs
```

```python
from diploidgraph import SequenceGraph, build_diploid
from diploidgraph.validate import brute_force_check

g = SequenceGraph.from_edges([("s", "a"), ("a", "t"), ("s", "b"), ("b", "t")])
gd, decomposition = build_diploid(g)
print(gd.num_vertices(), gd.num_het_vertices())   # 6 2
print(brute_force_check(g, gd)["exact"])           # True
```

## The construction

* **Homozygous copy.** One vertex `(v, v)` for every vertex `v` of the input and
  one edge `(u, u) -> (v, v)` for every input edge.
* **Heterozygous vertices.** For every P-node with terminals `(x, y)` and
  branches `B_1, ..., B_k`, one vertex `(B_i, B_j)` for every ordered pair
  `i != j`, with edges `(x, x) -> (B_i, B_j) -> (y, y)`. It means "haplotype 1
  crosses from x to y inside B_i and haplotype 2 inside B_j".
* **Non-series-parallel core.** Series and parallel reductions leave an
  irreducible core: a DAG whose edges each stand for a series-parallel
  subgraph. Two haplotypes that leave a core vertex along different core edges
  can meet again only at a core vertex. The core part of `G_D` follows such a
  divergence with pair states `[u | v]` (haplotype 1 at core vertex `u`,
  haplotype 2 at `v`), moving the haplotype that is earlier in a fixed
  topological order one core edge at a time, until both reach the same vertex.
  This makes the walk of every pair unique.
* **Size.** For series-parallel inputs, `|V(G_D)| = |V| + sum_P k_P (k_P - 1)`
  and `|E(G_D)| = |E| + 2 sum_P k_P (k_P - 1)`. A core with `n_c` vertices adds
  at most `n_c (n_c - 1)` pair states plus forks and steps bounded by core edges
  times core vertices. Nothing depends on the number of haplotypes.

## Correctness checks

A walk of `G_D` decodes to a set of ordered haplotype pairs. `G_D` is **exact**
if every ordered pair of source-to-sink paths of the input is decoded by exactly
one walk, and nothing else is decoded.

* `validate.brute_force_check` enumerates all pairs and all walks on small
  graphs and reports missing, spurious and duplicated pairs.
* `validate.weighted_walk_count` is a linear-time count that equals `N^2`
  (`N` = number of source-to-sink paths) for an exact `G_D`; it is used on large
  graphs where enumeration is impossible.

The test suite checks exactness by brute force on 40 random series-parallel
graphs and 60 random non-series-parallel graphs, and by the weighted count on
large graphs of both kinds. Tests run on every push (see badge).

## Examples

Six small worked examples (single bubble, three alleles, bubbles in series, a
nested bubble, a random nested graph, and a W-motif), each with input and
diploid graph figures and an exactness check: [examples/README.md](examples/README.md).

## Experiments

See [experiments/README.md](experiments/README.md) for:

* runtime, peak memory and diploid graph size vs. naive pairing on real
  gene/locus-level pangenome subgraphs (`real_graphs.py`, `plot_real.py`);
* which added vertices are used by observed haplotype pairs
  (`observed_support.py`, `explain_unsupported.py`);
* scaling on random series-parallel graphs (`synthetic_scaling.py`, `plot_synthetic.py`).

## Repository layout

```
src/diploidgraph/   package: seqgraph, gfa, decompose, diploid, validate, support, stats, export, cli
tests/              pytest suite
examples/           worked examples and figures
experiments/        scripts and figures for the paper's measurements
notebooks/          original prototype notebooks (superseded by the package)
```

## Citation

```bibtex
@article{jafarzadeh2025diploid,
  title   = {An Efficient Graph Algorithm for Diploid Local Ancestry Inference},
  author  = {Jafarzadeh, Nafiseh and Eizenga, Jordan M. and Paten, Benedict},
  journal = {bioRxiv},
  year    = {2025},
  doi     = {10.1101/2025.07.05.662656}
}
```
