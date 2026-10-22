# diploidgraph

[![tests](https://github.com/NafisehJ/spqr-diploid-genotyper/actions/workflows/tests.yml/badge.svg)](https://github.com/NafisehJ/spqr-diploid-genotyper/actions/workflows/tests.yml)

Compact **diploid sequence graphs** from two-terminal DAG sequence graphs.

A haploid sequence graph encodes haplotypes as source-to-sink paths. A diploid
sequence graph encodes **ordered pairs of haplotypes** as walks. The naive
construction creates one object per pair of haplotypes, which grows
quadratically in the number of haplotypes and exponentially in the number of
variant sites. This package keeps homozygous regions as a single copy of the
input graph and adds new vertices only where the two haplotypes can differ: at
the terminals of parallel structures (P-nodes) found by a series-parallel /
SPQR-tree decomposition.

Paper: N. Jafarzadeh, J. M. Eizenga, B. Paten. *An Efficient Graph Algorithm for
Diploid Local Ancestry Inference.* bioRxiv 2025.
[doi:10.1101/2025.07.05.662656](https://doi.org/10.1101/2025.07.05.662656)

## Status

| Component | State |
|---|---|
| Input checks (two-terminal DAG), GFA reader (S, L, P, W lines) | done |
| Series-parallel decomposition (linear-time series/parallel reductions) | done |
| Diploid construction, series-parallel graphs | done, verified exact |
| Non-series-parallel graphs (R-nodes, W-motifs, complexity graph C(G)) | in progress |
| Experiments: synthetic scaling, real pangenome loci | scripts ready |

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
diploidgraph graph.gfa --json stats.json
```

```python
from diploidgraph import SequenceGraph, build_diploid
from diploidgraph.validate import brute_force_check

g = SequenceGraph.from_edges([("s", "a"), ("a", "t"), ("s", "b"), ("b", "t")])
gd, decomposition = build_diploid(g)
print(gd.num_vertices(), gd.num_het_vertices())   # 6 2
print(brute_force_check(g, gd)["exact"])           # True
```

For an input that is not series-parallel, `build_diploid` raises
`NonSeriesParallelError`; `--allow-partial` (or `allow_partial=True`) builds the
series-parallel part only and reports `pair_coverage`, the fraction of haplotype
pairs represented.

## The construction

* **Homozygous copy.** One vertex `(v, v)` for every vertex `v` of the input and
  one edge `(u, u) -> (v, v)` for every input edge.
* **Heterozygous vertices.** For every P-node with terminals `(x, y)` and
  branches `B_1, ..., B_k`, one vertex `(B_i, B_j)` for every ordered pair
  `i != j`, with edges `(x, x) -> (B_i, B_j) -> (y, y)`. It means "haplotype 1
  crosses from x to y inside B_i and haplotype 2 inside B_j".
* **Size.** `|V(G_D)| = |V| + sum_P k_P (k_P - 1)` and
  `|E(G_D)| = |E| + 2 sum_P k_P (k_P - 1)`: linear in the input for bounded
  branching, independent of the number of haplotypes.

## Correctness checks

A walk of `G_D` decodes to a set of ordered haplotype pairs. `G_D` is **exact**
if every ordered pair of source-to-sink paths of the input is decoded by exactly
one walk, and nothing else is decoded.

* `validate.brute_force_check` enumerates all pairs and all walks on small
  graphs and reports missing, spurious and duplicated pairs.
* `validate.weighted_walk_count` is a linear-time count that equals `N^2`
  (`N` = number of source-to-sink paths) for an exact `G_D`; it is used on large
  graphs where enumeration is impossible.

The test suite checks exactness on random nested series-parallel graphs (brute
force) and on large ones (weighted count). Tests run on every push (see badge).

## Examples

`python examples/run_examples.py --png` writes the input graph, the diploid
graph and a summary table for each example to `examples/out/`.

| Example | |V| | s-t paths N | P-nodes | |V(G_D)| | het vertices | naive N(N-1) | exact |
|---|---|---|---|---|---|---|---|
| One bubble (SNP) | 4 | 2 | 1 | 6 | 2 | 2 | yes |
| Three alleles | 5 | 3 | 1 | 11 | 6 | 6 | yes |
| Two bubbles in series | 7 | 4 | 2 | 11 | 4 | 12 | yes |
| Nested bubble (insertion carrying a SNP) | 6 | 3 | 2 | 10 | 4 | 6 | yes |
| Random nested SP graph | 18 | 7 | 4 | 34 | 16 | 42 | yes |
| W-motif (non-SP) | 7 | 6 | 1 | - | - | 30 | pending |

Nested bubble: input (left) and diploid graph (right). Blue vertices are
homozygous; orange vertices are heterozygous pairs of branches, written with
nested alternatives as `(x|y)`.

<p>
<img src="examples/out/04_nested_input.png" height="160">
<img src="examples/out/04_nested_diploid.png" height="220">
</p>

## Experiments

See [experiments/README.md](experiments/README.md) for:

* runtime, peak memory and diploid graph size vs. naive pairing on random
  series-parallel graphs (`synthetic_scaling.py`, `plot_synthetic.py`);
* the same measurements on real gene/locus-level pangenome subgraphs (BRCA1,
  MHC; HPRC v2.0 Minigraph-Cactus) (`real_graphs.py`).

## Repository layout

```
src/diploidgraph/   package: seqgraph, gfa, decompose, diploid, validate, stats, export, cli
tests/              pytest suite
examples/           worked examples and figures
experiments/        scripts for the paper's measurements
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
