# diploidgraph

Compact diploid sequence graphs from two-terminal DAG sequence graphs.

Companion code for: N. Jafarzadeh, J. M. Eizenga, B. Paten,
*Diploid Sequence Graphs: Compact Representation of Haplotype Pairs via SPQR-Tree Decomposition* (RECOMB 2027 submission).

## Status

| Part | State |
|---|---|
| Two-terminal DAG checks, GFA reader | done |
| Series-parallel decomposition (P-nodes = hetero-positions) | done, linear-time reductions |
| Diploid construction, series-parallel case | done, verified exact |
| Non-series-parallel case (R-nodes, W-motifs, C(G)) | [TODO] |
| Experiment scripts (BRCA1, MHC) | [TODO] |

## Install and test

```bash
pip install -e ".[test]"
pytest
```

## Usage

```bash
diploidgraph BRCA1.gfa --allow-partial --json brca1_stats.json
```

Reports input size, number of P-nodes and R-nodes, diploid graph size,
naive baselines (all graph paths, and observed haplotypes), runtime per stage,
peak Python memory, and `pair_coverage`: the fraction of ordered haplotype pairs
represented by the diploid graph (1.0 means complete; below 1.0 means some pairs
diverge inside R-node regions, which are not extended yet).

## Definitions used in the code

* **Diploid graph G_D.** Homozygous vertex `("hom", v)` for each vertex v of G.
  For each P-node with terminals (x, y) and branches B_1..B_k, one heterozygous
  vertex `("het", p, i, j)` for each ordered pair i != j, placed between
  `("hom", x)` and `("hom", y)`. It stands for "haplotype 1 crosses x..y through
  B_i and haplotype 2 through B_j".
* **Size.** |V(G_D)| = |V| + sum_p k_p(k_p - 1), |E(G_D)| = |E| + 2 sum_p k_p(k_p - 1).
* **Exactness.** Every ordered pair of s-t paths of G is decoded by exactly one
  s-t walk of G_D, and nothing else is decoded.

## How correctness is checked

* `brute_force_check`: enumerates all ordered haplotype pairs and all walks of
  G_D on small graphs and compares them (missing / spurious / duplicated).
* `weighted_walk_count`: linear-time count that must equal N^2 (N = number of
  s-t paths) for an exact G_D; used on large graphs.

Tests: 40 random nested SP graphs (brute force) and large SP graphs (weighted count).
