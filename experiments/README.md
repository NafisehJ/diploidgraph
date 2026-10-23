# Experiments

All scripts write to `experiments/results/` (not tracked by git).

## 1. Synthetic scaling

Random two-terminal series-parallel DAGs of increasing size (several seeds per size).
Measures decomposition and construction time, peak Python memory (`tracemalloc`),
`|V(G_D)|`, `|E(G_D)|`, and the naive baseline `N(N-1)` (ordered pairs of
source-to-sink paths), reported as log10 because `N` grows exponentially.

```bash
python experiments/synthetic_scaling.py --out experiments/results/synthetic.csv
python experiments/plot_synthetic.py experiments/results/synthetic.csv --out experiments/results/synthetic.png
```

## 2. Real pangenome loci

Graphs: HPRC v2.0, Minigraph-Cactus, GRCh38-based, `ec1M` subgraphs from
<https://public.gi.ucsc.edu/cgl/ci/vg/vg-data/bakeoff/>.

```bash
# vg (bioconda build works on aarch64)
conda create -n vgenv -c conda-forge -c bioconda vg -y
conda activate vgenv

curl -L -O https://public.gi.ucsc.edu/cgl/ci/vg/vg-data/bakeoff/hprc-v2.0-mc-grch38.ec1M-BRCA1.gbz
curl -L -O https://public.gi.ucsc.edu/cgl/ci/vg/vg-data/bakeoff/hprc-v2.0-mc-grch38.ec1M-MHC.gbz
vg convert -f hprc-v2.0-mc-grch38.ec1M-BRCA1.gbz > BRCA1.gfa
vg convert -f hprc-v2.0-mc-grch38.ec1M-MHC.gbz > MHC.gfa
```

Then, in the environment where `diploidgraph` is installed:

```bash
python experiments/real_graphs.py BRCA1.gfa MHC.gfa --out experiments/results/real.csv
```

Reported per graph: input size, haplotype walks (total and distinct), whether
the graph is series-parallel, P-nodes and their maximum branching, time per
stage (read, decompose, build), peak memory, `|V(G_D)|`, `|E(G_D)|`,
new vertices by kind (P-node pairs, core forks, pair states, steps), the size of
the non-series-parallel core, the naive baselines (all graph paths, and observed
haplotypes), and `exactness_passed` (weighted walk count equals N^2).

## 3. Which new vertices are used by observed haplotype pairs

A new vertex of the diploid graph is *supported* if the walk of some ordered
pair of distinct observed haplotypes (full source-to-sink walks in the GFA)
passes through it; the others stand for pairs that arise only by recombination.
Reported per vertex kind (P-node pairs, core forks, pair states, steps).

```bash
python experiments/observed_support.py BRCA1.gfa --json brca1_support.json
python experiments/observed_support.py MHC.gfa --json mhc_support.json
```

The method never builds the diploid graph; it is checked against a brute-force
enumeration in `tests/test_support.py`.

## 4. Per-pair exactness check on real graphs

Builds the diploid graph and, for sampled ordered pairs of haplotypes, counts the
walks that decode to each pair (an exact graph gives exactly 1). Pairs come from
the observed haplotypes and from random source-to-sink paths, which include
recombinant pairs. Unlike the weighted count, this catches a missing pair that is
compensated by a duplicated one.

```bash
python experiments/pair_check.py BRCA1.gfa --pairs 500 --json brca1_pairs.json
```

## 5. Dependence on the vertex order

The P-node part of the diploid graph does not depend on the vertex order; the
core part does, because the haplotype that is earlier in a fixed topological
order moves first. This reports the size for the GFA order and for random orders.

```bash
python experiments/order_sensitivity.py BRCA1.gfa --orders 20 --json brca1_order.json
```

## 6. Diagnostics for the non-series-parallel core

* `core_blocks.py` splits the core at vertices that every haplotype passes and
  reports the block sizes (the pair states of the construction never span two blocks).
* `core_pairs.py` counts pairs of internally disjoint core paths between every
  pair of core vertices. It shows why the construction does not add one vertex per
  disjoint path pair: on BRCA1 it stops after more than 2,000,000 core paths.

## Everything at once

`bash reproduce.sh` (repository root) runs sections 2-5 for BRCA1 and MHC and
prints a summary table; `bash reproduce.sh --quick` runs BRCA1 only.
