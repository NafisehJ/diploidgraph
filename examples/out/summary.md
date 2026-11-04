| Example | Description | \|V\| | \|E\| | s-t paths N | P-nodes | **\|V(G_D)\| (ours)** | \|V\| of Cartesian product | exact |
|---|---|---|---|---|---|---|---|---|
| 01_bubble | One bubble (a SNP): s -> {a, b} -> t. | 4 | 4 | 2 | 1 | **6** | 14 | yes |
| 02_three_alleles | One bubble with three alleles. | 5 | 6 | 3 | 1 | **11** | 31 | yes |
| 03_two_bubbles | Two bubbles in series, and a large deletion s -> t that skips both. | 7 | 9 | 5 | 3 | **13** | 49 | yes |
| 04_nested | A bubble nested inside one allele of an outer bubble (an insertion carrying a SNP). | 6 | 7 | 3 | 2 | **10** | 36 | yes |
| 05_random_sp | Random nested series-parallel graph (seed 7). | 18 | 23 | 7 | 4 | **34** | 262 | yes |
| 06_w_motif | Smallest non-series-parallel graph (W-motif) in series with a bubble; core handled by pair states. | 7 | 9 | 6 | 1 | **23** | 31 | yes |
| 07_paper_figure1 | Worked example of the paper (Figure 1); non-series-parallel. | 6 | 8 | 4 | 1 | **22** | 46 | yes |
| 08_gfa | Synthetic locus read from a GFA file (examples/data/locus.gfa): SNPs, three alleles, a nested insertion, a W-motif and a large deletion. | 17 | 24 | 63 | 5 | **45** | 121 | yes |
| 09_brca1_window | BRCA1 window (HPRC v2.0, real data), cut with extract_window.py: a CA repeat. | 23 | 40 | 168 | 15 | **85** | 623 | yes |
