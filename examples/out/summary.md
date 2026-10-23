| Example | Description | \|V\| | \|E\| | s-t paths N | P-nodes | \|V(G_D)\| | new vertices | naive N(N-1) | exact |
|---|---|---|---|---|---|---|---|---|---|
| 01_bubble | One bubble (a SNP): s -> {a, b} -> t. | 4 | 4 | 2 | 1 | 6 | 2 | 2 | yes |
| 02_three_alleles | One bubble with three alleles. | 5 | 6 | 3 | 1 | 11 | 6 | 6 | yes |
| 03_two_bubbles | Two bubbles in series: pairs differ in either, both, or none. | 7 | 8 | 4 | 2 | 11 | 4 | 12 | yes |
| 04_nested | A bubble nested inside one allele of an outer bubble (an insertion carrying a SNP). | 6 | 7 | 3 | 2 | 10 | 4 | 6 | yes |
| 05_random_sp | Random nested series-parallel graph (seed 7). | 20 | 27 | 9 | 6 | 40 | 20 | 72 | yes |
| 06_w_motif | Smallest non-series-parallel graph (W-motif) in series with a bubble; core handled by pair states. | 7 | 9 | 6 | 1 | 23 | 16 | 30 | yes |
