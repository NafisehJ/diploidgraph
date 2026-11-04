# Worked examples

Eight small sequence graphs, from a single SNP to a non-series-parallel W-motif,
including the worked example of the paper and one graph read from a GFA file.
For each one: the input graph, its diploid graph, and an exactness check (every
ordered pair of source-to-sink paths of the input is decoded by exactly one walk
of the diploid graph, verified by brute force).

`python examples/run_examples.py --png` (run from the repository root) regenerates
the figures and the table below in `examples/out/`.

| Example | \|V\| | s-t paths N | P-nodes | **\|V(G_D)\| (ours)** | \|V\| of Cartesian product | exact |
|---|---|---|---|---|---|---|
| 1. One bubble (SNP) | 4 | 2 | 1 | **6** | 14 | yes |
| 2. Three alleles | 5 | 3 | 1 | **11** | 31 | yes |
| 3. Two bubbles in series and a large deletion | 7 | 5 | 3 | **13** | 49 | yes |
| 4. Nested bubble (insertion carrying a SNP) | 6 | 3 | 2 | **10** | 36 | yes |
| 5. Random nested series-parallel graph | 18 | 7 | 4 | **34** | 262 | yes |
| 6. Bubble followed by a W-motif (non-series-parallel) | 7 | 6 | 1 | **23** | 31 | yes |
| 7. Worked example of the paper (Figure 1) | 6 | 4 | 1 | **22** | 46 | yes |
| 8. Read from a GFA file (SNP and a 16 bp deletion) | 6 | 4 | 2 | **10** | 22 | yes |

**Cartesian product** is the Cartesian product with a fixed order of moves: the
same pair-state construction applied to the whole input, without decomposition.
It represents exactly the same haplotype pairs and is counted in the same unit
(vertices), so it shows the saving due to the decomposition.

**How to read the figures.** Blue ellipses `(v,v)` are the homozygous copy of the
input: both haplotypes at `v`. Orange boxes `(x, y)` are heterozygous vertices of a
P-node: haplotype 1 crosses the bubble through branch `x`, haplotype 2 through
branch `y` (`-` is a direct edge; nested alternatives are written `(b|c)`). In the
non-series-parallel example, a green box `(c, d)` is where the two haplotypes
separate (haplotype 1 goes to `c`, haplotype 2 to `d`), a grey diamond `[u | v]` is
a pair state (haplotype 1 at `u`, haplotype 2 at `v`), and an arrow labelled
`h1 → t` moves haplotype 1 to `t` while haplotype 2 waits, until both meet again.

## 1. One bubble (SNP)

Two alleles `a` and `b`. The diploid graph adds the two heterozygous vertices `(a, b)` and `(b, a)`.

<img src="out/01_bubble_input.png" height="110"> &nbsp;&nbsp; <img src="out/01_bubble_diploid.png" height="170">

## 2. Three alleles

A bubble with three branches adds one vertex per ordered pair of distinct branches: 3 x 2 = 6.

<img src="out/02_three_alleles_input.png" height="130"> &nbsp;&nbsp; <img src="out/02_three_alleles_diploid.png" height="260">

## 3. Two bubbles in series and a large deletion

Two bubbles in series, plus a direct edge `s -> t` that skips both: a large
deletion. The two haplotypes can differ at either bubble, at both, or at none,
and either one can carry the deletion. The deletion makes one more P-node with
two branches (the direct edge and the rest of the graph), so it adds only two
vertices: `(-, rest)` and `(rest, -)`. In all, 13 vertices against 49 for the
Cartesian product.

<img src="out/03_two_bubbles_input.png" height="110"><br>
<img src="out/03_two_bubbles_diploid.png" height="260">

## 4. Nested bubble (insertion carrying a SNP)

One allele of the outer bubble (`s -> t`) is an insertion `a (b|c) d` that itself
carries a SNP. The inner bubble is expanded inside the insertion; the outer bubble
pairs the whole insertion with the direct edge.

<img src="out/04_nested_input.png" height="130"> &nbsp;&nbsp; <img src="out/04_nested_diploid.png" height="210">

## 5. Random nested series-parallel graph

A larger random example with four nested bubbles.

<img src="out/05_random_sp_input.png" height="320"><br>
<img src="out/05_random_sp_diploid.png" width="100%">

## 6. Bubble followed by a W-motif (non-series-parallel)

A bubble `s -> {a, b} -> m`, then a W-motif with diagonal edge `c -> d`. The bubble
is handled as in example 1; inside the W-motif the two haplotypes are followed with
pair states until they meet again at `d` or `t`. The 6 single-haplotype moves are drawn as
labelled arrows; in the table they are counted as vertices (23 vertices in total).

<img src="out/06_w_motif_input.png" height="130"><br>
<img src="out/06_w_motif_diploid.png" width="100%">

## 7. Worked example of the paper (Figure 1)

The input graph of Figure 1 of the paper: edges `u -> u'`, `u -> v'`, `u' -> c`,
`u' -> d`, `c -> v'`, `d -> v'`, `u' -> v`, `v' -> v`. It is not series-parallel.
The diploid graph has 22 vertices and 32 edges: the homozygous copy (6), the
heterozygous vertices of the P-node between `u'` and `v'` (2), and the
divergence vertices and pair states of the core (14, counting the 6
single-haplotype moves as vertices, as in the table). Brute force confirms that
each of the N^2 = 16 ordered pairs of source-to-sink paths (identical pairs included)
is decoded by exactly one walk.

<img src="out/07_paper_figure1_input.png" height="130"><br>
<img src="out/07_paper_figure1_diploid.png" width="100%">

## 8. Read from a GFA file

A small GFA file, [data/locus.gfa](data/locus.gfa), in the format of a
pangenome graph: six segments, a SNP (segments `2` and `3`, `A` or `G`) and a
16 bp deletion (the link `4 -> 6` skips segment `5`), and five haplotype walks
(`W` lines).

```
S	1	ACGTTGCA
S	2	A
S	3	G
S	4	CCTAGGTA
S	5	TTGACCATGGCATTAG
S	6	GATC
L	1	+	2	+	0M
L	1	+	3	+	0M
L	2	+	4	+	0M
L	3	+	4	+	0M
L	4	+	5	+	0M
L	5	+	6	+	0M
L	4	+	6	+	0M
W	GRCh38	0	chr1	0	37	>1>2>4>5>6
W	HG001	1	chr1	0	21	>1>3>4>6
W	HG001	2	chr1	0	37	>1>2>4>5>6
W	HG002	1	chr1	0	37	>1>3>4>5>6
W	HG002	2	chr1	0	21	>1>2>4>6
```

Build the diploid graph, check it and report its size from the command line:

```bash
diploidgraph examples/data/locus.gfa --json stats.json
```

Main fields of the output (times and memory depend on the machine):

```
"n_vertices": 6,
"n_edges": 7,
"is_series_parallel": true,
"n_pnodes": 2,
"diploid_vertices": 10,
"diploid_edges": 15,
"exactness_check": "exact",
```

The same in Python:

```python
from diploidgraph import build_diploid
from diploidgraph.gfa import read_gfa

g = read_gfa("examples/data/locus.gfa").graph
gd, decomposition = build_diploid(g)
print(gd.num_vertices())   # 10
```

In the diploid graph, `(3, 2)` and `(2, 3)` are the two heterozygous states of
the SNP, and `(-, 5)` and `(5, -)` are the two states where exactly one
haplotype carries the deletion.

<img src="out/08_gfa_input.png" height="110"><br>
<img src="out/08_gfa_diploid.png" height="230">
