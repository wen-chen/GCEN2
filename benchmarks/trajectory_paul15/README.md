# GCEN2 trajectory networks vs an R implementation

How does co-expression change along a differentiation trajectory?
`gcen2 trajectory` orders the cells by pseudotime, builds a network in
overlapping windows and follows every edge across them. No published package
does this (tools for pseudotime data infer regulatory networks or lagged
correlations, a different question), so GCEN2 is compared with a direct R
implementation of the same definition (`02_reference.R`): WGCNA's `cor` and
`corPvalueStudent`, `p.adjust` and data.table, as an R user would write it.

**Data.** Paul et al., Cell 2015: 2,730 mouse bone marrow myeloid progenitor
cells (scanpy's `paul15`), counts per 10,000, the 2,000 most variable genes;
diffusion pseudotime rooted in the megakaryocyte-erythroid progenitors, as in
the scanpy DPT tutorial.

**Settings.** The same for both: windows of a fifth of the cells (546),
shifted by half a window (9 windows); metacells of 5 consecutive cells;
log2(x + 1); Pearson; |r| >= 0.5 and Benjamini-Hochberg FDR <= 0.001 over all
gene pairs; 8 threads.

## Run

Data: scanpy's copy of the paul15 data,
`https://exampledata.scverse.org/scanpy/paul15.h5` (sha256
`6161984f758dd464992edc23f1a8ab89b2081600130c6aa93a1ab12b3ada5bfb`);
01_prepare.py downloads it to `data/paul15/` if it is not there.

```bash
cd <data directory>
python .../01_prepare.py      # work/paul15_pt.h5ad, paul15_expr.tsv, paul15_pt.tsv
bash .../run_all.sh           # both, /usr/bin/time
python .../03_compare.py      # results/compare.tsv
```

## Results

`results/compare.tsv` (times: mean of two runs):

| | GCEN2 | R |
|---|---|---|
| time, peak memory | **1.2 s, 0.32 GB** | 6.1 s, 1.09 GB |
| edges per window | 349, 301, 608, 844, 1,216, 2,097, 2,693, 3,967, 9,931 | the same edges |
| correlations of the window networks | | within 5e-07 (GCEN2 writes 6 decimals) |
| edge dynamics (15,017 edges) | | r in every window within 1e-14; the same pattern, first and last window; trend within 4e-15 |

Patterns: 9,923 gained, 4,745 transient, 341 lost, 6 intermittent, 2 stable.
The networks grow towards the end of the trajectory, where the last window
holds the differentiated cells (pseudotime up to 1).
