# GCEN2 vs hdWGCNA on PBMC 10k

Cell-type co-expression networks from the 10x Genomics PBMC 10k v3 data set
(10,685 cells after QC), built by GCEN2 (`build_by_celltype`, as
`gcen2 celltype`) and by hdWGCNA 0.4.12 (Seurat 5.5.1, WGCNA 1.74, R 4.4.3)
from the same inputs:

* the same cells and cell-type labels (KMeans on 30 PCs, clusters labelled
  by marker genes; six types are compared: CD4 T 3,391, CD14 Mono 3,333,
  B 1,496, CD8 T 1,061, NK 616 and FCGR3A Mono 355 cells);
* the same PCA (passed to Seurat as the `pca` reduction) and gene lists
  (genes detected in >= 5% of the cell type's cells, as hdWGCNA's
  `fraction = 0.05`; 7,651-8,985 genes);
* the same metacell settings: k = 25, max_shared = 10, at most 1,000 per
  cell type, mean of the members, then log-normalised.

GCEN2 networks (`gcen2 celltype`'s defaults but the cutoff): Pearson,
r >= 0.3, Benjamini-Hochberg FDR <= 0.001 over all gene pairs, the p-values on
every gene pair's effective sample size (from the cells the metacells share;
the median gene has 30-327 effective samples for 44-1,000 metacells),
positive edges only (355-122,020 edges), modules by `gcen2.module.identify`
(similarity 0.5). hdWGCNA: `TestSoftPowers` and `ConstructNetwork` with
their defaults (signed network, automatic soft power: 4-8; dynamic tree cut,
minimum module size 50).

## Data

In a data directory of your choice (`<pbmc>`), `data/` holds:

| file | source | SHA-256 |
|---|---|---|
| `pbmc_10k_v3_filtered_feature_bc_matrix.h5` | 10x Genomics, `https://cf.10xgenomics.com/samples/cell-exp/3.0.0/pbmc_10k_v3/pbmc_10k_v3_filtered_feature_bc_matrix.h5` (CC BY 4.0) | `ebc5dedc938830e20f8e1aafb893b0a9e0bf88584f3ef2d00b232dd277a188af` |
| `string/9606.protein.info.v12.0.txt.gz` | STRING v12, `https://stringdb-downloads.org/download/protein.info.v12.0/9606.protein.info.v12.0.txt.gz` | `144de4b0d98c6a7dfde6ddc2591cf88657f27b989eadff4f501450c3ed1f0f1c` |
| `string/9606.protein.links.detailed.v12.0.txt.gz` | STRING v12, `https://stringdb-downloads.org/download/protein.links.detailed.v12.0/9606.protein.links.detailed.v12.0.txt.gz` | `28a2ea3b22c964e691c0f2f934c5e3e397a4b0a12c833d8787fa02b7847b91b6` |

`04_compare.py` writes the STRING gold standard to `work/string_gold.tsv`,
which the GTEx benchmark reuses (`GCEN2_STRING_GOLD`).

## Run

Environments: `../environment-r.yml` (R, hdWGCNA) and
`../requirements-python.txt`; `RSCRIPT` is that environment's `Rscript`
(default: `Rscript` on the `PATH`). All scripts run from the data directory.

```bash
cd <pbmc>
python -I .../01_preprocess.py              # QC, PCA, labels, gene lists -> work/
bash .../run_all.sh                         # both tools, 24 runs
python -I .../04_compare.py                 # results/*.tsv
python -I .../05_tom_modules.py             # results/tom_modules.tsv
python -I .../07_cell_level.py runs/hdwgcna  # results/cell_level.tsv
python -I .../08_seed_ari.py                 # results/seed_ari_all_genes.tsv
MAGIC_LIB=<dir> python -I .../09_strategies.py  # results/strategies.tsv
python -I .../10_effective_samples.py         # results/effective_samples.tsv
python -I .../11_diff_calibration.py          # results/diff_calibration.tsv
```

hdWGCNA's `library()` attaches enrichR, which contacts maayanlab.cloud and
fails when that times out; `02_hdwgcna.R` attaches enrichR first, with
retries. Nothing else needs the network.

**hdWGCNA metacell seeds.** `MetacellsByGroups` draws the metacell centres
with `sample()`, but subsets the Seurat object per group right before, and
subsetting a Seurat object resets R's RNG (after `set.seed(1)` or
`set.seed(2)`, `so[, cells]` leaves the same `runif()`). A seed set by the
caller therefore has no effect: hdWGCNA always builds the same metacells.
`02_hdwgcna.R` seeds `ConstructMetacells` through `trace()` so that the seed
runs below differ (seed 1 reproduces hdWGCNA's usual metacells).

The STRING v12 gold standard uses only the `experimental` and `database`
channels (score >= 700; 106,444 pairs), so STRING's own co-expression evidence
plays no part. 8 threads on both sides.

## Results

Time and peak memory (whole script, reading the data included; hdWGCNA's
includes `TestSoftPowers` and the TOM, which GCEN2 does not compute; GCEN2:
mean of two runs, `results/gcen2_time_gemm.txt`):

| cell type (cells, genes) | GCEN2 | hdWGCNA | faster |
|---|---|---|---|
| CD14 Mono (3,333, 8,985) | 7.3 s, 1.09 GB | 120 s, 6.0 GB | 16x |
| CD4 T (3,391, 7,651) | 6.3 s, 1.09 GB | 91 s, 4.3 GB | 14x |
| B (1,496, 7,664) | 3.1 s, 0.59 GB | 78 s, 4.5 GB | 25x |
| CD8 T (1,061, 7,784) | 3.0 s, 0.61 GB | 76 s, 4.2 GB | 26x |
| NK (616, 8,509) | 1.9 s, 0.38 GB | 80 s, 5.2 GB | 43x |
| FCGR3A Mono (355, 9,192) | 1.5 s, 0.34 GB | 82 s, 6.0 GB | 55x |

GCEN2 correlates the up to 1,000 metacells with the matrix-product kernel
(`kernel="auto"`, from 64 samples). With the per-pair kernel of GCEN v1 the
runs took 12.7, 9.9, 4.3, 3.9, 2.5 and 1.6 s (`03_gcen2.py ... pairwise`,
mean of two runs, `results/gcen2_time.txt`; measured on 2026-11-17, before the
per-pair effective samples, which add 0.1-0.3 s, and the FDR, which needs the
candidate edges of all pairs: up to 0.4 GB more memory, and 0.3 s less for
the small cell types, whose networks and modules are smaller).
Both kernels give the same edges (in the same order) and modules; the
correlations differ by at most 2e-14 and the p-values by at most 2e-11
(relative). FCGR3A Mono has fewer than 64 metacells and uses the per-pair
kernel either way.

Gene-pair rankings: fraction of the top 1,000 / 10,000 pairs that are STRING
interactions (all-pair baseline 0.15-0.18%):

| cell type | GCEN2 r | hdWGCNA r | hdWGCNA TOM |
|---|---|---|---|
| CD14 Mono | 52.4% / 21.1% | 54.0% / 21.2% | 23.2% / 13.8% |
| CD4 T | 71.9% / 21.6% | 63.6% / 18.7% | 88.2% / 17.3% |
| B | 26.0% / 13.8% | 27.4% / 14.5% | 32.2% / 14.0% |
| CD8 T | 77.9% / 24.5% | 78.7% / 25.2% | 81.7% / 24.0% |
| NK | 77.7% / 29.9% | 74.4% / 29.6% | 61.5% / 23.8% |
| FCGR3A Mono | 34.5% / 12.6% | 37.3% / 14.5% | 18.8% / 11.0% |

Agreement of the gene-gene correlations (correlation of all pairs; overlap of
the top 10,000 pairs):

| cell type | GCEN2 vs hdWGCNA | GCEN2, seed 42 vs 1 / 2 | hdWGCNA, seed 42 vs 1 / 2 |
|---|---|---|---|
| CD14 Mono | 0.938; 93% | 0.939 / 0.940; 93% | 0.938 / 0.940; 93% |
| CD4 T | 0.933; 89% | 0.939 / 0.938; 90-92% | 0.938 / 0.939; 91-92% |
| B | 0.926; 89% | 0.923 / 0.922; 89% | 0.920 / 0.923; 89% |

The two tools differ by as much as one tool with another metacell seed: the
metacell implementations agree up to sampling noise (smaller types: NK
0.83, FCGR3A Mono 0.73 between tools, with 130 and 44 metacells).

Modules (number; genes in modules; STRING enrichment of within-module pairs;
ARI of the modules of two disjoint halves of the cells, on the genes in
modules in both; ARI across metacell seeds):

| cell type | GCEN2 | hdWGCNA |
|---|---|---|
| CD14 Mono | 74; 1,130; 11x; 0.93; 0.96 | 4; 4,084; 2.2x; 0.83; 0.98 |
| CD4 T | 87; 893; 16x; 0.79; 0.98-0.99 | 2; 3,592; 1.3x; 0.79; 0.99 |
| B | 43; 359; 153x; 0.62; 0.93-0.96 | 5; 3,109; 1.6x; 0.29; 0.93 |
| CD8 T | 27; 277; 134x; 0.81 | 3; 3,950; 1.7x; 0.36 |
| NK | 32; 268; 208x; 0.89 | 13; 4,927; 1.6x; 0.10 |
| FCGR3A Mono | 9; 63; 17x; - | 29; 5,943; 1.7x; 0.02 |

GCEN2's modules are many small, tight groups covering fewer genes; hdWGCNA's
are a few large modules covering half or more of the genes. Larger modules
dilute within-module enrichment, so that figure does not rank the methods on
its own. hdWGCNA's modules are reproducible across halves for the two large
cell types (0.79-0.83) but not below ~1,500 cells (0.02-0.36), where its
automatic soft power also becomes unstable (14 for the FCGR3A Mono halves);
GCEN2's are 0.62-0.93 on the (fewer) genes they cover. The FCGR3A Mono
halves (14 and 17 metacells) have 9 and 18 edges, and only two genes are in
modules in both. Module ARI between the tools is 0.52-0.70.

The FDR is `gcen2 celltype`'s default since 2026-11-18. The per-pair
effective sample sizes are calibrated (null-gene benchmark below): about
0.001 of unrelated pairs pass p <= 0.001, so p alone keeps tens of thousands
of chance edges among a cell type's 30-42 million pairs. With p <= 0.001
alone (`--no-fdr`) the networks had 23,245-302,654 edges and the small
cell types 200-390 modules, mostly of chance edges, with split-half ARI 0.41
(B) and 0.45 (FCGR3A Mono). Before the per-pair sizes, GCEN2 tested every pair
on cells / k samples (`--effective-samples cells`), which is too strict for
most gene pairs: 52-119 modules covering 265-2,295 genes, split-half ARI
0.79-0.87.

**hdWGCNA's modules with GCEN2** (`gcen2 module -e metacells.tsv --preset
hdwgcna -t 8`, `05_tom_modules.py`). The preset is what `ConstructNetwork`
runs (`blockwiseConsensusModules` on one set: signed TOM, deepSplit 4, no PAM
stage, minKMEtoStay 0.2, minimum size 50, merge cut 0.2; its gene
reassignment compares lower-tail chi-square p-values and never moves a gene)
with the power `TestSoftPowers` picks (pickSoftThreshold on bicor, lowest
power above 3 with R^2 >= 0.8). On hdWGCNA's own metacell matrix
(`datExpr`) GCEN2 picks the same power and gives exactly hdWGCNA's modules
in all 24 runs (adjusted Rand index 1.0 over all genes), in 3.1-6.3 s and
0.6-1.0 GB for the whole command against 36-60 s for hdWGCNA's
`TestSoftPowers` + `ConstructNetwork` alone (`t_network`), 8.9-13x faster
(`results/tom_modules.tsv`, GCEN2 mean of two runs):

| cell type | metacells | power | modules | GCEN2 | hdWGCNA network |
|---|---|---|---|---|---|
| CD14 Mono | 1,000 | 4 | 4 | 6.3 s, 0.97 GB | 59 s |
| CD4 T | 1,000 | 6 | 2 | 4.6 s, 0.75 GB | 42 s |
| B | 428 | 8 | 5 | 3.7 s, 0.64 GB | 38 s |
| CD8 T | 346 | 5 | 3 | 3.5 s, 0.63 GB | 40 s |
| NK | 130 | 4 | 13 | 3.7 s, 0.69 GB | 48 s |
| FCGR3A Mono | 44 | 5 | 29 | 4.4 s, 0.77 GB | 47 s |

(Pearson instead of bicor for the power table would pick another power in 9
of the 24 runs.) On GCEN2's own metacells (same settings, sampled
independently) the same command reproduces hdWGCNA's modules as closely as
hdWGCNA reproduces itself with another metacell seed: ARI with hdWGCNA's
modules 0.68 / 0.80 / 0.74 for B / CD14 Mono / CD4 T, against 0.70-0.71 /
0.81 / 0.75-0.78 between hdWGCNA's seeds (same ARI over all genes; with all
three seeds of each tool, medians 0.77 between the tools, 0.77 between
GCEN2's seeds and 0.78 between hdWGCNA's, `08_seed_ari.py`,
`results/seed_ari_all_genes.tsv`); it falls
with the number of metacells as hdWGCNA's own reproducibility does (NK 0.43,
FCGR3A Mono 0.18).

Split-half correlation of the gene-gene correlations is the same for both
tools and falls with the cell count (CD4 T 0.41, CD14 Mono 0.37, CD8 T
0.17-0.18, B 0.12, NK 0.09, FCGR3A Mono 0.03-0.04): below ~1,000 cells the
gene-pair correlations of a cell type are mostly noise for either tool.

**Cell-level eigengenes and kME** (`07_cell_level.py`,
`results/cell_level.tsv`). After `ConstructNetwork`, hdWGCNA computes the
module eigengenes of the single cells (`ModuleEigengenes`) and the kME of
every gene on them (`ModuleConnectivity`), which also moves genes with a
negative kME for their module. `gcen2 celltype-modules` does the same by
default (`cell_eigengenes.tsv`, `cell_kme.tsv`; `--no-cell-eigengenes`,
`--no-reassign`). `02_hdwgcna.R` now runs both steps, and the comparison
gives GCEN2's `gcen2.single_cell.cell_eigengenes` hdWGCNA's modules and the
same cells:

| cell type | cells | modules | eigengenes, max \|Δ\| (max \|ME\|) | kME, max \|Δ\| | genes reassigned | final modules |
|---|---|---|---|---|---|---|
| CD14 Mono | 3,333 | 4 | 4.2e-13 (47) | 2.2e-14 | 14 | identical |
| CD4 T | 3,391 | 2 | 3.9e-13 (68) | 3.9e-14 | 82 | identical |
| B | 1,496 | 5 | 1.2e-13 (44) | 2.6e-14 | 44 | identical |
| CD8 T | 1,061 | 3 | 2.6e-13 (57) | 3.9e-14 | 110 | identical |
| NK | 616 | 13 | 1.2e-13 (36) | 9.7e-15 | 114 | identical |
| FCGR3A Mono | 355 | 29 | 3.7e-13 (26) | 2.2e-14 | 194 | identical |

The eigengenes use the code of the spatial workflow (`spot_eigengenes`,
already checked against hdWGCNA on Visium). The kME is accumulated over
blocks of cells, so a cell type of 80,000 cells needs no dense or float64
copy of its cells × genes matrix.

## Cells as samples: metacells, pseudo-replicates, single cells, MAGIC

`gcen2 celltype` makes metacells the samples of a cell type's network. It
also offers random pseudo-replicates (`--strategy pseudo_replicates`,
deprecated), and `build_dropout_aware(zero_handling="impute")` correlates
the cells after MAGIC imputation (warns). `09_strategies.py` builds the
network of each of the six cell types in six ways with GCEN2's own code
(Pearson, r >= 0.3, p <= 0.001, positive edges; the metacell settings above):

- **metacell**: the default: p-values on every pair's effective sample size
  from the metacells' overlap (`metacell_effective_samples`).
- **metacell, cells / k**: the same networks tested on cells / 25 samples
  (`--effective-samples cells`, the rule before 2026-11-18).
- **pseudo10** and **pseudo50**: the cells split at random into 10 (the
  default) or 50 groups. Library-size-normalised counts are summed and
  log-transformed. Summing raw counts makes the sums follow the groups'
  library sizes, and ribosomal protein pairs then rank first.
- **cells**: every log-normalised cell is a sample.
- **magic**: MAGIC 3.0.0, with its defaults, on the same cells.

False edges are measured with 200 null genes added to each cell type. Each
null gene takes a gene of the list and permutes its fraction of every cell's
counts across the cells. Its counts are then drawn as Poisson with mean
fraction × library size. A null gene therefore co-varies with nothing, and
its expected count is proportional to the library size.

A first version permuted the counts within 50 bins of library size. That
keeps E[count | library size], which for genes such as RPS2 is not
proportional to the library size: their share of the counts differs between
cell states with different depths. After normalisation those null genes
followed the states. In CD14 Mono this gave 0.19% of null pairs as metacell
edges, nearly all from three ribosomal genes (RPS2, RPS8, SLC25A6).

The table gives the fraction of the pairs with a null gene that become edges;
with valid p-values it is at most 0.0005. The ranking measure is the STRING
precision of the top 1,000 and top 10,000 pairs (as above). Reproducibility
is the overlap of the top 10,000 pairs of the two halves of the cells. Range
over the six cell types (`results/strategies.tsv`):

| | metacell | metacell, cells / k | pseudo10 | pseudo50 | cells | magic |
|---|---|---|---|---|---|---|
| edges | 23,282–304,059 | 4,191–323,932 | 15,352–24,312 | 21,413–29,523 | 1,488–4,597 | 9.3–17.2 million |
| null pairs that are edges | 1.2e-5–3.3e-4 | 0–1.9e-5 | 5.1–6.2e-4 | 5.7–6.8e-4 | 0–1.5e-5 | 0.12–0.27 |
| STRING precision, top 1,000 | 0.26–0.78 | (same) | 0.001–0.005 | 0.02–0.76 | 0.48–0.91 | 0.05–0.64 |
| STRING precision, top 10,000 | 0.13–0.30 | (same) | 0.002–0.003 | 0.007–0.18 | 0.20–0.36 | 0.04–0.19 |
| top 10,000 shared by the halves | 0.06–0.69 | (same) | 0.000–0.001 | 0.01–0.24 | 0.09–0.79 | 0.17–0.49 |

(All-pair STRING baseline 0.15–0.18%.)

- **Effective samples of metacells.** A gene unrelated to gene y and to
  the choice of metacells has metacell values whose covariance is the
  overlap, `C[a, b] = shared cells / k`. Its correlation with y therefore
  has variance `s(y) / M`, with `s(y) = y'Cy / y'y`. Across the null genes
  of CD14 Mono, the observed variance of the correlations matches
  `M / s(y)` closely; the other five cell types agree as well, from
  1,000 metacells down to FCGR3A Mono's 44 (`10_effective_samples.py`,
  `results/effective_samples.tsv`):

  | genes, by s(y) | effective samples, observed | `M / s(y)` |
  |---|---|---|
  | lowest half | 351 | 372 |
  | 50th–90th percentile | 260 | 274 |
  | 90th–99th percentile | 160 | 158 |
  | top 1% | 110 | 102 |

  The fixed cells / k (133) was thus too strict for 90% of the genes and
  too lenient for the smoothest 1% (cell-state genes). GCEN2 now gives every
  gene `M / s` and tests a pair on the smaller of its two values. The
  null-pair edge rate is 1.2e-5 to 3.3e-4, below the 5e-4 of valid
  p-values. With cells / k it was at most 1.9e-5. The median gene has 2.1
  to 2.9 times as many effective samples as cells / k, and the networks of
  the small cell types have up to 5.6 times as many edges (FCGR3A Mono).
  `tests/test_metacell.py` checks the calibration on simulated cells.
- **Differences between cell types** (`--diff`, Fisher's z). The test now
  uses the same per-gene effective samples: every shared edge on the smaller
  of its two genes' values in each cell type. Calibration: the two random
  halves of each cell type compared as two cell types, 2,000 random genes,
  all 2 million pairs (`11_diff_calibration.py`,
  `results/diff_calibration.tsv`). The halves differ only by chance, so 5%
  of the pairs should have p <= 0.05:

  | cell type (metacells per half) | per gene pair | cells / k |
  |---|---|---|
  | CD14 Mono (652, 651) | 4.9% | 0.66% |
  | CD4 T (732, 679) | 4.6% | 0.76% |
  | B (141, 137) | 4.2% | 0.36% |
  | CD8 T (108, 106) | 4.8% | 0.44% |
  | NK (37, 39) | 4.1% | 0.26% |
  | FCGR3A Mono (17, 14) | 4.8% | 0.36% |

  Per gene pair the rate is 3.7–5.5% in each fifth of the pairs by
  effective samples, from the smoothest genes to the noisiest; at
  p <= 0.01 it is 0.7–0.9% (nominal 1%). With cells / k it falls from
  0.95–2.9% for the smoothest fifth to 0.02–0.17% for the noisiest: the
  test missed most changes.
- **Pseudo-replicates.** Random groups of cells all have the same expected
  composition. They differ only by which cells they draw, so the
  correlations of their sums estimate the cell-level correlations from 10
  (or 50) samples. With 10, the ranking is at chance (STRING 1–3× the
  baseline) and the two halves share no top pairs. The 15,000–24,000 edges
  are about the number expected by chance among 30–42 million pairs at
  p <= 0.001, so nearly every edge is a false one.
  With 50 the top pairs carry some signal, but they are not reproducible.
- **MAGIC.** Diffusion makes every gene follow the cell-state graph, so
  unrelated genes become correlated. Among the pairs with a null gene,
  12–27% become edges. The networks have 9–17 million edges, a quarter to
  a half of all gene pairs, and every top-pair ranking is below that of the
  cells it imputes. The p-values of imputed cells, which count every cell as
  an independent sample, are not valid.
- **Single cells against metacells.** Plain Pearson on the single cells
  ranks the strongest pairs best (top 1,000: 0.48–0.91, against 0.26–0.78
  for metacells; metacells are ahead only for NK) and gives few false
  edges. At r >= 0.3, however, it keeps only 1,500–4,600 edges. Metacells
  raise the correlations of genes that follow cell states. This gives them
  6–95 times as many edges and, for every cell type, a more reproducible
  correlation structure: the correlation of all gene-pair correlations
  between the halves is 0.04–0.41, against 0.02–0.20 for the cells. These
  edges and that structure are what the modules are built from.

The network takes 0.1–1.4 s per cell type for every strategy except MAGIC
(5.1–12 s; 8 threads, `t_network`).
