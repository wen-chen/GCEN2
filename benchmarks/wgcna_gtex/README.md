# GCEN2 vs WGCNA on GTEx whole blood (bulk RNA-seq)

GTEx v10 whole blood, 803 samples: log2(TPM + 1) of the genes with median
TPM >= 1 (11,944 unique symbols); the 5,000 most variable for the main
comparison, 10,000 and all 11,944 for scaling. Both tools read the same
matrices.

* **GCEN2**: `gcen2 build -m pearson -c 0.4 -p 0.001 -f -s` (positive edges,
  FDR) and `gcen2 module` (similarity 0.5); `-c 0.8` on uncorrected data.
* **WGCNA 1.74**: `pickSoftThreshold` (signed; 12 when no power reaches the
  fit) and `blockwiseModules` (signed network and TOM, one block,
  minModuleSize 30, mergeCutHeight 0.25).

**Confounders.** The top expression PCs of GTEx carry technical variation
that correlates most genes (PC1-2: 61% of the variance; a third of all gene
pairs at r >= 0.5). As recommended for GTEx co-expression (Parsana et al.,
Genome Biol 2019), the main comparison uses the data with the top 5 PCs
regressed out (each split half corrected on its own). 5 PCs gave the best
STRING precision among 1, 3, 5, 10 and 20, so the absolute precisions on
corrected data are optimistic; both tools get the same input.

STRING v12 gold standard: experimental or database channel >= 700 (no
co-expression evidence). 8 threads, i7-13700.

## Data

GTEx Analysis v10, gene TPM of whole blood, from the GTEx Portal (open
access): `https://storage.googleapis.com/adult-gtex/bulk-gex/v10/rna-seq/tpms-by-tissue/gene_tpm_v10_whole_blood.gct.gz`
(SHA-256 `53212ab8f21417e1be5b5ea36ed0a86d89b8e1ecf1fb011162853812a789aaaa`),
saved as `data/gene_tpm_v10_whole_blood.gct.gz` in a data directory of your
choice (`<gtex>`). The STRING gold standard comes from the PBMC benchmark
(`../hdwgcna_pbmc`, `work/string_gold.tsv`): give its path in
`GCEN2_STRING_GOLD`, or copy it to `<gtex>/work/`.

## Run

Environments: `../environment-r.yml` (R, WGCNA) and
`../requirements-python.txt`; `RSCRIPT` is that environment's `Rscript`
(default: `Rscript` on the `PATH`). All scripts run from the data directory.

```bash
cd <gtex>
export GCEN2_STRING_GOLD=<pbmc>/work/string_gold.tsv
python -I .../01_prepare.py    # -> work/
bash .../run_all.sh            # both tools, 6 runs
python -I .../04_compare.py    # results/
python -I .../05_tom.py        # GCEN2's TOM (gcen2 build --tom) vs WGCNA's
python -I .../06_tom_modules.py  # gcen2 module -e (WGCNA-style modules) vs blockwiseModules
python -I .../07_tom_precision.py  # 1 - TOM in single vs double precision (results/tom_precision*.tsv)
```

## Results

**Correlations.** GCEN2's Pearson r equals `WGCNA::cor` to 2.8e-15 (all
44,850 pairs of 300 genes), so rankings by r are the same for both tools.

**STRING precision**, top 1,000 / 10,000 / 100,000 pairs (baseline 0.07%):

| ranking | corrected (5 PCs) | uncorrected |
|---|---|---|
| Pearson r (both tools) | 34.4% / 11.0% / 1.84% | 9.4% / 2.6% / 0.74% |
| GCEN2 Spearman r | 37.5% / 11.0% / 1.82% | 8.7% / 2.4% / 0.69% |
| WGCNA TOM | 41.0% / 13.1% / 1.93% | 0.0% / 0.04% / 0.07% |

On corrected data WGCNA's topological overlap ranks interacting pairs higher
than the correlation itself (+7 and +2 points at the top 1,000 / 10,000).
On uncorrected data `pickSoftThreshold` estimates soft power 1 and the
TOM falls to the baseline (its top pairs are hubs of the confounded block):
confounding hurts it more than plain correlation.

**GCEN2's TOM** (`gcen2 build --tom`, `05_tom.py`): with `--power auto`
GCEN2 picks the same soft power as `pickSoftThreshold` in all six runs (10,
1, 12, 9, 10, 12). On the 5,000 genes its TOM is within 9.9e-7 (relative) of
`blockwiseModules`' TOM over WGCNA's top 300,000 pairs, the top 10,000 pairs
are the same, and so is the STRING precision (corrected 41.0% / 13.1% /
1.93%, uncorrected 0.0% / 0.04% / 0.07%). All 12.5M pairs take 1.9 s. On
all 11,944 genes `--tom` adds about 4.4 s and 0.6 GB to the FDR build below
(7.0 s, 1.24 GB in all with the matrix-product correlation kernel, 15.4 s
with the per-pair one; `tom` lines of `results/time_gcen2_gemm.txt`; WGCNA
36 s, 3.9 GB).

**GCEN2's WGCNA-style modules** (`gcen2 module -e expr.tsv --signed
--min-size 30 --merge-cut 0.25`, `06_tom_modules.py`: TOM, average-linkage
clustering, dynamic tree cut, eigengene checks and merging as
`blockwiseModules`) give exactly WGCNA's modules in all six runs (adjusted
Rand index 1.0 over all genes, same modules and unassigned genes), with
WGCNA's whole run (pickSoftThreshold + blockwiseModules) for comparison
(`results/tom_modules.tsv`, mean of two runs):

| run | modules (genes) | GCEN2 | WGCNA |
|---|---|---|---|
| corrected, 5,000 genes | 19 (3,062) | 1.7 s, 0.40 GB | 17 s, 1.3 GB |
| split half A | 18 (2,129) | 1.4 s, 0.34 GB | 11 s, 0.94 GB |
| split half B | 20 (3,626) | 1.6 s, 0.34 GB | 12 s, 0.94 GB |
| uncorrected, 5,000 genes | 2 (4,714) | 1.8 s, 0.46 GB | 20 s, 1.1 GB |
| corrected, 10,000 genes | 27 (6,321) | 6.5 s, 1.1 GB | 29 s, 3.1 GB |
| corrected, 11,944 genes | 23 (4,431) | 9.1 s, 1.5 GB | 36 s, 3.9 GB |

The 1 - TOM matrix for modules is kept in double precision
(`07_tom_precision.py`, `results/tom_precision*.tsv`; corrected 5,000 genes).
Rounded to single precision, whose spacing near 1 is 6e-8, the
dissimilarities of weakly connected pairs become ties (every pair off by
1.5e-8, median); the average-linkage tree changes and 59 genes move to other
modules (tree cut ARI 0.97, final 0.98). A TOM computed in single precision
(as `build --tom`, errors up to 2.8e-7) or rounded to it, with 1 - TOM then in
double, gives the same modules: most pairs are off by far less (2e-10,
median). Memory is one genes x genes double matrix (correlations,
adjacency and 1 - TOM are built in place, and the clustering works in it)
plus the expression data; most of the time is the TOM (6.5 of 10.8 s for
11,944 genes when the AVX product kernel was added in commit 47b941f; 13 s
with SSE2 only).

**Modules** (5,000 genes, corrected):

| | GCEN2 | WGCNA |
|---|---|---|
| modules | 123 | 19 |
| genes in modules | 1,495 | 3,062 |
| largest module | 286 | 437 |
| within-module STRING enrichment | 30x | 7.6x |
| split-half ARI (401 vs 402 samples) | 0.74 | 0.79 |

Module ARI between the tools: 0.69 (1,253 genes in modules of both).
Uncorrected, both are dominated by one confounded block (GCEN2 largest
module 2,787 genes, WGCNA 2 modules covering 4,714 genes; ARI 0.06).

**Time and peak memory** (`gcen2 build` + `gcen2 module`, mean of two runs,
`results/time_gcen2_gemm.txt`; WGCNA: pickSoftThreshold + blockwiseModules):

| genes | GCEN2, FDR | GCEN2, no FDR | WGCNA |
|---|---|---|---|
| 5,000 | 0.8 s, 0.22 GB | 0.5 s, 0.13 GB | 17 s, 1.3 GB |
| 10,000 | 2.2 s, 0.48 GB | 1.0 s, 0.22 GB | 29 s, 3.1 GB |
| 11,944 | 2.8 s, 0.62 GB | 1.2 s, 0.25 GB | 36 s, 3.9 GB |

With 803 samples `build` uses the matrix-product correlation kernel
(`--kernel auto`, from 64 samples): the same edges and modules as the
per-pair kernel of GCEN v1 (`--kernel pairwise`), whose `build` took 1.8 s /
0.22 GB, 7.4 s / 0.47 GB and 11.1 s / 0.60 GB with FDR (1.4, 6.3 and 9.6 s
without; `pairwise` lines of `results/time_gcen2_gemm.txt`). On all 11,944
genes the two edge files differ in one p-value in the 7th printed digit
(9.501664e-33 vs 9.501665e-33) of 533,232 edges.

Same edges with and without FDR here (with 803 samples every pair at
r >= 0.4 has p far below 0.001). The FDR needs the p-value of every pair
with p <= 0.001 (28M of the 71M pairs on all genes, 16 bytes each while they
are ranked), hence its extra memory. Over the six runs GCEN2 is 13-21x faster
than WGCNA with 3.2-6.5x less memory with FDR (`results/time.tsv`); without
FDR 30-35x faster with 10-15x less memory.

Before the FDR path moved to C++ (commit ca5f8f5: numpy, every block
correlated twice) the FDR runs took 5.4 s / 0.45 GB, 23 s / 1.15 GB and
33 s / 1.46 GB; the edge files are byte-identical.
