# PyWGCNA vs WGCNA on GTEx whole blood

PyWGCNA (Rezaie et al., Bioinformatics 2023) reimplements the WGCNA
workflow in Python. Does it give WGCNA's modules? Its `findModules`
follows the WGCNA tutorial's step-by-step pipeline: `pickSoftThreshold`,
`adjacency`, `TOMsimilarity`, average-linkage clustering, `cutreeDynamic`
(hybrid, with the PAM stage, `pamRespectsDendro = FALSE`) and
`mergeCloseModules`. Its reference is therefore that same pipeline in WGCNA 1.74
(`02_wgcna_stepwise.R`). It is not `blockwiseModules`, which reassigns genes by
kME, respects the dendrogram in the PAM stage, and is what `gcen2 module -e`
reproduces (`../wgcna_gtex`).

Data and settings are those of `../wgcna_gtex`:
- the same six GTEx v10 whole-blood matrices: 5,000 / 10,000 / 11,944 genes,
  5 PCs regressed out or uncorrected, and the two split halves;
- a signed network and signed TOM, the power vector 1–10, 12–20,
  minModuleSize 30, merge cut 0.25 and deepSplit 2;
- 8 threads (BLAS threads for PyWGCNA), i7-13700, idle machine, mean of
  two runs.

Two of PyWGCNA's defaults differ from WGCNA's, and are set to WGCNA's here:
- `RsquaredCut` is 0.9 where WGCNA uses 0.85. With 0.9, PyWGCNA picks
  power 12 instead of 10 on the corrected 5,000 genes, from the same fit table.
- `TPMcutoff` (1, applied to the input) is switched off: the matrices are
  already filtered and log-transformed, and the PC-regressed values are
  centred near 0.

## Run

PyWGCNA 2.2.1 does not run with the library versions of the other
benchmarks:
- with pandas 3, `TOMsimilarity` writes into the read-only array that
  `DataFrame.to_numpy()` now returns;
- with scipy 1.18, `rankdata` returns floats, which `cutreeHybrid` uses as
  indices.

It runs in its own environment with the versions of its release, February
2025 (`requirements-pywgcna.txt`, Python 3.12).

```bash
cd <gtex>                                    # after ../wgcna_gtex (runs/wgcna/)
PYTHON=<pywgcna-env>/bin/python RSCRIPT=<Rscript with WGCNA> bash .../run_all.sh
GCEN2_STRING_GOLD=<pbmc>/work/string_gold.tsv python -I .../03_compare.py
<Rscript> .../04_treecut.R && <pywgcna-env>/bin/python -I .../04_treecut.py
```

## Results

**Up to the gene tree PyWGCNA agrees with WGCNA:**
- **Soft power:** the scale-free fit table is the same to the printed digits.
  With WGCNA's `RsquaredCut`, the power is the same in all six runs
  (10, 12, 9, 1, 10, 12).
- **TOM:** within 5.0e-9 (relative) of `blockwiseModules`' TOM over WGCNA's top
  300,000 pairs, with the same top 10,000 pairs (`results/tom.tsv`).
- **Gene tree:** the merge heights are those of WGCNA's `hclust` to 5.0e-9
  (`results/treecut.tsv`). This difference is the rounding of 1 − TOM to 8
  decimals that `findModules` applies before clustering.

**The modules are not WGCNA's** (`results/compare.tsv`). The ARI is over all
genes, with unassigned genes as one group:

| run | modules: PyWGCNA / WGCNA step by step | genes in modules | ARI with WGCNA step by step | ARI with blockwiseModules | smallest PyWGCNA module |
|---|---|---|---|---|---|
| corrected, 5,000 genes | 18 / 20 | 5,000 / 4,992 | 0.36 | 0.14 | 5 |
| split half A | 17 / 19 | 5,000 / 4,831 | 0.33 | 0.08 | 70 |
| split half B | 17 / 20 | 5,000 / 5,000 | 0.30 | 0.20 | 16 |
| uncorrected, 5,000 genes | 2 / 2 | 5,000 / 5,000 | 0.00 | 0.00 | 2 |
| corrected, 10,000 genes | 22 / 28 | 10,000 / 10,000 | 0.41 | 0.14 | 20 |
| corrected, 11,944 genes | 26 / 27 | 11,944 / 11,934 | 0.44 | 0.06 | 6 |

For comparison, `gcen2 module -e` gives `blockwiseModules`' modules in all
six runs (ARI 1.0, `../wgcna_gtex`). The two WGCNA pipelines also differ
from each other (ARI 0.07–0.85 between step by step and blockwiseModules), so
the step-by-step pipeline, not `blockwiseModules`, is the fair reference for
PyWGCNA.

**Where the modules diverge.** `04_treecut.py` runs PyWGCNA's
`cutreeHybrid` and WGCNA's `cutreeDynamic` on the same tree (corrected
5,000 genes, power 10):

- **The tree cut itself** (`pamStage = FALSE`):
  - PyWGCNA leaves 3,628 genes unassigned; WGCNA leaves 1,737.
  - On the 1,265 genes both assign, the clusters agree (ARI 0.97), so the
    branches found are largely the same, but PyWGCNA keeps fewer genes on each.
  - Over all genes the ARI is 0.19.
- **With the PAM stage** (as `findModules`): almost every gene is assigned
  (5 / 8 unassigned), from these different cores. The ARI is 0.34.
- **Renumbering:** `cutreeHybrid` renumbers its labels by module size
  (the "Value" column, which `findModules` uses). The unassigned genes
  (branch label 0) receive the label of a module (22 here), and the largest
  module receives 0.
  - `labels2colors` then gives every label a colour, so PyWGCNA never
    reports unassigned genes: "genes in modules" is always all genes.
  - The unassigned genes form a module of their own, below the minimum size
    of 30: the smallest modules in the table (5, 16, 2, 20 and 6 genes).

Within-module STRING enrichment (over the baseline of all pairs), corrected
5,000 genes:
- PyWGCNA 2.9×;
- WGCNA step by step 3.6×;
- blockwiseModules 7.6× (it keeps 3,062 genes in modules; the
  pipelines with `pamRespectsDendro = FALSE` assign almost all of them).

**Time and peak memory** (whole pipeline from the expression file):

| run | PyWGCNA | WGCNA step by step | WGCNA blockwiseModules | GCEN2 `module -e` |
|---|---|---|---|---|
| corrected, 5,000 genes | 27 s, 1.56 GB | 9.7 s, 1.24 GB | 17 s, 1.32 GB | 1.7 s, 0.40 GB |
| split half A | 22 s, 1.51 GB | 7.6 s, 1.14 GB | 11 s, 0.94 GB | 1.4 s, 0.34 GB |
| split half B | 19 s, 1.52 GB | 7.7 s, 1.15 GB | 12 s, 0.94 GB | 1.6 s, 0.34 GB |
| uncorrected, 5,000 genes | 22 s, 1.56 GB | 9.9 s, 1.17 GB | 20 s, 1.12 GB | 1.8 s, 0.46 GB |
| corrected, 10,000 genes | 121 s, 4.93 GB | 26 s, 3.08 GB | 29 s, 3.10 GB | 6.5 s, 1.08 GB |
| corrected, 11,944 genes | 197 s, 6.87 GB | 36 s, 4.04 GB | 36 s, 3.93 GB | 9.1 s, 1.45 GB |

- PyWGCNA's time includes drawing two figures (soft-threshold fit and
  eigengene tree), which `findModules` always draws.
- Most of the time is the tree cut and the module steps after the TOM: 190 s
  of `findModules` on 11,944 genes.
- On all 11,944 genes PyWGCNA takes 22× GCEN2's time and 4.7× its memory.

## Files

| file | contents |
|---|---|
| `01_pywgcna.py` | PyWGCNA on one matrix: `modules.tsv`, `info.tsv`, top TOM pairs |
| `02_wgcna_stepwise.R` | WGCNA's step-by-step pipeline with PyWGCNA's settings |
| `run_all.sh` | both on the six matrices, twice, with `/usr/bin/time` |
| `03_compare.py` | `results/compare.tsv` (power, modules, ARI, STRING, time), `results/tom.tsv` |
| `04_treecut.R`, `04_treecut.py` | the tree cut stage by stage: `results/treecut.tsv` |
| `requirements-pywgcna.txt` | PyWGCNA's environment |
