# GCEN2 vs hdWGCNA: module projection and preservation

Do the modules found in one data set exist in another? hdWGCNA answers with
two functions, and GCEN2 has a command for each:

| hdWGCNA | GCEN2 | what it gives |
|---|---|---|
| `ProjectModules` | `gcen2 project-modules` (`gcen2.single_cell.project_modules`) | the eigengene of every reference module on each query cell |
| `ModulePreservation` (WGCNA's `modulePreservation`) | `gcen2 module-preservation` (`gcen2.preservation.module_preservation`) | Z statistics of how well each module holds together in the query |

`ModulePreservation` compares each module's statistics in the query with the
same statistics for random gene sets of the same size (permutations). Its
Zsummary combines density (do the genes stay correlated?) and connectivity
(do the same genes stay central?). WGCNA reads Zsummary > 10 as strong
evidence of preservation and 2–10 as moderate.

**Data.** Human Lung Cell Atlas core, natural killer cells (see
[`../hdwgcna_consensus`](../hdwgcna_consensus/README.md)):
- reference: dataset Meyer_2019, 6,732 cells;
- query: dataset Krasnow_2020, 4,303 cells.

**Settings.** The same for both tools:
- reference: hdWGCNA's single-cell workflow (genes in at least 5% of the
  cells, metacells of 25 neighbours in the scANVI embedding,
  `ConstructNetwork` with its defaults, `ModuleEigengenes`,
  `ModuleConnectivity`);
- query: metacells built the same way, restricted to the reference genes;
- preservation: hdWGCNA's defaults, 500 permutations, seed 12345, unsigned
  network.

## Run

```bash
cd <hlca>/consensus     # after ../hdwgcna_consensus/01_prepare.py
python -I .../01_prepare.py nk.h5ad ref_meyer.h5ad Meyer_2019
bash .../run_all.sh   # each tool twice, then 03_compare.py
```

GCEN2's workflow:

```bash
gcen2 celltype-modules -i ref_meyer.h5ad --groupby cell_type --use-rep X_scanvi_emb \
  --n-pcs 30 -k 25 --max-shared 10 --write-metacells -o ref/ -t 8
gcen2 project-modules -i nk.h5ad -m ref/natural_killer_cell/modules.tsv \
  --groupby dataset --group Krasnow_2020 -o projected/
gcen2 module-preservation -r ref/natural_killer_cell/metacells.tsv \
  -m ref/natural_killer_cell/modules.tsv -q nk.h5ad --groupby dataset --group Krasnow_2020 \
  --use-rep X_scanvi_emb --n-pcs 30 -k 25 --max-shared 10 -o preservation/ -t 8
```

## Results

| | GCEN2 | hdWGCNA |
|---|---|---|
| reference genes | 4,898 | 4,898 (the same) |
| reference modules | 3 | 4 |
| whole workflow: time, peak memory (mean of 2) | **325 s, 2.54 GB** | **2,174 s, 3.87 GB** |
| projection | 1.6 s | 12.4 s |
| preservation, 500 permutations (on hdWGCNA's metacells) | 296 s | 2,077 s |

Source files: `results/compare.tsv`, `results/time.txt`,
`results/hdwgcna_timing.tsv`.

- **Projection.** GCEN2 is given hdWGCNA's reference modules, grey included
  as hdWGCNA projects it. All 5 eigengenes on the 4,303 query cells are
  within 1.6e-13 of hdWGCNA's.
- **Preservation.** GCEN2 is given hdWGCNA's own reference and query
  metacells and its modules. It produces the same modules (the 4, plus
  WGCNA's random "gold" set and grey) and all 20 observed and Z statistics:
  - observed statistics within 1.3e-8 relative, Z within 7.3e-8;
  - every Zsummary.pres the same to two decimals: blue 18.96, brown 48.43,
    turquoise 17.73, yellow 3.26 (grey 2.79, gold 12.66).

  How this is exact: WGCNA draws a random baseline module ("gold"),
  subsamples modules of more than 1,000 genes and permutes genes, all with
  R's `sample()`. GCEN2 makes the same draws with a port of R's random number
  generator. The remaining relative difference of about 1e-8 comes from
  WGCNA's correlation routine rounding differently in the last bit. The MAR
  statistic subtracts nearly equal sums and magnifies that.
- **Faster:** 7× on the preservation step. GCEN2 runs the permutations in
  parallel; their random draws are made first, in R's order, so the result
  does not depend on the number of threads.
- **Whole workflow** (each tool with its own metacells and modules,
  `results/zsummary.tsv`): GCEN2 finds 3 modules. Each is matched to the
  hdWGCNA module sharing the most genes:

  | GCEN2 module (genes) | hdWGCNA module | Jaccard | Zsummary.pres GCEN2 | Zsummary.pres hdWGCNA |
  |---|---|---|---|---|
  | ME0 (1,343) | turquoise | 0.79 | 20.6 | 17.7 |
  | ME1 (561) | blue | 0.88 | 19.1 | 19.0 |
  | ME2 (295) | brown | 0.25 | 21.7 | 48.4 |

  Both tools find all three modules strongly preserved in the query
  (Zsummary > 10). ME2 (ribosomal proteins) overlaps hdWGCNA's 75-gene
  brown module only partly, so its statistics differ. hdWGCNA's fourth
  module (yellow, 54 genes) is not preserved (Zsummary 3.3); GCEN2 has no
  counterpart for it.

The reference modules (`results/gcen2_ref_hub_genes.tsv`):
- ME0: GZMA, CYBA, RAC2, ACTB (cytotoxicity, cytoskeleton);
- ME1: NR4A2, ZNF331, MAFF, DUSP2 (immediate-early activation);
- ME2: ribosomal proteins.

## Files

| file | contents |
|---|---|
| `01_prepare.py` | the reference dataset's cells as their own `.h5ad` |
| `02_hdwgcna.R` | hdWGCNA: reference network, `ProjectModules`, query metacells, `ModulePreservation`; writes the modules, both metacell matrices, projected eigengenes, preservation tables and stage timings |
| `run_all.sh` | each tool twice, then `03_compare.py` |
| `03_compare.py` | the comparison (`results/compare.tsv`, `results/zsummary.tsv`) |
| `results/` | `compare.tsv`, `zsummary.tsv`, `time.txt`, hdWGCNA's stage timings and Z table, GCEN2's Z table, projected modules and reference hub genes |
