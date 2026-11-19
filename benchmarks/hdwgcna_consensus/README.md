# GCEN2 vs hdWGCNA: consensus network, harmonized eigengenes, DMEs and traits

hdWGCNA's workflow for data from several samples or studies, run on the same
cells by both tools:

1. **Consensus network.** Metacells are built within each cell type and
   dataset (`MetacellsByGroups(group.by = c(cell_type, dataset))`). Each
   dataset is one set (`SetMultiExpr`) and gets its own soft power
   (`TestSoftPowersConsensus`). `ConstructNetwork(consensus = TRUE)` then calls
   WGCNA's `blockwiseConsensusModules`: one TOM per set, the TOMs calibrated
   to the first set, their 0.3 quantile as the consensus TOM, then modules.
2. **Harmonized eigengenes.** `ModuleEigengenes(group.by.vars = donor_id)`
   runs Harmony on each module's principal components. `ModuleConnectivity`
   then computes the kME of every cell against the harmonized eigengenes, and
   genes with a negative kME for their module are reassigned.
3. **Differential module eigengenes and module-trait correlations.**
   `FindAllDMEs(group.by = dataset)` tests each module's harmonized
   eigengene between the cells of each dataset and all other cells (Seurat
   `FindMarkers`, Wilcoxon). `ModuleTraitCorrelation(traits = age, BMI, sex;
   group.by = dataset)` correlates the eigengenes with the donors' traits, for
   all cells and within each dataset.

GCEN2 runs both in one command:

```bash
gcen2 celltype-modules -i nk.h5ad --groupby cell_type --consensus-by dataset \
  --harmonize-by donor_id --use-rep X_scanvi_emb --n-pcs 30 -k 25 --max-shared 10 \
  --dme-by dataset --trait age --trait BMI --trait sex --trait-groupby dataset -o out/ -t 8
```

**Data.** Human Lung Cell Atlas core (Sikkema et al. 2023; see
[`../hdwgcna_hlca`](../hdwgcna_hlca/README.md)): the natural killer cells of
the 4 datasets with at least 790 of them, 16,006 cells from 38 donors. Traits:
age (`age_or_mean_of_age_range`), BMI (missing for 2,444 cells) and sex.

| dataset | cells |
|---|---|
| Banovich_Kropski_2020 | 4,179 |
| Krasnow_2020 | 4,303 |
| Lafyatis_Rojas_2019_10Xv2 | 792 |
| Meyer_2019 | 6,732 |

**Settings.** The same for both tools:
- genes with counts in at least 5% of the cells;
- metacells: 25 neighbours in the scANVI embedding, at most 10 shared, up to
  1,000 per dataset;
- hdWGCNA's consensus defaults: signed network, consensus quantile 0.3,
  deepSplit 4, minimum module size 50, merge cut 0.2.

## Run

```bash
cd <hlca>/consensus
python -I .../01_prepare.py ../raw/hlca_core.h5ad nk "natural killer cell" \
  Meyer_2019 Krasnow_2020 Banovich_Kropski_2020 Lafyatis_Rojas_2019_10Xv2
bash .../run_all.sh   # each tool twice, metacell seeds 1 and 2, then 03_compare.py
```

## Results

| | GCEN2 | hdWGCNA |
|---|---|---|
| genes | 5,080 | 5,080 (the same) |
| metacells | 3,361 | 3,359 |
| soft power per set (Banovich, Krasnow, Lafyatis, Meyer) | 8, 5, 6, 4 | 8, 5, 6, 4 |
| modules | 3 | 3 |
| time, peak memory (mean of 2) | **23.7 s, 2.51 GB** | **161 s, 8.13 GB** |

Source files: `results/compare.tsv`, `results/time.txt`,
`results/hdwgcna_timing.tsv`, `results/*_sets.tsv`.

- **Network step.** GCEN2 is given hdWGCNA's own metacells for each set
  (`gcen2.consensus.consensus_tom_modules` with `HDWGCNA_CONSENSUS`). It
  picks the same four powers and gives exactly hdWGCNA's modules (ARI 1.0),
  in 5.8 s against 58 s for `TestSoftPowersConsensus` + `ConstructNetwork`.
  The calibration draws its sample of TOM entries with R's random number
  generator and `sample()`. The consensus quantile uses WGCNA's own
  quickselect routine. Both are reproduced bit for bit.
- **Eigengenes on hdWGCNA's modules.**
  - Unharmonized eigengenes of every cell: within 1.0e-13 of hdWGCNA's.
  - Harmonized eigengenes: Pearson r 0.997–0.998 with hdWGCNA's, per module.
    That equals GCEN2's agreement with itself under another Harmony seed
    (r 0.997). Harmony is random (k-means seeding, the order of cell
    updates), so this is as close as two Harmony runs get. hdWGCNA's result
    does not change with R's seed only because Seurat's `RunPCA` sets the
    seed to 42 just before each Harmony run.
  - kME against the harmonized eigengenes: r 0.99994, at most 0.006 apart.
  - How much harmonizing changes the eigengenes here: r 0.965 between the
    harmonized and unharmonized eigengenes.
- **DMEs and module-trait correlations on hdWGCNA's harmonized
  eigengenes** (`gcen2.single_cell.find_all_dmes` and `module_trait_cor`):
  - DMEs: the same 12 rows in the same order and the same `pct.1` /
    `pct.2`. p-values, adjusted p-values and `avg_log2FC` agree within
    8e-14 (relative).
  - Correlations: all 45 (3 traits × 3 modules, for all cells and each of
    the 4 datasets) within 1.1e-14, with the same numbers of cells.
  - Correlation p-values and FDR agree within 1.2e-7 (relative) where
    hdWGCNA's are above 1e-12. hdWGCNA computes them as `2 * (1 - pt(...))`
    (`Hmisc::rcorr`), which loses digits below about 1e-12 and gives 0
    below about 1e-16 (12 of the 45 here). GCEN2 keeps them (for example
    1.5e-159 for BMI and ME1).
- **Whole workflow** (each tool with its own random metacells): ARI
  0.76–0.85 between GCEN2 and hdWGCNA over the 9 pairs of seeds. This is
  the same range as each tool against itself under another seed (GCEN2
  0.79–0.84, hdWGCNA 0.83–0.86).
- **Faster, less memory:** 6.8× faster with 3.2× less memory. hdWGCNA spends
  33 s on the metacells, 58 s on the network, 44 s on scaling, the
  harmonized eigengenes and kME, and 0.5 s on the DMEs and trait
  correlations.

The three consensus modules (`results/gcen2_hub_genes.tsv`):
- ME0: ACTB, RAC2, PFN1, CLIC1, CORO1A (cytoskeleton, cytotoxic machinery);
- ME1: NR4A2, CD69, DNAJA1, NFKBIA (activation and stress response);
- ME2: ribosomal proteins.

The DMEs (`results/gcen2_dmes.tsv`) and trait correlations
(`results/gcen2_trait_cor.tsv`) mostly show differences between datasets:
- ME1 (activation) is high in Meyer_2019 (log2 fold change 5.7) and low in
  the other three datasets.
- Over all cells, ME1 correlates with BMI (r 0.23) and with sex (r −0.24,
  higher in women). Each dataset has few donors, so these trait correlations
  are partly differences between datasets.

**A pitfall in hdWGCNA.** `TestSoftPowersConsensus` stores the power table,
and `ConstructNetwork` splits it with `dplyr::group_split`, which sorts the
sets alphabetically. It then pairs these powers with the sets of
`SetMultiExpr` by position. `SetMultiExpr` keeps the sets in order of first
appearance in the cells unless `multi_groups` is given. In this sample, the
default order gave each set another set's power: the powers 8, 5, 6, 4 went
to Krasnow, Meyer, Banovich, Lafyatis instead of Banovich, Krasnow,
Lafyatis, Meyer. The modules changed accordingly (ARI 0.84 against the
correct pairing). `02_hdwgcna.R` therefore passes the sets in sorted order.
GCEN2 always takes the sets in sorted order, so every set gets its own
power. With the mismatched powers passed explicitly, GCEN2 again gives
hdWGCNA's modules exactly (ARI 1.0).

## Files

| file | contents |
|---|---|
| `01_prepare.py` | one HLCA cell type over the given datasets → `.h5ad` and the R export |
| `02_hdwgcna.R` | hdWGCNA's consensus workflow with harmonized eigengenes, `FindAllDMEs` and `ModuleTraitCorrelation`; writes its metacell sets, modules before and after `ModuleConnectivity`, eigengenes, harmonized eigengenes, DMEs, trait correlations and stage timings |
| `run_all.sh` | each tool twice, then metacell seeds 1 and 2, then `03_compare.py` |
| `03_compare.py` | the comparison (`results/compare.tsv`) |
| `results/` | `compare.tsv`, `time.txt`, hdWGCNA's stage timings and sets, GCEN2's `celltypes.tsv`, sets, hub genes, DMEs and trait correlations |
