# GCEN2 vs hdWGCNA on spatial transcriptomics

WGCNA-style modules of two Visium samples and, without a Visium grid, one
Slide-seqV2 sample ([below](#data-without-a-visium-grid-slide-seqv2)) and one
Visium HD sample ([below](#visium-hd)). The
Visium modules are built by `gcen2 spatial-modules` and by hdWGCNA 0.4.12's
spatial workflow (Seurat 5.5.1, WGCNA 1.74, R 4.4.3; the
"hdWGCNA in spatial transcriptomics" tutorial), from the same counts and the
same groups:

| sample | spots | groups for the metaspots |
|---|---|---|
| DLPFC 151673 (human dorsolateral prefrontal cortex; Maynard et al. 2021, spatialLIBD) | 3,611 with a layer annotation | the manual cortical layers L1–L6, WM |
| 10x mouse brain, sagittal anterior (V1_Mouse_Brain_Sagittal_Anterior, Space Ranger 1.1.0) | 2,695 | k-means (12, seed 0) of 30 PCs of log-normalized HVGs; 3 groups have < 50 spots and are skipped by both tools |

Both tools follow the same steps:
1. Genes are those with counts in at least `round(0.05 × spots)` spots (`SetupForWGCNA(gene_select = "fraction")`).
2. Metaspots are built within each group (`MetaspotsByGroups`: each centre spot is summed with its six grid neighbours) and then log-normalized (`NormalizeMetacells`).
3. Modules are found on all metaspots: soft power by `TestSoftPowers`, then `ConstructNetwork` with hdWGCNA's defaults.
4. Module eigengenes are computed on the spots (`ModuleEigengenes`).
5. kME is computed on the spots, and genes with a negative kME for their module are reassigned (`ModuleConnectivity`).

hdWGCNA's `MetaspotsByGroups` calls `GetAssayData(slot = )`, which SeuratObject 5 made defunct. `02_hdwgcna.R` patches only that input check. The aggregation itself already uses `LayerData`.

## Run

```bash
cd <spatial>      # a data directory of your choice; raw/ holds the downloads (below)
python -I .../01_prepare.py dlpfc raw/dlpfc_151673 work/dlpfc
python -I .../01_prepare.py tenx raw/mouse_brain_anterior work/mouse_anterior
bash .../run_all.sh        # both tools twice per sample, then 04_compare.py
```

Downloads:
- **DLPFC**: `https://spatial-dlpfc.s3.us-east-2.amazonaws.com/h5/151673_filtered_feature_bc_matrix.h5`. From `github.com/LieberInstitute/HumanPilot`: `10X/151673/{tissue_positions_list.txt, scalefactors_json.json, tissue_lowres_image.png}` and `10X/barcode_level_layer_map.tsv`. SHA-256: matrix `216e8010e1c313925afa0f174207da7c8947497c960e81b07edbd16176460413`; positions `c1ce0eafc855982bb2a9773f147ce6c5aef68dfd06600caa8a7e7a0e6f7cc5e8`; scale factors `91c4beb9ca8145855a6a6df0a00eda6541931fc12d82d09cf8f01c9c15179dbc`; image `2ca692ceb00b205bd2f75182935a018753b08c63fbba7e142fd154db22547fcf`; layer map `b5d81f0e61810bb74d8a02495eda8f84e0eb9335211dcdd03513d747036e5ca9`.
- **Mouse brain**: `https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/V1_Mouse_Brain_Sagittal_Anterior/V1_Mouse_Brain_Sagittal_Anterior_{filtered_feature_bc_matrix.h5,spatial.tar.gz}`. SHA-256: matrix `56078d8d6fe6c13de248fdb1c518b691cdef78fb00021b659786b4a47c6656d5`; spatial `5f41a803e2bd69fa4dfca6abc8fa2d4e0d76aeb6c72d7038a5fdcf9cc50a36f8`.

Machine: Intel Core i7-13700, 123 GB RAM, Ubuntu 26.04.1; 8 threads for both tools. Time and peak memory cover the whole command, reading the data included, with `--plot`. hdWGCNA's time also includes attaching enrichR, which contacts maayanlab.cloud.

## Results

| | DLPFC 151673 | mouse brain anterior |
|---|---|---|
| genes / metaspots | 9,317 / 895 | 12,861 / 604 |
| soft power, modules | 7, 2 | 9, 11 |
| same genes, metaspots (names and order) | yes | yes |
| datExpr, max \|Δ\| | 5e-10 (TSV precision; 6e-15 in memory) | 5e-10 |
| modules after ConstructNetwork, ARI | 1.0 | 1.0 |
| genes reassigned by spot kME | 0 | 1 (Sec62) |
| final modules, ARI | 1.0 | 1.0 |
| spot eigengenes, max \|Δ\| | 5e-9 (TSV precision; 1e-13 in memory) | 5e-9 |
| spot kME, max \|Δ\| | 5e-11 | 5e-11 |
| **GCEN2** time, peak memory (mean of 2) | **9.0 s, 1.28 GB** | **18.9 s, 1.94 GB** |
| **hdWGCNA** time, peak memory (mean of 2) | **113 s, 5.85 GB** | **140 s, 10.6 GB** |
| faster / less memory | 13× / 4.6× | 7.4× / 5.5× |

Both tools ran twice, alternating (`results/time.txt`: GCEN2 8.1 and 9.6 s
on DLPFC, 18.1 and 17.4 s on the mouse brain; hdWGCNA 116 and 111 s, 140
and 141 s). hdWGCNA's own stage timers (`results/*_hdwgcna_timing.tsv`)
give the following for DLPFC and the mouse brain respectively:
- metaspots: 7 s and 10 s;
- soft power plus network: 63 s and 83 s;
- eigengenes plus kME on the spots: 25 s and 29 s.

In GCEN2 most of the time goes to the TOM (C++).

The spot eigengenes recover the tissue structure:
- **DLPFC**: ME1, whose hub genes include CRYAB, is the white matter (r = 0.86 with the WM indicator). ME0, with SNAP25, VSNL1 and NRGN, is the grey matter.
- **Mouse brain**: the modules separate the striatum, the cortex, the olfactory bulb layers and smaller structures (`results/*_spatial_eigengenes.png`).

## Data without a Visium grid: Slide-seqV2

hdWGCNA's `MetaspotsByGroups` needs the Visium grid. For other platforms
(Slide-seq beads, Visium HD bins, Stereo-seq, single cells of imaging
platforms) `gcen2 spatial-modules` builds metaspots from the coordinates
alone (`--metaspots knn`, chosen automatically when `obs` has no hexagonal
grid): each centre with its k − 1 nearest spots of its group, centres in a
random order (`--seed`) kept when they share at most `--max-shared` spots
(default round(2k/7)) with every metaspot kept before. hdWGCNA users do the
same with `MetacellsByGroups` on the spatial coordinates, so that is the
comparison here (`02_hdwgcna_knn.R`).

**Sample**: Slide-seqV2 mouse hippocampus, puck 200115_08 (Stickels et al.
2021): 41,786 beads of 10 µm. squidpy's `slideseqv2` file has log-normalized
values of 4,000 highly variable genes; `01_prepare_slideseq.py` recovers
their integer counts exactly, a median of 87 per bead. Both tools use the
same settings:
- all beads in one group;
- k = 100 beads per metaspot (about the tissue area a Visium metaspot of 7
  spots captures), at most 29 shared, summed;
- genes with counts in at least 1% of the beads (1,770 genes; 5% keeps 604).

```bash
cd <spatial>
python -I .../01_prepare_slideseq.py raw/slideseqv2/slideseqv2.h5ad work/slideseqv2
bash .../run_slideseq.sh   # both tools twice, then metacell seeds 1 and 2; 05_compare_slideseq.py
```

Download: `https://ndownloader.figshare.com/files/28242783` (squidpy's
`slideseqv2`; SHA-256 `224a0f2b3d4f6c0ac1583c93bb6bfa910a986df82302c4a927d61883b8b63d8b`).

| | GCEN2 | hdWGCNA |
|---|---|---|
| genes | 1,770 | 1,770 (the same) |
| metaspots / metacells | 602 | 598 |
| soft power, modules | 7, 6 | 6, 5 |
| time, peak memory (mean of 2) | **7.2 s, 0.65 GB** | **220 s, 2.97 GB** |

- **Network step**: `gcen2 module -e --preset hdwgcna` on hdWGCNA's own
  metacell matrix gives exactly its modules (ARI 1.0), in 0.72 s against
  14 s for `TestSoftPowers` + `ConstructNetwork`.
- **Whole workflow** (each tool with its own random metacells): GCEN2's
  modules agree with hdWGCNA's (ARI 0.66–0.85 over the 9 pairs of seeds) as
  well as each tool agrees with itself under another seed (GCEN2 0.64–0.87,
  hdWGCNA 0.68–0.72). Five of GCEN2's six bead eigengenes have a counterpart
  in hdWGCNA's (|r| 0.97–1.00). The sixth (ME5: Gfap, Slc1a3, Gja1,
  astrocytes) is a module hdWGCNA did not find with its default seed
  (best |r| 0.38).
- **Faster / less memory**: 31× / 4.5×. hdWGCNA spends 156 s building
  the metacells, 14 s on the network and 37 s on eigengenes and kME
  (`results/slideseqv2_hdwgcna_timing.tsv`).

The modules map the structures of the section
(`results/slideseqv2_spatial_eigengenes.png`):
- ME0 (Nrgn, Neurod2): the CA and dentate gyrus cell layers;
- ME1 (Shox2, Plcb4, Lef1): the thalamus;
- ME2 (Plp1, Mbp, Mobp): white matter;
- ME3 (Syt1, Snap25): neuronal layers;
- ME4 (Pou4f1, Gpr151): the medial habenula;
- ME5 (Gfap, Gja1): astrocytes.

### knn metaspots on the Visium samples

How much does the metaspot method change the result? `06_knn_vs_grid.py`
runs `--metaspots knn` (k = 7, from the pixel coordinates; seeds 42, 1, 2)
on the two Visium samples and compares the result with the grid run above
(`results/knn_vs_grid.tsv`; ARI over all genes):

| | DLPFC 151673 | mouse brain anterior |
|---|---|---|
| metaspots, grid / knn | 895 / 699–705 | 604 / 476–487 |
| ARI, knn vs grid | 0.75–0.76 | 0.54–0.59 |
| ARI, knn seed vs seed | 0.78 | 0.62–0.74 |
| ARI, grid on 95% of the spots vs grid on all | 0.87–0.88 | 0.58–0.64 |
| spot eigengenes: lowest best \|r\| with a knn eigengene | 1.00 | 0.94 |

Every grid module has a knn module with the same spatial pattern (|r| ≥
0.94). On the gene level, knn agrees with the grid about as well as knn
agrees with itself under another seed, and a little less than the grid
method agrees with itself when 5% of the spots are left out. The gene
assignment of the 11 mouse-brain modules is the least stable in all three
comparisons. The grid method remains the default for Visium.

## Visium HD

Visium HD bins lie on a square grid, so `spatial-modules --metaspots auto`
takes knn metaspots for them, and on the hdWGCNA side this is again
`MetacellsByGroups` on the spatial coordinates (`02_hdwgcna_knn.R`).

**Sample**: 10x Visium HD mouse brain (FFPE, coronal; Space Ranger 3.0.0),
the 16 µm bins: 98,917 bins under the tissue, a median of 667 counts per bin.
Groups: k-means (12, seed 0) of 30 PCs of log-normalized HVGs, as for the
Visium mouse brain (448–31,019 bins per group). Both tools use the same
settings:
- k = 64 bins per metaspot (16,384 µm², about the area of the 7 spots of a
  Visium metaspot), at most 18 shared, summed;
- genes with counts in at least 5% of the bins (4,316 genes).

```bash
cd <spatial>
tar xzf raw/visium_hd_mouse_brain/binned_outputs.tar.gz -C raw/visium_hd_mouse_brain \
  binned_outputs/square_016um/filtered_feature_bc_matrix.h5 \
  binned_outputs/square_016um/spatial/tissue_positions.parquet
python -I .../01_prepare_visium_hd.py raw/visium_hd_mouse_brain/binned_outputs/square_016um work/visium_hd
bash .../run_visium_hd.sh   # both tools twice, then metacell seeds 1 and 2; 05_compare_slideseq.py
```

Download: `https://cf.10xgenomics.com/samples/spatial-exp/3.0.0/Visium_HD_Mouse_Brain/Visium_HD_Mouse_Brain_binned_outputs.tar.gz`
(4.6 GB; SHA-256 `0f3a8a5cc796230413dd3c66f886936a9b8152f6407b74d9bd9a2c9ce41dc076`). The parquet positions need pyarrow; without it in the GCEN2
environment, convert them once to `spatial/tissue_positions.csv` with any
Python that has it, and `01_prepare_visium_hd.py` reads the CSV.

| | GCEN2 | hdWGCNA |
|---|---|---|
| genes | 4,316 | 4,316 (the same) |
| metaspots / metacells | 1,978 | 1,986 |
| soft power, modules | 6, 6 | 6, 6 |
| time, peak memory (mean of 2) | **19.2 s, 2.23 GB** | **570 s, 18.8 GB** |

(`results/visium_hd_compare.tsv`, `results/time_visium_hd.txt`)

- **Network step**: `gcen2 module -e --preset hdwgcna` on hdWGCNA's own
  metacell matrix gives exactly its modules (ARI 1.0), in 2.9 s against
  24 s for `TestSoftPowers` + `ConstructNetwork`.
- **Whole workflow** (each tool with its own random metacells): ARI
  0.915–0.945 between the tools over the 9 pairs of seeds, the same as each
  tool against itself under another seed (GCEN2 0.925–0.944, hdWGCNA
  0.917–0.943). Every GCEN2 bin eigengene has a counterpart in hdWGCNA's
  (|r| 0.98–1.00).
- **Faster / less memory**: 30× / 8.4×. hdWGCNA spends 235 s building the
  metacells, 24 s on the network and 286 s on the bin eigengenes and kME
  (`results/visium_hd_hdwgcna_timing.tsv`). GCEN2 computes the eigengenes
  and kME of the bins in row blocks of the sparse counts, without a dense
  bins × genes matrix.

The modules (`results/visium_hd_spatial_eigengenes.png`,
`results/visium_hd_hub_genes.tsv`):
- ME0 (Stxbp1, Syp, Snap25): neuronal grey matter;
- ME1 (Cldn11, Mag, Cnp): white matter (oligodendrocytes);
- ME2 (Nrgn, Phyhip, Chn1): cortex and the hippocampal CA layers;
- ME3 (Atp2a2, Kcnc2, Pdp1): the thalamus;
- ME4 (Gprasp1, Maged1, Ndn): hypothalamus and ventral areas;
- ME5 (Apoe, Aldoc, Atp1a2, Mt1): astrocytes, strongest at the pial surface.

## `gcen2 spatial-build` on DLPFC

`spatial-build` does not reproduce a published tool. It builds correlation
networks of the spots themselves, without metaspots or modules. The
command offers three modes:

- **filtered:** keeps the edges between genes with significant positive
  spatial autocorrelation: Moran's I > 0 and FDR <= 0.05 (`--moran-fdr`)
  in a one-sided test with the moments of I under randomisation (Cliff and
  Ord). Before this test was added, the rule was I > 0 alone
  (`--moran-fdr 1`, row "I > 0").
- **weighted:** smooths every spot with its neighbours before correlating.
- **regions:** builds one network per region.

`07_spatial_build.py` runs the three modes and the plain network on
DLPFC 151673:

- 3,611 spots and 9,317 genes;
- Spearman, p <= 0.05 and 6 neighbours, the command's defaults;
- |r| >= 0.6, the default cutoff, and also 0.3;
- regions: the seven manually annotated layers.

As in the PBMC benchmark, 200 null genes are added. Each takes a real gene,
permutes its fractions of the spots' counts across the spots, and draws
Poisson counts of fraction × library size. Of these, 150 pass the 5%
detection filter. The edges are
compared with STRING (the PBMC gold standard; all-pair baseline 0.13%).
Results are in `results/spatial_build.tsv` and
`results/spatial_build_moran.tsv`:

| mode | edges (genes), \|r\| >= 0.6 | STRING precision | edges (genes), \|r\| >= 0.3 | STRING precision | null edges at 0.3 |
|---|---|---|---|---|---|
| plain | 39 (14) | 0.74 | 2,930 (208) | 0.020 | 0 |
| filtered (FDR <= 0.05) | 39 (14) | 0.74 | 2,929 (206) | 0.020 | 0 |
| I > 0 | 39 (14) | 0.74 | 2,930 (208) | 0.020 | 0 |
| weighted | 593 (145) | 0.088 | 285,890 (2,569) | 0.007 | 0 |
| regions, pooled | 65 (64) | 0.45 | 26,997 (7,566) | 0.003 | 726 |

- **Plain network.** At the default cutoff the network is small and
  precise. Its 14 genes are three white-matter genes (MBP, PLP1, TF), nine
  mitochondrial genes and two secretoglobins (SCGB1D2, SCGB2A2), among the
  genes with the highest Moran's I.
- **filtered.** It removes almost nothing here: 1 edge at 0.3. With
  I > 0 alone it removes nothing. 82% of the genes have I > 0, including
  every gene with an edge, but so do 53% of the null genes: I > 0 is not a
  test. The test keeps 37% of the genes and 1.3% of the null genes, close
  to its nominal rate. The co-expressed genes of this section are strongly
  spatial, so the filter matters only for networks with genes that have no
  spatial pattern.
- **weighted.** Smoothing raises the correlations: 15 times as many edges
  at 0.6 and 100 times as many at 0.3, but with STRING precision 3–8 times
  lower. Smoothing gives every gene a spatial pattern: after it, the
  median I of the null genes is 0.61 (0.64 for real genes). No null gene
  gets an edge at either cutoff, but the p-values count every spot as an
  independent sample, which smoothed spots are not. (A first version of the
  null genes permuted counts within bins of library size. They inherited
  the spatial pattern of the library size and got 710 edges at 0.3.)
- **regions.** The layer networks at 0.6 are small and precise (precision
  0.46–0.78, L6 has no STRING pair). At 0.3, the thin layers with
  218–273 spots (L1, L2, L4) give most of the edges, at precision 0.3–1.9%.

In short, the plain network and its filtered version nearly agree on this sample.
The weighted mode trades precision for edges. Modules of spatial data are
better found with `spatial-modules`, which reproduces hdWGCNA.

## Files

| file | contents |
|---|---|
| `01_prepare.py` | Visium → AnnData (counts, `array_row`/`array_col`, `region`, `obsm["spatial"]`) and the R export |
| `02_hdwgcna.R` | hdWGCNA's spatial workflow; writes datExpr, modules before / after ModuleConnectivity (with kME), spot MEs, stage timings |
| `run_all.sh` | both tools twice per sample, plus untimed gcen2 runs for the comparison (`--write-metaspots`, `--no-reassign`) |
| `04_compare.py` | the comparison table (`results/compare.tsv`) |
| `01_prepare_slideseq.py` | Slide-seqV2 → AnnData with the recovered counts, `obsm["spatial"]`, and the R export |
| `02_hdwgcna_knn.R` | hdWGCNA's workflow with `MetacellsByGroups` on the spatial coordinates (no grid) |
| `run_slideseq.sh` | Slide-seqV2: both tools twice, then metacell seeds 1 and 2 |
| `05_compare_slideseq.py` | the knn comparison: `python 05_compare_slideseq.py RUNS_DIR [SAMPLE] [TIME_FILE]` (`results/slideseqv2_compare.tsv`, `results/visium_hd_compare.tsv`) |
| `01_prepare_visium_hd.py` | Visium HD bins (one bin size) → AnnData with k-means regions, and the R export |
| `run_visium_hd.sh` | Visium HD: both tools twice, then metacell seeds 1 and 2 |
| `06_knn_vs_grid.py` | knn against grid metaspots on the Visium samples (`results/knn_vs_grid.tsv`) |
| `07_spatial_build.py` | `gcen2 spatial-build` on DLPFC: the three modes against the plain network, STRING and null genes (`results/spatial_build.tsv`, `results/spatial_build_moran.tsv`; `STRING_GOLD` = the PBMC benchmark's `work/string_gold.tsv`) |
| `results/` | `compare.tsv`, `time.txt`, hdWGCNA stage timings, gcen2 summaries, hub genes and spatial eigengene maps; the `slideseqv2_*` and `visium_hd_*` files, `time_slideseqv2.txt`, `time_visium_hd.txt`, `knn_vs_grid.tsv` and the `spatial_build*` files |
