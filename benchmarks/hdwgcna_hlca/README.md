# GCEN2 vs hdWGCNA on an atlas: the Human Lung Cell Atlas (585k cells)

`gcen2 celltype-modules` against hdWGCNA 0.4.12's single-cell workflow
(Seurat 5.5.1, WGCNA 1.74, R 4.4.3) on the core of the Human Lung Cell Atlas
(Sikkema et al. 2023, Nat Med; CELLxGENE dataset
`066943a2-fdac-4b29-b348-40cede398e4e`): 584,944 cells, 27,402 genes,
1.14 billion non-zero counts, 50 cell types. The settings are those of
`../hdwgcna_pbmc`, on two sizes: the whole atlas and a random subset of
50,000 cells (seed 0).

| step | both tools |
|---|---|
| genes | per cell type, genes with counts in ≥ 5% of its cells (`SetupForWGCNA(gene_select = "fraction")`) |
| metacells | per cell type, k = 25 nearest cells in the atlas' integrated scANVI embedding (`X_scanvi_emb`, 30 dims; hdWGCNA's `reduction`), at most 10 shared, up to 1,000 metacells, averaged, then log-normalized |
| network | `TestSoftPowers` (signed, bicor) + `ConstructNetwork` defaults (GCEN2: `--preset hdwgcna`) |

hdWGCNA runs on the 10 largest cell types (21,285–80,113 cells each,
430,202 in all; `results/types.txt`): its workflow takes 0.7–6 minutes per
type. GCEN2 runs on the same 10 types and, separately, on all 50.

R holds a sparse matrix of at most 2^31 − 1 non-zero values. The atlas has
fewer (1.14 billion), but `02_hdwgcna.R` still reads each cell type on its own
(the HDF5 file for R is sorted by cell type, so a type is one slice), as one
would for a larger atlas. Each type then goes through the whole workflow.

## Run

```bash
cd <hlca>      # a data directory of your choice; raw/hlca_core.h5ad: the CELLxGENE download below
python -I .../01_prepare.py raw/hlca_core.h5ad work 50000
bash .../run_all.sh     # both sizes; then 04_compare.py
```

Download: `https://datasets.cellxgene.cziscience.com/688185ad-11c2-4172-a53a-f4f1f4076860.h5ad`
(5.9 GB; counts in `raw.X`; SHA-256 `1cbdde1e513a31bdb6a3f70296e2d0612760d2af4032e239c94596ea87b44007`).

For the metacell seeds (below), `celltype-modules ... --seed 1` and `--seed 2` on
the same 10 types, into `runs/gcen2/<size>_top10_seed<s>`, and
`02_hdwgcna.R work/hlca_<size>.h5 work/types.txt runs/hdwgcna/<size>_seed<s> <s>`
for hdWGCNA, then `04_compare.py` again.

## Results

Per cell type (`results/compare.tsv`):

| | 50,000 cells | all 584,944 cells |
|---|---|---|
| same genes as hdWGCNA | 10 of 10 types | 10 of 10 |
| network step: `gcen2 module -e --preset hdwgcna` on hdWGCNA's metacells, ARI with hdWGCNA | 1.0 in 10 of 10 | 1.0 in 10 of 10 |
| network step time per type, GCEN2 / hdWGCNA | 1.9–8.2 s / 24–65 s | 1.9–8.5 s / 25–61 s |
| whole workflow (each tool's own metacells), ARI GCEN2 vs hdWGCNA, median (range) | 0.65 (0.15–0.90) | 0.63 (0.14–0.82) |
| GCEN2 alone, metacell seed vs seed, median (range) | 0.69 (0.47–0.84) | 0.61 (0.26–0.79) |
| hdWGCNA alone, metacell seed vs seed, median (range) | 0.64 (0.15–0.81) | 0.60 (0.12–0.79) |

ARI is over all genes, unassigned genes as one group. Seeds are 42, 1 and 2
for each tool: the tool-vs-tool row has all 9 seed pairs per type, the
seed-vs-seed rows the 3 pairs of one tool (with hdWGCNA's seed 42 only, the
tool-vs-tool row was 0.64 (0.37–0.90) and 0.61 (0.32–0.82)).
- Given the same metacells, GCEN2 gives exactly hdWGCNA's modules in all
  20 runs.
- With their own metacells, the two tools agree with each other as well as
  each agrees with itself under another seed (medians 0.60–0.69). The
  difference comes from which cells are drawn into metacells, not from the
  method; the low values (club cells at 50,000, CD8 T cells on the whole
  atlas) are low for hdWGCNA's own seeds too.

**Time and peak memory** (`results/time.tsv`; GCEN2 mean of 2 runs, hdWGCNA
1 run; 8 threads; both include the cell-level eigengenes and kME):

| | 50,000 cells | all 584,944 cells |
|---|---|---|
| GCEN2, 10 types | 75 s, 2.5 GB | 211 s, 11.9 GB |
| GCEN2, 10 types, `--backed` | 71 s, 1.7 GB | 206 s, 4.6 GB |
| hdWGCNA, 10 types | 821 s, 8.1 GB | 2,952 s, 43.9 GB |
| GCEN2, all cell types (34 / 49 with ≥ 100 cells) | 172 s, 2.5 GB | 466 s, 11.9 GB |
| GCEN2, all cell types, `--backed` | 176 s, 1.7 GB | 461 s, 4.6 GB |

On the 10 types GCEN2 is 11x (50,000 cells) and 14x (whole atlas) faster,
with 3.3x and 3.7x less memory (9.5x with `--backed`). On the whole atlas
most of GCEN2's memory is the count matrix itself (1.14 billion values as
float32 CSR, about 9 GB). With `--backed` the matrix stays on disk and each
cell type's rows are read when it is reached; the outputs are identical and
the time the same. The rows of a cell type are read once and shared by the
gene selection, the metacells and the cell level, and the cell-level
eigengenes and kME work through blocks of rows. Cell types are processed one
at a time: the metacells of one type are built when it is reached, and its
results are written and released before the next. Per cell type, hdWGCNA
spends 24–65 s on `TestSoftPowers` + `ConstructNetwork`, 8–42 s (50,000
cells) and 52–268 s (whole atlas) on the cell-level eigengenes and kME. Its
metacells take 7–14 s on the 50,000-cell subset but 25–267 s on the whole
atlas, where they grow with the cell type (267 s for the 78,816 alveolar
macrophages; `results/hdwgcna_timing_*.tsv`).

## Files

| file | contents |
|---|---|
| `01_prepare.py` | CELLxGENE h5ad → counts AnnData per size (`X_scanvi_emb`, `cell_type`) and an HDF5 file for R sorted by cell type; the 10 largest types |
| `02_hdwgcna.R` | hdWGCNA per cell type (reads one type's slice at a time) |
| `run_all.sh` | both sizes: GCEN2 on 10 and on all types (twice), hdWGCNA on 10 (once) |
| `04_compare.py` | genes, network-step and whole-workflow ARI, seed ARI, times (`results/compare.tsv`, `results/time.tsv`) |
| `results/` | the tables above, hdWGCNA's per-type timings and seed runs (`time_hdwgcna_seeds.txt`), GCEN2's per-type summaries (`gcen2_celltypes_*.tsv`) |
