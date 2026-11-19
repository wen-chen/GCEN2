# Changelog

All notable changes to GCEN2. The format follows
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/); versions follow
[PEP 440](https://peps.python.org/pep-0440/).

## [Unreleased] — 2.0.0.dev0

GCEN2 is a rewrite of [GCEN v1](https://github.com/wen-chen/GCEN) (C++
command-line tools) as a Python package with a C++ extension. It keeps every
v1 command and file format and adds WGCNA / hdWGCNA module analysis, single
cell, spatial, multi-omics and trajectory analyses, publication figures and
reports. Development history: `docs/development_log.md` and `git log` (earlier plans in
`docs/plans/archive/`); tests: `docs/test_report.md`; current state in
`docs/status.md`.

### Added

**GCEN v1 commands** (`gcen2 normalize`, `filter`, `build`, `module`,
`annotate`, `enrich`, `rwr`, `merge`, `extract`, `shuffle`, `stats`):
the same inputs and, where results are meant to agree, the same results as
GCEN v1 0.6.3 (identical edges and modules; see `docs/migration.rst`).

**Network construction**
- Pearson, Spearman, dropout-aware Pearson (`pearson_nonzero`,
  `pearson_pairwise_complete`); p-values, Benjamini-Hochberg FDR; edges
  streamed to the file block by block.
- A matrix-product (GEMM) correlation kernel for many samples
  (`--kernel gemm`) besides v1's pairwise kernel.
- WGCNA topological overlap as an extra edge score (`build --tom`).
- AnnData (`.h5ad`) input for `build` and `module -e`.

**WGCNA / hdWGCNA modules**
- `gcen2 module -e`: WGCNA's blockwiseModules (soft power by
  pickSoftThreshold, TOM, average-linkage tree, dynamic tree cut, eigengene
  merging); `--preset hdwgcna` for hdWGCNA's ConstructNetwork; eigengenes,
  kME, hub genes, module-trait correlations, gene significance, soft-power
  table, dendrogram and hub TOM for figures.
- `gcen2 consensus-modules`: WGCNA consensus networks over several sets.
- `gcen2.preservation` / `gcen2 module-preservation`: WGCNA modulePreservation
  (Zsummary, median rank) with R's random number generator, so the
  permutations are those of R.

**Single cell**
- kNN metacells (hdWGCNA MetacellsByGroups); `gcen2 celltype` (one network
  per cell type, metacell or pseudo-replicate samples) and
  `gcen2 celltype-modules` (the hdWGCNA workflow per cell type: metacells,
  modules, cell eigengenes and kME, reassignment of genes with negative kME,
  hub genes; `--backed` for atlases, one cell type at a time).
- Consensus networks over samples or studies (`--consensus-by`), Harmony
  integrated eigengenes (`--harmonize-by`, `gcen2.harmony`), differential
  module eigengenes (`--dme-by`) and module-trait correlations (`--trait`).
- `gcen2 project-modules` (hdWGCNA ProjectModules).
- `metacell_members.tsv`: the cells of every metacell.

**Spatial**
- `gcen2 spatial-modules`: the hdWGCNA spatial workflow (Visium grid
  metaspots, or kNN metaspots for Slide-seq, Visium HD and other platforms),
  spot eigengenes, region correlations.
- `gcen2 spatial-build`: networks filtered by spatial autocorrelation,
  on spatially smoothed expression, or per region; `--moran` writes every
  gene's Moran's I (`gcen2.spatial.spatial_autocorrelation`).

**Other analyses**
- `gcen2 trajectory`: networks in sliding windows along pseudotime, edge and
  degree dynamics.
- `gcen2 fuse` (weighted fusion of network layers, e.g. co-expression, PPI,
  regulation) and `gcen2 rwr-multiplex` (random walk on a multiplex network).
- `gcen2 diff`: differential co-expression (Fisher's z) between two networks,
  with Benjamini-Hochberg FDR (`fdr` column; `-f` calls changed edges by FDR)
  and Spearman's variance 1.06 / (n − 3) (`-m spearman`), as DGCA.

**Figures and reports**
- `gcen2.plot` / `gcen2 plot` (43 commands): figures for every analysis in one
  style (journal widths, editable text in PDF / SVG, one colour per module in
  every figure, WGCNA colour names or Okabe-Ito); PDF, PNG and SVG.
- `gcen2 report`: HTML or PDF report with vector figures and numbered
  captions; WGCNA sections, cell-type overview and metacells (`--coords`),
  merged / extracted networks, Moran's I; `--figures-dir` saves every figure.
- `gcen2.plotting`: the earlier plotting API, kept for compatibility.

**Python API**: every command is a function (`gcen2.network.build`,
`gcen2.module.identify`, `gcen2.tom_modules`, `gcen2.single_cell`,
`gcen2.spatial`, ...); `Network.to_adjacency()` returns the weighted
adjacency matrix (`weight=`, `absolute=`, `genes=`).

**Data, tests, packaging**
- `data/example/`: the GCEN v1 sample data, a PBMC 10k subset and a Visium
  mouse brain section (27 MB); `tools/make_example_data.py` rebuilds them.
- GO, association and gene-list files may be gzip-compressed.
- Golden-output regression tests (`tests/golden/`): every command on fixed
  inputs, outputs frozen; tests for every command, plot command and public
  function.
- Benchmarks against GCEN v1, WGCNA, hdWGCNA, DGCA and RandomWalkRestartMH
  (`benchmarks/`), reproducible
  elsewhere: environments (`benchmarks/environment-r.yml`,
  `requirements-python.txt`), downloads with checksums in each README, no
  fixed paths in the scripts.
- Supported: Python 3.9 to 3.14; numpy ≥ 1.24, pandas ≥ 2.0, scipy ≥ 1.10,
  anndata ≥ 0.10, matplotlib ≥ 3.7, typer ≥ 0.9 (the oldest versions tested).
- `LICENSE` (GPL-3.0-or-later), `CITATION.cff`, `CONTRIBUTING.md`,
  `CODE_OF_CONDUCT.md` (Contributor Covenant 2.1), issue and pull request
  templates; `.readthedocs.yaml` for the documentation.
- The source distribution includes the example data, so its tests run
  (verified: built from the sdist in a clean environment, 912 passed, 1 skipped).
- `docs/development_log.md` (development by phase) and `docs/test_report.md`.
- Benchmarks against DGCA 1.0.3 (`diff`, `benchmarks/dgca_gtex`),
  RandomWalkRestartMH (`rwr-multiplex`) and R implementations of the fusion
  (`benchmarks/multilayer_string`) and of the networks along pseudotime
  (`benchmarks/trajectory_paul15`); `tests/data/rwrmh/` holds
  RandomWalkRestartMH reference scores.
- A biological case study on GTEx whole blood (`benchmarks/case_gtex_blood`):
  GO enrichment, hub genes and donor-trait correlations of the WGCNA-style
  modules (checked in two halves of the donors), and cross-validated RWR
  function prediction against neighbour voting, a degree-preserving shuffle
  and node degree.
- Thread scaling and parameter sensitivity (`benchmarks/scaling`): `-t` 1–16
  for `build`, `module -e`, `diff` and `celltype-modules`; the correlation
  cutoff, GCEN module similarity and WGCNA module parameters against STRING.
- PyWGCNA 2.2.1 against WGCNA (`benchmarks/pywgcna_gtex`): the same soft
  power, TOM and gene tree, but different modules (ARI 0.30–0.44 with WGCNA's
  step-by-step pipeline on the corrected data), traced to its dynamic tree cut.
- Cells as correlation samples (`benchmarks/hdwgcna_pbmc/09_strategies.py`):
  metacells against random pseudo-replicates, single cells and MAGIC
  imputation, with null genes for false edges, STRING and split halves.
- `spatial-build` on DLPFC (`benchmarks/hdwgcna_spatial/07_spatial_build.py`):
  the filtered, weighted and regions modes against the plain network.

### Changed

- `spatial-build -m filtered` keeps the edges between genes with significant
  positive spatial autocorrelation: Moran's I > 0 and a Benjamini–Hochberg
  FDR <= `--moran-fdr` (default 0.05) of a one-sided test under
  randomisation, instead of I > 0 alone (`--moran-fdr 1`). The `--moran`
  table gains `pval` and `fdr` columns.
- `gcen2 celltype` (`build_by_celltype`) tests every gene pair of a metacell
  network on its own effective sample size: each gene gets M / s(y), with
  s(y) = y'Cy / y'y and C the metacells' shared cells / k
  (`metacell_effective_samples`); a pair takes the smaller of its two genes'
  values. The former fixed n_cells // k (`--effective-samples cells`) was
  too strict for most genes and too lenient for the smoothest cell-state
  genes; the new sizes match null-gene correlations on PBMCs
  (`benchmarks/hdwgcna_pbmc/10_effective_samples.py`). Networks change:
  1.0–5.6 times as many edges on PBMCs, most in the small cell types.
- `gcen2 celltype` (`build_by_celltype(fdr=True)`) applies the
  Benjamini-Hochberg FDR over all gene pairs by default, `-p` applying to it
  (`--no-fdr` for the former behaviour). With calibrated p-values, p <= 0.001
  alone kept about 0.1% of a cell type's 30–40 million pairs as chance edges,
  which made up most networks and modules of small cell types. PBMC networks:
  355–122,020 edges instead of 23,245–302,654.
- Log-normalisation (`spatial-build`, metacells of `celltype` /
  `celltype-modules`, cell eigengenes) computes x · 10,000 / library with one
  rounding, so mathematically equal ratios give equal values. Spearman
  networks of normalised data no longer depend on the CPU: on CI runners with
  AVX-512, NumPy's `log1p` differs in the last bit and reordered near-equal
  values (correlations then differed in the sixth decimal). Values change in
  the last bits.
- C++ builds with Apple's libc++ (no floating-point `std::from_chars`) and
  with MSVC (a missing `<iterator>`; `1.0 / 0.0` replaced by
  `numeric_limits::infinity()`).
- Text files are read and written as UTF-8 on every platform (Windows
  defaulted to its ANSI code page): the network report of `plot network
  --interactive` (pyvis wrote it in the default encoding), tree and seed files
  of `gcen2 plot`, and the developer tools.
- Output files have `\n` line endings on every platform: every `to_csv` call
  passes `lineterminator="\n"` and text files are opened with `newline="\n"`
  (pandas and Python wrote `\r\n` on Windows, so files differed from Linux
  and from the C++ writers).
- k nearest neighbours (spatial graphs of `spatial-build` and Moran's I, kNN
  metaspots, metacells) break ties of equal distance by index
  (`gcen2._knn.knn_indices`). `cKDTree` orders equidistant points differently
  on different platforms, and grids of spots or bins have many, so macOS gave
  other neighbours than Linux. Results on grids can change where ties were
  broken differently.
- `differential.compare` without `n_a` / `n_b` takes each network's p-value
  sample size from its metadata: the per-gene effective samples of a metacell
  network, else `p_value_samples`, else `n_samples`. Before, a metacell
  network from `build_by_celltype` was tested on its number of metacells
  (the quick start's example did this).
- The migration guide covers GCEN v1's four remaining utilities: Python
  functions for `data_stat`, `generate_expr_matrix_from_rsem` and
  `generate_expr_matrix_from_stringtie` that give v1's values (tested against
  v1 0.6.3 in `tests/test_migration_recipes.py`, which runs the guide's own
  code), and what `calculate_accuracy` computed (no counterpart).
- Reading a header-less GCEN v1 CSV matrix whose first gene name is a quoted
  field holding a comma no longer takes that line for a header (the header
  check now splits the line with quotes, as the reader does). The migration
  guide now says that GCEN2 has no `csv_to_tsv` / `tsv_to_csv` commands and
  converts in Python (`gcen2.io.read_csv` + `write_tsv`).
- `gcen2 celltype --diff` tests every shared edge on the effective sample
  sizes of the networks' own p-values (the smaller of its two genes' in each
  cell type) instead of n_cells // k: on random halves of PBMC cell types
  4.1–4.9% of the pairs have p <= 0.05, against 0.26–0.76% before
  (`benchmarks/hdwgcna_pbmc/11_diff_calibration.py`). `differential.compare`
  and `compare_files` accept mappings gene → samples for `n_a` / `n_b`.
  `celltypes.tsv`'s `effective_samples` and the diff summary's `n_a` / `n_b`
  are now the median gene's.
- `gcen2 celltype --strategy pseudo_replicates` is deprecated (FutureWarning):
  random groups of cells rank gene pairs at chance level.
- `build_dropout_aware(zero_handling="impute")` warns that the p-values of
  MAGIC-imputed cells are not valid.

- `gcen2 annotate` and `gcen2 enrich` (and `gcen2.annotation.enrich` /
  `annotate_go` / `annotate_kegg`) add an `fdr` column: the Benjamini-Hochberg
  adjusted p-value of each term over every term annotated in the background
  (as R's `p.adjust(method = "BH")`), computed exactly in the same C++ pass
  as the counts and p-values; `-p` still filters on the p-value. The other
  columns are unchanged (golden outputs regenerated for the new column).

- The command line is a package (`gcen2/cli/`, one module per group of
  commands) and `gcen2 --help` lists the commands by group (GCEN v1
  commands, WGCNA modules, single cell, spatial, trajectory, multi-omics,
  report, plot). Every command reports errors the same way (`ERROR: ...`,
  exit code 1, no traceback; e.g. a missing input file of `spatial-build`
  or `trajectory`) and prints warnings as `WARNING: ...`.
- The options `gcen2 --threads / --seed / --verbose` before the command are
  gone: they were never used (each command has its own `-t` and `--seed`,
  as the documentation says).
- A missing input file of a `gcen2.plot` function raises `GCENIOError` (as
  everywhere else; it was `GCENInputError`).
- `single_cell.build_by_celltype(max_genes_per_chunk=...)` no longer splits
  the genes (the argument is accepted and ignored): the correlation kernel
  already works through the gene pairs in blocks, and with `fdr=True` the
  split adjusted the p-values within each block instead of over all pairs.
- Eigen is no longer bundled (the kernels stopped using it) and OpenMP is no
  longer linked (the kernels use `std::thread`): clone without `--recursive`,
  no third-party C++ library besides pybind11.
- `gcen2 rwr-multiplex` / `multilayer.rwr_multiplex` column-normalise the
  whole supra-adjacency matrix by default, as RandomWalkRestartMH (the
  RWR-M authors' implementation), and give its scores (`--transition supra`;
  weighted layers divided by their largest weight, as it does). The previous
  per-layer normalisation, where every gene jumps with probability `delta`,
  is `--transition layer` (unchanged results). Default scores change.
- `gcen2 diff` / `differential.compare` write an `fdr` column (Benjamini-
  Hochberg over the tested pairs) before `category`; `-f / --fdr` calls
  `shared_changed` by FDR, `-m spearman` uses Spearman's variance
  1.06 / (n − 3). `gcen2 celltype --diff` uses that variance when both
  networks are Spearman networks (the default), so fewer edges are called
  changed.

### Changed (from GCEN v1)

Deliberate differences; details and how to reproduce v1 output in
`docs/migration.rst`.
- `build -f` also drops edges with FDR > p cutoff (v1 only adds the column).
- RWR: seed vector normalised (values are probabilities), seeds not listed,
  `>=` cutoff, `|r|` as edge weight; ties ranked by gene name.
- Enrichment: one-sided hypergeometric test (v1: two-sided Fisher), fold
  enrichment column; v1's p-values are reproduced from GCEN2's counts.
- Output files have header lines; `module` writes one `module gene` pair per
  line; `merge` writes the direction as its own column.
- `filter` uses the sample standard deviation (ddof = 1).

### Performance

C++ kernels for correlation (pairwise and GEMM), TOM, FDR, module
identification, enrichment, RWR, network statistics, edge shuffling, merge,
extract, diff and fusion; most commands read and write files without pandas.
Against GCEN v1, WGCNA and hdWGCNA GCEN2 was 1.1× to 89× faster, against DGCA
(Pearson, Spearman) 2.6× and 185× and RandomWalkRestartMH 58×
(`benchmarks/README.md`, measured on an idle machine on 2026-11-17).

- `gcen2 diff -t` (default 2): the edge files are parsed, sorted and written
  by several threads and the FDR is sorted in parallel (12.5M pairs: 13.6 s
  → 5.5 s with 8 threads); the output does not depend on the thread count.
  `celltype --diff` uses the command's `-t`.

### Fixed

- GCEN v1's `merge` appended the direction to the correlation string.
- On Python 3.9 no `gcen2` command ran (annotations that typer evaluates).
- With typer < 0.12, list options left out (e.g. `--group`) selected nothing.
- R's random number generator failed with numpy 1.x.
- TIFF / EPS figures failed with matplotlib < 3.9.
- Example and test `.h5ad` files were unreadable for anndata < 0.11.
- RWR ranked tied genes differently with numpy 1 and 2.
- `gcen2.plotting.network_interactive` failed with matplotlib ≥ 3.9 and wrote
  pyvis's `lib/` into the working directory.
- A missing output directory raised a traceback; it is now created.
- Clear errors for a missing report template and for weasyprint without its
  system libraries.
- The BLAS thread limit (`-t`) missed scipy's own OpenBLAS when it was first
  loaded inside the limit: `gcen2 module -e` computed the eigengenes on every
  logical CPU (GTEx, 5,000 genes: 3.3 → 1.7 s; the same modules).
- `gcen2 plot enrich-dot` and `enrich-heatmap` labelled the colour bar
  "−log10 p" when the tables are coloured by an adjusted column (`fdr`,
  `p_adj`); the label now names the column used, as `enrich-bar` does.

Found while benchmarking: hdWGCNA pairs the soft powers with the wrong sets in
a consensus network unless the sets come in sorted order
(`benchmarks/hdwgcna_consensus/README.md`).
