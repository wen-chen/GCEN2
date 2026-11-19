# GCEN2

**A Python Toolkit for Gene Co-Expression Network Analysis with Single-Cell Resolution and Interactive Visualization**

---

**Software Requirements Specification**

| | |
|---|---|
| Document Version | 2.2 |
| Date | November 2026 (2.1: November 2026; 2.0: March 2026) |
| Author | Wen Harold Chen |
| Affiliation | Hunan University of Chinese Medicine |
| Contact | chenwen@biochen.org |

---

## Table of Contents

- [1. Project Overview](#1-project-overview)
  - [1.1 Background](#11-background)
  - [1.2 Project Goals](#12-project-goals)
  - [1.3 Target Users](#13-target-users)
  - [1.4 Document Conventions](#14-document-conventions)
- [2. System Architecture](#2-system-architecture)
  - [2.1 Architecture Overview](#21-architecture-overview)
  - [2.2 C++ vs. Python Division of Labor](#22-c-vs-python-division-of-labor)
- [3. Functional Requirements](#3-functional-requirements)
  - [3.1 Core Pipeline (Retained from GCEN v1)](#31-core-pipeline-retained-from-gcen-v1)
  - [3.2 New Feature: Single-Cell Support](#32-new-feature-single-cell-support)
  - [3.3 New Feature: Visualization](#33-new-feature-visualization)
  - [3.4 New Feature: Research Report Generation](#34-new-feature-research-report-generation)
  - [3.5 New Feature: Differential Co-Expression Analysis](#35-new-feature-differential-co-expression-analysis)
  - [3.6 New Feature: Trajectory Co-Expression Networks](#36-new-feature-trajectory-co-expression-networks)
  - [3.7 New Feature: Multi-Omics Network Integration](#37-new-feature-multi-omics-network-integration)
  - [3.8 New Feature: Spatial Transcriptomics Support](#38-new-feature-spatial-transcriptomics-support)
  - [3.9 New Feature: LLM-Assisted Interpretation](#39-new-feature-llm-assisted-interpretation)
  - [3.10 New Feature: WGCNA Module Analysis](#310-new-feature-wgcna-module-analysis)
  - [3.11 Example Data and Reproducibility](#311-example-data-and-reproducibility)
- [4. API Design](#4-api-design)
  - [4.1 Python API Overview](#41-python-api-overview)
  - [4.2 Examples and CLI](#42-examples-and-cli)
- [5. Non-Functional Requirements](#5-non-functional-requirements)
  - [5.1 Performance](#51-performance)
  - [5.2 Compatibility](#52-compatibility)
  - [5.3 Dependencies](#53-dependencies)
  - [5.4 Testing](#54-testing)
  - [5.5 Documentation](#55-documentation)
  - [5.6 Licensing](#56-licensing)
- [6. Status](#6-status)
- [7. Known Issues from GCEN v1 (Status)](#7-known-issues-from-gcen-v1-status)

---

## 1. Project Overview

### 1.1 Background

GCEN (Gene Co-Expression Network) is an open-source command-line toolkit for gene co-expression network analysis and lncRNA annotation, first published in *Current Issues in Molecular Biology* (2022). Written entirely in C++, GCEN provides high-performance computation but lacks Python ecosystem integration, visualization capabilities, and single-cell transcriptomics support. These limitations increasingly restrict its applicability as the field moves toward single-cell resolution and Python-based bioinformatics workflows.

### 1.2 Project Goals

GCEN2 is a major upgrade that aims to address the following objectives:

1. Redesign the software architecture with Python as the primary user interface while retaining C++ for performance-critical computation.
2. Add comprehensive visualization and automated research report generation capabilities.
3. Provide dedicated support for single-cell RNA-seq data analysis.
3a. Reproduce the module analyses of WGCNA (bulk) and hdWGCNA (single cell, spatial) with the same results and much less time and memory.
4. Introduce new analytical features including differential co-expression analysis, trajectory-based dynamic networks, multi-omics integration, and LLM-assisted biological interpretation.
5. Maintain backward compatibility with existing GCEN v1 workflows through a CLI compatibility layer.
6. Fix known bugs from GCEN v1 (network_merge direction column, deviation ddof, CSV parsing).

### 1.3 Target Users

- Bioinformatics researchers performing gene co-expression network analysis
- Biologists working with bulk or single-cell RNA-seq data who need functional annotation
- Computational biologists building multi-omics analysis pipelines
- Graduate students and postdocs requiring accessible, code-free network analysis tools

### 1.4 Document Conventions

| Priority | Description |
|---|---|
| P0 (Must Have) | Core functionality required for initial release |
| P1 (Should Have) | Important features planned for initial release |
| P2 (Nice to Have) | Features for subsequent releases |
| ✅ Done | Implemented and tested |
| 🔲 Pending | Not yet implemented |

---

## 2. System Architecture

### 2.1 Architecture Overview

GCEN2 adopts a layered architecture with Python as the orchestration layer and C++ as the computation engine. The architecture consists of four layers: the User Interface Layer (Python API + CLI), the Application Logic Layer (Python modules), the Computation Engine Layer (C++ with pybind11 bindings), and the Data Layer (file I/O and format conversion).

### 2.2 C++ vs. Python Division of Labor

| Component | Language | Rationale |
|---|---|---|
| Pearson/Spearman correlation | C++17 (pybind11) | O(n²) pairwise computation, multithreaded |
| Hypergeometric p-value | C++17 (pybind11) | Called millions of times in annotation |
| Jaccard similarity + DFS | C++17 (pybind11) | O(n²) similarity matrix, memory intensive |
| TOM, average-linkage clustering | C++17 (pybind11) | O(n³) matrix products on genes × genes |
| RWR | C++17 (pybind11) | Sparse matrix-vector iterations |
| Network files: build, statistics, merge, extract, shuffle, diff, fusion | C++17 (pybind11) | Read and write large edge lists without pandas |
| Data normalization | Python (numpy/pandas) | I/O bound, pandas more flexible |
| Data filtering | Python (numpy/pandas) | Simple statistics |
| File I/O and format conversion | Python (pandas/AnnData) | Native support for h5ad, CSV, TSV |
| Visualization | Python (matplotlib/pyvis) | Rich plotting ecosystem |
| Report generation | Python (Jinja2) | Template-based, flexible output formats |

The files of the repository are described in the Code Guide, the build in the
Implementation Document.

---

## 3. Functional Requirements

### 3.1 Core Pipeline (Retained from GCEN v1)

All existing GCEN functionality is preserved and exposed through both the Python API and the CLI compatibility layer.

#### 3.1.1 Data Preprocessing [P0] ✅

- Normalization methods: Upper Quartile (upqt), Median, DESeq, TMM, Housekeeping Gene (hkg)
- Filtering by expression mean, standard deviation, or percentile
- Log transformation: log (natural), log2, log10 with configurable pseudo-count
- Sample standard deviation (`ddof=1`) used throughout, fixing GCEN v1's population std deviation

#### 3.1.2 Network Construction [P0] ✅

- Pearson and Spearman correlation with multithreaded C++ computation (pairwise kernel as GCEN v1, or a matrix-product kernel for many samples)
- P-value (two-sided t-test) and Benjamini-Hochberg FDR correction
- Signed/unsigned network options
- Network merging (with direction as a separate column, fixing GCEN v1 bug), extraction, shuffling, and statistics
- Degree distribution, density, heterogeneity, clustering coefficient, transitivity, connected components

#### 3.1.3 Module Identification [P0] ✅

- Jaccard similarity-based module detection with DFS clustering
- Configurable similarity cutoff and multithreaded C++ computation
- Singleton filtering (modules with < 2 genes excluded)

#### 3.1.4 Functional Annotation [P0] ✅

- GO annotation with OBO file parsing and DAG-based term propagation
- KEGG annotation via K2KO mapping files
- Enrichment analysis using hypergeometric test (C++ kernel for batch calls)
- Standalone enrichment for arbitrary gene lists against a background

#### 3.1.5 Random Walk with Restart [P0] ✅

- Gene function prediction via network propagation
- Support for directed/undirected and weighted/unweighted networks
- Seed gene weighting; Python fallback when C++ extension is unavailable

### 3.2 New Feature: Single-Cell Support

#### 3.2.1 Data Format Adaptation [P0] ✅

Native support for AnnData (`.h5ad`) as the primary single-cell data format. Integration with scanpy workflows is a core design requirement. Sparse matrices are handled at the Python boundary before conversion to dense numpy arrays for the C++ kernels. `network.build` / `build_to_file` accept an AnnData directly (`layer`, `use_highly_variable`), and `gcen2 build` reads `.h5ad` with `--layer`, `--hvg` and `--groupby`/`--group` cell subsets; the sparse matrix is densified only after subsetting, straight into genes × cells order.

#### 3.2.2 Metacell / Pseudo-bulk Aggregation [P0] ✅

Aggregation of single-cell data into pseudo-bulk samples by cell type or any `obs` key. Supports `sum` and `mean` aggregation. Groups with fewer than `min_cells` cells are excluded with a warning.

kNN metacells (`gcen2.single_cell.metacells`, after hdWGCNA `MetacellsByGroups`): within each group (one or more `obs` keys), each cell's `k` nearest neighbours in PCA space (`X_pca` or computed) form a candidate metacell; candidates sharing more than `max_shared` cells with a chosen one are rejected; members are averaged (or summed). Returns an AnnData that `network.build` accepts directly.

#### 3.2.3 Cell-Type-Specific Networks [P0] ✅

Build separate co-expression networks for each cell type. By default (`strategy="metacell"`) the correlation samples are the kNN metacells of each cell type, log-normalised when the data are counts. Because metacells overlap and each averages `k` cells, p-values and FDR are computed on the effective number of independent samples rather than the number of metacells; counting every metacell yields thousands of spurious edges on simulated data. The effective number is computed per gene pair from the cells the metacells share (`metacell_effective_samples`: M / s(y), the smaller of the two genes'), which matches the variance of null-gene correlations on PBMCs; the former fixed `n_cells // k` (`--effective-samples cells`) was too strict for most genes and too lenient for the smoothest cell-state genes. Because these p-values are calibrated, `p <= 0.001` alone would keep about 0.1% of the 30–40 million gene pairs of a cell type as chance edges; the Benjamini–Hochberg FDR over all pairs is therefore on by default (`fdr=True`; `--no-fdr`). CLI: `gcen2 celltype` writes one network per cell type plus a summary table and, with `--diff`, the differential co-expression of every pair of cell types (Fisher's z; every edge on the effective sample sizes of the networks' own p-values, i.e. the smaller of its two genes' in each cell type; calibrated on random halves of PBMC cell types, `benchmarks/hdwgcna_pbmc/11_diff_calibration.py`). The former strategy, cells randomly partitioned into `pseudo_replicates` (default 10) summed samples, remains available as `strategy="pseudo_replicates"` but is deprecated (FutureWarning): on simulated module data it recovers ~11% of the true top edges versus ~100% for metacells, and on PBMCs 10 groups rank gene pairs at chance level (`benchmarks/hdwgcna_pbmc/09_strategies.py`).

#### 3.2.4 Dropout Handling [P1] ✅

Per-gene-pair masked Pearson, computed in the C++ kernel (`network.build` methods `pearson_nonzero` and `pearson_pairwise_complete`, p-values on the number of cells used), exposed by `single_cell.build_dropout_aware(zero_handling=...)`:

- `both_nonzero` (default): only the cells where both genes are non-zero (zeros treated as missing values). On simulated data with random technical dropout (60–80%) it recovers 87–99% of the true top edges versus 21–63% for plain Pearson; when zeros mostly reflect low expression (UMI-like sampling) it does worse than plain Pearson (1–17% vs 76–99%).
- `pairwise_complete`: only zero-zero cells ignored. Kept for completeness and warns: the remaining cells with one gene at zero make unrelated genes anti-correlated (mean r ≈ −0.43 on the simulation), so it fails in every scenario tested.
- `impute`: MAGIC imputation (requires `magic-impute`). Warns: imputed cells are not independent samples, so the p-values are not valid; on PBMCs 12–27% of the pairs with a permuted null gene became edges.

Before this change `pairwise_complete` was documented but not applied (it built a plain network).

#### 3.2.5 Scalability Optimization [P1] ✅

Gene pairs are correlated in blocks and written block by block (`build_to_file`, `gcen2 build`), so the memory of a network does not grow with the number of pairs; highly variable gene subsetting (`use_highly_variable`, `--hvg`); atlases are read one cell type at a time from disk (`celltype-modules --backed`). WGCNA modules need one genes × genes matrix in memory (as WGCNA's single block).

### 3.3 New Feature: Visualization

#### 3.3.0 Publication Figures (`gcen2.plot`) [P0] ✅

One figure module for every analysis (43 `gcen2 plot` commands, also `gcen2.plot.<name>` in Python): WGCNA / hdWGCNA module figures (soft power, dendrogram, eigengene network, module-trait, kME, hub networks, module UMAP, TOM heatmap), single-cell and spatial figures (eigengenes on embeddings and tissue, dot / violin / radar plots, DMEs, metacells, preservation, Moran's I), network, enrichment, RWR, differential, trajectory, fusion, merge and extract figures. One style: journal column widths (89 / 120 / 183 mm), 5-7 pt text kept editable in PDF / SVG, one colour per module in every figure (WGCNA colour names or Okabe-Ito); PDF, PNG and SVG. The subsections below describe the earlier `gcen2.plotting` API, kept for compatibility.

#### 3.3.1 Network Visualization [P0] ✅

- Interactive HTML visualization (pyvis): force-directed layout, node color by module or degree, edge width by correlation, tooltip with gene name and degree
- Static publication-quality figures (networkx + matplotlib): spring layout, color-coded modules

#### 3.3.2 Enrichment Visualization [P0] ✅

- Dot plot: x-axis = enrichment ratio, y-axis = term name, size = gene count, color = −log₁₀(p-value)
- Bar plot: horizontal bars ranked by p-value with color by ontology namespace

#### 3.3.3 Expression Visualization [P1] ✅

- Box plots and violin plots for expression distribution QC (per-sample)
- Heatmaps (seaborn clustermap) for gene expression patterns within identified modules

#### 3.3.4 Degree Distribution [P1] ✅

Degree distribution histogram with power-law fit overlay.

#### 3.3.5 Differential Network Visualization [P1] ✅

Edge-colored matplotlib figure showing condition-specific (A_only, B_only) and shared rewired/stable edges.

### 3.4 New Feature: Research Report Generation

#### 3.4.1 Automated Report Pipeline [P1] ✅

Jinja2-based HTML report (with optional PDF export via weasyprint) integrating:
1. Network overview: node/edge counts, density, degree distribution plot
2. Merged and extracted networks; spatial autocorrelation (Moran's I)
3. Module results: module count, size distribution, top modules; WGCNA-style module sections (soft power, dendrogram, eigengenes, kME, hub networks, traits, DMEs, consensus sets), cell-type overview and metacells
4. Functional annotation: enriched terms plot + summary table
5. RWR prediction: ranked gene table

Figures are `gcen2.plot` figures embedded as vector SVG with numbered captions; `--figures-dir` also saves each as PDF, PNG and SVG.

#### 3.4.2 PDF Export [P1] ✅

Optional PDF generation via weasyprint (optional dependency). Automatically activated when `output` ends with `.pdf`.

### 3.5 New Feature: Differential Co-Expression Analysis

#### 3.5.1 Condition-Specific Network Comparison [P1] ✅

Fisher's Z-transformation test for each gene pair present in both networks (variance 1.06/(n−3) for Spearman correlations), with Benjamini-Hochberg FDR over all tested pairs. Edge categories: `A_only`, `B_only`, `shared_changed`, `shared_stable` (changed: p, or optionally the FDR, at or below the cutoff).

**Implemented:** `gcen2/differential.py`, `gcen2 diff` (`-m spearman`, `-f`). Verified against DGCA 1.0.3 on 12.5 million gene pairs, Pearson and Spearman: the same z, p and FDR (`benchmarks/dgca_gtex`).

### 3.6 New Feature: Trajectory Co-Expression Networks

#### 3.6.1 Pseudotime-Based Dynamic Networks [P2] ✅

For single-cell data with pseudotime information, construct co-expression networks along the trajectory using sliding windows. Window size, step size, and correlation thresholds are user-configurable.

**Implemented:** `gcen2/trajectory.py` — `build_trajectory(data, pseudotime, window_size, step, metacell_size, ...)` takes AnnData (an `obs` pseudotime column) or a genes × samples DataFrame (e.g. a bulk time course). Windows hold equal numbers of cells (count or fraction) and overlap by default (step = half a window); optional metacells average consecutive cells within a window. The result holds one `Network` per window and an edge-dynamics table: the correlation of every edge in every window, the windows where it passes the cutoffs, and a pattern (`stable`, `gained`, `lost`, `transient`, `intermittent`) with the change and trend along pseudotime. CLI: `gcen2 trajectory`; figure: `gcen2 plot trajectory`.

### 3.7 New Feature: Multi-Omics Network Integration

#### 3.7.1 Heterogeneous Network Fusion [P2] ✅

Extend `network.merge()` to support heterogeneous network types (co-expression, PPI, regulatory). Weighted fusion and multi-layer RWR propagation are the primary integration strategies.

**Implemented:** `gcen2/multilayer.py` (a separate module; `network.merge()` keeps its GCEN v1 semantics of summing signed correlations). Layers are `Network`s, DataFrames or edge-list files (GCEN networks, STRING-style scored PPI files, two-column regulatory lists). `fuse()` scales each layer's weights (`max`, `rank` or `none`) and scores every edge by the weighted mean over layers, recording which layers support it; `min_layers` / `min_score` filter the result. `rwr_multiplex()` runs RWR-M (Valdeolivas et al., *Bioinformatics* 2019) on the supra-adjacency matrix, by default normalised as in the authors' RandomWalkRestartMH (the same scores; `transition="layer"` gives every gene the jump probability `delta`), using the C++ RWR kernel, and aggregates the per-layer scores per gene (geometric mean by default). With one layer it reproduces `rwr.predict`. CLI: `gcen2 fuse`, `gcen2 rwr-multiplex`.

### 3.8 New Feature: Spatial Transcriptomics Support

#### 3.8.1 Spatial Co-Expression Networks [P1] ✅

Support AnnData input with `obsm['spatial']` 2D coordinates (Visium, MERSCOPE, Xenium). Three analysis modes:

1. **Spatial Filtering** (`build_spatial_filtered`) — build standard co-expression network, then retain only edges where both genes have significant positive spatial autocorrelation: Moran's I > 0 and Benjamini–Hochberg FDR ≤ `moran_fdr` (default 0.05) of a one-sided test with the moments of I under randomisation (Cliff and Ord); `moran_fdr=1` keeps the former rule, I > 0.
2. **Spatial-Weighted Correlation** (`build_spatial_weighted`) — smooth expression by Gaussian-weighted neighbourhood averaging before building a standard co-expression network; bandwidth auto-estimated from median neighbour distance.
3. **Region-Based Networks** (`build_spatial_regions`) — partition tissue into spatial regions (manual `obs` column annotation or automatic KMeans clustering on coordinates), build one network per region (analogous to `build_by_celltype`).

Neighbor graph: KNN or radius cutoff; mutually exclusive. Exposed as `gcen2.spatial` module and `gcen2 spatial-build` CLI command; `--moran` writes every gene's Moran's I (`spatial_autocorrelation`, also after smoothing in weighted mode). Spatial modules (the hdWGCNA workflow): see 3.10.3.

**Input:** `AnnData` with `obsm['spatial']` (shape `n_spots × 2`).
**Output:** `Network` (filtered/weighted) or `dict[str, Network]` (regions).
**Dependencies:** `anndata`, `numpy`, `scipy` (already in `single-cell` extra).
**Implemented:** `gcen2/spatial/` — `neighbors.py`, `core.py`, `modules.py`, `__init__.py`.

### 3.9 New Feature: LLM-Assisted Interpretation

#### 3.9.1 Automated Biological Interpretation [P2] 🔲

Pass enrichment results, key module information, and network statistics to a Large Language Model via API (Claude API, OpenAI API, or local models). Generated text is included in the research report with clear AI attribution.

### 3.10 New Feature: WGCNA Module Analysis

#### 3.10.1 WGCNA Modules [P1] ✅

`gcen2 module -e` (`gcen2.tom_modules`) reproduces WGCNA's `blockwiseModules`: soft-thresholding power by `pickSoftThreshold` (scale-free fit), signed or unsigned adjacency, topological overlap, average-linkage clustering of 1 − TOM, dynamic tree cut, merging of modules with similar eigengenes; module eigengenes, kME, hub genes, eigengene-trait correlations and gene significance. `--preset hdwgcna` reproduces hdWGCNA's `ConstructNetwork` (bicor power, signed, deepSplit 4). Outputs for figures: soft-power table, dendrogram, hub TOM. Verified: the same modules as WGCNA 1.74 and hdWGCNA 0.4.12 (ARI 1.0).

#### 3.10.2 Consensus Modules and Module Preservation [P1] ✅

- `gcen2 consensus-modules`: WGCNA consensus modules (one TOM per data set, calibrated, combined by a quantile).
- `gcen2 module-preservation`: WGCNA `modulePreservation` (Zsummary, median rank) with R's random number generator, so that the permutations are those of R; the same statistics as WGCNA (relative difference < 7e-8).

#### 3.10.3 Single-Cell and Spatial Workflows (hdWGCNA) [P1] ✅

- `gcen2 celltype-modules`: per cell type, kNN metacells, modules, cell eigengenes and kME, reassignment of genes with negative kME, hub genes; consensus over samples (`--consensus-by`), Harmony-integrated eigengenes (`--harmonize-by`), differential module eigengenes (`--dme-by`), module-trait correlations (`--trait`); `--backed` for atlases. Writes the members of every metacell.
- `gcen2 project-modules`: hdWGCNA `ProjectModules`.
- `gcen2 spatial-modules`: hdWGCNA's spatial workflow on Visium grid metaspots or kNN metaspots (Slide-seq, Visium HD), spot eigengenes and region correlations.

### 3.11 Example Data and Reproducibility

#### 3.11.1 Example Data [P0] ✅

The repository holds example data for every workflow (`data/example/`, 27 MB): the GCEN v1 sample data with GO and KEGG annotations, a PBMC 10k subset (1,800 cells, three cell types) and a Visium mouse brain section, with sources and licences; GO and annotation files may be gzip-compressed.

#### 3.11.2 Deterministic Results [P0] ✅

The same input gives the same output with any number of threads and every supported numpy / pandas version (ties ordered explicitly, e.g. RWR by gene name). Figures and reports are reproducible (fixed SVG ids, no creation dates).

---

## 4. API Design

### 4.1 Python API Overview

The Python API follows scanpy-style conventions. Key design principles:
- AnnData as the primary data container for single-cell workflows
- pandas DataFrame as the primary format for bulk workflows
- Sensible defaults matching the GCEN v1 recommended pipeline
- Threads and random seeds as arguments of each function and options of each command

### 4.2 Examples and CLI

Examples of the Python API and the command line are in the quickstart and the
tutorials (`docs/`); every command and option is in `docs/cli.rst`. The GCEN v1
commands keep v1's flags (`docs/migration.rst`). GCEN2 adds `module -e`,
`consensus-modules`, `module-preservation`, `celltype`, `celltype-modules`,
`project-modules`, `spatial-modules`, `spatial-build`, `trajectory`, `fuse`,
`rwr-multiplex`, `merge`, `extract`, `shuffle`, `diff` and `plot` (43 figure
commands). Commands that compute in parallel take `-t` / `--threads`; random
steps take `--seed` (there are no global options).

---

## 5. Non-Functional Requirements

### 5.1 Performance

| Metric | Requirement | Notes |
|---|---|---|
| Network build (10k genes) | < 5 minutes on 8-core machine | C++ multithreaded kernel |
| Single-cell (50k cells) | < 30 minutes end-to-end | With pseudo-bulk aggregation |
| Memory (10k genes) | < 4 GB peak RAM | Dense matrix in C++ |
| Report generation | < 60 seconds | HTML output |

### 5.2 Compatibility

- Python >= 3.9; tested on 3.9, 3.10, 3.11, 3.12, 3.13 and 3.14, and on 3.9 with the oldest supported dependencies
- C++ compiler supporting **C++17** (GCC >= 7, Clang >= 5, MSVC >= 2017; built and tested with GCC 15 so far)
- Supported platforms (project goal): Linux x86-64 (AMD64), Windows x86-64 (AMD64) and macOS on Apple Silicon (arm64); built from source and tested on all three by the CI (2026-11-18: Windows Server 2022 and macOS 15 with Python 3.12, Linux with Python 3.9-3.14 and the oldest dependencies). Other platforms are not targets.
- Distribution (project goal): source code (GitHub, https://github.com/wen-chen/GCEN2), PyPI (pre-built wheels for the three platforms; `pip install gcen2`) and Bioconda (Linux and macOS; Bioconda does not build for Windows; `conda install -c conda-forge -c bioconda gcen2`). Website: https://www.biochen.org/gcen2/

### 5.3 Dependencies

| Package | Purpose | Version |
|---|---|---|
| numpy | Numerical computation | >= 1.24 |
| pandas | Data manipulation and I/O | >= 2.0 |
| scipy | Statistical functions | >= 1.10 |
| threadpoolctl | Limit BLAS threads | >= 3.0 |
| typer | Command-line interface | >= 0.9 |
| pybind11 | C++/Python binding | >= 2.10 |
| scikit-build-core | Build system | >= 0.5 |
| anndata | Single-cell data container (optional) | >= 0.10 |
| scanpy | Single-cell preprocessing (optional) | >= 1.9 |
| matplotlib | Static plotting (optional) | >= 3.7 |
| seaborn | Statistical visualization (optional) | >= 0.12 |
| pyvis | Interactive network visualization (optional) | >= 0.3 |
| networkx | Network analysis (optional) | >= 3.1 |
| Jinja2 | Report template engine (optional) | >= 3.1 |
| weasyprint | PDF report export (optional) | any |

The lower bounds are the oldest versions the test suite was run with.

### 5.4 Testing

- Every C++ kernel compared with its Python reference implementation
- A test for every command, every `gcen2 plot` command and every public function; a new command fails the suite until it has one
- Golden-output regression tests (`tests/golden/`): the baseline for refactoring and performance work
- Reference tests against GCEN v1, WGCNA, hdWGCNA and RandomWalkRestartMH results; integration tests on the example data
- Benchmarks against GCEN v1, WGCNA, hdWGCNA, DGCA and RandomWalkRestartMH (`benchmarks/`)
- Continuous integration on Linux (Python 3.9-3.14, oldest dependencies), macOS and Windows

How the tests are organised: Implementation Document §4; results: `docs/test_report.md`.

### 5.5 Documentation

- API reference generated by Sphinx with autodoc ✅
- CLI reference generated from the commands ✅
- Quickstart and tutorials (bulk, single cell, spatial, trajectories and multi-omics, figures) on the example data ✅
- Developer guide (`CONTRIBUTING.md`) ✅
- Migration guide from GCEN v1 to GCEN2 ✅
- `CHANGELOG.md`, `CITATION.cff` ✅

### 5.6 Licensing

GCEN2 is released under the GNU General Public License, version 3 or later (GPL-3.0-or-later), consistent with GCEN v1 (`LICENSE`). The example data keep their own licences. Repository: `github.com/wen-chen/GCEN2`.

---

## 6. Status

| Requirements | Status |
|---|---|
| 3.1 Core pipeline (GCEN v1) | ✅ same results as GCEN v1 0.6.3 |
| 3.2 Single cell; 3.10 WGCNA / hdWGCNA modules | ✅ same results as WGCNA 1.74 and hdWGCNA 0.4.12 |
| 3.3 Visualization; 3.4 Reports | ✅ |
| 3.5 Differential; 3.6 Trajectory; 3.7 Multi-omics | ✅ same results as DGCA 1.0.3, RandomWalkRestartMH and R implementations of fusion and trajectory networks |
| 3.8 Spatial | ✅ same results as hdWGCNA 0.4.12 |
| 3.9 LLM-assisted interpretation (P2) | 🔲 not started |
| 3.11 Example data, reproducibility | ✅ |
| 5.1 Performance | ✅ benchmarks re-timed on an idle machine after the feature freeze (2026-11-17) |
| 5.2 Compatibility | Linux (Python 3.9-3.14, oldest dependencies), Windows and macOS (Apple Silicon) built and tested by the CI; compiler minimums not verified |

What is left: `docs/status.md`; how the work went: `docs/development_log.md`.

---

## 7. Known Issues from GCEN v1 (Status)

| Issue | Status | Resolution |
|---|---|---|
| `network_merge.cpp` appends string "1"/"-1" to numeric correlation values | ✅ Fixed | `network.merge()` outputs direction as a separate column |
| `deviation()` computes population std (ddof=0) instead of sample std (ddof=1) | ✅ Fixed | `filter_genes()` uses `ddof=1` throughout |
| Naive CSV parsing (character replacement, no quoted field support) | ✅ Fixed | no `csv_to_tsv` / `tsv_to_csv` commands; `gcen2.io.read_csv` (pandas) reads quoted fields, `write_tsv` / `write_csv` convert |
| Inconsistent error handling (`network_extract`, `network_shuffle` no exit on failure) | ✅ Fixed | Unified `GCENIOError` raised on all file open failures |
| C++ extension not importable after editable install | ✅ Fixed | scikit-build-core's redirect editable mode (`pip install --no-build-isolation -e .`) |

---

*End of Document*
