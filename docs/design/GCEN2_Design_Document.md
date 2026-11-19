# GCEN2

**A Python Toolkit for Gene Co-Expression Network Analysis with Single-Cell Resolution and Interactive Visualization**

---

**Software Design Document**

| | |
|---|---|
| Document Version | 2.2 |
| Date | November 2026 (2.1: November 2026; 2.0: March 2026) |
| Author | Wen Harold Chen |
| Affiliation | Hunan University of Chinese Medicine |
| Contact | chenwen@biochen.org |

---

## Table of Contents

- [1. Introduction](#1-introduction)
- [2. System Architecture](#2-system-architecture)
- [3. Data Models](#3-data-models)
- [4. Module Design](#4-module-design)
- [5. Algorithm Details](#5-algorithm-details)
- [6. Build, Tests and Migration](#6-build-tests-and-migration)

---

## 1. Introduction

### 1.1 Purpose

This Software Design Document (SDD) describes the technical architecture, data models, module interfaces, algorithms, and implementation strategies for GCEN2. It serves as the primary reference for developers during implementation and translates the requirements defined in the GCEN2 Software Requirements Specification (SRS v2.2) into a concrete engineering blueprint. The build, the coding rules, the test system and the release are in the Implementation Document; the files of the repository in the Code Guide.

### 1.2 Scope

This document covers the internal design of all GCEN2 components: the C++ computation engine with pybind11 bindings, the Python application layer, the visualization subsystem, the single-cell analysis module, the report generation pipeline, and the CLI compatibility layer. External system integrations (LLM APIs, KEGG REST API) are described at the interface level.

### 1.3 Design Principles

- **Performance where it matters:** C++ handles only O(n²) or higher computational kernels; everything else is Python for maintainability.
- **GIL safety at boundaries:** pybind11 bindings release the GIL only after all Python buffer/array access is complete. `py::call_guard<py::gil_scoped_release>` must NOT be used on lambdas that access `py::buffer_info` — use scoped `py::gil_scoped_release release` inside the lambda body instead.
- **Results independent of the machine:** no floating-point contraction, fixed summation orders, per-thread results merged and sorted, explicit tie-breaks; the same output with any number of threads and on any CPU.
- **C++17:** the kernels use C++17 (`std::string_view`, `std::to_chars` where the library has it, class template argument deduction).
- **Scanpy-native:** AnnData is the canonical single-cell data container.
- **Fail fast, fail loud:** Input validation at the Python boundary; C++ kernels assume validated data and throw `std::runtime_error` for unexpected conditions.
- **Reproducibility:** All random operations accept a seed parameter; default pipelines are deterministic.

### 1.4 Reference Documents

- GCEN2 Software Requirements Specification (SRS), v2.2, November 2026; Implementation Document, v3.0; Code Guide
- Chen et al. *GCEN: An Easy-to-Use Toolkit for Gene Co-Expression Network Analysis and lncRNAs Annotation.* Curr Issues Mol Biol 2022; 44(4): 1479-1487. doi:10.3390/cimb44040100
- Langfelder & Horvath. *WGCNA: an R package for weighted correlation network analysis.* BMC Bioinformatics 2008; 9: 559. doi:10.1186/1471-2105-9-559
- Morabito et al. *hdWGCNA identifies co-expression networks in high-dimensional transcriptomics data.* Cell Rep Methods 2023; 3(6): 100498. doi:10.1016/j.crmeth.2023.100498
- GCEN v1 source code: `github.com/wen-chen/GCEN`

---

## 2. System Architecture

### 2.1 Layered Architecture

GCEN2 is organized into four distinct layers. Each layer communicates only with its immediate neighbors:

| Layer | Name | Technology | Responsibility |
|---|---|---|---|
| L4 | User Interface | CLI (typer), Jupyter | Command parsing, parameter validation, output formatting |
| L3 | Application Logic | Python 3.9+ | Workflow orchestration, data transformation, I/O |
| L2 | Computation Engine | C++17 + pybind11 | Correlation, enrichment, module detection, RWR |
| L1 | Data Layer | pandas, AnnData, h5py | File I/O, format conversion, sparse matrix management |

### 2.2 Package Dependency Graph

Internal module dependencies follow a strict acyclic graph (dependencies point downward only):

```
gcen2.cli → gcen2.report → gcen2.plot (figures; gcen2.plotting: earlier API)
   │
   ├── gcen2.single_cell (metacells, celltype networks / modules, DMEs, projection)
   │        ├── gcen2.consensus ─┐
   │        └── gcen2.harmony    │
   ├── gcen2.spatial ────────────┤
   ├── gcen2.preservation ── gcen2._rrng (R's RNG)
   │                             ↓
   │                     gcen2.tom_modules (WGCNA modules)
   ├── gcen2.trajectory, gcen2.multilayer, gcen2.differential
   │        ↓
   ├── gcen2.annotation, gcen2.rwr
   │        ↓
   └── gcen2.module ← gcen2.network ← gcen2.preprocessing
                ↓
          gcen2.core (_gcen_core C++ extension)
                       ├── .correlation  (pairwise / GEMM kernels, network files)
                       ├── .enrichment
                       ├── .module       (Jaccard modules, TOM, hierarchical clustering)
                       ├── .network      (edge files: stats, merge, extract, shuffle, diff, fusion)
                       └── .rwr
```

`gcen2.plot` depends on no analysis module: it reads result files or
objects, so figures can be redrawn without recomputing.

All modules also depend on `gcen2.models`, `gcen2.exceptions`, `gcen2.config`, and `gcen2.io` as cross-cutting utilities.

`import gcen2` is lazy (PEP 562 `__getattr__` in `gcen2/__init__.py`):
submodules and the top-level helpers (`read_tsv`, `GCENConfig`, ...) are
imported on first access, so a command only loads the libraries it uses.
`gcen2.models`, `gcen2.io`, `gcen2.network`, `gcen2.module`, `gcen2.rwr` and
`gcen2.annotation` import pandas and scipy inside the functions that need
them. With the compiled extension, `gcen2 build`, `module`, `rwr` and
`annotate` run without pandas or scipy: about 31 MB and ~1 s of imports
instead of 84 MB and ~5 s on WSL; `gcen2 build` on a plain matrix does not
import numpy either (`build_from_file`, about 22 MB before reading the data)
(`tests/test_io.py::test_build_path_does_not_import_pandas`,
`tests/test_annotation.py::TestPandasFreeAnnotate::test_cli_paths_do_not_import_pandas`).

Text inputs (OBO, associations, gene lists, module files, expression
matrices) are read with `_matrix_io.iter_lines` (binary, 1 MB buffer,
decode per line): iterating a text-mode file is ~4x slower on 9p file
systems such as WSL `/mnt/*` (0.94 s vs 0.26 s for the 34 MB GO OBO file).

### 2.3 Thread Safety Model

C++ computation kernels release the GIL inside the computation loop (after all Python object access is done) and re-acquire it before returning results. Thread-local result accumulation with mutex-protected merge is used to avoid data races. The thread count is given to each function (`threads`) and command (`-t`; default 2, matching GCEN v1).

**Critical rule:** Never use `py::call_guard<py::gil_scoped_release>()` on a function/lambda that accesses `py::buffer_info` or any Python-managed memory. Always extract raw pointers first, then release the GIL inside a scoped block.

### 2.4 Error Handling Strategy

Input validation occurs at the Python boundary (L3). C++ kernels throw `std::runtime_error`; pybind11 translates these to Python `RuntimeError`. A unified exception hierarchy:

```python
class GCENError(Exception): pass
class GCENInputError(GCENError): pass   # bad input file or parameter
class GCENComputeError(GCENError): pass # numerical failure in C++
class GCENIOError(GCENError): pass      # file read/write failure
```

---

## 3. Data Models

### 3.1 Core Data Structures

All data containers are defined in `gcen2/models.py` as Python dataclasses wrapping standard objects.

#### 3.1.1 ExpressionMatrix

| Attribute | Type | Description |
|---|---|---|
| `data` | `pd.DataFrame \| AnnData` | Underlying expression data |
| `gene_names` | `np.ndarray[str]` | Gene identifiers |
| `sample_names` | `np.ndarray[str]` | Sample/cell identifiers |
| `is_single_cell` | `bool` | AnnData-backed flag |
| `is_normalized` | `bool` | Set after normalization |
| `normalization_method` | `str \| None` | Method used |

Note: Most Python API functions accept `pd.DataFrame | AnnData` directly (via `_unwrap()` helper) without requiring an ExpressionMatrix wrapper, for scanpy-style ergonomics.

#### 3.1.2 Network

| Attribute | Type | Description |
|---|---|---|
| `edges` | `pd.DataFrame` | Columns: `gene_a`, `gene_b`, `correlation`, `p_value` (+ `fdr` if computed) |
| `node_count` | `int` | Unique genes |
| `edge_count` | `int` | Edges |
| `is_signed` | `bool` | Signed network: negative correlations are excluded (GCEN v1 `-s`) |
| `method` | `str` | `'pearson'` or `'spearman'` |
| `metadata` | `dict` | Build parameters, source file, etc. |

Key methods: `to_networkx()`, `to_adjacency(weight="correlation", absolute=False, genes=None)` → weighted symmetric `scipy.sparse.csr_array` (sign kept; `absolute=True` for |r|, `weight=None` for 0/1), `to_edgelist(path)` (v1 number formatting for standard edge columns), `subnetwork(genes)`, `stats()` → `NetworkStats`.

#### 3.1.3 NetworkStats

| Attribute | Type |
|---|---|
| `node_count`, `edge_count` | `int` |
| `density`, `average_degree`, `max_degree` | `float` |
| `connected_components`, `giant_component_size` | `int` |
| `heterogeneity`, `clustering_coefficient`, `transitivity` | `float` |

Computed via `NetworkStats.compute(network, threads=2)` classmethod. With the
compiled extension, degrees, triangles per node and connected components come
from `_gcen_core.network.graph_stats` on the adjacency of
`module._unique_adjacency` (nodes in sorted order); the scipy implementation
(`_compute_sparse`, triangles from row-chunked `A @ A`) is the fallback.
`NetworkStats.from_adjacency(indptr, indices)` computes the same from a CSR
adjacency, as `network.stats_from_file(path)` does from
`module.read_file_adjacency` without pandas or scipy (`gcen2 stats`).
`NetworkStats.with_degrees(indptr, indices)` also returns the node degrees
(networkx `G.degree` values: a self-loop adds 2), for
`network.summary_from_file`.

#### 3.1.4 ModuleSet

| Attribute | Type | Description |
|---|---|---|
| `modules` | `list[set[str]]` | Gene sets (each set = one module) |
| `module_count` | `int` | `len(modules)` |
| `source_network` | `Network \| None` | Reference network |
| `similarity_cutoff` | `float` | Jaccard cutoff used |

#### 3.1.5 EnrichmentResult

| Attribute | Type | Description |
|---|---|---|
| `results` | `pd.DataFrame` | Columns: `term_id`, `name`, `namespace`, `enrichment`, `study_count`, `study_n`, `pop_count`, `pop_n`, `p_value` |
| `ontology_type` | `str` | `'GO'` or `'KEGG'` |
| `gene_or_module` | `str` | Query identifier |
| `p_cutoff` | `float` | P-value cutoff |

#### 3.1.6 RWRResult

| Attribute | Type | Description |
|---|---|---|
| `ranked_genes` | `pd.DataFrame` | Columns: `gene`, `probability`; sorted descending |
| `seed_genes` | `list[str]` | Seed genes |
| `gamma` | `float` | Restart probability |
| `steps` | `int` | Iterations to convergence |

#### 3.1.7 DiffCoExResult

| Attribute | Type | Description |
|---|---|---|
| `results` | `pd.DataFrame` | Columns: `gene_a`, `gene_b`, `cor_a`, `cor_b`, `z_score`, `p_value`, `category` |
| `network_a_name`, `network_b_name` | `str` | Labels for the two networks |
| `method` | `str` | `'fisher_z'` |
| `p_cutoff` | `float` | P-value cutoff |

### 3.2 File Formats

| Format | Description | Read/Write | Module |
|---|---|---|---|
| TSV | Tab-separated (GCEN v1 native) | R/W | `gcen2.io` |
| CSV | Comma-separated (quoted fields supported) | R/W | `gcen2.io` |
| h5ad | AnnData HDF5 for single-cell | R/W | `gcen2.io` |
| OBO | Gene Ontology OBO format | R | `gcen2.annotation` |
| HTML | Interactive network or research report | W | `gcen2.plotting`, `gcen2.report` |
| PDF | Research report (via weasyprint) | W | `gcen2.report` |

### 3.3 Configuration

`GCENConfig` dataclass in `gcen2/config.py` (public, kept for compatibility). Module-level `default_config` instance. The library functions and commands do not read it: threads and random seeds are arguments of each function (`threads`, `seed`) and options of each command (`-t`, `--seed`).

| Parameter | Default | Env Var | Description |
|---|---|---|---|
| `n_threads` | 2 | `GCEN2_THREADS` | Default C++ thread count |
| `random_seed` | 42 | — | Reproducibility seed |
| `verbosity` | 0 | — | 0=silent, 1=progress, 2=debug |
| `temp_dir` | None | — | Temp directory (auto if None) |
| `figure_dpi` | 300 | — | DPI for static figures |
| `figure_format` | `'pdf'` | — | `pdf`, `png`, or `svg` |

---

## 4. Module Design

### 4.1 gcen2.core — C++ Computation Engine

The C++ source is compiled into a single shared library `_gcen_core` (a pybind11 extension module). It exposes four submodules, each registered in `src/bindings/module.cpp` and implemented in `src/bindings/*.cpp`.

#### 4.1.1 `_gcen_core.correlation`

**`pairwise_correlation(matrix, method='spearman', n_threads=2, min_abs_cor=0, max_p=inf, positive_only=False, row_begin=0, row_end=-1) → (indices, correlations, p_values)`**

- `matrix`: `np.ndarray` shape `(n_genes, n_samples)`, `float64` (any layout; forced to C order)
- `min_abs_cor` / `max_p` / `positive_only`: in-kernel filters, so memory scales with retained edges
- `row_begin` / `row_end`: only pairs `(i, j)` with `row_begin <= i < row_end`, `j > i` (used for block streaming)
- Returns: `(indices: int32 array (n_edges, 2) sorted by (i, j), correlations: float64 array, p_values: float64 array)`
- GIL: released via `py::gil_scoped_release release` inside the lambda body, after `buf.request()` and pointer extraction
- Threading: work units of up to 256 pairs (a row, or a column range of a long row) are handed out dynamically, so even a block of a few long rows keeps every thread busy; `build_network_file` runs all blocks on one `WorkerPool` (threads started once)

**Implementation note:** `py::call_guard<py::gil_scoped_release>()` must NOT be used here — it caused a segfault because it released the GIL before `matrix.request()` was called.

#### 4.1.2 `_gcen_core.enrichment`

**`hypergeometric_pvalue(pop_n, pop_count, study_n, study_count) → float`**

**`batch_hypergeometric(pop_n, pop_counts, study_n, study_counts) → np.ndarray[float64]`**

- Log-factorial table precomputed once per batch call (size: `pop_n + 1`)
- All results clamped to `[0.0, 1.0]` to handle floating-point precision accumulation
- `pop_counts`, `study_counts`: `np.ndarray` dtype `int32`
- Two-sided Fisher test as in GCEN v1; kept for reference (not used by `gcen2.annotation`)

**`hypergeom_sf(k, M, n, N, n_threads=1) → np.ndarray[float64]`**

One-sided upper tail `P(X >= k)`, `X ~ Hypergeom(M, n, N)`, elementwise over
equal-length float arrays holding integers; multithreaded, GIL released. Same
algorithm as the numpy fallback `annotation._hypergeom_sf_numpy`: each tail is
summed away from the mode with the pmf ratio recurrence (log-factorial table
from `lgamma`), so results agree with it and with `scipy.stats.hypergeom.sf`
to ~1e-10 relative.

**`format_annotation_rows(term_prefix, terms, enrichment, study_count, study_n, pop_count, pop_n, p_value, bounds, n_threads=1) → list[bytes]`**

Formats annotation result rows (one `bytes` per query group
`bounds[g]:bounds[g+1]`) byte-identically to `DataFrame.to_csv(sep="\t")`:
`term_prefix[term]` holds the pre-quoted `term_id name namespace` fields;
floats are written as Python `repr()` (`std::to_chars` shortest digits, or a
round-trip `snprintf` search where `to_chars` for `double` is unavailable).
`format_py_float(x)` exposes the float formatting for tests.

#### 4.1.3 `_gcen_core.module`

**`jaccard_modules(indptr, indices, similarity_cutoff, n_threads) → labels`**

- Symmetric CSR adjacency (int64 `indptr`, int32 `indices`; a self-loop is a diagonal entry), from `module._unique_adjacency` or from a file by `read_edge_adjacency(path, skip_header, with_counts, descending)` (C++ reader, nodes sorted by name)
- Jaccard similarity of the linked pairs in parallel (per-thread union-find), then the connected components of the pairs above the cutoff; one int32 label per node (-1: in no module)

**Jaccard formula:** `J(i,j) = (|N(i)∩N(j)| + a(i,j)) / (|N(i)| + |N(j)| - |N(i)∩N(j)| - a(i,j))` where `a(i,j) = 1` if edge `(i,j)` exists.

**Important:** For a 2-node graph with a single edge, `J = 1.0` (the edge itself contributes `a_b = 1`). This is correct algorithm behavior — both nodes form a module regardless of the cutoff (as long as cutoff < 1.0).

#### 4.1.4 `_gcen_core.rwr`

**`random_walk_restart(indptr, indices, data, P0, gamma=0.5, max_iter=1000000, tol=1e-10) → (PT, steps)`**

- `indptr` (int64), `indices` (int32), `data` (float64): CSR of the column-normalised adjacency `W`, shape `(n, n)`
- `P0`: seed vector, `np.ndarray` shape `(n,)`, `float64`
- Iterates `PT ← (1-γ) W PT + γ P0` until the L1 change is below `tol` (GCEN v1 rule); GIL released

**`random_walk_restart_counts(indptr, indices, counts, P0, gamma=0.5, max_iter=1000000, tol=1e-10) → (PT, steps)`**

- Same walk on the unnormalised adjacency with int32 entry multiplicities `counts`; `W[i, j] = (1 / col_sum[j]) * counts` is formed inside the iteration (count-1 entries use a per-step `inv[j] * PT[j]`), so no float64 copy of the matrix is stored and the result is bit-identical to `random_walk_restart` with `data = (1.0 / col_sums)[indices] * counts`

#### 4.1.5 `_gcen_core.network`

**`graph_stats(indptr, indices, n_threads=2) → dict`**

- Input: symmetric CSR adjacency, rows sorted and unique, self-loop = diagonal entry
- Returns `degree_no_loop`, `self_loop`, `triangles` (int64 arrays, zero-copy) and `edge_count`, `components`, `giant_component`
- Triangles: each edge is oriented from the lower to the higher (degree, id) rank, so every triangle is found once from its lowest-ranked node; per node the out-neighbours are marked and scanned in the out-lists of those neighbours. Rows are handed out dynamically, one count array per thread. Components by iterative DFS.
- Full 33k-gene sample (3.6M edges): 2.4 s at 4 threads; `gcen2 stats` 843 MB / 51 s → 84 MB / 6.2 s, identical values.

**`EdgeSwapper(a, b, base)`**, `.run(i, j) → applied`, `.shuffle(state, n_swaps) → (applied, state)`, `.b()`

- Degree-preserving swaps for `network.shuffle`: `run` takes candidate edge pairs `(i[t], j[t])` drawn by numpy; `shuffle` draws them itself from numpy's stream reproduced in C++ (`Pcg64`, below), per batch of 2^20 as `rng.integers(0, n, m)` for `i` then for `j`. Either way a seed gives the same network as the pure-Python loop
- Edges stored as `{int32 a, int32 b}` records (one cache line per lookup; 3.6M edges = 29 MB, mostly in L3)
- Edge keys `min * base + max` in an open-addressing set (linear probing, backward-shift deletion, load factor ≤ 1/4), `uint32` keys when `base ≤ 65535`, else `uint64`; the table asks for transparent huge pages (`huge_buffer.hpp`, `madvise(MADV_HUGEPAGE)`). The edge records of the candidates 16 ahead and the four key slots of the candidate 8 ahead (two tested, two erased if the swap is applied) are prefetched
- Acceptance depends only on key equality, so the layout changes give identical swaps. 3.6M-edge sample: 96 → 56 ns per attempt (unprefetched erase slots 96 → 72, edge records 72 → 57); a denser table (load ≤ 1/2) is slower (71–77 ns)
- 198k edges × 100: 20M swaps in 3 s (Python loop 40 s); `gcen2 shuffle` 47 s → 9.2 s on WSL, identical output

**`Pcg64`** (C++ only, `pcg64.hpp`; used by `EdgeSwapper.shuffle` and `shuffle_edge_file`)

- numpy's PCG64 (XSL-RR 128/64: `state = state * M + inc`, output `rotr(hi ^ lo, hi >> 58)`) started from `bit_generator.state`, so seeding (`SeedSequence`) stays in numpy; `next32` returns the low half of a 64-bit output and buffers the high half (`has_uint32` / `uinteger`)
- `integers(high, n)` = `Generator.integers(0, high, n)` (int64): `random_bounded_uint64_fill` with Lemire's rejection on 32-bit draws for `high - 1 ≤ 2^32 - 1`, on 64-bit draws above; the final state equals numpy's
- 128-bit product via `unsigned __int128` or, without it (MSVC), 32-bit limbs (checked against `__int128` on 20M products)
- `pcg64.hpp`; checked against numpy 2.5 for highs across the 32/64-bit boundaries, with and without a buffered half

**`shuffle_edge_file(path, skip_header, n_columns, v1_layout, header, output, swap_factor, state) → dict`**

- `network.shuffle` from a file to a file (`edge_table.cpp`): genes interned, the value fields of each row kept as text, `EdgeSwapper.shuffle`, then rows written as `Network.to_edgelist` would: the v1 layout (`gene_a gene_b correlation p_value [fdr]`) through `format_edge_lines` (`%f` / `%e`), any other layout as `DataFrame.to_csv` of the dtypes pandas infers per column (int64 → integer; float64, also int with missing values → `repr()`, NaN empty; object → the unquoted text, quoted when needed)
- Returns `ok = False` (with `reason`, nothing written) when pandas' reading is not certain: bool-like or whitespace-padded values, `inf` / `nan` spellings, numeric-looking tokens outside `[+-]digits[.digits][e±digits]`, more than 17 significant digits, rows longer than the header, missing gene names; the caller then uses the pandas pipeline

**`extract_edge_file(path, skip_header, n_columns, v1_layout, header, output, genes) → dict`**

- `network.extract` from a file to a file with the same reader and writer as `shuffle_edge_file`: only the rows whose first or second gene (unquoted) is in `genes` are kept and written
- The column dtypes and the checks above come from every row, as `read_network` reads the whole file (an int column with a missing value outside the extracted rows is still written as float); a missing gene name only matters in a kept row
- Replaces `filter_edge_lines` (C++ line filter, then pandas parsed the kept lines), whose dtypes came from the kept lines only and could differ from `read_network` in that case
- Fields are unquoted as pandas' C tokenizer reads them (`edge_text.hpp`, shared by all the file functions): a field starting with `"` is quoted, `""` inside it is a `"`, and after the closing quote the rest of the field is literal (`"q"t"` → `qt"`). Lines are split at tabs first, so a quoted field holding a tab or a line break is not supported

**`merge_edge_files(files, cor_cutoff) → (bytes, node_count, edge_count)`**

- `network.merge` on files: `files` holds (path, skip_header, correlation column) tuples
- Correlations of each ordered pair summed with the Kahan compensation of pandas' `group_sum` (missing values skipped, an infinite compensation reset to 0); pairs sorted by name; sums written as Python `repr()` (`py_repr.hpp`, shared with the enrichment formatter); gene names holding `"`, a tab or a line break quoted with `"` doubled, as `to_csv` writes them
- pandas ≥ 3.1 sums with `np.add.reduceat` instead when the input rows are already in pair order with no missing value; that can differ in the last bit when a pair has three or more values. Real merges are not in pair order, so this does not come up there.
- test_1 + test_2 (550k edges): `gcen2 merge` 240 MB / 1.0 s → 135 MB / 0.45 s on native Linux (WSL: 192 MB / 6.7 s → 111 MB / 2.7 s), byte-identical output

**`diff_edge_files(a, b) → DiffTable`**

- `differential.compare` on two files, each a (path, skip_header, correlation column) tuple. Rows become (gene id, gene id, correlation) records (an NA or missing correlation is NaN; a missing gene name is an error), re-keyed as `(min rank, max rank)` of the sorted names, stable-sorted with the last row of each edge kept, then merge-joined
- `DiffTable.cor_a()` / `cor_b()` / `flags()` (bit 0: in A, bit 1: in B) / `lo()` / `hi()` (indices into `names()`, the sorted gene names) are read-only numpy views of the C++ vectors; the Fisher z test of the shared edges runs in Python (`_fisher_z_arrays`, the same numpy/scipy calls on the same arrays as `compare`, hence the same bits)
- `DiffTable.write(path, z, p, pval)` writes the table as `DataFrame.to_csv(sep="\t")` does (floats as `repr()`, NaN empty, names quoted when needed) in 4 MB slices and returns the category counts

### 4.2 gcen2.preprocessing

**`normalize(expr, method='tmm', housekeeping_genes=None) → DataFrame | ExpressionMatrix`**

Routes to one of five implementations:

| Method | Algorithm |
|---|---|
| `upqt` | Per-sample upper quartile of nonzero values; scale so mean upper quartile is preserved |
| `median` | Per-sample median of nonzero values; scale to mean median |
| `deseq` | Log-space geometric mean per gene; per-sample ratio; median of ratios as size factor |
| `tmm` | Weighted trimmed mean of M-values; reference sample = closest upper quartile to mean |
| `hkg` | Per-sample median of housekeeping genes; scale to mean housekeeping median |

**`filter_genes(expr, cutoff_mean=0.0, cutoff_sd=0.0, percent_mean=1.0, percent_sd=1.0) → DataFrame`**

Uses `ddof=1` (sample standard deviation), fixing GCEN v1's population std bug.

**`log_transform(expr, base='log2', pseudo=1.0) → DataFrame`**

Supports `log`, `log2`, `log10` with pseudo-count offset.

### 4.3 gcen2.network

**`build(expr, method='spearman', cor=0.1, pval=0.001, signed=False, fdr=False, log=None, all_edges=False, threads=2) → Network`**

1. Validate shape, numeric dtype, sample count ≥ 3
2. Apply log transform if `log` is not None
3. Convert to `float64` row-major numpy array
4. Call `_gcen_core.correlation.network_edges()` (`src/network_stream.cpp`, one
   block): correlations, `cor` / `pval` / sign filters and, if `fdr=True`, the
   Benjamini-Hochberg values with `n_total = n*(n-1)//2` (bit-identical to
   `_benjamini_hochberg`)
5. Return `Network`

**`build_to_file(expr, path, ..., max_pairs_per_block=500_000) → {"node_count", "edge_count"}`**

Same parameters and byte-identical output as `build(...).to_edgelist(path)`, but
rows of the pair matrix are processed in blocks of ≤ `max_pairs_per_block` pairs
and written immediately (peak memory: one block; 1.26 GB when built in memory),
all in C++ (`_gcen_core.correlation.build_network_array`, `stream_network` in
`src/network_stream.cpp`). With `fdr=True` one pass keeps the p-values of all
tests with p ≤ `pval` (8 bytes each; ranked with a parallel sort, then BH as
numpy computes it: `p * m / rank`, cumulative minimum from the right, capped at
1, looked up by binary search) and the edges that pass the other filters (24
bytes each, up to 2^24; beyond that the blocks are correlated again to write
them). The kernels skip the p-value (incomplete beta) of pairs below the |r|
implied by `p ≤ pval` at the full sample count (`min_abs_cor_for_p`, with a
margin), which leaves the filtered pairs unchanged. GTEx whole blood, 11,944
genes x 803 samples, 8 threads: 11.6 s / 0.64 GB with FDR (numpy path: 33 s /
1.5 GB), 10.1 s / 0.19 GB without.
Lines are formatted in C++ by `_gcen_core.correlation.EdgeFormatter` with GCEN
v1 number formatting (correlation `%f`, p-value/FDR `%e`, 6 decimals); the
edge lines are byte-identical to v1 `network_build` output (all 12.6M lines of
`-a` on 5k genes, all 3.9M edges on the full sample).
`expr` may also be a plain genes x samples array with `gene_names=`, as read by
`gcen2._matrix_io.read_matrix` (numpy only, same header detection and values
as `io.read_tsv`, bit-identical on the sample data).

**`build_from_file(input_path, path, ..., header=None) → {"node_count", "edge_count"}`**

Used by `gcen2 build` (`gcen2/_build_file.py`, re-exported as
`network.build_from_file`); same file as `read_matrix` + `build_to_file`, but
without importing numpy (≈17 MB, 1–2 s of start-up on WSL). For plain matrix
files, `_gcen_core.correlation.build_network_file` reads the matrix in C++
(UTF-8, tab-separated, plain decimal values parsed with `std::from_chars`,
short rows padded with NaN, header rules of `read_matrix`), ranks it once
(`prepare_rows`; the numpy path re-ranks per block), computes each row block
with `correlate_rows`, applies the same final mask and writes with the shared
formatter (`src/edge_format.hpp`, also behind `EdgeFormatter`). Any other
input (NA strings, quotes or spaces around numbers, `inf`/`nan`, invalid
UTF-8, fewer than 3 samples, header/row width mismatch, an output that
cannot be opened) returns `None` before anything is written, and the numpy
path then handles or reports it exactly as before; `--log` always uses the
numpy path (FDR is handled by both, with the same `stream_network`). Full 33k-gene sample at 4 threads: 34 MB peak (numpy path
56 MB, with pandas 112 MB; v1 13 MB); 5k genes 29 MB. Lines are formatted and
written in slices of 256k edges on the same `WorkerPool` as the correlations
(no thread start-up per block). The C++ path therefore defaults to
`max_pairs_per_block=500_000` (`MAX_PAIRS`, now also the default of
`build_to_file`, which streams in C++ as well): `-a` on 5k genes (12.6M lines) peaks at 116 MB /
6.1 s (181 MB at 2M, numpy path 555 MB / 10.8 s; 86 MB at 200k), while the
full sample at `-c 0.9` takes the same time at 2 threads and +0.2 s (11.4 s)
at 4 threads. Output byte-identical
at all benchmark sizes (`tests/test_build_file.py` covers options, header and
line-ending variants, ragged rows and every fallback case).

**Correlation kernel (`kernel="auto"`; `gcen2 build --kernel auto|pairwise|gemm`, also `spatial-build`)**

Pearson / Spearman edges come from one of two kernels, chosen in
`make_correlate` (`src/network_stream.cpp`) for every entry point
(`network_edges`, `build_network_array`, `build_network_file`):
`"pairwise"` is GCEN v1's two-pass sum per gene pair (`correlate_rows`);
`"gemm"` standardises each row once (`prepare_gemm_rows`: z = (x − mean) /
sqrt(Σ(x − mean)²), zeros for a constant row, NaN for a non-finite one) and
takes r = z_i · z_j from `gemm_nt` tiles of 96 rows × 256 columns, handed
out dynamically to the `WorkerPool` (`correlate_rows_gemm`; r clipped to
[−1, 1], then the same cut-off, sign and p-value filter and (i, j) order).
`gemm_nt`'s fixed summation order makes the result independent of threads,
tiling and instruction set (AVX = baseline, `tests/test_gemm_kernel.py`).
`"auto"` uses `"gemm"` from `GEMM_MIN_SAMPLES` = 64 samples
(`kGemmMinSamples`), so small bulk data (and GCEN v1's sample data, 15
samples) keep v1's exact sums. r agrees with the per-pair kernel to 1e-15–
1e-14; an edge exactly on a cut-off can therefore differ. Speed (4,000
genes, 8 threads, random data): equal at 16 samples, 1.8× at 64, 3.8× at
128, 10–12× at 512–2,048. Measured: GTEx 11,944 genes × 803 samples FDR
build 13.3 → 2.7 s (same 533,232 edges; one p-value differs in the 7th
printed digit); DLPFC Visium `spatial-build` (9,317 genes × 3,611 spots)
39 → 3.3 s (identical files); PBMC metacell networks 1.3–1.7× faster
overall (same edges and modules). The extra memory is one standardised
copy of the matrix while it is prepared. The dropout-aware methods ignore
the kernel.

**TOM edge score (`tom=True`, `power="auto"`, `tom_min=0.0`; `gcen2 build --tom [--power N|auto] [--tom-min X]`)**

Optional `tom` column with WGCNA's topological overlap (`src/tom.cpp`,
`network.TOM_DOC`), for `pearson` / `spearman`. Adjacency
`((1 + r) / 2)^power` with `signed=True`, `|r|^power` otherwise (WGCNA
`TOMType = networkType`, `TOMDenom = "min"`); `TOM = (Σ_u a_iu a_uj + a_ij) /
(min(k_i, k_j) + 1 − a_ij)`. `power="auto"` reproduces WGCNA 1.74
`pickSoftThreshold` (default powers 1–10, 12–20; `scaleFreeFitIndex` with R's
`cut` / `tapply` / `hist` rules, plain R² > 0.85; otherwise the WGCNA FAQ power
for the number of samples, with a warning); `network.pick_soft_threshold()`
returns the table. Computation: r of all pairs by a double-precision
product of standardised rows (64-row tiles), stored as a float n×n matrix
(4 bytes × genes²: 0.57 GB at 11,944 genes) and turned into the adjacency
in place; `l_ij` from single-precision products of fixed 128×256 tiles,
computed only for tiles holding emitted edges (`stream_network` scores each
block before `tom_min` filters it). All products go through `gemm_nt`
(`src/tom_gemm.{hpp,cpp}`): a fixed summation order (4 interleaved lane sums
for double, 8 for float, each term a separately rounded product, then a
pairwise reduction) implemented by an AVX kernel (`tom_gemm_avx.cpp`,
compiled with `-mavx`, chosen at run time when the CPU has AVX) and a
baseline kernel (SSE2, also on MSVC x64; NEON on AArch64; plain C++
otherwise or with `-DGCEN2_GEMM_SCALAR`; `GCEN2_SIMD=baseline` forces it at
run time), all without FP contraction, so results are bit-identical across
instruction sets, compilers and tilings. Checked: GCC 15 and Clang 19 builds
(SSE2, AVX, `-march=native`, plain C++) and an AArch64 build (NEON and plain
C++, run under qemu) give the same output hashes (with FP contraction, as
AArch64 GCC does by default, they differ; all C++ is therefore compiled
with `-ffp-contract=off`, which also keeps the RWR kernels bit-identical to
numpy on AArch64 and with `GCEN2_NATIVE=ON`); the whole extension
cross-compiles for AArch64 through CMake and the test suite passes under
qemu (566 tests; those that start a subprocess cannot run there). CMake adds the AVX file only for an x86 target: `CMAKE_OSX_ARCHITECTURES`
must then be `x86_64` alone (cibuildwheel's arm64 cross build on an Intel
Mac and universal2 builds get the baseline kernel), and the compiler must
accept `-mavx`. 6×2 register blocks; each pass keeps the 256-term slices of a
128×128 panel in cache (looping over the whole row length for every
register block made the kernel memory-bound beyond 4 threads). AVX is 2×
the SSE2 Eigen product it replaced (Eigen is no longer a dependency); the TOM of 11,944 genes takes 6.5 s
instead of 13 s on 8 threads. Tiles run on the `WorkerPool`, so values do
not depend on threads or block size. Genes with a constant or non-finite row
get adjacency 0. The power and fit table go to `Network.metadata["tom"]` /
the `"tom"` key of `build_to_file` / `build_from_file`. Against WGCNA:
pickSoftThreshold R² within 1e-8 and TOM within 7e-7 (relative) on
`tests/data/wgcna_tom`; on GTEx whole blood (5,000 genes) TOM within 9.9e-7 of
`blockwiseModules`' TOM, same powers on all six benchmark runs, same STRING
precision. 11,944 genes × 803 samples, 8 threads: `--tom` adds 4.9 s and
0.6 GB to the FDR build (11.6 s → 16.6 s, 1.24 GB).

**`_benjamini_hochberg(p_values, n_total)`**

Uses `n_total` (all possible gene pairs) as the denominator, not just the edges passing the initial filter. This matches GCEN v1 behavior and is critical for correct FDR estimation.

**`merge(networks, cor_cutoff=0.5) → Network`**

Per ordered pair `(gene_a, gene_b)` (a pair and its reverse are separate rows): sum of the correlations over all networks, kept when `|sum| >= cor_cutoff`, rows sorted by pair. The sign is a separate `direction` column (`1`, `-1`, or `0`), fixing the GCEN v1 bug where direction was appended to the numeric correlation value.

**`merge_files(paths, output, cor_cutoff=0.5) → dict`**

Same file as `merge([read_network(p) for p in paths]).to_edgelist(output)`, computed by `_gcen_core.network.merge_edge_files` without pandas. Falls back to that call without the extension, or when the gene columns are not the first two. Used by `gcen2 merge`.

**`extract(network, genes) → Network`** / **`extract_file(path, genes, output) → dict`**

Edges with at least one end in `genes`. `extract_file` writes the same file as `extract(read_network(path), genes).to_edgelist(output)`, computed by `_gcen_core.network.extract_edge_file` without pandas; it falls back to that pipeline without the extension, when the gene columns are not the first two, or when the extension declines the file. Used by `gcen2 extract`, 3.6M-edge sample on native Linux: 20 genes 527 MB / 1.5 s (pandas) → 64 MB / 0.55 s; 2,000 genes (574k edges kept) 502 MB / 1.7 s → 147 MB / 0.55 s. Most of the time is reading the 195 MB file.

`read_network` reads the two gene columns as strings. Numeric gene IDs (e.g. Entrez) therefore keep their text and match string gene lists. Before, pandas read them as integers, and `extract` found no edges for them.

**`shuffle(network, swap_factor=100, seed=42) → Network`**

Degree-preserving edge swap randomization. `swap_factor` × edge_count swap
attempts, with candidates drawn from `numpy.random.default_rng(seed)` in batches
of 2^20: by `_gcen_core.network.EdgeSwapper.shuffle` with the stream
reproduced in C++, or by numpy for the pure-Python loop without the
extension; same result.

**`shuffle_file(path, output, swap_factor=100, seed=42) → dict`**

Same file as `shuffle(read_network(path), swap_factor, seed).to_edgelist(output)`, computed by `_gcen_core.network.shuffle_edge_file` without pandas. Falls back to that pipeline without the extension, when the gene columns are not the first two (or names repeat), or when the extension declines the file. Used by `gcen2 shuffle`, native Linux, byte-identical output:

| Input, swap factor | before | after |
|---|---|---|
| 3.6M-edge sample, 100 | 38.6 s / 1.44 GB | 22.7 s / 416 MB |
| 3.6M-edge sample, 5 | 5.0 s / 1.44 GB | 2.5 s / 414 MB |
| test_1 (198k edges, 6 columns), 100 | 1.4 s / 170 MB | 1.2 s / 99 MB |

**`stats(network) → NetworkStats`**

Delegates to `NetworkStats.compute(network)`.

**`stats_from_file(path, threads=2) → NetworkStats`**

Same values as `stats(read_network(path))`, read by `module.read_file_adjacency`
without pandas or scipy. Used by `gcen2 stats` (`-t` threads).

**`summary_from_file(path, threads=2) → NetworkSummary`**

`NetworkSummary(stats, degrees, method, is_signed)`: what `report.generate` needs of a network, without its edge table. Statistics and degrees come from one `read_file_adjacency` + `graph_stats` pass; `method` / `is_signed` are the `Network` defaults, as for a network read by `read_network`. Without the extension it returns `read_network(path)`. Used by `gcen2 report`: 3.6M-edge sample 13–16 s / 1.71 GB → 2.5 s / 112 MB on native Linux, byte-identical HTML.

### 4.4 gcen2.module

**`identify(network, similarity=0.5, threads=2) → ModuleSet`**

1. Validate `0 < similarity < 1`
2. Build the symmetric 0/1 CSR adjacency (`_unique_adjacency`): gene codes from
   categorical columns (`read_network(..., genes_only=True)` in the CLI), int32
   indices, bool data, duplicates merged by `tocsr()`
3. Call `_gcen_core.module.jaccard_modules(indptr, indices, cutoff, n_threads)`
   (scipy `_jaccard_modules` fallback if C++ unavailable)
4. Filter out singletons (sets with < 2 genes)
5. Return `ModuleSet`

**`modules_from_file(path, similarity=0.5, threads=2) → list[set[str]]`**

Same modules and order as `identify(read_network(path)).modules`, used by
`gcen2 module`. `_gcen_core.module.read_edge_adjacency(path, skip_header)`
reads the first two columns in C++ (4 MB chunks, hash map of names, names
sorted in byte order = pandas category order) and returns the deduplicated
CSR directly; the header is detected in Python
(`_matrix_io.network_has_header`, same keys as `read_network`). No pandas or
scipy is imported (the CLI also writes the module file itself), which saves
~4 s of start-up on WSL. Falls back to `identify` without the extension.

Memory: the reader keeps edge endpoints as uint16 node ids while there are
fewer than 65,536 nodes (int32 after that) and hands its vectors to numpy
without copying; `jaccard_modules` unites linked pairs in a per-thread
union-find (n ints each) merged at the end, instead of storing every linked
pair (the components, and so the labels, are the same).

Full 33k-gene sample (3.9M edges): 9.3 s at 4 threads, 85 MB peak (was
101 MB; with pandas and scipy: 14.8 s / 271 MB; v1: 300 s / 749 MB). 5k
genes: 1.1 s (was 5.3 s).

### 4.4.1 gcen2.tom_modules (WGCNA modules)

**`tom_modules(expr, method="pearson", signed=False, power="auto", min_size=20, deep_split=2, merge_cut=0.15, detect_cut=0.995, log=None, threads=2, sample_names=None) → TomModules`**

WGCNA 1.74 `blockwiseModules` on one block (all genes in one TOM), used by
`gcen2 module -e expr.tsv|data.h5ad` (`--signed`, `--power`, `--min-size`,
`--deep-split`, `--merge-cut`; for .h5ad `--layer`, `--hvg`, `--groupby`,
`--group` as `build`); defaults are blockwiseModules' (TOMType =
networkType). Steps:

1. `_gcen_core.module.tom_dissimilarity`: 1 − TOM of all genes in **double**
   precision (`TomCore<double>` in `src/tom.cpp`, same tiles and formula as
   the float `tom` edge column). Average-linkage trees are sensitive to
   ~1e-8 in the dissimilarities: on GTEx (5,000 genes) 1 − TOM rounded to
   float32 (spacing 6e-8 near 1) turns near-ties into ties and moves 59 genes
   (ARI 0.98), while a TOM computed in or rounded to float32, with 1 − TOM
   then in double, gives the same modules
   (`benchmarks/wgcna_gtex/07_tom_precision.py`).
2. `hclust_average_inplace` (`src/hclust.cpp`): fastcluster's NN-chain
   average linkage (its search order, `<` tie rule, `s·a + t·b` update,
   stable sort by height, union-find relabelling to R's merge convention).
   Identical merges and bit-identical heights to `fastcluster::hclust` on the
   reference dissimilarities. The upper triangle of the n×n matrix is the
   working storage and is restored from the lower one, so no condensed copy.
3. `cutree_hybrid`: line-by-line port of `dynamicTreeCut::cutreeHybrid`
   1.63-1 (branch bookkeeping, core scatter / gap / split-height criteria,
   PAM stage respecting the dendrogram, size-ordered labels). Where R would
   stop (PAM nearest label 0) the objects stay unassigned.
4. `wgcna_post_process`: eigengenes (first principal component of the
   scaled module genes, from the smaller Gram matrix, aligned with the
   average), kME core check (`> 0.5` for `min_size/3` genes), `kME < 0.3`
   removal, reassignment by `corPvalueFisher` ratio `< 1e-6` (to the
   eigengene's column number, as WGCNA writes it).
5. `merge_close_modules`: eigengene 1 − r, average linkage, cut at
   `merge_cut`, iterated; labels 1..K by size.

**hdWGCNA preset** (`tom_modules(expr, **HDWGCNA)`, `--preset hdwgcna`):
hdWGCNA 0.4.12 `ConstructNetwork` is `blockwiseConsensusModules` on one set;
for one set the consensus and calibration steps are the identity, and the
differences from blockwiseModules are its defaults (deepSplit 4,
`pamStage = FALSE`, minKMEtoStay 0.2, minModuleSize 50, mergeCutHeight 0.2)
and its reassignment, which compares lower-tail `pchisq` log p-values so
that a better module never passes the threshold (`reassign_threshold = 0`).
Its power comes from `TestSoftPowers` (pickSoftThreshold with
`corFnc = bicor`, powers 1..10, 12..30, lowest power above 3 with
SFT.R.sq >= 0.8): `power_method="bicor"` (C++ `tom_rows("bicor")`: x −
median weighted by (1 − u²)², u = (x − median) / (9 MAD), MAD 0 falls back
to x − mean, correlated without re-centring; bicor within 3e-16 of WGCNA),
`power_r2`, `min_power`, `powers`. Same power and modules as hdWGCNA in all
24 PBMC runs on hdWGCNA's datExpr (`benchmarks/hdwgcna_pbmc/05_tom_modules.py`).

Downstream, as in the WGCNA tutorial (identical to WGCNA within 1e-14 on
the fixture, `make_downstream_reference.R`): `TomModules.eigengenes` /
`eigengene_frame()` (`moduleEigengenes` of the merged labels, ME0 for the
unassigned genes, `var_explained` = mean squared correlation of the
eigengene with its genes), `kme()` (`signedKME`, `corPvalueStudent`; NaN for
left-out genes) and `trait_cor(traits)` (`corAndPvalue(MEs, traits, use =
"p")`: pairwise complete samples matched by name, long table module, trait,
cor, p, n). CLI: `--me`, `--kme`, `--traits` + `--trait-cor` write them as
TSV with the modules numbered as in the `-o` file (`ME<k>` is module k,
WGCNA's ME0 is `ME_unassigned`).

Genes with a constant or non-finite row are left out (label 0). numpy's BLAS
is limited to `threads` (threadpoolctl, a dependency): with one thread per
logical CPU the small eigengene products were 20× slower on a 24-thread CPU.
Memory: one n×n double matrix: `tom_dissimilarity` builds the
correlations, the adjacency and then 1 − TOM in the output array
(`TomCore::dissimilarity_inplace`: bands of ~n/16 rows, whose 128×256 l
tiles go to a buffer before the band's upper triangle is overwritten; the
lower triangle keeps the adjacency until it is mirrored at the end), and
the clustering works in place. Against WGCNA: every step identical on
`tests/data/wgcna_tom` (`make_modules_reference.R`), and the final modules
are the same partition as `blockwiseModules` on all six GTEx whole-blood
runs (`benchmarks/wgcna_gtex/06_tom_modules.py`): 5,000 genes 1.7 s /
0.40 GB for the whole command (WGCNA 17 s / 1.3 GB), 11,944 genes 9.1 s /
1.5 GB (WGCNA 36 s / 3.9 GB).

Outputs for figures (module functions of a result `res`):
`soft_power_table(res)` (the fit of every power tried), `dendrogram_dict(res)`
(merge matrix, heights, unmerged and merged labels), `hub_tom(res, n_hubs=25)`
(TOM between each module's genes and its top-kME genes) and
`gene_significance(res, traits)`. CLI: `--soft-power`, `--dendrogram`,
`--tom-hubs`, `--gs`.

### 4.4.2 gcen2.consensus (WGCNA consensus modules)

**`consensus_tom_modules(sets, gene_names=None, set_names=None, signed=False, power="auto", consensus_quantile=0.0, calibration="single quantile", calibration_quantile=0.95, ...) → TomModules`** (CLI `gcen2 consensus-modules`; `celltype-modules --consensus-by`)

WGCNA 1.74 `blockwiseConsensusModules` on one block: the TOM of each set
with its own soft power; calibration ("single quantile": each set's TOM
raised to the power that maps its 0.95 quantile onto the first set's, the
quantile taken on entries drawn by R's `sample()` after `set.seed`, as
`consensusTOM`); per gene pair the `consensus_quantile` quantile over the
sets (WGCNA's `rowQuantileC`, C++ `row_quantile`); then average linkage,
dynamic tree cut and eigengene merging on the consensus dissimilarity, with
consensus eigengene dissimilarities. `HDWGCNA_CONSENSUS` holds hdWGCNA's
`ConstructNetwork(consensus = TRUE)` defaults (signed, quantile 0.3).
Found while validating: hdWGCNA pairs each set's power with the sets in
sorted order, which mismatches them unless the sets come sorted
(`benchmarks/hdwgcna_consensus/README.md`).

### 4.4.3 gcen2.preservation (WGCNA module preservation)

**`module_preservation(ref, test, modules, n_permutations=500, seed=12345, network_type="unsigned", max_gold_module_size=1000, max_module_size=1000, threads=2) → (Z table, observed statistics)`** (CLI `gcen2 module-preservation`)

WGCNA 1.74 `modulePreservation` for one reference and one test set, as
hdWGCNA's `ModulePreservation` calls it: a random "gold" module, modules
larger than `max_module_size` subsampled; quality and preservation
statistics (propVarExplained, meanSignAwareKME, separability,
meanSignAwareCorDat, meanAdj, meanMAR, cor.kIM, cor.kME, cor.cor,
cor.kMEall), Z-scores against `n_permutations` permutations of the module
labels, Zsummary and median rank. Every random step uses R's Mersenne
Twister and `sample()` (`gcen2._rrng.RRandom`, the same draws as R after
`set.seed(seed)`), so the gold module and the permutations are WGCNA's own:
the same statistics as WGCNA (relative difference < 7e-8), 6.8x faster on
the HLCA benchmark (`benchmarks/hdwgcna_preservation`).

### 4.5 gcen2.annotation

**`load_go(obo_file) → GOOntology`**

Parses OBO format. `GOOntology.propagate()` traverses the DAG upward, adding ancestor terms to each gene's annotation set. `GOTerm` preserves `alt_id` mappings and `is_obsolete` flags.

**`load_association(assoc_file, ontology=None) → dict[str, set[str]]`**

Supports single-column-per-term and multi-term-per-row formats. GO propagation applied automatically when ontology is provided.

**`annotate_go(network_or_modules, obo, assoc, pval=0.05, threads=2) → dict[str, EnrichmentResult]`**

For network-based: each unannotated gene is tested by enriching its neighbor set against background. For module-based: each module is tested as a unit (GCEN v1 rules, `build_annotation_queries`). Study counts for a chunk of queries come from one sparse product (queries × genes) @ (genes × terms) (`_iter_enrichment`); one-sided p-values from `_gcen_core.enrichment.hypergeom_sf()`.

**`annotate_to_dir(network_or_modules, out_dir, assoc, go=None, kegg=None, pval=0.05, threads=2, chunk_size=500) → int`**

Streaming version used by `gcen2 annotate`: writes `<gene>.go` / `module_<i>.go` (+ `module_<i>.txt`) chunk by chunk with `format_annotation_rows`, so memory does not grow with the number of results. Files are byte-identical to `to_csv` of the `annotate_go` / `annotate_kegg` results. The CLI passes a `NetworkAdjacency` (`NetworkAdjacency.from_file`, the C++ reader of §4.4); `build_annotation_queries` gives the same queries and background for it as for the `Network`. Inside `annotate_to_dir`, gene terms and queries are integer CSR arrays (`_GeneTerms`, `_Queries`) rather than dicts of string sets: GO propagation codes the ancestor closure of each listed term once and writes each gene's propagated term codes straight into the CSR, and `_Queries.from_adjacency` builds the network queries from the adjacency arrays. Counts, p-values and the cutoff are fused in `_gcen_core.enrichment.significant_counts` (per-thread dense term counter; p from the same routine as `hypergeom_sf`, so bit-identical), so only significant entries are materialised; `study_counts` + `hypergeom_sf` (or scipy without the extension) remain for other callers and for inputs with `study_n > pop_n`. `GOOntology.propagate` computes the ancestor closure once per distinct term. Full 33k-gene sample (3.9M edges, 3.15M result rows): 9.2 s / 234 MB at 4 threads (v1: 952 s / 284 MB), previously 10.0 s / 478 MB and before that 65 s / 1.65 GB; 5k genes 2.9 s / 161 MB.

**`enrich(query_genes, background_genes, ontology, assoc, pval=0.05) → EnrichmentResult`**

General-purpose enrichment for an arbitrary gene list.

**Python fallback:** without the C++ extension, p-values come from the vectorised numpy `_hypergeom_sf_numpy()` and files are written by pandas (same bytes).

### 4.6 gcen2.rwr

**`predict(network, seed_genes, gamma=0.5, prob_cutoff=0.01, directed=False, weighted_network=False, weighted_genes=False, gene_weights=None, max_iter=1_000_000, tol=1e-10) → RWRResult`**

1. Build adjacency matrix `W` (edge weights if `weighted_network=True`, else 1.0); column-normalize
2. Build seed vector `P0` (gene weights if `weighted_genes=True`, else uniform); normalize to sum=1
3. Call `_gcen_core.rwr.random_walk_restart(W, P0, gamma, max_iter, tol)` (or `_rwr_python()` fallback)
4. Filter by `prob_cutoff`, sort descending, wrap in `RWRResult`

**`predict_from_file(path, seed_genes, gamma=0.5, prob_cutoff=0.01) → (genes, probabilities)`**

Used by `gcen2 rwr` (undirected, unweighted); no pandas or scipy. The C++
reader returns the CSR with entry multiplicities (`with_counts`, parallel
edges summed as in `Network.sparse_adjacency`) and each row in descending
column order (`descending`), the order scipy's `W @ diags(1 / col_sums)`
leaves them in, so the floating-point sums and the output are bit-identical
to `predict()`. Ties are ordered as pandas `sort_values` (`_rank_desc`).
The reader hands its vectors to numpy without copying, and the walk uses
`random_walk_restart_counts`, so the matrix is never stored as float64.
Full sample: 102 MB peak (was 218 MB; 616 MB before `predict_from_file`).

### 4.7 gcen2.single_cell

**`pseudobulk(adata, groupby, layer=None, method='sum', min_cells=10) → pd.DataFrame`**

Returns genes×groups DataFrame. Groups with < `min_cells` cells are excluded.

**`build_by_celltype(adata, celltype_key, use_highly_variable=False, min_cells=10, strategy="metacell", k=25, max_shared=None, max_metacells=1000, use_rep=None, n_pcs=30, layer=None, normalize=None, pseudo_replicates=10, seed=42, effective="overlap", fdr=True, **build_kwargs) → dict[str, Network]`** (CLI `gcen2 celltype`)

For each cell type, a correlation network (`gcen2.network.build`) on its kNN
metacells (`strategy="metacell"`, default; p-values on every pair's
effective number of samples, below; `effective="cells"`: `n_cells // k`), or on `pseudo_replicates` random groups of its
cells (`strategy="pseudo_replicates"`; deprecated, warns: random groups rank
gene pairs at chance level). `fdr=True` by default (Benjamini–Hochberg over
all gene pairs, `pval` applying to it; `gcen2 celltype --no-fdr`), unlike
`network.build`: with calibrated p-values, `pval` alone keeps about
`pval` × 30–40 million pairs of a cell type as chance edges.

**`metacell_effective_samples(expr, members, k, method, log) → Series`** —
metacells sharing cells are correlated: a gene unrelated to gene y and to the
choice of metacells has metacell values with covariance
`C[a, b] = |cells(a) ∩ cells(b)| / k`, so its correlation with y has variance
`s(y) / M`, `s(y) = y'Cy / y'y` (y centred; ranks for Spearman). Every gene
gets `M / s`, clipped to [3, M]; a pair is tested on the smaller of its two
genes' values (`_effective_pvalues`, `t` with `n − 2` degrees of freedom;
FDR over all pairs as before). `C` is `O Oᵀ / k` with `O` the sparse
metacell × cell membership. With even, complete coverage a perfectly smooth
gene gets `n_cells / k`, the former fixed value; noisy genes get up to M.
Calibrated on PBMC null genes and in `tests/test_metacell.py`.

**`write_celltype_networks(networks, output_dir, diff=False, pval=0.05, threads=2)`** —
`--diff` tests the shared edges of every pair of cell types on the sample
sizes of the networks' own p-values: `gene_effective_samples` per edge
(smaller gene, in each cell type) when both networks have them, else
`effective_samples()` (`n_cells // k`; pairs below 4 skipped). Calibrated
on two random halves of each PBMC cell type (`11_diff_calibration.py`).

**`modules_by_celltype(adata, celltype_key, groups=None, fraction=0.05, use_highly_variable=False, k=25, max_shared=None, max_metacells=1000, min_cells=100, use_rep=None, n_pcs=30, layer=None, normalize=None, seed=42, threads=2, preset="hdwgcna", cell_level=True, reassign=True, **tom_kwargs) → dict[str, TomModules]`** (`gcen2/single_cell/modules.py`; CLI `gcen2 celltype-modules`)

The hdWGCNA workflow per cell type: genes with a value > 0 in at least
`round(fraction * n)` of the type's n cells (`fraction_genes`, R's round,
half to even; hdWGCNA `SetupForWGCNA(gene_select = "fraction")` per group), metacells of all types
chosen in one pass on all genes (`metacell._metacell_plan`: the members, one random
generator over the types, as one `metacells(groupby=celltype_key, groups=...)` call; so the
normalisation and neighbour space do not depend on the gene selection, as in
hdWGCNA, which normalises metacells before subsetting genes), then
`tom_modules` on the type's genes × metacells with the `preset` options
(`HDWGCNA` or blockwiseModules defaults) updated by `tom_kwargs`.
`TomModules.metadata` keeps `cell_type`, `n_cells`, `n_metacells`,
`n_genes` and the metacell matrix. `iter_modules_by_celltype` (same
arguments) yields `(cell type, TomModules)` one type at a time and builds a
type's metacell profiles (`metacell._metacell_profiles`, row by row, so the
same values as all at once) only when it is reached; `modules_by_celltype`
is `dict(iter_modules_by_celltype(...))`. The CLI writes and releases each
type before the next, so memory holds the expression plus one type's
metacells and TOM: on the 50,000-cell HLCA subset (34 types, 27,402 genes)
2.4 GB, against 9.3 GB when all types' metacells (dense, up to 1,000 per type × 27,402 genes)
and results were kept. `write_celltype_modules(results, out,
n_hubs=10, write_metacells=False)` writes `<type>/modules.tsv`,
`eigengenes.tsv`, `kme.tsv`, `hub_genes.tsv` (`TomModules.hub_genes`: top
kME per module, |kME| for unsigned networks, as hdWGCNA `GetHubGenes`),
optionally `metacells.tsv`, `cell_eigengenes.tsv` / `cell_kme.tsv`
(below), and `celltypes.tsv` (with `reassigned`). The result's own
eigengenes and kME (`eigengenes.tsv`, `kme.tsv`, hub genes) are those of the
metacells.

**Cell level** (`cell_level=True`, `--cell-eigengenes`, default), as
hdWGCNA's `ModuleEigengenes` and `ModuleConnectivity` on the type's cells:
`cell_eigengenes(X, gene_idx, labels, levels, normalize=True)` works through
blocks of rows (about 10^7 stored values each). Per module it fills a dense
cells × module-genes block with the LogNormalized values (x · 10^4 /
library, library over all genes; log1p) and computes the eigengene with
`spatial.module_eigengene_scores` (ScaleData clipped at 10, PC1 score
without re-centring, signed by the mean: the code of `spot_eigengenes`,
always in column-major layout because BLAS rounds differently by layout, so
the same bits as on spots). kME is the Pearson correlation of every gene
with every eigengene, accumulated over the row blocks (Σ x (y − ȳ), Σ x,
Σ x²), so an 80,000-cell type needs neither a dense nor a float64 copy of
its cells × genes matrix.

**Backed input** (`--backed`; `metacell._rows`): with an AnnData opened as
`read_h5ad(path, backed="r")` the expression stays on disk. Each cell type's
rows are read once (`X[cells]` of the CSR dataset) and shared by the gene
selection (`_positive_counts`, blocks of stored values), the metacell
profiles (members are re-indexed into those rows) and the cell level.
`_is_counts` reads the matrix in blocks of 20,000 rows. A CSC X is refused.
The outputs are identical to in-memory input. With `reassign=True` (default, ≥ 2 modules)
`spatial.reassign_negative_kme` applies `ReassignModules`; eigengenes are not
recomputed (nor in hdWGCNA). `metadata`: `cell_eigengenes` (cells ×
`ME<k>`), `cell_kme` (genes × `kME<k>`), `reassigned`. PBMC 10k: the same gene lists and
soft powers as hdWGCNA for all six benchmark types, modules as close to
hdWGCNA's as hdWGCNA across metacell seeds (ARI 0.71-0.79 for the four
large types); all eight types in 44 s / 2.0 GB.

**Consensus, Harmony, DMEs and traits** (`modules_by_celltype` arguments):
`consensus_by` (an `obs` column): the metacells of each value of the column
are a set and the modules are `consensus_tom_modules` of the sets;
`harmonize_by`: harmonized cell eigengenes, as hdWGCNA's
`ModuleEigengenes(group.by.vars = ...)` (per module, PCA of the scaled module
genes, `gcen2.harmony.harmony_integrate` of the `npcs` components, the
harmonized eigengene the first corrected component signed by the average
expression); `dme_by`: `single_cell.dme.find_all_dmes` of the cell
eigengenes per value of the column (hdWGCNA `FindAllDMEs`: Seurat
`FindMarkers` with Wilcoxon tests on the eigengenes clipped at 0);
`traits` / `trait_groupby`: `module_trait_cor` (hdWGCNA
`ModuleTraitCorrelation`). Written as `cell_harmonized_eigengenes.tsv`,
`dmes.tsv`, `trait_cor.tsv`, `sets.tsv` and `set_eigengenes.tsv`; every
cell type also gets `metacell_members.tsv` (metacell, cell, center).

**`project_modules(adata, modules, min_overlap=0.5, layer=None, cells=None, harmonize_by=None, npcs=50, seed=0, ...)`** (`single_cell/projection.py`; CLI `gcen2 project-modules`): hdWGCNA `ProjectModules`. Modules with at least `min_overlap` of their genes in the query keep those genes; their eigengenes on the query cells are computed as `cell_eigengenes` does (optionally harmonized). The projected eigengenes are within 2e-13 of hdWGCNA's.

**`build_dropout_aware(expr, method='spearman', zero_handling='pairwise_complete', **build_kwargs) → Network`**

- `pairwise_complete`: compute per-gene-pair correlation using only samples where both genes are nonzero
- `impute`: apply MAGIC before standard correlation (requires `magic-impute` package); warns that the p-values of imputed cells are not valid

### 4.8 gcen2.differential

**`compare(network_a, network_b, method='fisher_z', pval=0.05, n_a=None, n_b=None, fdr=False, correlation='pearson') → DiffCoExResult`**

For each undirected gene pair in the union of both networks (names ordered
within the pair; for a pair listed more than once the last row wins), apply
Fisher's Z-transformation test to pairs present in both (§5.6), and adjust
their p-values by Benjamini-Hochberg (`fdr` column). Categories: `A_only`,
`B_only` (z, p and FDR empty), `shared_changed` (p, or the FDR with
`fdr=True`, ≤ `pval`), `shared_stable`. `n_a` / `n_b` may be mappings gene →
effective samples (`_pair_samples`: each edge on the smaller of its two
genes' values; missing gene: `GCENInputError`; fewer than 4: z and p NaN,
`shared_stable`), as `write_celltype_networks` passes for metacell networks.
Without `n_a` / `n_b`, each network's metadata gives them (`_default_samples`:
`gene_effective_samples`, else `p_value_samples`, else `n_samples`, else 30),
so a metacell network is never tested on its number of metacells. Vectorised: genes become integer codes in sorted
name order, each pair the key `lo * n + hi`, and the union is a sorted merge
of the two key arrays, so rows come out sorted by `gene_a, gene_b` (the
row-by-row version iterated a set, in an order that changed between runs);
values are identical. 198k + 349k edges: 45 s / 465 MB → 11.7 s / 274 MB.

**`compare_files(path_a, path_b, output, pval=0.05, n_a=None, n_b=None, fdr=False, correlation='pearson', threads=2) → dict`**

Writes the same file as `compare(read_network(a), read_network(b), ...).results.to_csv(output, sep="\t", index=False)` and returns the count of each category present. The extension reads, joins and writes the edges (`_gcen_core.network.diff_edge_files`); only the Fisher z test of the shared edges runs in numpy/scipy (in chunks on `threads` threads), the FDR in C++ (`bh_adjust`), and pandas is not imported. Falls back to that pipeline without the extension or when the gene columns are not the first two. `n_a` / `n_b` default to 30 (`read_network` gives no `n_samples`); for mappings gene → samples the per-edge sizes come from `DiffTable.names()` / `lo()` / `hi()`. Used by `gcen2 diff`, native Linux, byte-identical output:

| A vs B | before | after |
|---|---|---|
| test_1 vs test_2 (548k edges) | 1.7 s / 336 MB | 0.76 s / 159 MB |

The files are read in 64 MB blocks of whole lines, each parsed by `threads`
threads with their own name tables (merged in file order, so the result does
not depend on the thread count); the edges are sorted by a stable parallel
sort and the table is formatted by `threads` threads in blocks written in
order. Two all-pairs networks of 5,000 genes (12.5M pairs each): 13.6 s with
1 thread, 5.5 s with 8, the same file.
| 3.6M-edge sample vs test_1 | 8.2 s / 1.15 GB | 1.9 s / 237 MB |
| 3.6M-edge sample vs itself | 11.9 s / 1.82 GB | 3.1 s / 561 MB |

**`_fisher_z_arrays(r_a, r_b, n_a, n_b) → (z_scores, p_values)`**

```
z_a = 0.5 * ln((1 + r_a) / (1 - r_a))
z_b = 0.5 * ln((1 + r_b) / (1 - r_b))
z_test = (z_a - z_b) / sqrt(1/(n_a-3) + 1/(n_b-3))
```

Correlations are clipped to ±0.9999 first; p = 2 · `norm.sf(|z_test|)`.
`_fisher_z_arrays` applies the same operations element-wise.

### 4.9 gcen2.plot and gcen2.plotting

#### 4.9.1 gcen2.plot (publication figures)

A package of 52 functions behind 43 `gcen2 plot` commands, one figure per
analysis step. Layout:

| File | Content |
|---|---|
| `_style.py` | The style (`rc_params`, `style()` context), widths 89 / 120 / 183 mm, module colours (`module_color`, `module_label`: WGCNA's `standardColors` by module number, or Okabe-Ito), `set_defaults(palette, formats)`, `headless_backend()` (Agg when no backend was chosen and the session is not interactive: matplotlib's automatic choice probes `DISPLAY` and can hang for minutes on a stale X server) |
| `_save.py` | `save_figure(fig, output, formats)`: a stem gives PDF, PNG (300 dpi) and SVG, a file name one format; no creation dates and fixed SVG ids, so the same data give the same files |
| `_io.py` | Readers for the result files (tables, modules, eigengenes, kME) |
| `_common.py` | `canvas`, `finish`, dendrogram drawing from R-convention merges, `repel_labels` (overlap-free labels kept as point offsets), `pack_components` (network components laid out one by one and packed), `plain_log_axis`, `draw_venn`, `draw_upset` |
| `wgcna.py`, `single_cell.py`, `network.py`, `enrichment.py`, `function.py`, `diff.py`, `trajectory.py`, `multiomics.py`, `qc.py` | The figures by analysis |

Conventions: every function takes result files (paths) or the objects
themselves, draws inside `style()`, returns the `Figure` and saves it when
`output` is given; text is 5-7 pt and stays text in PDF / SVG
(`pdf.fonttype 42`, `svg.fonttype none`); dense point layers are
rasterized inside the vector figure. The CLI declares the figure commands
(`gcen2/cli/plot.py`, `_plot`) and imports the plotting code only when a
command runs, so `gcen2 --help` and the analysis commands do not import
matplotlib.

#### 4.9.2 gcen2.plotting (earlier API)

Kept for compatibility. All functions follow: accept a domain object → return `matplotlib.figure.Figure` (for static plots) or write HTML file (for interactive). Optional `output` path triggers file save. In Jupyter, figures display inline.

| Function | Input | Output Type |
|---|---|---|
| `network_interactive(network, output, color_by, modules, annotation)` | Network | HTML file (pyvis) |
| `network_static(network, output, modules, figsize)` | Network | Figure |
| `dotplot(enrichment, top_n, output, title, figsize)` | EnrichmentResult | Figure |
| `barplot(enrichment, top_n, output, title, figsize)` | EnrichmentResult | Figure |
| `expression_boxplot(expr, output, figsize)` | DataFrame | Figure |
| `expression_violin(expr, output, figsize)` | DataFrame | Figure |
| `heatmap(expr, modules, output, figsize)` | DataFrame + ModuleSet | Figure |
| `degree_distribution(network, output, figsize)` | Network | Figure |
| `degree_histogram(degrees, output, figsize)` | node degrees (any order) | Figure (the same as `degree_distribution`) |
| `diff_network(diff, output, figsize)` | DiffCoExResult | Figure |
| `eigengene_heatmap(eigengenes, output, figsize)` | samples × eigengenes DataFrame | Figure (WGCNA `plotEigengeneNetworks` heatmap) |
| `module_trait_heatmap(trait_cor, output, figsize)` | long module-trait table | Figure (r and p per cell) |
| `spatial_eigengenes(coords, eigengenes, output, ncols, point_size, flip_y)` | spots × 2 coordinates + spots × eigengenes DataFrame | Figure: one spatial map per module, colour scale symmetric around 0 (99th percentile of \|ME\|) |

**Interactive network (pyvis):**
- Node size: log-scaled degree
- Node color: module membership (categorical) or degree (sequential)
- Edge width: absolute correlation coefficient
- Tooltip: gene name, degree, top 3 GO terms (if annotation provided)
- Physics: Barnes-Hut for > 500 nodes; force-atlas2 otherwise
- Output: one self-contained HTML file (`cdn_resources="in_line"`; pyvis's default wrote a `lib/` directory into the working directory)

### 4.10 gcen2.report

**`generate(network, modules, enrichment, rwr, expr, output, template) → None`**

1. Build context dict via `_build_context()`: assemble statistics and draw the figures with `gcen2.plot`, embedded as vector SVG with numbered captions (`_Figures`; `figures_dir` also saves each as PDF, PNG and SVG). `network` may be a `Network` or a `network.NetworkSummary` (statistics and degree histogram from it, no `to_networkx()`); both give the same report. Further sections: `celltypes` (`celltype_overview(dir, coords, basis)`: cell types and metacells on an embedding), `merge` and `extract` (figures of the inputs of `gcen2 merge` / `extract`), `moran` (Moran's I)
2. Render default or user-provided Jinja2 template
3. Write HTML. If `output.endswith('.pdf')`, convert via weasyprint

**WGCNA-style module sections** (`generate(..., wgcna=...)`): `wgcna_section(modules, eigengenes, kme, hubs, trait_cor, title, n_hubs)` reads the `gcen2 module -e` outputs (paths or DataFrames) into a section: per module its size and hub genes (from `--hubs`, else the top `n_hubs` by kME from `--kme`, else the first genes), the eigengene correlation heatmap (`plotting.eigengene_heatmap`: r between eigengenes ordered by average linkage on 1 − r, as WGCNA `plotEigengeneNetworks`) and the module-trait heatmap (`plotting.module_trait_heatmap`: r and Student p per cell, as the WGCNA tutorial's `labeledHeatmap`) with the ten strongest associations. `celltype_module_sections(dir)` gives one section per cell type of a `gcen2 celltype-modules` directory (cells, metacells, genes and power as a note). `ME_unassigned` is left out of both heatmaps.

Default template sections: Data Quality Control, Network Overview, Merged and Extracted Networks, Spatial Autocorrelation, Module Results, WGCNA-style Modules (cell-type overview first), Functional Annotation, RWR Predictions, Methods Summary. A missing template is a `GCENIOError`; PDF output without weasyprint, or without its system libraries, a `GCENInputError` with instructions.

### 4.11 gcen2.cli

Typer-based CLI (`gcen2.cli:app`, entry point `gcen2`), one module per group of commands in `gcen2/cli/` (`network`, `modules`, `single_cell`, `spatial`, `trajectory`, `multilayer`, `report`, `plot`; shared helpers in `_common`). There are no global options: commands that compute in parallel take `-t / --threads`, random steps `--seed`. `docs/cli.rst` (generated by `tools/make_cli_docs.py`) lists every option.

Every command function is wrapped (`wrap_command`, applied to all registered commands): errors of GCEN2 (`GCENError`) are printed as `ERROR: <message>` with exit code 1 and user warnings as `WARNING: <message>`, so the commands only raise; the parent directories of its outputs (the `output` option and options whose help starts with "write" / "also write") are created first, and list options left out get `None` (typer < 0.12 passes an empty tuple). Command annotations use `Optional[X]`, not `X | None`: typer evaluates them, and `X | None` fails on Python 3.9 (`tests/test_api_gaps.py` checks this).

| Subcommand | Python API | Key Options |
|---|---|---|
| `normalize` | `preprocessing.normalize()` | `-i`, `-o`, `-m` method, `-l` log |
| `filter` | `preprocessing.filter_genes()` | `-i`, `-o`, `-p` percent_mean, `--percent-sd`, `--cutoff-mean`, `--cutoff-sd` |
| `build` | `network.build_from_file()` | `-i`, `-o`, `-m`, `-c`, `-p`, `-f` fdr, `-s` signed, `-a` all edges, `-l` log, `-t`, `--max-pairs`, `--tom` / `--power` / `--tom-min` (TOM column), `--kernel auto\|pairwise\|gemm` |
| `module` | `module.identify()`; `tom_modules.tom_modules()` with `-e` | `-i`, `-o`, `-s` similarity, `-t`; `-e` expression matrix or .h5ad (WGCNA TOM modules) with `-m`, `--signed`, `--power`, `--min-size`, `--deep-split`, `--merge-cut`, `-l`, `--me`, `--kme`, `--hubs`/`--n-hubs`, `--traits`/`--trait-cor`, `--layer`/`--hvg`/`--groupby`/`--group`, `--preset wgcna\|hdwgcna` |
| `annotate` | `annotation.annotate_to_dir()` | `-g` obo, `-k` K2ko, `-a` assoc, `-n` network, `-m` modules, `-o`, `-p`, `-t` |
| `rwr` | `rwr.predict()` | `-n` network, `-g` genes, `-o`, `--gamma`, `--prob-cutoff` |
| `enrich` | `annotation.enrich()` | `-e` query, `-b` background, `-g`, `-a`, `-o`, `-p` |
| `stats` | `network.stats_from_file()` | `-i`, `-t` |
| `report` | `network.summary_from_file()` + `report.generate()` (`MPLBACKEND` defaults to `Agg`: no display probe, ~1.5 s over SSH X forwarding); `report.wgcna_section()` / `celltype_module_sections()` / `celltype_overview()` | `-n`, `-m` modules, `-o`; WGCNA sections: `--me`, `--kme`, `--hubs`, `--trait-cor`, `--soft-power`, `--dendrogram`, `--tom-hubs` (describe `-m`), `--celltype-modules DIR` with `--coords` / `--basis`, `--n-hubs`; `--merge-input` / `--merged`, `--extract-genes`, `--moran`; `--figures-dir` |
| `merge` | `network.merge_files()` | `-i` (multiple), `-o`, `-c` cor_cutoff |
| `extract` | `network.extract_file()` | `-i` network, `-g` genes, `-o` |
| `shuffle` | `network.shuffle_file()` | `-i`, `-o`, `--swap-factor`, `--seed` |
| `diff` | `differential.compare_files()` | `-A`, `-B` networks, `-o`, `-p`, `--n-a`, `--n-b`, `-m`, `-f`, `-t` |
| `spatial-build` | `spatial.build_spatial_filtered/weighted/regions()`; `spatial.spatial_autocorrelation()` with `--moran` | `-i` h5ad, `-o`, `-m` mode, `--n-neighbors`, `--radius`, `--bandwidth`, `--region-key`, `--n-regions`, `-c`, `-p`, `--method`, `--layer`, `--fraction`, `--normalize/--no-normalize`, `-t`, `--kernel`, `--moran` |
| `spatial-modules` | `spatial.spatial_modules()` + `write_spatial_modules()` (+ `plotting.spatial_eigengenes()` with `--plot`, `MPLBACKEND` defaulting to `Agg`) | `-i` h5ad, `-o` dir, `--groupby`, `--fraction`, `--hvg`, `--genes`, `--min-spots`, `--mode sum\|average`, `--metaspots auto\|grid\|knn`, `--row-key`, `--col-key`, `--spatial-key`, `-k`, `--max-shared`, `--seed`, `--layer`, `--preset hdwgcna\|wgcna`, `--power`, `--min-size`, `--deep-split`, `--merge-cut`, `--n-hubs`, `--spot-eigengenes/--no-spot-eigengenes`, `--reassign/--no-reassign`, `--plot`, `--write-metaspots`, `-t` |
| `celltype` | `single_cell.build_by_celltype()` + `write_celltype_networks()` | `-i` h5ad, `--groupby`, `-o` dir, `--group`, `--strategy`, `-k`, `--max-shared`, `--max-metacells`, `--min-cells`, `--use-rep`, `--layer`, `--hvg`, `-m`, `-c`, `-p`, `--fdr/--no-fdr` (default on), `--effective-samples`, `--diff` |
| `celltype-modules` | `single_cell.iter_modules_by_celltype()` + `write_celltype_modules()` (one cell type at a time) | `-i` h5ad, `--groupby`, `-o` dir, `--group`, `--fraction`, `--hvg`, `-k`, `--max-shared`, `--max-metacells`, `--min-cells`, `--use-rep`, `--layer`, `--preset hdwgcna\|wgcna`, `--power`, `--min-size`, `--deep-split`, `--merge-cut`, `--n-hubs`, `--write-metacells`, `--write-tom`, `--seed`, `--cell-eigengenes/--no-cell-eigengenes`, `--reassign/--no-reassign`, `--backed`, `--consensus-by`, `--consensus-quantile`, `--harmonize-by`, `--npcs`, `--dme-by`, `--trait`, `--trait-groupby` |
| `consensus-modules` | `consensus.consensus_tom_modules()` | `-e` (one per set), `-o` dir, `--set-name`, `--preset`, `--signed`, `--power`, `--consensus-quantile`, `--min-size`, `--deep-split`, `--merge-cut`, `-l`, `--n-hubs`, `-t` |
| `project-modules` | `single_cell.project_modules()` | `-i` h5ad, `-m` modules, `-o` dir, `--groupby`, `--group`, `--min-overlap`, `--harmonize-by`, `--layer`, `--npcs`, `--seed` |
| `module-preservation` | `preservation.module_preservation()` | `-r` reference expression, `-m` modules, `-q` test (TSV, or .h5ad: metacells built from it), `-o` dir, `--permutations`, `--seed`, `--network-type`, metacell options for an .h5ad query, `-t` |
| `trajectory` | `trajectory.build_trajectory()` | `-i` h5ad or TSV, `--pseudotime`, `-o` dir, `-w`, `-s`, `--metacell`, `-m`, `-c`, `-p`, `-f`, `--signed`, `-l`, `-t`, `--plot` |
| `fuse` | `multilayer.fuse_to_file()` | `-i` (one per layer), `-o`, `--name`, `-w`, `--normalize`, `--min-layers`, `--min-score`, `-t` |
| `rwr-multiplex` | `multilayer.rwr_multiplex()` | `-n` (one per layer), `-g`, `-o`, `--name`, `--gamma`, `--delta`, `--weighted`, `--aggregate`, `--prob-cutoff`, `--transition` |
| `plot <figure>` | `gcen2.plot.<function>` | 43 figure commands; all take `-o`, `--size`, `--format`, `--palette` |

### 4.12 gcen2.spatial

Sub-package with three spatial co-expression network methods (AnnData with `obsm['spatial']`, 2D coordinates; `scipy.spatial.KDTree` and `scipy.sparse`, the correlations in C++ through `network.build`) and WGCNA-style spatial modules (`gcen2.spatial.modules`, the hdWGCNA spatial workflow; the TOM and clustering in C++ through `tom_modules`).

#### gcen2.spatial.neighbors

**`gcen2._knn.knn_indices(points, k, workers=1)`** — the k nearest points of every point, itself first, then by squared distance and, for equal distances, by index; where the k-th and (k + 1)-th are equally far, all points within that distance are gathered first. Used by `build_neighbor_graph`, kNN metaspots and metacells, so the neighbours are the same on every platform (cKDTree's order of equidistant points is not; grids of spots or bins have many).

**`build_neighbor_graph(coords, n_neighbors=None, radius=None) → scipy.sparse.csr_matrix`**

Exactly one of `n_neighbors` or `radius` must be specified.

- **KNN mode:** query `n_neighbors + 1` nearest points (including self), drop self, compute Gaussian weights `exp(-d² / 2σ²)` where `σ = median(nonzero distances)`; symmetrize via `W.maximum(W.T)`.
- **Radius mode:** `KDTree.query_pairs(radius)` returns symmetric pairs; weights are 1.0.

Input validation: `n_neighbors >= 1`, `radius >= 0`. The resulting matrix is a CSR sparse matrix of shape `(n_spots, n_spots)` with non-negative entries.

#### gcen2.spatial.core

**`_morans_i(X, W)`** — Per-gene Moran's I spatial autocorrelation.

```
I = (n / W_sum) * (z' W z) / (z' z)
```
where `z = X - mean(X)` (column-wise). Returns shape `(n_genes,)`.

**`_get_expression(adata, layer, fraction=0.0, normalize=None)`** — the genes with a value > 0 in at least `round(fraction * n_spots)` spots (`single_cell.fraction_genes`; 0: all genes), Seurat LogNormalize of integer counts (`normalize=None`: automatic), dense `(n_spots, n_genes)` of those genes only, and their names. The three build functions take `fraction` (default 0.05) and `normalize` (default None) and use it, so the Moran's I, the smoothing and the correlations work on log-normalized expression of expressed genes.

**`_morans_test(X, W)`** — Moran's I, z and the one-sided p-value for positive autocorrelation, from the mean and variance of I under randomisation (Cliff and Ord 1981; any distribution of the values, asymmetric `W`). `_bh(p)`: Benjamini–Hochberg.

**`build_spatial_filtered(..., moran_fdr=0.05)`** — Standard network → Moran's I test per gene → keep edges where both genes have `I > 0` and FDR ≤ `moran_fdr` (`moran_fdr=1`: `I > 0` only).

**`spatial_autocorrelation(adata, layer=None, fraction=0.05, normalize=None, n_neighbors=None, radius=None, smoothed=False, bandwidth=None) → DataFrame`** — every gene's Moran's I (`gene`, `morans_i`, `pval`, `fdr`, `mean`, `detected`; with `smoothed` also `morans_i_smoothed`, the I after the smoothing of `build_spatial_weighted`), on the genes and neighbour graph of the build functions; rows by decreasing I. `gcen2 spatial-build --moran` writes it.

**`build_spatial_weighted`** — Gaussian neighbourhood smoothing (`_smooth`, shared with `spatial_autocorrelation`):
1. Row-normalize `W + I` (adding self-loops)
2. Smooth: `X_smooth = W_norm @ X`
3. Build standard network on smoothed expression

Optional `bandwidth` parameter overrides the auto-estimated Gaussian σ.

**`build_spatial_regions`** — Per-region networks:
- Manual: split by `adata.obs[region_key]`
- Auto: KMeans on `obsm['spatial']` coordinates into `n_regions` clusters
- Regions with < 3 spots emit `UserWarning` and are skipped
- Returns `dict[region_name → Network]`

#### gcen2.spatial.modules

hdWGCNA's spatial workflow (`SetupForWGCNA(gene_select = "fraction")` →
`MetaspotsByGroups` → `NormalizeMetacells` → `SetDatExpr` on all metaspots →
`TestSoftPowers` / `ConstructNetwork` → `ModuleEigengenes` on the spots).

**`metaspots(adata, groupby=None, row_key="array_row", col_key="array_col", min_spots=50, layer=None, mode="sum", method="auto", spatial_key="spatial", k=7, max_shared=None, seed=42, threads=1) → AnnData`** — `method="grid"`: hdWGCNA `ConstructMetaspots` per group (groups sorted by name; groups with fewer than `min_spots` spots, or with fewer than two centre columns, skipped with a warning). On the Visium grid coordinates of the group: centre columns are every 4th column from the smallest (positions 4, 8, … of `min..max`; the middle one of three, the first of one or two); centre rows the rows of the same parity as the centre columns. Members of a centre: the group's spots at the smallest non-zero grid distance in the group (the four diagonal neighbours on a Visium hexagon), then the spots of the same row in the columns two positions before and after in the group's sorted unique columns, in the group's spot order; centres with fewer than two members are dropped. As in hdWGCNA a spot can be counted twice when the group has no diagonal neighbours. Counts are summed (or averaged) through a sparse aggregation matrix. Returns metaspots × genes (sparse if the input is), index the centre barcodes, `obs` the centre's obs plus `n_spots`, `uns["metaspot_members"]` the member indices (centre first), `uns["metaspot_params"]` the method and its parameters.

`method="knn"` is for platforms without a Visium grid (Slide-seq beads, Visium HD bins, Stereo-seq, single cells of imaging platforms) and needs only `obsm[spatial_key]` (2 or 3 columns, any unit). hdWGCNA has no metaspot method for these data; its users build metacells in the spatial coordinates with `MetacellsByGroups`, and this is the same idea: within each group (at least `max(min_spots, k)` spots), every spot's `k` nearest spots of the group (itself first; `scipy.spatial.cKDTree`); centres tried in a random order (`numpy.random.default_rng(seed)`, one generator over the groups in sorted order) and kept when the neighbourhood shares at most `max_shared` spots with every neighbourhood kept so far (`gcen2.single_cell.metacell._select`, the greedy choice of `metacells`). Defaults: `k = 7` (a centre and its ring, the size of a Visium hexagon) and `max_shared = round(2k / 7)` (2 for k = 7, the overlap of neighbouring grid metaspots). `k` should grow as the units shrink, so that a metaspot holds about as many counts as a Visium one: about 100 for 10 µm Slide-seq beads. Members are the k nearest spots whatever their distance; there is no radius cut-off.

`method="auto"` (default) is the grid when `obs` has integer `row_key` / `col_key` columns whose sum has a single parity (the hexagonal Visium layout) and knn otherwise, including the square grid of Visium HD bins.

**`log_normalize(X, columns=None, scale=1e4)`** — Seurat LogNormalize of the rows (library size over all columns), dense output of `columns` only. Computed as `x * scale / library` with one rounding (`x * scale` is exact for counts; metacells likewise from the members' sums, cell eigengenes the same way), so mathematically equal ratios give equal values. Then the Spearman ranks of the result, and the networks built from them, do not depend on the last bit of the platform's `log1p`, which differs between CPUs (AVX-512 runners reordered near-equal values; checked by perturbing `log1p` by ±1 ulp as a function of its input: the golden outputs do not change).

**`spot_eigengenes(data, labels, levels) → (n_spots, n_modules)`** — hdWGCNA `ComputeModuleEigengene` on the spots: Seurat `ScaleData` of the module's genes (z-scores with the n − 1 SD, clipped above at 10), then `RunPCA`'s first component score `u[:, 0] * s[0]` of that matrix without re-centring (from the top eigenvector of the smaller Gram matrix, `scipy.linalg.eigh(subset_by_index)`), signed to correlate positively with the module's average log-normalized expression.

**`spatial_modules(adata, groupby=None, fraction=0.05, use_highly_variable=False, genes=None, min_spots=50, mode="sum", row_key, col_key, layer=None, threads=2, preset="hdwgcna", spot_level=True, reassign=True, method="auto", spatial_key="spatial", k=7, max_shared=None, seed=42, **tom_kwargs) → TomModules`** — raw counts required. Metaspots by `metaspots(method, spatial_key, k, max_shared, seed)`. Genes: `fraction_genes` over all spots (and `highly_variable`), or the list `genes`. Metaspots, LogNormalize, then `tom_modules(genes × metaspots, method="pearson", **HDWGCNA | tom_kwargs)`. `metadata`: `n_spots`, `n_metaspots`, `n_genes`, `metaspot_params`, `groups`, `metaspots` (genes × metaspots DataFrame) and, with `spot_level`, `spot_eigengenes` (spots × `ME<k>`, modules 1..K) computed on the log-normalized spots and `spot_kme` (genes × `kME<k>`: Pearson correlation of each gene's log-normalized spot expression with the spot eigengenes, hdWGCNA `ModuleConnectivity`, which uses `corSparse`). With `reassign=True` (default) and at least two modules, `reassign_negative_kme` then applies hdWGCNA's `ReassignModules`: a gene of module k > 0 with a negative spot kME for k moves to the module of its largest spot kME (first on ties), or to unassigned if all are negative; `metadata["reassigned"]` lists those genes. Eigengenes are not recomputed afterwards (as in hdWGCNA).

**`group_trait_cor(spot_me, groups)`** — `cor_pairwise` of the spot eigengenes with one 0/1 indicator per group: `module, trait, cor, p, n`.

**`write_spatial_modules(res, output_dir, groups=None, n_hubs=10, write_metaspots=False)`** — `write_module_files` (modules, metaspot eigengenes, kME, hub genes), `spot_eigengenes.tsv` (columns renamed `ME0..ME<K-1>` as the module numbers of `modules.tsv`), `spot_kme.tsv` (gene, module, `kME0..`), `group_cor.tsv` (with `groups`), `metaspots.tsv` (optional), `summary.tsv` (n_spots, n_metaspots, n_genes, power, modules, genes_in_modules, reassigned, groups, metaspots = grid or knn).

**Validation** (`benchmarks/hdwgcna_spatial`): on DLPFC 151673 (layers) and the 10x mouse brain sagittal anterior sample (k-means regions), the same genes, metaspots (names and order) and datExpr as hdWGCNA 0.4.12, the same modules after ConstructNetwork and after ModuleConnectivity's reassignment (ARI 1.0 for both), and the same spot eigengenes and spot kME. On Slide-seqV2 (41,786 beads, no grid; k = 100, the same 1,770 genes), `module -e --preset hdwgcna` on the metacells of hdWGCNA's `MetacellsByGroups` in the spatial coordinates gives its modules exactly (ARI 1.0); with each tool's own random metacells the ARI between the tools (0.66–0.85) matches each tool's seed-to-seed ARI (GCEN2 0.64–0.87, hdWGCNA 0.68–0.72), in 7.2 s against 220 s. On the Visium samples knn metaspots (k = 7) give spot eigengenes matching the grid ones (|r| ≥ 0.95) and gene-level ARI with the grid modules close to knn's own seed-to-seed ARI. `tests/data/hdwgcna_spatial/make_reference.R` pins the metaspot members, normalized metaspots and spot eigengenes on a synthetic grid with holes and a skipped small group.

### 4.13 gcen2.trajectory

**`build_trajectory(data, pseudotime, window_size=0.2, step=None, metacell_size=1, layer=None, use_highly_variable=False, method="spearman", cor=0.5, ...) → TrajectoryNetworks`** (CLI `gcen2 trajectory`)

Cells (an AnnData with a pseudotime `obs` column) or samples of a time
course (genes × samples DataFrame) are ordered by pseudotime and cut into
overlapping windows of equal size (a count or a fraction; step half a window
by default); optional metacells average consecutive cells within a window.
One `network.build` per window; every edge seen in any window is followed
across all windows (its correlation, the windows where it passes the
cutoffs, a pattern `stable`, `gained`, `lost`, `transient` or
`intermittent`, the change and trend along pseudotime), and gene degrees
per window.

### 4.14 gcen2.multilayer

A layer is any undirected gene network (co-expression, PPI, regulatory):
a `Network`, DataFrame or edge-list file, weighted by `|correlation|`, a
`weight` column or 1. **`fuse(layers, weights=None, normalize="max", min_layers=1, min_score=0.0) → Network`**: per edge the weighted mean of
its normalised layer weights (absent = 0), with the layers supporting it;
`fuse_to_file` reads, merges, scores and writes in C++
(`_gcen_core.network`). **`rwr_multiplex(layers, seed_genes, gamma=0.5, delta=0.5, weighted=False, aggregate="geometric", prob_cutoff=0.0, transition="supra", ...) → RWRResult`**: RWR-M
(Valdeolivas et al., Bioinformatics 2019) on the supra-adjacency matrix
(layers scaled by `1 - delta`, `delta / (L - 1)` between the copies of a
gene) with the C++ RWR kernel. `transition="supra"` column-normalises the
whole matrix, as the authors' RandomWalkRestartMH (same scores,
`benchmarks/multilayer_string`); `"layer"` normalises each layer first, so
every gene jumps with probability `delta`. Per-gene scores aggregated over
layers; ties ranked by gene name. One layer gives `rwr.predict`.

### 4.15 gcen2.harmony

**`harmony_integrate(embedding, batches, theta=2.0, sigma=0.1, nclust=None, max_iter=10, ...) → ndarray`**: Harmony (Korsunsky et al. 2019)
as R's harmony 2.0.5 `RunHarmony` with its defaults, in double precision:
cosine-normalized cells, `min(round(N / 30), 100)` soft clusters seeded by
harmony's k-means, rounds of diversity-penalized soft clustering and the
mixture-of-experts ridge correction. Used for hdWGCNA's harmonized module
eigengenes (`celltype-modules --harmonize-by`, `project-modules
--harmonize-by`): r >= 0.997 with hdWGCNA's on the HLCA benchmark.

---

## 5. Algorithm Details

### 5.1 Pairwise Correlation

**Pearson (C++):**
Compute per-gene mean, then standard covariation formula. Two-sided p-value via the incomplete beta function (Lentz continued fraction algorithm, same as GCEN v1 `pearson.hpp`).

**Spearman (C++):**
Pre-rank each gene's expression values (average rank for ties, same as GCEN v1 `base.hpp::get_rank()`), then apply Pearson on the rank vectors.

**Two-section threading:**
Gene pair space `[0, n_genes)` partitioned across `n_threads` threads. Each thread `n` handles:
- Section I: `[n*n_genes/(2*n_threads), (n+1)*n_genes/(2*n_threads))`
- Section II: `[(2*n_threads-n-1)*n_genes/(2*n_threads), (2*n_threads-n)*n_genes/(2*n_threads))`

This interleaved assignment ensures uniform workload (each thread handles equal-length row ranges from both ends of the gene list).

### 5.2 FDR Correction (Benjamini-Hochberg)

The total test count `n_total = n_genes * (n_genes - 1) / 2` (all possible pairs), not just filtered edges. This matches GCEN v1 and is required for correct FDR estimation. Implementation: sort p-values descending, compute `fdr[i] = p[i] * n_total / rank`, then apply monotone correction forward to ensure non-increasing values.

### 5.3 Module Identification (Jaccard + DFS)

**Phase 1 — Similarity (parallel):**

For each node pair `(i, j)`:
```
J(i,j) = (|N(i) ∩ N(j)| + a(i,j)) / (|N(i)| + |N(j)| - |N(i) ∩ N(j)| - a(i,j))
```
where `a(i,j) = 1` if edge `(i,j)` exists in the original network, else 0. Pairs with `J > cutoff` are added to the similarity graph.

**Phase 2 — DFS (single-threaded):**
Each connected component of the similarity graph is a module. DFS starts from each unvisited node.

**Edge case:** A 2-node graph `{G1, G2}` connected by a single edge gives `J(G1,G2) = 1.0`, forming a module regardless of cutoff (as long as cutoff < 1.0). This is correct algorithm behavior.

### 5.4 Random Walk with Restart

Given column-normalized adjacency matrix `W`, seed vector `P0`, restart probability `γ`:

```
P_{t+1} = (1 - γ) · W · P_t + γ · P_0
```

Convergence: `||P_{t+1} - P_t||₁ < tol`. Default `tol = 1e-10`, `max_iter = 1,000,000`.

**Column normalization:** `W[col] /= sum(W[col])`. For directed networks, out-edges of each node are normalized. For undirected, the symmetric adjacency matrix is used as-is and normalized column-wise.

### 5.5 Hypergeometric P-value (Fisher's Exact Test)

Given a 2×2 contingency table `[[a, b], [c, d]]`:
```
P = sum over x of C(a+b,x)*C(c+d,a+c-x) / C(n, a+c)   for all x with P(x) <= P(observed)
```

Computed via log-factorial table to avoid numerical overflow. Log-factorials accumulated as `log_fac[i] = log_fac[i-1] + log(i)`. The result is clamped to `[0.0, 1.0]` to handle floating-point precision accumulation in the sum.

### 5.6 Differential Co-Expression (Fisher's Z Test)

For a gene pair with correlations `r_a` (condition A, `n_a` samples) and `r_b` (condition B, `n_b` samples):

```
z_a = 0.5 * ln((1 + r_a) / (1 - r_a))
z_b = 0.5 * ln((1 + r_b) / (1 - r_b))
z_test = (z_a - z_b) / sqrt(v/(n_a-3) + v/(n_b-3))
```

`v = 1` for Pearson and `v = 1.06` for Spearman correlations (Fieller et
al., Biometrika 1957; as DGCA); `r` is clipped to ±0.9999 first. `z_test`
follows the standard normal under H₀; the two-sided p-values of all tested
pairs are adjusted by Benjamini-Hochberg (as R's `p.adjust(p, "BH")`).
Verified against DGCA 1.0.3 (`benchmarks/dgca_gtex`).

### 5.7 Pseudo-bulk Strategy

**With biological replicates:** cells grouped by `sample × celltype`, summed to produce pseudo-bulk profiles.

**Without biological replicates (single sample):** cells randomly partitioned into `pseudo_replicates` (default 10) artificial samples. Partitioning is seeded (`seed` parameter) for reproducibility. Deprecated: random groups have the same expected composition, so their sums estimate cell-level correlations from a few samples (on PBMCs, chance-level rankings).

---

## 6. Build, Tests and Migration

- **Build, dependency groups, CI, release:** Implementation Document §1, §4, §5.
- **Test design** (golden outputs, reference results, two implementations
  compared, completeness checks): Implementation Document §4; results in
  `docs/test_report.md`.
- **Migration from GCEN v1** (command mapping, verified equivalence,
  deliberate differences): `docs/migration.rst`.

---

*End of Document*
