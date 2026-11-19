# GCEN2 benchmarks

Each comparison has its own directory with scripts, raw results and a
README with the details:

| directory | against | data |
|---|---|---|
| this one (`compare_v1_v2.py`) | GCEN v1 0.6.3, command by command | the v1 sample matrix (32,940 genes × 15 samples) and subsets |
| [`wgcna_gtex/`](wgcna_gtex/README.md) | WGCNA 1.74 | GTEx v10 whole blood, 803 samples |
| [`pywgcna_gtex/`](pywgcna_gtex/README.md) | PyWGCNA 2.2.1 against WGCNA 1.74 (step by step and blockwiseModules) and GCEN2 | GTEx v10 whole blood, the six matrices of `wgcna_gtex/` |
| [`hdwgcna_pbmc/`](hdwgcna_pbmc/README.md) | hdWGCNA 0.4.12, single-cell workflow | 10x PBMC 10k, six cell types |
| [`hdwgcna_hlca/`](hdwgcna_hlca/README.md) | hdWGCNA 0.4.12, single-cell workflow | Human Lung Cell Atlas core, 584,944 cells, and a 50,000-cell subset |
| [`hdwgcna_consensus/`](hdwgcna_consensus/README.md) | hdWGCNA 0.4.12, consensus network, harmonized eigengenes, DMEs and module-trait correlation | HLCA natural killer cells over 4 datasets, 38 donors |
| [`hdwgcna_preservation/`](hdwgcna_preservation/README.md) | hdWGCNA 0.4.12, ProjectModules and ModulePreservation | HLCA natural killer cells, one dataset onto another |
| [`hdwgcna_spatial/`](hdwgcna_spatial/README.md) | hdWGCNA 0.4.12, spatial workflow | Visium DLPFC 151673 and 10x mouse brain anterior; Slide-seqV2 (no grid); Visium HD mouse brain, 16 µm bins |
| [`dgca_gtex/`](dgca_gtex/README.md) | DGCA 1.0.3, differential co-expression | GTEx v10 whole blood, males vs females |
| [`multilayer_string/`](multilayer_string/README.md) | RandomWalkRestartMH 1.13.2 (multiplex RWR); an R implementation of the network fusion | GTEx co-expression, STRING v12 and TRRUST v2 layers |
| [`trajectory_paul15/`](trajectory_paul15/README.md) | an R implementation (WGCNA functions) of the networks along pseudotime | Paul et al. 2015 myeloid differentiation, 2,730 cells |
| [`scaling/`](scaling/README.md) | no reference tool: thread scaling (`-t` 1–16) of `build`, `module -e`, `diff`, `celltype-modules`; parameter sensitivity of networks and modules | GTEx v10 whole blood; HLCA alveolar macrophages |
| [`case_gtex_blood/`](case_gtex_blood/README.md) | no reference tool: a biological case study (module GO enrichment and donor traits, cross-validated RWR function prediction) | GTEx v10 whole blood, 803 samples; Gene Ontology 2026-08-05 |

## Reproducing

- **Environments**: `environment-r.yml` (R 4.4.3 with WGCNA 1.74, Seurat 5.5.1,
  harmony 2.0.5 and the other dependencies of hdWGCNA; hdWGCNA 0.4.12,
  DGCA 1.0.3 and RandomWalkRestartMH 1.13.2 installed as the file shows) and `requirements-python.txt` (the Python packages,
  with GCEN2 installed from this repository). The `run_*.sh` scripts use
  `RSCRIPT` (default: `Rscript` on the `PATH`).
- **Data**: each benchmark's README lists its downloads (URLs; checksums where
  the source does not version its files). Scripts run from a data directory
  of your choice (`data/` for the downloads, `work/` and `runs/` are written
  there); nothing refers to a fixed path.
- **Order**: the GTEx benchmark reuses the STRING gold standard written by the
  PBMC benchmark (`GCEN2_STRING_GOLD`); the DGCA and multi-layer benchmarks
  start from the GTEx matrices (`wgcna_gtex/01_prepare.py`) and the
  multi-layer one from the STRING download of the PBMC benchmark; the
  consensus and preservation benchmarks start from the HLCA download.
- **Re-timing**: `RETIMING_PLAN.md` (how the times below were measured, and
  its record table).

## Overview

All on one machine: Intel Core i7-13700 (16 cores / 24 threads), 123 GB RAM,
Ubuntu 26.04.1, g++ 15.2.0. All times were measured on 2026-11-17 on an
otherwise idle machine, one benchmark at a time, both tools in the same
session: GCEN2 at commit 91bf932, except `gcen2 module -e` (the WGCNA-style
modules of GTEx and the network-step rows), measured again at faed7ac, which
limits scipy's BLAS threads as well (at 91bf932 its eigengenes ran on all 24
threads). Each benchmark's README gives the number of runs
([`RETIMING_PLAN.md`](RETIMING_PLAN.md), record table).

Time and peak memory (`/usr/bin/time`) cover the whole command or script,
reading the input included: 4 threads against GCEN v1, 8 against the
other tools. A range spans the six PBMC cell types. "Faster"
and "less memory" are the other tool's figure divided by GCEN2's.

The table is generated from the committed results by `overview.py` (also
written to `overview.tsv`); `tests/test_benchmark_overview.py` checks it is up
to date.

<!-- overview:begin -->
| vs | data | task | GCEN2 | other tool | faster | less memory | results |
|---|---|---|---|---|---|---|---|
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | build | 3.1 s, 52 MB | 19.1 s, 13 MB | 6.1× | 4.0× more | identical edges (3,892,110), max \|Δr\| = 0 |
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | module | 2.5 s, 101 MB | 112 s, 748 MB | 44× | 7.4× | identical modules (1,657) |
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | rwr | 0.78 s, 120 MB | 55.0 s, 7.54 GB | 70× | 64× | max \|Δp\| = 2.5e-09, top 100 identical |
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | annotate | 2.8 s, 252 MB | 213 s, 281 MB | 75× | 1.1× | runtime only (one- vs two-sided test) |
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | stats | 1.6 s, 102 MB | 140 s, 578 MB | 89× | 5.7× | same statistics |
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | extract | 0.74 s, 58 MB | 0.81 s, 3 MB | 1.1× | 21× more | identical edges |
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | merge | 8.0 s, 1.23 GB | 12.2 s, 1.14 GB | 1.5× | 1.1× more | same edge count (9,820,818); sums differ by design |
| GCEN v1 0.6.3 | 32,940 genes × 15 samples, 3.9M edges | shuffle | 3.8 s, 387 MB | 126 s, 567 MB | 33× | 1.5× | both keep every degree; different algorithms |
| WGCNA 1.74 | GTEx whole blood, 11,944 genes × 803 samples | WGCNA modules (`module -e`) | 9.1 s, 1.45 GB | 36.1 s, 3.93 GB | 4.0× | 2.7× | same modules in 6 of 6 runs (ARI 1.0) |
| WGCNA 1.74 | GTEx whole blood, 11,944 genes × 803 samples | correlation network + modules, FDR (`build` + `module`) | 2.8 s, 633 MB | 36.1 s, 3.93 GB | 13× | 6.4× | different method (r and FDR, not TOM) |
| hdWGCNA 0.4.12 | PBMC 10k, 6 cell types, 44–1,000 metacells | hdWGCNA modules (`module -e --preset hdwgcna`) | 3.5–6.3 s, 645–993 MB | 38.4–58.6 s | 9.1–13× | – | same modules in 24 of 24 runs (ARI 1.0); hdWGCNA time: soft power + network only |
| hdWGCNA 0.4.12 | PBMC 10k, 6 cell types, 355–3,391 cells | correlation network + modules per cell type (`03_gcen2.py`) | 1.5–7.3 s, 352 MB–1.09 GB | 75.9–120 s, 4.21–6.02 GB | 14–55× | 4.0–18× | different method (r, no TOM); whole script vs hdWGCNA's |
| hdWGCNA 0.4.12 | Human Lung Cell Atlas core (584,944 cells), the 10 largest cell types (430,202 cells) | single-cell workflow (`celltype-modules`) | 211 s, 11.9 GB | 2,952 s, 43.9 GB | 14× | 3.7× | same genes in 10 of 10 types; same modules on hdWGCNA's metacells in 10 of 10 (ARI 1.0) |
| hdWGCNA 0.4.12 | HLCA natural killer cells, 16,006 cells in 4 datasets, 5,080 genes | consensus network, harmonized MEs, DMEs, traits (`celltype-modules --consensus-by --harmonize-by --dme-by --trait`) | 23.7 s, 2.51 GB | 161 s, 8.13 GB | 6.8× | 3.2× | same modules on hdWGCNA's metacells (ARI 1.0); harmonized MEs r >= 0.997; same DMEs and trait correlations; own metacells: ARI 0.76–0.85, as seed to seed |
| hdWGCNA 0.4.12 | HLCA natural killer cells, Meyer_2019 (6,732 cells) onto Krasnow_2020 (4,303 cells), 4,898 genes | reference modules, projection and preservation, 500 permutations (`celltype-modules`, `project-modules`, `module-preservation`) | 325 s, 2.54 GB | 2,174 s, 3.87 GB | 6.7× | 1.5× | projected MEs within 2e-13; on hdWGCNA's metacells the same preservation statistics (Z within 7e-08 relative) |
| hdWGCNA 0.4.12 | Visium DLPFC 151673, 3,611 spots, 9,317 genes | spatial workflow (`spatial-modules`) | 9.0 s, 1.28 GB | 113 s, 5.85 GB | 13× | 4.6× | same modules (ARI 1.0), spot MEs within 5e-09 (TSV precision) |
| hdWGCNA 0.4.12 | Visium mouse brain anterior, 2,695 spots, 12,861 genes | spatial workflow (`spatial-modules`) | 18.9 s, 1.94 GB | 140 s, 10.6 GB | 7.4× | 5.5× | same modules (ARI 1.0), spot MEs within 5e-09 (TSV precision) |
| hdWGCNA 0.4.12 | Slide-seqV2 hippocampus (no grid), 41,786 beads, 1,770 genes | spatial workflow, knn metaspots (`spatial-modules -k 100`) | 7.2 s, 670 MB | 220 s, 2.97 GB | 31× | 4.5× | same modules on hdWGCNA's metacells (ARI 1.0); own metacells: ARI 0.66–0.85, as seed to seed |
| hdWGCNA 0.4.12 | Visium HD mouse brain, 16 µm bins (square grid), 98,917 bins, 4,316 genes | spatial workflow, knn metaspots (`spatial-modules -k 64`) | 19.2 s, 2.23 GB | 570 s, 18.8 GB | 30× | 8.4× | same modules on hdWGCNA's metacells (ARI 1.0); own metacells: ARI 0.92–0.94, as seed to seed |
| DGCA 1.0.3 | GTEx whole blood, 5,000 genes, 541 male vs 262 female | differential co-expression of every pair, Pearson (`build` × 2 + `diff`) | 9.0 s, 1.26 GB | 23.2 s, 4.03 GB | 2.6× | 3.2× | same test on 12.5M pairs: r within 5e-07 (6 decimals), z within 9e-15 on the same r; 21 pairs differ at p = 0.05, 0 at FDR = 0.05 |
| DGCA 1.0.3 | GTEx whole blood, 5,000 genes, 541 male vs 262 female | differential co-expression of every pair, Spearman (`build` × 2 + `diff`) | 8.6 s, 1.27 GB | 1,582 s, 4.33 GB | 185× | 3.4× | same test on 12.5M pairs: r within 5e-07 (6 decimals), z within 7e-15 on the same r; 19 pairs differ at p = 0.05, 4 at FDR = 0.05 |
| RandomWalkRestartMH 1.13.2 | GTEx co-expression + STRING + TRRUST, 18,082 genes, 3 layers | multiplex RWR (`rwr-multiplex`) | 0.49 s, 179 MB | 28.3 s, 383 MB | 58× | 2.1× | scores within 1e-12, top 100 identical |
| R (data.table) | GTEx co-expression + STRING + TRRUST, 18,082 genes, 3 layers | network fusion (`fuse`) | 0.42 s, 194 MB | 0.57 s, 153 MB | 1.4× | 1.3× more | same 376,750 edges, scores within 6e-16 |
| R (WGCNA) | paul15 myeloid differentiation, 2,730 cells, 2,000 genes, 9 windows | networks along pseudotime (`trajectory`) | 1.2 s, 329 MB | 6.1 s, 1.09 GB | 5.2× | 3.4× | same edges in every window, same patterns of 15,017 edges |

Where GCEN2 is behind GCEN v1 on the gene subsets (the fixed cost is the Python + numpy start-up, about 0.1 s and 60 MB):

- `build`: more memory on 4 of 4 (GCEN2 46–52 MB, v1 8–13 MB)
- `module`: more memory on 1 of 4 (GCEN2 58 MB, v1 27 MB)
- `stats`: more memory on 1 of 4 (GCEN2 58 MB, v1 22 MB)
- `extract`: slower on 3 of 4 subsets (0.16 s vs 0.03 s on 5,014 genes); more memory on 4 of 4 (GCEN2 57–58 MB, v1 3 MB)
- `merge`: slower on 1 of 4 subsets (0.27 s vs 0.18 s on 5,014 genes); more memory on 4 of 4 (GCEN2 92 MB–1.23 GB, v1 29 MB–1.14 GB)
- `shuffle`: more memory on 3 of 4 (GCEN2 76–252 MB, v1 21–216 MB)
<!-- overview:end -->

What the comparisons show:
- Where both tools compute the same thing, GCEN2 gives the same result: the
  same edges and modules as GCEN v1 (RWR probabilities within 2.5e-9), the
  same modules as WGCNA and hdWGCNA (adjusted Rand index 1.0 in every run),
  DGCA's differential correlation test, RandomWalkRestartMH's multiplex RWR
  scores, and the R implementations of the fusion and of the networks along
  pseudotime. Rows marked "different method" compare GCEN2's own correlation
  workflow with the other tool's for scale only.
- GCEN2 is faster in every comparison on whole data sets. The gain is
  smallest where the work is mostly reading and writing files (`extract`,
  `merge` against GCEN v1; `diff` of two all-pairs networks against DGCA) or
  a single join (`fuse` against data.table).
- The cost is memory on small inputs: every GCEN2 command starts a Python
  interpreter with numpy (about 60 MB), so it can use more memory than v1's
  C++ binaries (the list under the table).

## GCEN v1 vs GCEN2

`compare_v1_v2.py` runs GCEN v1 and GCEN2 through their command-line tools on
subsets of the sample expression matrix. For each step it records wall time,
peak memory and how closely the results agree.

### Build GCEN v1 0.6.3

The source package of the release holds the sample data as well
(`sample_data/`: the expression matrix, GO files, seed and gene lists; GCEN v1
reads them uncompressed):

```bash
curl -LO https://github.com/wen-chen/GCEN/releases/download/v0.6.3/gcen-0.6.3-source.tar.gz
tar xzf gcen-0.6.3-source.tar.gz && cd gcen-0.6.3-source
# newer GCC needs <cstdint>/<limits> for the bundled robin_hood.h
make CCFLAGS="-std=c++11 -O3 -Wall -static-libstdc++ -D LINUX -include cstdint -include limits"
```

### Run

```bash
pip install --no-build-isolation -e .         # GCEN2, see ../CONTRIBUTING.md
python benchmarks/compare_v1_v2.py --v1 /path/to/gcen-0.6.3-source \
    --data /path/to/gcen-0.6.3-source/sample_data \
    --sizes 5000 10000 20000 0 --threads 4 --mem-limit-gb 32 --out /fast/disk/bench
```

`--data` defaults to `data/sample` of this repository, a local copy of that
directory (not committed).

Size `0` means all genes. Seed genes for RWR are always kept in each subset.
`results.json` has one row per step and tool, including the agreement metrics;
`summary.tsv` has the time and memory table.

### Committed results

`results.json` / `summary.tsv` were measured on 2026-11-17 with the command
above (GCEN2 commit 91bf932), on native Linux with the data, the virtualenv
and the outputs on a local file system. The all-genes rows are in the
overview above; `summary.tsv` has every subset.

| | |
|---|---|
| CPU | Intel Core i7-13700 (16 cores / 24 threads, 30 MB L3); `--threads 4` |
| RAM | 123 GB; `--mem-limit-gb 32` per process |
| OS | Ubuntu 26.04.1 LTS, kernel 7.0.0-43, ext4 on NVMe |
| Compiler | g++ 15.2.0 (both tools) |
| GCEN v1 | 0.6.3, built as above |
| GCEN2 | 2.0.0.dev0 (commit 91bf932), Python 3.14.10, numpy 2.5.4, pandas 3.1.1, scipy 1.18.1 |

The 15-sample matrix stays on the per-pair correlation kernel (the
matrix-product one starts at 64 samples), so `build` gives v1's edges bit for
bit. Against the previous native run (commit a2c3822, 2026-11-15) the GCEN2
times moved by at most 7%; GCEN v1's by a few percent, except `build` on
20,003 genes (7.1 -> 8.5 s) and `stats` on all genes (149 -> 140 s). Older results (WSL2, 4 threads, 7.8 GB
RAM, virtualenv on `/mnt/d`, 6 GB limit) are in the git history; the v1 RWR on
all genes failed there.

### Notes

- Each process runs under an address-space limit (`--mem-limit-gb`, default 6;
  use 32 for the full sample): the dense v1 RWR needs 8·n² bytes (7.7 GB peak
  on all genes).
- Every v2 row includes the CLI start-up time, reported separately as the
  `startup` row (`gcen2 --help`: 0.1 s and 30 MB here, 1.6 s on WSL with the
  virtualenv under `/mnt/*`). `import gcen2` is lazy; with the compiled
  extension `build`, `module`, `rwr`, `annotate`, `stats`, `extract`,
  `merge` and `shuffle` import neither pandas nor scipy (`diff` imports scipy
  for the Fisher z test but not pandas). On small networks `extract` is still
  slower than v1's plain line copy (0.15 s vs 0.02 s at 5,000 genes): that is
  the Python start-up and the typer / numpy imports.
- `annotate` is compared on runtime only. Queries, background and counts are
  the same as v1 (checked in `tests/test_annotation.py::TestV1AnnotateFixture`),
  but GCEN2 uses a one-sided test and v1 a two-sided one, so the sets of
  significant terms differ.
- `stats`, `extract`, `merge` and `shuffle` read the same input for both tools:
  the v1 network of the `build` step (GCEN2 reads v1 files; v1 cannot read the
  GCEN2 header line). `extract` uses `data/sample/gene_list.txt`.
- `merge` combines that network with a second one, built by GCEN2 at
  `cor - 0.05` (0.85) with its header removed; this build is not timed. It is
  compared on runtime and edge count only: v1 merges unordered pairs and sums
  positive and negative correlations apart (appending the sign to the value),
  GCEN2 sums ordered pairs into one value with a separate `direction` column.
- `shuffle` uses `--swap-factor 10` (v1 `-s 10`) for both. The algorithms
  differ: v1 picks random nodes and a random neighbour of each and repeats
  until `10 × edges` swaps *succeed*; GCEN2 picks random edge pairs and makes
  `10 × edges` *attempts* (about 2/3 succeed on the 3.6M-edge sample network). The agreement
  checks that both keep every node degree and reports the fraction of edges
  changed (v1 0.77, GCEN2 0.79 on all genes). At its default `-s 100`, v1
  would need ten times as many successful swaps (roughly 20 minutes on all
  genes, extrapolated, not measured).
