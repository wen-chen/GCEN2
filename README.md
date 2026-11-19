# GCEN2

[![platforms](https://img.shields.io/badge/platforms-Linux%20%7C%20Windows%20%7C%20macOS-yellowgreen)](https://www.biochen.org/gcen2/install.html)
[![Python](https://img.shields.io/badge/python-3.9%E2%80%933.14-blue)](https://www.biochen.org/gcen2/install.html)
[![license GPL-3.0-or-later](https://img.shields.io/badge/license-GPL--3.0--or--later-blue)](LICENSE)
[![CI](https://github.com/wen-chen/GCEN2/actions/workflows/ci.yml/badge.svg)](https://github.com/wen-chen/GCEN2/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/gcen2)](https://pypi.org/project/gcen2/)
[![Bioconda](https://img.shields.io/conda/vn/bioconda/gcen2?label=bioconda)](https://anaconda.org/bioconda/gcen2)
[![stars](https://img.shields.io/github/stars/wen-chen/GCEN2?style=social)](https://github.com/wen-chen/GCEN2/stargazers)
[![GitHub issues](https://img.shields.io/github/issues/wen-chen/GCEN2)](https://github.com/wen-chen/GCEN2/issues)

**A Python Toolkit for Gene Co-Expression Network Analysis with Single-Cell Resolution and Interactive Visualization**

GCEN2 is a major upgrade of [GCEN v1](https://github.com/wen-chen/GCEN): a Python package
and command-line tool (`gcen2`) with C++ kernels for speed. It does

- **everything GCEN v1 does**, with the same commands, files and results: normalization,
  filtering, correlation networks, modules, GO / KEGG annotation, enrichment, random walk
  with restart, merge, extract, shuffle, statistics;
- **differential co-expression** between two conditions (DGCA's test);
- **WGCNA module analysis**: soft power, TOM, dynamic tree cut, eigengenes, kME, hub genes,
  traits, consensus modules over several data sets, module preservation;
- **the hdWGCNA single-cell and spatial workflows**: metacells, modules per cell type,
  cell and spot eigengenes, Harmony, differential module eigengenes, module projection;
  Visium, Slide-seq and Visium HD;
- **networks along a trajectory and multi-omics networks** (fusion, multiplex random walk);
- **publication figures and reports** for every analysis.

Against GCEN v1, WGCNA, hdWGCNA, DGCA and RandomWalkRestartMH it gives the same results,
faster and mostly with less memory (see [Benchmarks](#benchmarks)).

[![The GCEN2 workflow](https://www.biochen.org/gcen2/static/img/workflow.png)](https://www.biochen.org/gcen2/)

*The GCEN2 workflow: four kinds of data, from input to modules and their comparison, each
step labelled with its `gcen2` command.*

- Website: <https://www.biochen.org/gcen2/>
- Source code: <https://github.com/wen-chen/GCEN2>
- GCEN v1: <https://www.biochen.org/gcen/> and <https://github.com/wen-chen/GCEN>

## Installation

GCEN2 is distributed in three ways, for Python 3.9 to 3.14:

| | Linux (x86-64 / AMD64) | Windows (x86-64 / AMD64) | macOS (Apple Silicon) |
|---|---|---|---|
| PyPI (pre-built wheels) | ✓ | ✓ | ✓ |
| Bioconda | ✓ | — (Bioconda does not build for Windows) | ✓ |
| Source code | ✓ | ✓ | ✓ |

PyPI:

```bash
pip install "gcen2[all]"        # or gcen2 (core), gcen2[single-cell], [visualization], [report]
```

Bioconda (Linux, macOS):

```bash
conda install -c conda-forge -c bioconda gcen2
```

From source (a C++17 compiler and CMake are needed; `pip` fetches CMake if it is missing):

```bash
git clone https://github.com/wen-chen/GCEN2
cd GCEN2
pip install ".[all]"
```

Other platforms (for example Linux on ARM or macOS on Intel) are not supported
targets; installing from source may work but is not tested.

## Quick start

The GCEN v1 pipeline on the GCEN v1 sample data
([`data/example/bulk/`](data/example/README.md); run in an empty directory):

```bash
EX=/path/to/GCEN2/data/example/bulk
gcen2 normalize -i $EX/gene_expr.tsv -o norm.tsv -m tmm
gcen2 filter    -i norm.tsv -o filt.tsv -p 0.75
gcen2 build     -i filt.tsv -o net.network -m spearman -c 0.8 -p 0.001 -f -t 8
gcen2 module    -i net.network -o modules.txt -s 0.5 -t 8
gcen2 annotate  -g $EX/go-basic.obo.gz -a $EX/gene_go.assoc.gz -n net.network -o annotation/
gcen2 rwr       -n net.network -g $EX/rwr_seed_genes.list -o rwr.tsv
gcen2 report    -n net.network -m modules.txt -o report.html
```

WGCNA modules of the most variable genes, a figure, and the hdWGCNA workflow on the
example single cells:

```bash
gcen2 filter -i norm.tsv -o var.tsv -p 0.75 --percent-sd 0.2
gcen2 module -e var.tsv -o wgcna/modules.tsv -l log2 --signed --kme wgcna/kme.tsv -t 8
gcen2 plot kme -i wgcna/kme.tsv -o fig/kme                  # PDF, PNG and SVG
gcen2 celltype-modules -i $EX/../single_cell/pbmc_subset.h5ad --groupby cell_type \
      -o modules/ -t 8
```

The same in Python:

```python
import gcen2

expr = gcen2.read_tsv("filt.tsv")
net = gcen2.network.build(expr, method="spearman", cor=0.8, pval=0.001, fdr=True, threads=8)
modules = gcen2.module.identify(net, similarity=0.5)
wgcna = gcen2.tom_modules.tom_modules(gcen2.read_tsv("var.tsv"), log="log2", signed=True,
                                      threads=8)
```

## Documentation

| | |
|---|---|
| [Quickstart](docs/quickstart.rst) | installation, first analysis |
| Tutorials on the example data | [bulk](docs/tutorials/bulk.rst), [single cell](docs/tutorials/single_cell.rst), [spatial](docs/tutorials/spatial.rst), [trajectories and multi-omics](docs/tutorials/trajectory_multiomics.rst), [figures and reports](docs/tutorials/figures.rst) |
| [Command-line reference](docs/cli.rst) | every command and option |
| [Python API](docs/api.rst) | every module and function |
| [Migration from GCEN v1](docs/migration.rst) | command mapping, verified equivalence, deliberate differences |
| [Example data](data/example/README.md) | bulk, single-cell and spatial data, with sources and licences |
| [For developers](docs/README.md) | requirements, design and implementation documents, code guide, status, test report |

Build the documentation site with `pip install -r docs/requirements.txt` and
`sphinx-build -b html docs docs/_build/html`.

## Benchmarks

On one 16-core machine (data, scripts and the full table in
[`benchmarks/`](benchmarks/README.md); all times measured on an idle machine):

| Compared with | Data | Results | Speed |
|---|---|---|---|
| GCEN v1 0.6.3 | v1 sample data (32,940 genes) | same edges and modules | 6–89× faster (merge, extract 1.1–1.5×) |
| WGCNA 1.74 | GTEx whole blood (11,944 genes × 803 samples) | same modules in all six runs | 4.0×, 2.7× less memory |
| hdWGCNA 0.4.12 | PBMC 10k | same modules in all 24 runs | 9.1–13× (network step) |
| hdWGCNA 0.4.12 | Human Lung Cell Atlas (584,944 cells) | same genes; same modules on the same metacells | 14×; all 49 cell types with ≥ 100 cells in about 8 min |
| hdWGCNA 0.4.12 | consensus and Harmony; projection and preservation | same powers, modules, DMEs, preservation statistics | about 7× |
| hdWGCNA 0.4.12 | Visium, Slide-seqV2, Visium HD | same metaspots, modules and spot eigengenes | 7.4–31× |
| DGCA 1.0.3 | GTEx whole blood, males vs females (12.5M gene pairs) | same differential correlation test and FDR, Pearson and Spearman | 2.6× (Pearson), 185× (Spearman) |
| RandomWalkRestartMH 1.13.2 | co-expression + STRING + TRRUST layers | same multiplex RWR scores | 58× |
| R implementations (no published package) | the same layers; paul15 myeloid differentiation | same fused network; same networks along pseudotime | 1.4×; 5.2× |

## Contributing

See [`CONTRIBUTING.md`](CONTRIBUTING.md) (setup, tests, rules) and the
[code of conduct](CODE_OF_CONDUCT.md). Changes are listed in [`CHANGELOG.md`](CHANGELOG.md).

## Citation

Please cite GCEN2 as described in [`CITATION.cff`](CITATION.cff) (GCEN v1: Chen W. et al.,
*Curr. Issues Mol. Biol.* 2022, doi:10.3390/cimb44040100), and WGCNA / hdWGCNA when you use
their methods.

## License

GNU General Public License, version 3 or later (GPL-3.0-or-later); see [LICENSE](LICENSE).
The example data keep their own licences ([`data/example/README.md`](data/example/README.md)).
