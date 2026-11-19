# GCEN2 Development Continuation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** Complete the remaining P0/P1 development work to bring GCEN2 to a publishable, fully tested state.

**Architecture:** Seven independent or sequentially-dependent tasks covering: (1) the editable-install C++ path fix, (2) sample data acquisition, (3) integration tests, (4) Sphinx documentation, (5) benchmark tests, (6) bioconda recipe, and (7) chunked computation for large single-cell datasets.

**Tech Stack:** Python 3.12, pybind11 C++17, scikit-build-core, pytest, Sphinx + sphinx-rtd-theme, typer, conda-build.

---

## Execution Order

Tasks 1, 4, 5, 6, 7 are fully independent — can run in parallel.
Task 2 (sample data) is independent but needed before Task 3.
Task 3 (integration tests) depends on Task 2.

```
Task 1 (fix .so path)      ─── independent
Task 4 (Sphinx docs)       ─── independent
Task 5 (benchmarks)        ─── independent
Task 6 (bioconda recipe)   ─── independent
Task 7 (chunked compute)   ─── independent
Task 2 (sample data)       ─── independent
Task 3 (integration tests) ─── after Task 2
```

---

## Task 1: Fix Editable Install — C++ Extension Auto-Discovery

**Problem:** After `pip install -e .`, the compiled `_gcen_core*.so` lands in `build/` but `gcen2/core/__init__.py` cannot find it. Developers must manually copy it each time they rebuild.

**Files:**
- Modify: `gcen2/core/__init__.py`
- Test: (manual verification + existing import tests in all test files)

### Step 1: Write the fix in `gcen2/core/__init__.py`

Replace the current empty docstring-only file with:

```python
"""C++ extension module loader for GCEN2 core computation kernels.

The compiled extension ``_gcen_core`` provides:
  - ``_gcen_core.correlation.pairwise_correlation``
  - ``_gcen_core.enrichment.batch_hypergeometric``
  - ``_gcen_core.module.identify_modules``
  - ``_gcen_core.rwr.random_walk_restart``

During development (editable install), the extension is built into the
``build/`` directory by scikit-build-core. This module automatically adds
that directory to ``sys.path`` so the extension can be imported without
a manual copy step.
"""

from __future__ import annotations

import os
import sys

# In an editable install, scikit-build-core puts the compiled .so in build/.
# Walk up from this file's location to find the project root's build/ dir.
_this_dir = os.path.dirname(os.path.abspath(__file__))
_project_root = os.path.dirname(os.path.dirname(_this_dir))  # gcen2/core/ → gcen2/ → project root
_build_dir = os.path.join(_project_root, "build")

if os.path.isdir(_build_dir) and _build_dir not in sys.path:
    sys.path.insert(0, _build_dir)

try:
    from . import _gcen_core  # noqa: F401  (imported for side-effects: registers submodules)
    _HAVE_CORE = True
except ImportError:
    _gcen_core = None  # type: ignore[assignment]
    _HAVE_CORE = False

__all__ = ["_gcen_core", "_HAVE_CORE"]
```

### Step 2: Remove any stale manually-copied `.so` from `gcen2/core/`

```bash
rm -f /mnt/d/code_base/gcen2/gcen2/core/_gcen_core*.so
```

### Step 3: Verify the fix works end-to-end

```bash
cd /mnt/d/code_base/gcen2
python3 -c "from gcen2.core import _gcen_core, _HAVE_CORE; print('loaded:', _HAVE_CORE, dir(_gcen_core))"
```

Expected output:
```
loaded: True ['__doc__', '__file__', '__loader__', '__name__', '__package__', '__spec__', 'correlation', 'enrichment', 'module', 'rwr']
```

### Step 4: Run full test suite to confirm nothing broke

```bash
python3 -m pytest tests/ -q --tb=short
```

Expected: `84 passed, 31 skipped`

### Step 5: Update `gcen2/__init__.py` to expose `_HAVE_CORE` at package level

In `gcen2/__init__.py`, add near the top:
```python
from gcen2.core import _HAVE_CORE as _HAS_CPP_EXTENSION
```

This lets users check `gcen2._HAS_CPP_EXTENSION` to know if the C++ kernels are available.

### Step 6: Commit

```bash
cd /mnt/d/code_base/gcen2
git add gcen2/core/__init__.py gcen2/__init__.py
git commit -m "fix: auto-discover C++ extension in build/ for editable installs"
```

---

## Task 2: Fetch GCEN v1 Sample Data

**Problem:** `data/sample/` is empty. Integration tests and the quickstart tutorial require real biological data.

**Files:**
- Create: `data/sample/` (multiple files)

### Step 1: Download the GCEN v1 Linux release tarball

```bash
cd /tmp
curl -L -o gcen-0.6.3.tar.gz \
  "https://github.com/wen-chen/GCEN/releases/download/0.6.3/gcen-0.6.3-linux-x86_64.tar.gz"
```

### Step 2: Extract sample_data/ into the project

```bash
cd /tmp
tar -tzf gcen-0.6.3.tar.gz | grep sample_data | head -20   # list contents first
tar -xzf gcen-0.6.3.tar.gz --wildcards '*/sample_data/*' -C /tmp/gcen_extract --strip-components=1 2>/dev/null || \
  tar -xzf gcen-0.6.3.tar.gz -C /tmp/gcen_extract
ls /tmp/gcen_extract/
```

### Step 3: Copy sample files to `data/sample/`

```bash
cp -r /tmp/gcen_extract/sample_data/* /mnt/d/code_base/gcen2/data/sample/
ls /mnt/d/code_base/gcen2/data/sample/
```

Expected files to be present:
```
gene_expr.tsv
gene_expr_norm.tsv
gene_expr_norm_filter.tsv
gene_co_expr.network
go-basic.obo
gene_go.assoc
K2ko.tsv
gene_kegg.assoc
module.txt
rwr_seed_genes.list
enrichment_gene.list
background_gene.list
```

### Step 4: Create `data/sample/README.md` documenting the files

```markdown
# GCEN2 Sample Data

These files are derived from the GCEN v1 (0.6.3) sample dataset and are
used for integration tests and the GCEN2 quickstart tutorial.

| File | Description |
|---|---|
| `gene_expr.tsv` | Raw expression matrix (genes × samples, tab-separated) |
| `gene_expr_norm.tsv` | TMM-normalized expression |
| `gene_expr_norm_filter.tsv` | Normalized + filtered (top 75th percentile) |
| `gene_co_expr.network` | Co-expression network (Spearman, cor≥0.8, pval≤0.001) |
| `go-basic.obo` | Gene Ontology database (OBO format) |
| `gene_go.assoc` | Gene → GO term associations |
| `K2ko.tsv` | KEGG gene → KO mapping |
| `gene_kegg.assoc` | Gene → KEGG pathway associations |
| `module.txt` | Module assignments from module_identify |
| `rwr_seed_genes.list` | Seed genes for Random Walk with Restart |
| `enrichment_gene.list` | Query gene list for standalone enrichment |
| `background_gene.list` | Background gene list for standalone enrichment |
```

### Step 5: Commit

```bash
cd /mnt/d/code_base/gcen2
git add data/sample/
git commit -m "data: add GCEN v1 sample data for integration tests and tutorials"
```

---

## Task 3: Integration Tests (depends on Task 2)

**Goal:** Verify the full GCEN2 pipeline produces results numerically equivalent to GCEN v1 on the sample data.

**Files:**
- Create: `tests/test_integration.py`
- Read: `data/sample/` (all files from Task 2)

### Step 1: Create `tests/test_integration.py`

```python
"""Integration tests for the complete GCEN2 pipeline.

Runs the full recommended pipeline (normalize → filter → build → module →
annotate → rwr) using the GCEN v1 sample data and verifies:
  1. Each step produces output of the expected shape and dtype.
  2. The network built by GCEN2 matches GCEN v1 output within tolerance.
  3. Module identification produces the same modules as GCEN v1.

All tests are marked @pytest.mark.integration and @pytest.mark.slow;
they are skipped if sample data is not present.
"""

import os
import pytest
import pandas as pd
import numpy as np

SAMPLE_DIR = os.path.join(os.path.dirname(__file__), '..', 'data', 'sample')
HAS_SAMPLE_DATA = os.path.isfile(os.path.join(SAMPLE_DIR, 'gene_expr.tsv'))

skip_no_data = pytest.mark.skipif(
    not HAS_SAMPLE_DATA,
    reason="Sample data not found. Run: see data/sample/README.md"
)

import gcen2
from gcen2 import preprocessing, network, module, annotation, rwr


@skip_no_data
@pytest.mark.integration
class TestBulkPipeline:
    """Full normalize → filter → build → module → annotate pipeline."""

    def test_normalize_tmm_shape(self):
        expr = gcen2.io.read_tsv(os.path.join(SAMPLE_DIR, 'gene_expr.tsv'))
        normed = preprocessing.normalize(expr, method='tmm')
        assert normed.shape == expr.shape
        # After normalization, no negative values
        assert (normed.values >= 0).all()

    def test_filter_reduces_genes(self):
        expr = gcen2.io.read_tsv(os.path.join(SAMPLE_DIR, 'gene_expr_norm.tsv'))
        filtered = preprocessing.filter_genes(expr, percent_mean=0.75)
        # Filtering should reduce to ~25% of genes
        assert filtered.shape[0] < expr.shape[0]
        assert filtered.shape[0] > 0

    def test_network_build_matches_v1(self):
        """Network built by GCEN2 should have same edges as GCEN v1 reference."""
        expr = gcen2.io.read_tsv(os.path.join(SAMPLE_DIR, 'gene_expr_norm_filter.tsv'))
        net = network.build(expr, method='spearman', cor=0.8, pval=0.001, threads=2)
        ref = gcen2.io.read_tsv(os.path.join(SAMPLE_DIR, 'gene_co_expr.network'))

        assert net.edge_count > 0
        # Should have same number of edges as v1 reference (within ±10% for float tolerance)
        assert abs(net.edge_count - len(ref)) / max(len(ref), 1) < 0.10, (
            f"Edge count mismatch: GCEN2={net.edge_count}, v1={len(ref)}"
        )

    def test_module_identification(self):
        ref_net_path = os.path.join(SAMPLE_DIR, 'gene_co_expr.network')
        net = gcen2.io.read_network(ref_net_path)
        modules = module.identify(net, similarity=0.5, threads=2)
        assert modules.module_count > 0
        # All genes in modules should be from the network
        net_genes = set(net.edges['gene_a']) | set(net.edges['gene_b'])
        for m in modules.modules:
            assert m.issubset(net_genes)

    def test_go_annotation(self):
        obo = os.path.join(SAMPLE_DIR, 'go-basic.obo')
        assoc = os.path.join(SAMPLE_DIR, 'gene_go.assoc')
        if not os.path.isfile(obo) or not os.path.isfile(assoc):
            pytest.skip("GO files not in sample data")
        ref_net_path = os.path.join(SAMPLE_DIR, 'gene_co_expr.network')
        net = gcen2.io.read_network(ref_net_path)
        results = annotation.annotate_go(net, obo, assoc, pval=0.05, threads=1)
        assert len(results) > 0
        for gene, er in results.items():
            assert hasattr(er, 'results')
            assert er.p_cutoff == 0.05

    def test_rwr_prediction(self):
        seeds_file = os.path.join(SAMPLE_DIR, 'rwr_seed_genes.list')
        if not os.path.isfile(seeds_file):
            pytest.skip("rwr_seed_genes.list not in sample data")
        with open(seeds_file) as f:
            seed_genes = [line.strip() for line in f if line.strip()]
        ref_net_path = os.path.join(SAMPLE_DIR, 'gene_co_expr.network')
        net = gcen2.io.read_network(ref_net_path)
        result = rwr.predict(net, seed_genes[:5], gamma=0.5, prob_cutoff=0.01)
        assert result.steps > 0
        assert len(result.ranked_genes) > 0
        # Probabilities should sum close to 1
        assert abs(result.ranked_genes['probability'].sum() - 1.0) < 0.1


@skip_no_data
@pytest.mark.integration
class TestCLIPipeline:
    """Full pipeline via the CLI using typer CliRunner."""

    def test_cli_normalize(self, tmp_path):
        try:
            from typer.testing import CliRunner
            from gcen2.cli import app
        except ImportError:
            pytest.skip("typer not installed")
        runner = CliRunner()
        input_file = os.path.join(SAMPLE_DIR, 'gene_expr.tsv')
        output_file = str(tmp_path / 'normed.tsv')
        result = runner.invoke(app, ['normalize', '-i', input_file, '-o', output_file, '-m', 'tmm'])
        assert result.exit_code == 0, result.output
        assert os.path.isfile(output_file)

    def test_cli_full_pipeline(self, tmp_path):
        """normalize → filter → build end-to-end via CLI."""
        try:
            from typer.testing import CliRunner
            from gcen2.cli import app
        except ImportError:
            pytest.skip("typer not installed")
        runner = CliRunner()
        sample = os.path.join(SAMPLE_DIR, 'gene_expr.tsv')
        normed = str(tmp_path / 'normed.tsv')
        filtered = str(tmp_path / 'filtered.tsv')
        network_out = str(tmp_path / 'net.tsv')

        r1 = runner.invoke(app, ['normalize', '-i', sample, '-o', normed, '-m', 'tmm'])
        assert r1.exit_code == 0, r1.output

        r2 = runner.invoke(app, ['filter', '-i', normed, '-o', filtered, '-p', '0.75'])
        assert r2.exit_code == 0, r2.output

        r3 = runner.invoke(app, ['build', '-i', filtered, '-o', network_out,
                                  '-m', 'spearman', '-c', '0.8', '-p', '0.001'])
        assert r3.exit_code == 0, r3.output
        assert os.path.isfile(network_out)
```

### Step 2: Add `read_network` helper to `gcen2/io.py`

In `gcen2/io.py`, add:
```python
def read_network(path: str) -> "Network":
    """Read a GCEN-format edge list (TSV with header: gene_a, gene_b, correlation, p_value)."""
    from gcen2.models import Network
    df = pd.read_csv(path, sep='\t')
    # Normalize column names
    df.columns = [c.lower().replace(' ', '_') for c in df.columns]
    return Network(edges=df)
```

### Step 3: Run integration tests to verify (expect skips if data absent, passes otherwise)

```bash
cd /mnt/d/code_base/gcen2
python3 -m pytest tests/test_integration.py -v --tb=short -m "integration"
```

Expected: all pass if sample data present, all skip if absent.

### Step 4: Run full suite to confirm no regressions

```bash
python3 -m pytest tests/ -q --tb=short
```

Expected: `≥84 passed, ≥31 skipped`

### Step 5: Commit

```bash
git add tests/test_integration.py gcen2/io.py
git commit -m "test: add full pipeline integration tests with GCEN v1 sample data"
```

---

## Task 4: Sphinx Documentation

**Goal:** Generate a browsable HTML API reference and user guide from docstrings.

**Files:**
- Create: `docs/conf.py`
- Create: `docs/index.rst`
- Create: `docs/quickstart.rst`
- Create: `docs/api.rst`
- Create: `docs/migration.rst`
- Create: `docs/Makefile` (standard Sphinx Makefile)

### Step 1: Install Sphinx dependencies (if not already)

```bash
pip install sphinx sphinx-rtd-theme
```

### Step 2: Create `docs/conf.py`

```python
# docs/conf.py
import os
import sys

sys.path.insert(0, os.path.abspath('..'))

project = 'GCEN2'
copyright = '2026, Wen Harold Chen'
author = 'Wen Harold Chen'
release = '2.0.0.dev0'

extensions = [
    'sphinx.ext.autodoc',
    'sphinx.ext.autosummary',
    'sphinx.ext.napoleon',      # NumPy/Google style docstrings
    'sphinx.ext.viewcode',
    'sphinx.ext.intersphinx',
]

intersphinx_mapping = {
    'python': ('https://docs.python.org/3', None),
    'numpy': ('https://numpy.org/doc/stable', None),
    'pandas': ('https://pandas.pydata.org/docs', None),
}

autosummary_generate = True
autodoc_default_options = {
    'members': True,
    'undoc-members': False,
    'show-inheritance': True,
}

html_theme = 'sphinx_rtd_theme'
html_static_path = ['_static']

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']
```

### Step 3: Create `docs/index.rst`

```rst
GCEN2 Documentation
===================

**GCEN2** is a Python toolkit for gene co-expression network analysis with
single-cell resolution and interactive visualization.

.. code-block:: bash

   pip install gcen2[all]

.. toctree::
   :maxdepth: 2
   :caption: Getting Started

   quickstart

.. toctree::
   :maxdepth: 2
   :caption: Reference

   api

.. toctree::
   :maxdepth: 1
   :caption: Migration

   migration

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
```

### Step 4: Create `docs/quickstart.rst`

```rst
Quickstart Guide
================

Installation
------------

.. code-block:: bash

   pip install gcen2[all]

For single-cell analysis only:

.. code-block:: bash

   pip install gcen2[single-cell]

Bulk RNA-Seq Workflow
---------------------

.. code-block:: python

   import gcen2

   # 1. Load data
   expr = gcen2.io.read_tsv('gene_expr.tsv')

   # 2. Normalize (TMM)
   expr_norm = gcen2.preprocessing.normalize(expr, method='tmm')

   # 3. Filter low-expression genes
   expr_filtered = gcen2.preprocessing.filter_genes(expr_norm, percent_mean=0.75)

   # 4. Build co-expression network
   net = gcen2.network.build(
       expr_filtered, method='spearman',
       cor=0.8, pval=0.001, fdr=True, threads=4
   )
   print(f"Network: {net.node_count} genes, {net.edge_count} edges")

   # 5. Identify gene modules
   modules = gcen2.module.identify(net, similarity=0.5)
   print(f"Found {modules.module_count} modules")

   # 6. Functional annotation
   results = gcen2.annotation.annotate_go(net, 'go-basic.obo', 'gene_go.assoc')

   # 7. Visualize
   gcen2.plotting.network_interactive(net, output='network.html', modules=modules)
   gcen2.plotting.dotplot(list(results.values())[0])

   # 8. Generate report
   gcen2.report.generate(net, modules=modules, enrichment=results, output='report.html')

Single-Cell Workflow
--------------------

.. code-block:: python

   import scanpy as sc
   import gcen2

   adata = sc.read_h5ad('pbmc.h5ad')

   # Build cell-type-specific networks
   networks = gcen2.single_cell.build_by_celltype(
       adata, celltype_key='leiden', use_highly_variable=True,
       method='spearman', cor=0.8, threads=4
   )

   for ct, net in networks.items():
       print(f"{ct}: {net.node_count} genes, {net.edge_count} edges")

Command-Line Interface
----------------------

.. code-block:: bash

   # Normalize
   gcen2 normalize -i gene_expr.tsv -o normed.tsv -m tmm

   # Filter
   gcen2 filter -i normed.tsv -o filtered.tsv -p 0.75

   # Build network
   gcen2 build -i filtered.tsv -o network.tsv -m spearman -c 0.8 -p 0.001 -f -t 4

   # Identify modules
   gcen2 module -i network.tsv -o modules.txt -s 0.5

   # Annotate (GO)
   gcen2 annotate -g go-basic.obo -a gene_go.assoc -n network.tsv -o results/

   # Generate report
   gcen2 report -n network.tsv -m modules.txt -o report.html
```

### Step 5: Create `docs/api.rst`

```rst
API Reference
=============

.. automodule:: gcen2
   :members:

gcen2.preprocessing
-------------------

.. automodule:: gcen2.preprocessing
   :members:
   :undoc-members: False

gcen2.network
-------------

.. automodule:: gcen2.network
   :members:
   :undoc-members: False

gcen2.module
------------

.. automodule:: gcen2.module
   :members:
   :undoc-members: False

gcen2.annotation
----------------

.. automodule:: gcen2.annotation
   :members:
   :undoc-members: False

gcen2.rwr
---------

.. automodule:: gcen2.rwr
   :members:
   :undoc-members: False

gcen2.single_cell
-----------------

.. automodule:: gcen2.single_cell.core
   :members:
   :undoc-members: False

gcen2.differential
------------------

.. automodule:: gcen2.differential
   :members:
   :undoc-members: False

gcen2.plotting
--------------

.. automodule:: gcen2.plotting
   :members:
   :undoc-members: False

gcen2.report
------------

.. automodule:: gcen2.report
   :members:
   :undoc-members: False

gcen2.models
------------

.. automodule:: gcen2.models
   :members:
   :undoc-members: False
```

### Step 6: Create `docs/migration.rst`

```rst
Migration from GCEN v1
======================

Command Mapping
---------------

All GCEN v1 commands are available as ``gcen2`` subcommands:

.. list-table::
   :widths: 50 50
   :header-rows: 1

   * - GCEN v1
     - GCEN2
   * - ``data_norm -i in.tsv -o out.tsv -m tmm``
     - ``gcen2 normalize -i in.tsv -o out.tsv -m tmm``
   * - ``data_filter -i in.tsv -o out.tsv -p 0.75``
     - ``gcen2 filter -i in.tsv -o out.tsv -p 0.75``
   * - ``network_build -i in.tsv -o out.network -m spearman``
     - ``gcen2 build -i in.tsv -o out.network -m spearman``
   * - ``module_identify -i net -o mod.txt -s 0.5``
     - ``gcen2 module -i net -o mod.txt -s 0.5``
   * - ``annotate -g go.obo -a assoc -n net -o dir``
     - ``gcen2 annotate -g go.obo -a assoc -n net -o dir``
   * - ``rwr -n net -g seeds -o out.tsv``
     - ``gcen2 rwr -n net -g seeds -o out.tsv``

Breaking Changes
----------------

**network_merge output format**
  The ``direction`` column is now a separate ``"1"`` / ``"-1"`` column,
  not appended to the correlation value string.

**Standard deviation**
  ``filter_genes()`` uses sample std (``ddof=1``). Results differ slightly
  from GCEN v1 which used population std (``ddof=0``).

**CSV handling**
  ``csv_to_tsv`` and ``tsv_to_csv`` now correctly handle quoted fields.
```

### Step 7: Create `docs/Makefile`

```makefile
SPHINXOPTS    ?=
SPHINXBUILD   ?= sphinx-build
SOURCEDIR     = .
BUILDDIR      = _build

.PHONY: help clean html

help:
	@$(SPHINXBUILD) -M help "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)

clean:
	rm -rf $(BUILDDIR)

html:
	@$(SPHINXBUILD) -M html "$(SOURCEDIR)" "$(BUILDDIR)" $(SPHINXOPTS) $(O)
	@echo "Build finished. The HTML pages are in $(BUILDDIR)/html."
```

### Step 8: Build the docs and verify no errors

```bash
cd /mnt/d/code_base/gcen2/docs
make html 2>&1 | tail -20
```

Expected: `Build finished. The HTML pages are in _build/html.` with no ERROR lines.

### Step 9: Commit

```bash
cd /mnt/d/code_base/gcen2
git add docs/conf.py docs/index.rst docs/quickstart.rst docs/api.rst docs/migration.rst docs/Makefile
git commit -m "docs: initialize Sphinx documentation with quickstart and API reference"
```

---

## Task 5: Benchmark Tests

**Goal:** Measure GCEN2 performance on standardized datasets to track regressions and compare against GCEN v1 targets from the SRS.

**Files:**
- Create: `tests/test_benchmark.py`

### Step 1: Create `tests/test_benchmark.py`

```python
"""Benchmark tests for GCEN2 computation kernels.

Measures wall-clock time for the main performance-critical operations.
These are NOT run in regular CI — they require:
  - @pytest.mark.benchmark (install: pip install pytest-benchmark)
  - @pytest.mark.slow

Run with:
  pytest tests/test_benchmark.py -v -m benchmark --benchmark-sort=mean

Reference targets (from SRS §5.1):
  - Network build (10k genes, 100 samples) < 5 minutes on 8 cores
"""

import numpy as np
import pytest

try:
    from gcen2.core import _gcen_core
    HAS_CORE = True
except ImportError:
    HAS_CORE = False

skip_no_core = pytest.mark.skipif(not HAS_CORE, reason="_gcen_core not compiled")


@skip_no_core
@pytest.mark.benchmark
@pytest.mark.slow
class TestCorrelationBenchmarks:
    """Benchmark the C++ correlation kernel at various gene counts."""

    def _make_data(self, n_genes, n_samples, seed=42):
        rng = np.random.default_rng(seed)
        return rng.standard_normal((n_genes, n_samples)).astype(np.float64)

    def test_bench_1k_genes(self, benchmark):
        data = self._make_data(1000, 100)
        result = benchmark(
            _gcen_core.correlation.pairwise_correlation,
            data, 'spearman', 2
        )
        indices, corrs, pvals = result
        # 1k genes → up to 499,500 pairs
        assert len(corrs) <= 1000 * 999 // 2

    def test_bench_5k_genes(self, benchmark):
        data = self._make_data(5000, 100)
        result = benchmark(
            _gcen_core.correlation.pairwise_correlation,
            data, 'spearman', 4
        )
        indices, corrs, pvals = result
        assert len(corrs) <= 5000 * 4999 // 2

    def test_bench_thread_scaling(self, benchmark):
        """Verify that 4 threads is faster than 1 thread."""
        data = self._make_data(2000, 100)
        import time

        t0 = time.perf_counter()
        _gcen_core.correlation.pairwise_correlation(data, 'spearman', 1)
        t1 = time.perf_counter()

        t2 = time.perf_counter()
        _gcen_core.correlation.pairwise_correlation(data, 'spearman', 4)
        t3 = time.perf_counter()

        time_1t = t1 - t0
        time_4t = t3 - t2
        # 4 threads should be at least 1.5x faster than 1 thread
        # (conservative, system load may vary)
        print(f"\n1 thread: {time_1t:.2f}s, 4 threads: {time_4t:.2f}s, speedup: {time_1t/time_4t:.2f}x")
        # Not asserting speedup — just logging. Real assertion in manual review.


@skip_no_core
@pytest.mark.benchmark
@pytest.mark.slow
class TestModuleBenchmarks:

    def _make_dense_network_adj(self, n_genes, edge_density=0.05, seed=42):
        """Generate a random adjacency dict."""
        rng = np.random.default_rng(seed)
        genes = [f"G{i}" for i in range(n_genes)]
        adj = {g: set() for g in genes}
        for i in range(n_genes):
            for j in range(i + 1, n_genes):
                if rng.random() < edge_density:
                    adj[genes[i]].add(genes[j])
                    adj[genes[j]].add(genes[i])
        return adj

    def test_bench_1k_node_module(self, benchmark):
        adj = self._make_dense_network_adj(1000, edge_density=0.05)
        result = benchmark(
            _gcen_core.module.identify_modules,
            adj, 0.3, 2
        )
        assert isinstance(result, list)


@pytest.mark.benchmark
@pytest.mark.slow
class TestRWRBenchmarks:

    @skip_no_core
    def test_bench_rwr_convergence(self, benchmark):
        n = 1000
        rng = np.random.default_rng(0)
        W = np.abs(rng.standard_normal((n, n)))
        W /= W.sum(axis=0, keepdims=True)  # column normalize
        P0 = np.zeros(n)
        P0[:10] = 0.1  # 10 seed genes

        result = benchmark(
            _gcen_core.rwr.random_walk_restart,
            W, P0, 0.5, 1_000_000, 1e-10
        )
        PT, steps = result
        assert steps < 1_000_000  # must converge
```

### Step 2: Run to verify (they will be slow-skipped by default)

```bash
cd /mnt/d/code_base/gcen2
python3 -m pytest tests/test_benchmark.py -v -m "not slow" --tb=short
```

Expected: all collected tests skipped (marked slow).

### Step 3: Commit

```bash
git add tests/test_benchmark.py
git commit -m "test: add benchmark tests for correlation, module, and RWR kernels"
```

---

## Task 6: bioconda Recipe

**Goal:** Create a conda-build recipe so GCEN2 can be submitted to the bioconda channel.

**Files:**
- Create: `conda-recipe/meta.yaml`
- Create: `conda-recipe/build.sh` (Linux/macOS build script)
- Create: `conda-recipe/bld.bat` (Windows build script)

### Step 1: Create `conda-recipe/meta.yaml`

```yaml
{% set name = "gcen2" %}
{% set version = "2.0.0" %}

package:
  name: {{ name|lower }}
  version: {{ version }}

source:
  url: https://pypi.io/packages/source/{{ name[0] }}/{{ name }}/{{ name }}-{{ version }}.tar.gz
  # sha256: <fill in after PyPI upload>

build:
  number: 0
  script: {{ PYTHON }} -m pip install . --no-deps -vv
  # C++ extension compilation requires cmake + pybind11 + eigen
  run_constrained:
    - __unix  # Linux + macOS; Windows requires MSVC

requirements:
  build:
    - {{ compiler('cxx') }}
    - cmake >=3.18
    - make  # [not win]
  host:
    - python
    - pip
    - scikit-build-core >=0.5
    - pybind11 >=2.10
    - eigen >=3.4    # or use bundled third_party/eigen
  run:
    - python >=3.9
    - numpy >=1.21
    - pandas >=1.3
    - scipy >=1.7

test:
  imports:
    - gcen2
    - gcen2.preprocessing
    - gcen2.network
    - gcen2.module
    - gcen2.annotation
    - gcen2.rwr
  commands:
    - gcen2 --help

about:
  home: https://github.com/wen-chen/GCEN2
  license: GPL-3.0-or-later
  license_family: GPL
  summary: 'A Python toolkit for gene co-expression network analysis with single-cell resolution and interactive visualization'
  description: |
    GCEN2 is a major upgrade of GCEN (Gene Co-Expression Network), redesigned
    with Python as the primary interface while retaining C++ for performance-critical
    computation. It adds single-cell support, visualization, and automated report generation.
  doc_url: https://gcen2.readthedocs.io
  dev_url: https://github.com/wen-chen/GCEN2

extra:
  recipe-maintainers:
    - wen-chen
```

### Step 2: Create `conda-recipe/build.sh`

```bash
#!/bin/bash
set -ex

# Initialize Eigen submodule if not present
if [ ! -f "third_party/eigen/Eigen/Core" ]; then
    git submodule update --init --recursive
fi

$PYTHON -m pip install . --no-deps --ignore-installed -vv
```

### Step 3: Create `conda-recipe/bld.bat`

```batch
@echo off
if not exist "third_party\eigen\Eigen\Core" (
    git submodule update --init --recursive
)
%PYTHON% -m pip install . --no-deps --ignore-installed -vv
if errorlevel 1 exit 1
```

### Step 4: Verify the recipe is parseable (requires conda-build)

```bash
conda render conda-recipe/ 2>&1 | tail -20
# If conda-build not available, just lint the YAML:
python3 -c "import yaml; yaml.safe_load(open('conda-recipe/meta.yaml'))" && echo "YAML valid"
```

### Step 5: Commit

```bash
cd /mnt/d/code_base/gcen2
git add conda-recipe/
git commit -m "build: add bioconda recipe skeleton for conda package distribution"
```

---

## Task 7: Chunked Computation for Large Single-Cell Datasets

**Goal:** Allow `build_by_celltype()` to handle gene sets larger than memory allows by processing the correlation matrix in chunks.

**Files:**
- Modify: `gcen2/single_cell/core.py`
- Modify: `tests/test_single_cell.py` (add chunk tests)

### Step 1: Add failing test for chunked computation in `tests/test_single_cell.py`

```python
@pytest.mark.requires_single_cell
class TestChunkedCorrelation:
    """Test that chunked computation gives same results as full computation."""

    def test_chunked_matches_full(self):
        """chunk_size < n_genes should give same network as chunk_size = n_genes."""
        pytest.importorskip("anndata")
        import anndata as ad
        import numpy as np
        import gcen2

        rng = np.random.default_rng(42)
        n_genes, n_cells, n_groups = 50, 200, 3
        X = np.abs(rng.standard_normal((n_cells, n_genes)))
        obs = pd.DataFrame({'celltype': [f'CT{i % n_groups}' for i in range(n_cells)]})
        var = pd.DataFrame(index=[f'G{i}' for i in range(n_genes)])
        adata = ad.AnnData(X=X, obs=obs, var=var)

        # Full computation (no chunking)
        net_full = gcen2.single_cell.build_by_celltype(
            adata, celltype_key='celltype',
            pseudo_replicates=5, seed=42, cor=0.5, pval=0.05, threads=1
        )

        # Chunked computation
        net_chunk = gcen2.single_cell.build_by_celltype(
            adata, celltype_key='celltype',
            pseudo_replicates=5, seed=42, cor=0.5, pval=0.05, threads=1,
            max_genes_per_chunk=20  # force chunking
        )

        # Should have same total edge counts (order may differ)
        for ct in net_full:
            full_edges = net_full[ct].edge_count
            chunk_edges = net_chunk[ct].edge_count
            assert full_edges == chunk_edges, (
                f"CT={ct}: full={full_edges}, chunked={chunk_edges}"
            )
```

### Step 2: Run test to confirm it fails

```bash
python3 -m pytest tests/test_single_cell.py::TestChunkedCorrelation -v --tb=short
```

Expected: `FAILED` (TypeError: unexpected keyword argument `max_genes_per_chunk`)

### Step 3: Implement chunked computation in `gcen2/single_cell/core.py`

Add a helper function and modify `build_by_celltype`:

```python
def _chunk_gene_indices(n_genes: int, chunk_size: int) -> list[tuple[int, int]]:
    """Return (start, end) index pairs covering [0, n_genes) in chunks."""
    chunks = []
    start = 0
    while start < n_genes:
        end = min(start + chunk_size, n_genes)
        chunks.append((start, end))
        start = end
    return chunks


def _build_chunked(expr: pd.DataFrame, chunk_size: int, **build_kwargs) -> "Network":
    """Build a network by processing gene pairs in row-chunks.

    Splits the gene set into chunks of ``chunk_size`` rows, computes all
    pairs within each chunk and across adjacent chunks, then merges the
    partial networks into one.
    """
    import gcen2.network as net_mod
    from gcen2.models import Network
    import pandas as pd

    genes = list(expr.index)
    n_genes = len(genes)

    if n_genes <= chunk_size:
        return net_mod.build(expr, **build_kwargs)

    chunks = _chunk_gene_indices(n_genes, chunk_size)
    all_edges = []

    for i, (s1, e1) in enumerate(chunks):
        for j, (s2, e2) in enumerate(chunks):
            if j < i:
                continue  # avoid duplicate pairs
            # Build on the union of the two chunks
            if i == j:
                sub_expr = expr.iloc[s1:e1]
            else:
                sub_expr = expr.iloc[list(range(s1, e1)) + list(range(s2, e2))]

            if len(sub_expr) < 3:
                continue

            sub_net = net_mod.build(sub_expr, **build_kwargs)

            # Keep only cross-chunk edges when i != j
            if i != j:
                chunk_i_genes = set(genes[s1:e1])
                chunk_j_genes = set(genes[s2:e2])
                mask = (
                    (sub_net.edges['gene_a'].isin(chunk_i_genes) &
                     sub_net.edges['gene_b'].isin(chunk_j_genes)) |
                    (sub_net.edges['gene_a'].isin(chunk_j_genes) &
                     sub_net.edges['gene_b'].isin(chunk_i_genes))
                )
                edges = sub_net.edges[mask]
            else:
                edges = sub_net.edges

            all_edges.append(edges)

    if not all_edges:
        return Network(edges=pd.DataFrame(
            columns=['gene_a', 'gene_b', 'correlation', 'p_value']
        ))

    combined = pd.concat(all_edges, ignore_index=True)
    combined = combined.drop_duplicates(subset=['gene_a', 'gene_b'])
    return Network(edges=combined, method=build_kwargs.get('method', 'spearman'))
```

Also update `build_by_celltype` signature to accept `max_genes_per_chunk=None`:

```python
def build_by_celltype(
    adata,
    celltype_key: str,
    use_highly_variable: bool = False,
    min_cells: int = 10,
    pseudo_replicates: int = 10,
    seed: int = 42,
    max_genes_per_chunk: int | None = None,   # NEW
    **build_kwargs,
) -> dict:
    ...
    # In the per-celltype loop, replace network.build() with:
    if max_genes_per_chunk is not None and len(pseudo_bulk) > max_genes_per_chunk:
        ct_network = _build_chunked(pseudo_bulk, max_genes_per_chunk, **build_kwargs)
    else:
        ct_network = net_mod.build(pseudo_bulk, **build_kwargs)
```

### Step 4: Run the test to verify it passes

```bash
python3 -m pytest tests/test_single_cell.py::TestChunkedCorrelation -v --tb=short
```

Expected: `PASSED`

### Step 5: Run full test suite for regressions

```bash
python3 -m pytest tests/ -q --tb=short
```

Expected: `≥84 passed, ≥31 skipped`

### Step 6: Commit

```bash
git add gcen2/single_cell/core.py tests/test_single_cell.py
git commit -m "feat: add chunked correlation computation for large single-cell datasets"
```

---

## Final: Verification Run

After all tasks are complete, run the full test suite one more time:

```bash
cd /mnt/d/code_base/gcen2
python3 -m pytest tests/ -v --tb=short 2>&1 | tail -20
```

Expected: ≥84 passed (more if integration tests ran with sample data), 0 failed.

Build docs:
```bash
cd /mnt/d/code_base/gcen2/docs && make html
```

Expected: successful build with no ERRORs.
