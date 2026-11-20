Quickstart Guide
================

Installation
------------

GCEN2 is distributed as source code (https://github.com/wen-chen/GCEN2), on
PyPI and on Bioconda. Supported platforms are Linux x86-64 (AMD64), Windows
x86-64 (AMD64) and macOS on Apple Silicon. PyPI has pre-built wheels for all
three; Bioconda covers Linux and macOS (it does not build for Windows):

.. code-block:: bash

   conda install -c conda-forge -c bioconda gcen2

Install GCEN2 with all optional dependencies from PyPI:

.. code-block:: bash

   pip install gcen2[all]

For specific feature sets:

.. code-block:: bash

   pip install gcen2                  # core: C++ kernels, Python API, gcen2 command
   pip install gcen2[single-cell]     # + AnnData/scanpy support
   pip install gcen2[visualization]   # + matplotlib/pyvis plots
   pip install gcen2[report]          # + Jinja2/weasyprint reports
   pip install gcen2[all]             # everything

Tested with Python 3.9 to 3.14; the oldest versions tested are numpy 1.24,
pandas 2.0, scipy 1.10, anndata 0.10, matplotlib 3.7 and typer 0.9.

Building from Source
--------------------

.. code-block:: bash

   git clone https://github.com/wen-chen/GCEN2
   cd GCEN2
   pip install -e ".[dev]" --no-build-isolation

Requirements: Python ≥ 3.9 and a C++17 compiler. Tested with GCC 11, 14 and 15,
Clang 19, Apple Clang 17 and MSVC 19.44 (Visual Studio 2022); older C++17
compilers may work but are not tested.

Example data
------------

The repository's ``data/example/`` holds the GCEN v1 sample data (``bulk/``), a PBMC
subset (``single_cell/pbmc_subset.h5ad``) and a Visium mouse brain section
(``spatial/mouse_brain_anterior_subset.h5ad``); ``data/example/README.md`` lists the
commands for each. The examples below use the files of ``data/example/bulk/``
(``go-basic.obo.gz`` and ``gene_go.assoc.gz`` are read as they are).

Bulk RNA-Seq Workflow
---------------------

.. code-block:: python

   import gcen2

   # 1. Load data (TSV: genes × samples)
   expr = gcen2.io.read_tsv('gene_expr.tsv')

   # 2. Normalize using TMM
   expr_norm = gcen2.preprocessing.normalize(expr, method='tmm')

   # 3. Filter low-expression genes (keep top 75th percentile)
   expr_filt = gcen2.preprocessing.filter_genes(expr_norm, percent_mean=0.75)

   # 4. Build co-expression network (Spearman, with FDR)
   net = gcen2.network.build(
       expr_filt, method='spearman',
       cor=0.8, pval=0.001, fdr=True, threads=4
   )
   print(f"Network: {net.node_count} genes, {net.edge_count} edges")

   # 5. Identify gene modules (Jaccard + DFS)
   modules = gcen2.module.identify(net, similarity=0.5, threads=4)
   print(f"Modules: {modules.module_count}")

   # 6. GO enrichment annotation
   results = gcen2.annotation.annotate_go(
       net, 'go-basic.obo', 'gene_go.assoc', pval=0.05
   )

   # 7. Visualize
   gcen2.plotting.network_interactive(net, output='network.html', modules=modules)
   gcen2.plotting.dotplot(list(results.values())[0], output='dotplot.pdf')

   # 8. Gene function prediction via RWR
   rwr_result = gcen2.rwr.predict(net, seed_genes=['GENE_A', 'GENE_B'], gamma=0.5)
   print(rwr_result.ranked_genes.head(10))

   # 9. Generate HTML report
   gcen2.report.generate(
       net, modules=modules, enrichment=results,
       rwr=rwr_result, output='report.html'
   )

Single-Cell Workflow
--------------------

.. code-block:: python

   import scanpy as sc
   import gcen2

   # Load AnnData (e.g., PBMC dataset)
   adata = sc.read_h5ad('pbmc.h5ad')

   # Standard scanpy preprocessing
   sc.pp.normalize_total(adata)
   sc.pp.log1p(adata)
   sc.pp.highly_variable_genes(adata, n_top_genes=3000)

   # Build cell-type-specific co-expression networks; within each cell type,
   # similar cells are first averaged into kNN metacells (strategy='metacell')
   networks = gcen2.single_cell.build_by_celltype(
       adata,
       celltype_key='leiden',           # obs column with cell type labels
       use_highly_variable=True,        # restrict to HVGs
       method='spearman',
       cor=0.8,
       threads=4
   )

   for ct, net in networks.items():
       print(f"{ct}: {net.node_count} genes, {net.edge_count} edges")

   # Differential co-expression between two cell types (Fisher's z on the
   # networks' effective sample sizes, as `gcen2 celltype --diff`)
   diff = gcen2.differential.compare(
       networks['T_cell'], networks['B_cell']
   )
   gcen2.plotting.diff_network(diff, output='diff_network.pdf')

Command-Line Interface
----------------------

All GCEN v1 commands are available as ``gcen2`` subcommands:

.. code-block:: bash

   # Normalize expression data
   gcen2 normalize -i gene_expr.tsv -o normed.tsv -m tmm

   # Filter low-expression genes
   gcen2 filter -i normed.tsv -o filtered.tsv -p 0.75

   # Build co-expression network
   gcen2 build -i filtered.tsv -o network.tsv -m spearman -c 0.8 -p 0.001 -f -t 4

   # Identify gene modules
   gcen2 module -i network.tsv -o modules.txt -s 0.5 -t 4

   # GO annotation
   gcen2 annotate -g go-basic.obo -a gene_go.assoc -n network.tsv -o results/

   # Gene function prediction
   gcen2 rwr -n network.tsv -g seed_genes.list -o rwr_result.tsv

   # Generate report
   gcen2 report -n network.tsv -m modules.txt -o report.html

   # Network statistics
   gcen2 stats -i network.tsv
