Single cell (hdWGCNA workflow)
==============================

This tutorial runs hdWGCNA's workflow with GCEN2 on
``data/example/single_cell/pbmc_subset.h5ad``: 1,800 PBMC cells of three
cell types (CD14 monocytes, B and NK cells), raw counts, with a PCA and a
UMAP. The whole tutorial takes about 3 minutes.

.. code-block:: bash

   SC=/path/to/GCEN2/data/example/single_cell/pbmc_subset.h5ad

GCEN2 reads AnnData (``.h5ad``) files; ``X`` (or ``--layer``) must hold the
counts of each cell. The neighbours for the metacells are found in
``obsm["X_pca"]`` (``--use-rep``), or in a PCA that GCEN2 computes.

Modules per cell type
---------------------

.. code-block:: bash

   gcen2 celltype-modules -i $SC --groupby cell_type -o modules/ -t 8 \
         --write-metacells --write-tom --dme-by cluster --trait total_counts

For every cell type this:

1. aggregates its cells into metacells (each a cell and its 24 nearest
   neighbours, sharing at most 10 cells with another metacell; hdWGCNA's
   ``MetacellsByGroups``), and normalises them;
2. keeps the genes detected in at least 5 % of the cells;
3. finds the modules as hdWGCNA's ``ConstructNetwork`` does (signed network,
   soft power by the bicor scale-free fit, deepSplit 4);
4. computes the eigengene of every module in every cell and each gene's kME
   (``ModuleEigengenes``, ``ModuleConnectivity``), and moves genes with a
   negative kME out of their module;
5. here also: differential module eigengenes between the clusters of the cell
   type (``--dme-by``, hdWGCNA ``FindAllDMEs``) and the correlation of the
   eigengenes with a cell trait (``--trait``).

.. code-block:: text

   CD14 Mono: 800 cells, 216 metacells, 6736 genes, power 4, 9 modules (1915 genes; 31 reassigned)
   NK: 400 cells, 67 metacells, 6412 genes, power 5, 15 modules (3442 genes; 61 reassigned)
   B: 600 cells, 101 metacells, 6269 genes, power 8, 22 modules (3369 genes; 76 reassigned)

``modules/`` holds ``celltypes.tsv`` (one line per cell type) and a directory
per cell type with ``modules.tsv``, ``eigengenes.tsv`` and ``kme.tsv`` (of the
metacells), ``hub_genes.tsv``, ``cell_eigengenes.tsv`` and ``cell_kme.tsv``,
``dmes.tsv``, ``trait_cor.tsv``, ``metacell_members.tsv`` (the cells of each
metacell), ``soft_power.tsv``, ``dendrogram.json`` and, with the options above,
``metacells.tsv`` and ``tom_hubs.tsv``.

More options:

- ``--consensus-by sample``: one consensus network over the samples (or
  studies) of each cell type (hdWGCNA ``ConstructNetwork(consensus = TRUE)``);
- ``--harmonize-by sample``: also Harmony-integrated eigengenes;
- ``--backed``: for an atlas, read one cell type at a time from disk;
- ``--preset wgcna``: WGCNA's settings instead of hdWGCNA's.

Figures
-------

.. code-block:: bash

   gcen2 plot celltype-overview -i modules/celltypes.tsv -o fig/overview
   gcen2 plot metacells -i modules/ -c $SC -o fig/metacells
   gcen2 plot soft-power -i modules/B/soft_power.tsv -o fig/b_soft_power
   gcen2 plot dendrogram -i modules/B/dendrogram.json -o fig/b_dendrogram
   gcen2 plot me-embedding -i modules/B/cell_eigengenes.tsv -c $SC -o fig/b_eigengenes
   gcen2 plot me-dotplot -i modules/CD14_Mono/cell_eigengenes.tsv -g $SC --groupby cluster -o fig/mono_dotplot
   gcen2 plot dme-volcano -i modules/CD14_Mono/dmes.tsv -o fig/mono_dmes
   gcen2 plot module-trait -i modules/B/trait_cor.tsv -o fig/b_traits
   gcen2 plot hub-network -t modules/B/tom_hubs.tsv -k modules/B/kme.tsv -m ME1 -o fig/b_hub_me1

``me-embedding`` colours the cells of the UMAP (``--basis``, default
``X_umap``) by each module's eigengene; ``metacells`` shows where the metacells
lie among the cells. All figures in one report, also saved one by one:

.. code-block:: bash

   gcen2 report --celltype-modules modules/ --coords $SC -o report.html --figures-dir fig/report

Modules in another cell type or data set
----------------------------------------

Project the B-cell modules onto the NK cells (hdWGCNA ``ProjectModules``),
and test whether they are preserved there (WGCNA ``modulePreservation``, as
hdWGCNA runs it; 100 permutations here, 500 by default):

.. code-block:: bash

   gcen2 project-modules -i $SC -m modules/B/modules.tsv -o projected/ --groupby cell_type --group NK
   gcen2 module-preservation -r modules/B/metacells.tsv -m modules/B/modules.tsv \
         -q $SC --groupby cell_type --group NK -o preservation/ --permutations 100 -t 8
   gcen2 plot preservation -i preservation/preservation_Z.tsv -o fig/preservation

``projected/projected_eigengenes.tsv`` holds the eigengenes of the B-cell
modules in each NK cell. The preservation statistics (``preservation_Z.tsv``):
a Zsummary above 10 means a module is strongly preserved, between 2 and 10
moderately, below 2 not. Most B-cell modules are not preserved in NK cells;
ME12 (Zsummary 31.6) is.

Correlation networks per cell type
----------------------------------

GCEN's correlation networks, one per cell type (on metacells, with p-values
from the effective number of samples and the FDR over all gene pairs, ``-p``
applying to it; ``--no-fdr`` for the uncorrected p-values), and the
differences between them:

.. code-block:: bash

   gcen2 celltype -i $SC --groupby cell_type -o networks/ -c 0.5 --diff -t 8
   gcen2 plot celltype-networks -i networks/ -o fig/networks

``networks/diff/`` compares every pair of cell-type networks (Fisher's z, each
edge on the effective samples of its two genes in each cell type):
edges specific to one type, and shared edges whose correlation changed.

The same in Python
------------------

.. code-block:: python

   import anndata as ad
   import gcen2

   adata = ad.read_h5ad("data/example/single_cell/pbmc_subset.h5ad")
   results = gcen2.single_cell.modules_by_celltype(adata, "cell_type", threads=8)
   b = results["B"]
   print(b.power, b.modules()[:2])          # modules as gene lists
   me = b.metadata["cell_eigengenes"]       # cells x modules
   gcen2.plot.me_embedding(me, adata, output="fig/b_eigengenes")
