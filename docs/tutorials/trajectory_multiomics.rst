Trajectories and multi-omics
============================

Two analyses beyond single networks: how co-expression changes along a
trajectory, and how several networks (co-expression, protein interactions,
regulation) are combined. ``benchmarks/`` compares them with
RandomWalkRestartMH (the multiplex random walk) and with R implementations of
the trajectory networks and the fusion, which have no published package;
the results agree (``benchmarks/README.md``).

Networks along pseudotime
-------------------------

``gcen2 trajectory`` cuts the cells, ordered by pseudotime, into overlapping
windows of equal size, builds a network in each window and follows every edge
across them. It needs a pseudotime: an ``obs`` column of an ``.h5ad`` (e.g.
scanpy's ``dpt_pseudotime``), or a two-column TSV of sample and pseudotime
for an expression TSV (a bulk time course).

The example cells (``data/example/single_cell/pbmc_subset.h5ad``) are not a
developmental trajectory, but a diffusion pseudotime shows the commands:

.. code-block:: python

   import anndata as ad
   import numpy as np
   import scanpy as sc

   adata = ad.read_h5ad("data/example/single_cell/pbmc_subset.h5ad")
   x = adata.copy()
   sc.pp.normalize_total(x); sc.pp.log1p(x); sc.pp.pca(x, n_comps=20); sc.pp.neighbors(x)
   x.uns["iroot"] = int(np.argmin(x.obsm["X_pca"][:, 0]))
   sc.tl.diffmap(x); sc.tl.dpt(x)
   adata.obs["dpt_pseudotime"] = x.obs["dpt_pseudotime"]
   adata.write_h5ad("cells_pt.h5ad")

.. code-block:: bash

   gcen2 trajectory -i cells_pt.h5ad --pseudotime dpt_pseudotime -o trajectory/ \
       -w 0.2 --metacell 5 -l log2 -c 0.7 -t 8 --plot

``-w 0.2``: each window holds a fifth of the cells; windows overlap by half
(``-s``). ``--metacell 5`` averages five consecutive cells first, which
removes much of the dropout noise. Cells without a finite pseudotime are left
out with a warning. ``trajectory/`` holds one network per window
(``window_<i>.network``), ``edge_dynamics.tsv`` (the correlation of every edge
in every window and its pattern: ``stable``, ``gained``, ``lost``,
``transient`` or ``intermittent``), ``degree_dynamics.tsv`` and, with
``--plot``, the overview figure (also ``gcen2 plot trajectory -i
trajectory/``).

Combining networks
------------------

``gcen2 fuse`` scales the edge weights of each layer (``--normalize max``,
``rank`` or ``none``) and scores every edge by the weighted mean over the
layers, recording which layers support it; ``gcen2 rwr-multiplex`` ranks genes
by a random walk across all layers at once (RWR-M, Valdeolivas et al. 2019;
the same scores as their RandomWalkRestartMH, ``--transition layer`` for a
jump probability ``--delta`` that does not depend on the degree).
A layer is any edge-list file: a GCEN network, a STRING-style scored
interaction file, or a two-column list of regulator-target pairs.

The example data hold no protein-interaction layer, so the example combines
the Spearman and Pearson networks of the bulk data
(``data/example/bulk``); with your own data, replace the second layer by,
e.g., ``string_ppi.tsv``:

.. code-block:: bash

   B=/path/to/GCEN2/data/example/bulk
   gcen2 build -i $B/gene_expr.tsv -o spearman.network -m spearman -c 0.9 -t 8
   gcen2 build -i $B/gene_expr.tsv -o pearson.network -m pearson -c 0.9 -t 8
   gcen2 fuse -i spearman.network -i pearson.network --name spearman --name pearson \
       --min-layers 2 -o fused.tsv -t 8
   gcen2 rwr-multiplex -n spearman.network -n pearson.network \
       -g $B/rwr_seed_genes.list -o rwr_multiplex.tsv
   gcen2 plot fusion -i fused.tsv -o fig/fusion

``--min-layers 2`` keeps the edges found in both layers; ``-w`` gives each
layer a weight (default 1). ``rwr_multiplex.tsv`` has the combined score of
every gene and its score in each layer.

The same in Python
------------------

.. code-block:: python

   import anndata as ad
   import gcen2

   traj = gcen2.trajectory.build_trajectory(
       ad.read_h5ad("cells_pt.h5ad"), "dpt_pseudotime", window_size=0.2,
       metacell_size=5, log="log2", cor=0.7, threads=8)
   traj.dynamics.query("pattern == 'gained'")       # edges that appear late
   traj.write("trajectory/")

   layers = {"spearman": "spearman.network", "pearson": "pearson.network"}
   fused = gcen2.multilayer.fuse(layers, min_layers=2)
   seeds = gcen2.read_gene_list("data/example/bulk/rwr_seed_genes.list")
   ranked = gcen2.multilayer.rwr_multiplex(layers, seeds).ranked_genes
