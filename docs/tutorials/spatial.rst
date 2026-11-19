Spatial transcriptomics
=======================

This tutorial uses ``data/example/spatial/mouse_brain_anterior_subset.h5ad``:
a 10x Visium section of the mouse brain (2,695 spots, the 2,000 genes detected
in the most spots, raw counts), with the Visium grid (``obs["array_row"]``,
``obs["array_col"]``), spot coordinates (``obsm["spatial"]``) and 12 regions
(``obs["region"]``, clusters of the expression). It takes under a minute.

.. code-block:: bash

   SP=/path/to/GCEN2/data/example/spatial/mouse_brain_anterior_subset.h5ad

Spatial modules (hdWGCNA workflow)
----------------------------------

.. code-block:: bash

   gcen2 spatial-modules -i $SP --groupby region -o modules/ -t 8 --plot --write-tom

Like hdWGCNA's ``MetaspotsByGroups``, neighbouring spots of the same region are
summed into metaspots on the Visium grid (here 604 metaspots; regions with
fewer than 50 spots are skipped, ``--min-spots``). The modules are found on
the metaspots (8 modules of 1,190 genes), and each module's eigengene is
computed in every spot:

.. code-block:: text

   2695 spots, 604 grid metaspots, 2000 genes, power 10: 8 modules (1190 genes; ...) → modules/

``modules/`` holds the modules, eigengenes, kME and hub genes (as ``gcen2
module -e`` writes them), ``spot_eigengenes.tsv`` and ``spot_kme.tsv``,
``group_cor.tsv`` (eigengenes against the regions) and, with ``--plot``,
``spatial_eigengenes.png``.

Without a Visium grid (Slide-seq, Visium HD bins, other platforms), use
``--metaspots knn -k 100``: each metaspot is a spot and its nearest
neighbours in ``obsm["spatial"]``.

.. code-block:: bash

   gcen2 plot me-embedding -i modules/spot_eigengenes.tsv -c $SP --basis spatial --spatial -o fig/spatial_me
   gcen2 plot module-trait -i modules/group_cor.tsv --annotate stars -o fig/regions
   gcen2 plot kme -i modules/kme.tsv -o fig/kme
   gcen2 plot hub-network-all -t modules/tom_hubs.tsv -k modules/kme.tsv -o fig/hubs

In the maps the modules follow the anatomy: one marks the olfactory bulb,
another the striatum, others the fibre tracts and the cortex. ``--spatial``
draws the spots with equal axes and the image's y axis.

Networks with spatial information
---------------------------------

``gcen2 spatial-build`` builds GCEN correlation networks in three ways:

.. code-block:: bash

   # keep edges between genes with significant positive spatial autocorrelation
   # (one-sided Moran's I test, FDR <= 0.05; --moran-fdr 1: I > 0 only)
   gcen2 spatial-build -i $SP -m filtered -o filtered.network -c 0.6 -t 8 --moran moran.tsv
   # correlations of spatially smoothed expression (Gaussian-weighted neighbours)
   gcen2 spatial-build -i $SP -m weighted -o smoothed.network -c 0.6 -t 8 --moran moran_smoothed.tsv
   # one network per region
   gcen2 spatial-build -i $SP -m regions --region-key region -o regions/ -c 0.6 -t 8
   gcen2 plot moran -i moran_smoothed.tsv -o fig/moran

``--moran`` writes each gene's Moran's I (and, in weighted mode, its I after
smoothing), with its p-value and FDR: here 1,992 of 2,000 genes are
significantly autocorrelated (1,999 have I > 0), and the
smoothing raises the autocorrelation of every gene, which is why the smoothed
network has many more edges (27,558) than the filtered one (482).

Per-region networks need enough spots: a region of a few dozen spots gives
noisy correlations, and its network is larger, not more reliable (here the
31-spot region R05 gives the largest network). Prefer regions with hundreds of
spots, or ``spatial-modules``, which pools neighbouring spots first.

The same in Python
------------------

.. code-block:: python

   import anndata as ad
   import gcen2

   adata = ad.read_h5ad("data/example/spatial/mouse_brain_anterior_subset.h5ad")
   res = gcen2.spatial.spatial_modules(adata, groupby="region", threads=8)
   moran = gcen2.spatial.spatial_autocorrelation(adata, n_neighbors=6)
   print(res.power, len(res.modules()), moran.head())
