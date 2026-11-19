Figures and reports
===================

Every analysis has publication figures: ``gcen2 plot <figure>`` on the files a
command wrote, or the same function in Python (``gcen2.plot``). The
:doc:`../cli` lists all 43 figure commands; the tutorials show them in use.

One style
---------

- **Widths of journal columns**: ``--size single`` (89 mm), ``onehalf``
  (120 mm), ``double`` (183 mm) or a width in mm; the height follows the
  content.
- **Text** 5 to 7 pt, sans-serif (Arial or Helvetica when installed), kept
  as text in PDF and SVG so that it can be edited.
- **One colour per module in every figure**: WGCNA's colours, labelled with
  both names (``ME0 (turquoise)``, ``ME1 (blue)``, ...); ``--palette
  okabe-ito`` for the colour-blind-safe Okabe-Ito palette (and labels
  ``ME0``, ``ME1``). Unassigned genes are grey.
- **Formats**: PDF, PNG (300 dpi) and SVG by default; ``--format`` (repeatable)
  for others (``tiff``, ``eps``). Outputs are reproducible: the same data give
  the same files.

.. code-block:: bash

   gcen2 plot kme -i modules/B/kme.tsv -o fig/kme                    # fig/kme.pdf, .png, .svg
   gcen2 plot kme -i modules/B/kme.tsv -o fig/kme.png                # one PNG
   gcen2 plot kme -i modules/B/kme.tsv -o fig/kme --size single --format pdf --palette okabe-ito

The figures by analysis
-----------------------

.. list-table::
   :header-rows: 1
   :widths: 25 75

   * - Analysis
     - ``gcen2 plot ...``
   * - WGCNA / hdWGCNA modules
     - ``soft-power``, ``dendrogram``, ``eigengene-network``, ``module-trait``, ``kme``,
       ``mm-gs``, ``hub-network``, ``hub-network-all``, ``module-umap``, ``tom-heatmap``,
       ``sample-tree``, ``module-sizes``
   * - Single cell, spatial
     - ``me-embedding`` (spot maps with ``--spatial``, projected eigengenes), ``me-dotplot``,
       ``me-violin``, ``me-radar``, ``dme-volcano``, ``dme-lollipop``, ``celltype-overview``,
       ``metacells``, ``set-eigengenes``, ``preservation``, ``moran``
   * - Expression QC
     - ``expression-qc``, ``mean-variance``, ``sample-correlation``
   * - Networks
     - ``degree-distribution``, ``edge-weights``, ``network-view``, ``module-size-hist``,
       ``celltype-networks``, ``null-distribution``, ``merge-overlap``, ``extract-overview``
   * - Function
     - ``enrich-bar``, ``enrich-dot``, ``enrich-heatmap``, ``rwr-rank``, ``rwr-network``
   * - Differential, trajectory, multi-omics
     - ``diff-overview``, ``diff-network``, ``trajectory``, ``fusion``

Some figures need files that a command writes only when asked: the hub-network
and module-UMAP figures read ``tom_hubs.tsv`` (``celltype-modules
--write-tom``, ``module -e --tom-hubs``); ``soft_power.tsv`` and
``dendrogram.json`` come from ``celltype-modules``, ``spatial-modules`` and
``consensus-modules`` (and ``module -e --soft-power / --dendrogram``);
``metacells`` reads ``metacell_members.tsv`` (``celltype-modules``), ``moran``
the table of ``spatial-build --moran``; ``merge-overlap`` and
``extract-overview`` read the inputs and output of ``merge`` and ``extract``.

In Python
---------

Each function takes the result files (paths) or the objects themselves,
returns the matplotlib figure, and saves it when ``output`` is given:

.. code-block:: python

   import gcen2

   fig = gcen2.plot.kme("modules/B/kme.tsv")              # a matplotlib Figure
   fig.axes[0].set_title("B cells, ME0")                   # change it as usual
   gcen2.plot.save_figure(fig, "fig/kme")                  # PDF, PNG and SVG

   gcen2.plot.set_defaults(palette="okabe-ito", formats=("pdf",))
   gcen2.plot.dendrogram("modules/B/dendrogram.json", output="fig/dendrogram")

``hub_network``, ``mm_gs`` and ``module_sizes`` also draw into a given axes
(``ax=``), to combine panels in one figure.

Reports
-------

``gcen2 report`` writes an HTML report (or PDF, with ``-o report.pdf`` and
weasyprint installed) with the results of the files it is given, every figure
as a vector image with a numbered caption:

.. code-block:: bash

   # correlation network and modules
   gcen2 report -n net.network -m module.txt -o report.html
   # WGCNA modules (gcen2 module -e outputs)
   gcen2 report -m modules.tsv --me eigengenes.tsv --kme kme.tsv --hubs hubs.tsv \
                --soft-power soft_power.tsv --dendrogram dendrogram.json -o report.html
   # celltype-modules, with the metacells on the UMAP
   gcen2 report --celltype-modules modules/ --coords cells.h5ad -o report.html
   # merged and extracted networks, Moran's I
   gcen2 report -n net.network --merge-input a.network --merge-input b.network \
                --merged ab.network --extract-genes genes.txt --moran moran.tsv -o report.html

``--figures-dir fig/`` also saves every report figure as PDF, PNG and SVG,
named by section, ready for a manuscript.
