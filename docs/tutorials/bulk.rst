Bulk RNA-seq
============

This tutorial runs the GCEN workflow and a WGCNA module analysis on the GCEN
v1 sample data in ``data/example/bulk`` (32,940 genes x 15 samples, with GO
and KEGG annotations). Every command runs in seconds (the WGCNA modules
in about 3 seconds and 0.5 GB of memory on 8 threads). Work in an empty directory and set
``EX`` to the example data:

.. code-block:: bash

   EX=/path/to/GCEN2/data/example/bulk

Correlation network
-------------------

Normalise (TMM), keep the 75 % of genes with the highest mean, and build a
Spearman network (``|r| >= 0.8``, ``p <= 0.001``, FDR-adjusted, 8 threads):

.. code-block:: bash

   gcen2 normalize -i $EX/gene_expr.tsv -o norm.tsv -m tmm
   gcen2 filter    -i norm.tsv -o filt.tsv -p 0.75
   gcen2 build     -i filt.tsv -o net.network -m spearman -c 0.8 -p 0.001 -f -t 8
   gcen2 stats     -i net.network

``build`` reports 10,253 genes and 21,022 edges; ``stats`` prints the density,
degrees, connected components and clustering coefficient. The network file
has one edge per line (``gene_a gene_b correlation p_value fdr``), with the
numbers written as GCEN v1 writes them. Other correlation measures:
``-m pearson``; for many samples ``--kernel gemm`` computes the correlations as
a matrix product (the same edges, much faster).

Modules and function
--------------------

.. code-block:: bash

   gcen2 module   -i net.network -o module.txt -s 0.5 -t 8
   gcen2 annotate -g $EX/go-basic.obo.gz -a $EX/gene_go.assoc.gz -m module.txt -o annot_go/ -t 8
   gcen2 annotate -k $EX/K2ko.tsv -a $EX/gene_kegg.assoc -m module.txt -o annot_kegg/ -t 8

``module`` groups genes whose network neighbourhoods overlap (Jaccard
similarity >= 0.5): 2,013 modules. ``annotate -m`` tests each module for GO
(or KEGG) terms against all module genes and writes ``module_<i>.go`` (or
``.kegg``) and the genes of the module (``module_<i>.txt``). With ``-n
net.network`` instead of ``-m``, ``annotate`` predicts the function of every
unannotated gene from its annotated neighbours (guilt by association). The GO
files may be gzip-compressed, as here.

Candidate genes for a function, by a random walk with restart from known
genes, and a standalone enrichment test:

.. code-block:: bash

   gcen2 rwr     -n net.network -g $EX/rwr_seed_genes.list -o rwr.tsv
   gcen2 enrich  -e $EX/enrichment_gene.list -b $EX/background_gene.list \
                 -g $EX/go-basic.obo.gz -a $EX/gene_go.assoc.gz -o enrich.tsv
   gcen2 extract -i net.network -g $EX/gene_list.txt -o sub.network

``rwr.tsv`` ranks the genes by their probability (the seeds themselves are not
listed); ``extract`` keeps the edges of the listed genes (62 genes, 52 edges).

Figures and report:

.. code-block:: bash

   gcen2 plot degree-distribution -i net.network -o fig/degree
   gcen2 plot network-view -i net.network -m module.txt -o fig/network
   gcen2 plot rwr-rank -i rwr.tsv -g $EX/rwr_seed_genes.list -o fig/rwr
   gcen2 plot enrich-bar -i enrich.tsv -o fig/enrich
   gcen2 report -n net.network -m module.txt -o report.html

Each figure is written as PDF, PNG and SVG (``fig/degree.pdf``, ...).

WGCNA modules
-------------

``gcen2 module -e`` runs WGCNA's module detection on an expression matrix:
the soft-thresholding power, the topological overlap (TOM), average-linkage
clustering, the dynamic tree cut and merging of modules with similar
eigengenes (``blockwiseModules`` defaults). Here on the genes with the highest
variation (``--percent-sd 0.2`` keeps 6,588 genes), log2-transformed, as a
signed network:

.. code-block:: bash

   gcen2 filter -i norm.tsv -o var.tsv -p 0.75 --percent-sd 0.2
   gcen2 module -e var.tsv -o wgcna/modules.tsv -l log2 --signed -t 8 \
          --me wgcna/eigengenes.tsv --kme wgcna/kme.tsv --hubs wgcna/hubs.tsv \
          --soft-power wgcna/soft_power.tsv --dendrogram wgcna/dendrogram.json \
          --tom-hubs wgcna/tom_hubs.tsv

The command prints the scale-free fit of every power tried. With only 15
samples no power reaches the fit (R² > 0.85), and GCEN2 takes the power the
WGCNA FAQ recommends for that number of samples (18) and says so; WGCNA itself
recommends at least 20 samples. The result: 32 modules (``ME0`` to ``ME31``,
module 0 being the largest), their eigengenes, every gene's module membership
(kME) and the hub genes of each module. With ``--traits traits.tsv
--trait-cor trait_cor.tsv`` (samples as rows, numeric columns) the eigengenes
are also correlated with the traits; ``--gs`` writes the gene significance.
``--preset hdwgcna`` uses hdWGCNA's settings instead of WGCNA's.

.. code-block:: bash

   gcen2 plot soft-power -i wgcna/soft_power.tsv -o fig/soft_power
   gcen2 plot dendrogram -i wgcna/dendrogram.json -o fig/dendrogram
   gcen2 plot eigengene-network -i wgcna/eigengenes.tsv -o fig/eigengenes
   gcen2 plot kme -i wgcna/kme.tsv -o fig/kme
   gcen2 plot hub-network-all -t wgcna/tom_hubs.tsv -k wgcna/kme.tsv -o fig/hubs
   gcen2 report -m wgcna/modules.tsv --me wgcna/eigengenes.tsv --kme wgcna/kme.tsv \
          --hubs wgcna/hubs.tsv --soft-power wgcna/soft_power.tsv \
          --dendrogram wgcna/dendrogram.json --tom-hubs wgcna/tom_hubs.tsv \
          -o wgcna_report.html

Modules found in several data sets (e.g. two tissues) come from
``gcen2 consensus-modules -e set1.tsv -e set2.tsv -o consensus/``; whether a
module of one data set is preserved in another, from ``gcen2
module-preservation`` (see :doc:`single_cell`).

The same in Python
------------------

.. code-block:: python

   import gcen2

   expr = gcen2.io.read_tsv("data/example/bulk/gene_expr.tsv")
   norm = gcen2.preprocessing.normalize(expr, method="tmm")
   filt = gcen2.preprocessing.filter_genes(norm, percent_mean=0.75)
   net = gcen2.network.build(filt, method="spearman", cor=0.8, pval=0.001, fdr=True,
                             threads=8)
   modules = gcen2.module.identify(net, similarity=0.5, threads=8)
   go = gcen2.annotation.annotate_go(modules, "data/example/bulk/go-basic.obo.gz",
                                     "data/example/bulk/gene_go.assoc.gz")
   gcen2.plot.degree_distribution(net, output="fig/degree")
