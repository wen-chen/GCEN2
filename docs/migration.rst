Migration from GCEN v1
======================

GCEN2 reads every GCEN v1 input file format unchanged: header-less expression
matrices, ``.network`` files, module files and seed-gene lists. Each v1 tool
has a ``gcen2`` sub-command, and most flags keep their v1 names (``build``
supports ``-m -c -p -s -f -a -l -t``). Exceptions: ``rwr`` uses ``--gamma`` /
``--prob-cutoff`` instead of ``-r`` / ``-p``. Like v1, ``gcen2 build`` writes
edges while computing (block-wise), so memory stays small for large networks
(``network.build_to_file``; ``network.build`` returns the network in memory).
Differences in results are listed under *Behaviour Differences* below.

Command Mapping
---------------

.. list-table::
   :widths: 48 52
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
   * - ``network_stat -i net``
     - ``gcen2 stats -i net``
   * - ``enrich -e genes -b bg -g go.obo -a assoc -o out``
     - ``gcen2 enrich -e genes -b bg -g go.obo -a assoc -o out``
   * - ``csv_to_tsv -i in.csv -o out.tsv`` / ``tsv_to_csv -i in.tsv -o out.csv``
     - no command; in Python ``gcen2.io.write_tsv(gcen2.io.read_csv("in.csv").data, "out.tsv")``
       (and ``read_tsv`` / ``write_csv``), see *CSV files* below
   * - ``data_stat -i in.tsv``
     - no command; ``data_stat()`` below
   * - ``generate_expr_matrix_from_rsem -i samples.txt -o out.tsv [-t]``
     - no command; ``matrix_from_rsem()`` below
   * - ``generate_expr_matrix_from_stringtie -i samples.txt -o out.tsv [-t]``
     - no command; ``matrix_from_stringtie()`` below
   * - ``calculate_accuracy -g go-basic.obo -a assoc -n net``
     - no counterpart, see below

GCEN v1 Utilities Without a Command
-----------------------------------

Four tools of GCEN v1's ``util/`` directory have no ``gcen2`` command. The
first three take a few lines of Python; each function below gives the same
values as the v1 tool (checked against GCEN v1 0.6.3 in
``tests/test_migration_recipes.py``).

``data_stat`` prints the mean (± standard deviation) of the genes' means and of
their standard deviations, with the 5th to 50th percentiles, to choose the
cutoffs of ``filter``. As v1, it uses the population standard deviation and
takes percentiles as the element at position ``int(n * q)`` of the sorted values:

.. code-block:: python

   import numpy as np
   import gcen2

   def data_stat(path):
       """GCEN v1 data_stat: the genes' means and standard deviations."""
       x = gcen2.io.read_tsv(path).data.to_numpy(dtype=float)
       for name, v in (("mean", x.mean(axis=1)), ("sd", x.std(axis=1))):
           s = np.sort(v)
           pct = "; ".join(f"{q:.0%}: {s[int(len(s) * q)]:g}"
                           for q in (0.05, 0.1, 0.2, 0.25, 0.33, 0.5))
           print(f"gene {name}: {v.mean():g} +- {v.std():g} ({pct})")

``generate_expr_matrix_from_rsem`` and ``generate_expr_matrix_from_stringtie``
build a genes × samples matrix from a tab-separated list of sample IDs and
result files (RSEM ``genes.results``; StringTie GTF files, whose transcripts'
values are summed per gene). FPKM by default, TPM with ``value="TPM"`` (v1:
``-t``). Write the result with ``.to_csv("out.tsv", sep="\t")``; v1 writes
the genes in no particular order and with fewer digits:

.. code-block:: python

   import re
   import pandas as pd

   def _samples(sample_list):
       return pd.read_csv(sample_list, sep="\t", header=None, index_col=0).iloc[:, 0]

   def matrix_from_rsem(sample_list, value="FPKM"):
       """GCEN v1 generate_expr_matrix_from_rsem."""
       return pd.DataFrame({s: pd.read_csv(p, sep="\t", index_col=0)[value]
                            for s, p in _samples(sample_list).items()})

   def matrix_from_stringtie(sample_list, value="FPKM"):
       """GCEN v1 generate_expr_matrix_from_stringtie."""
       def gene_values(gtf):
           rows = []
           for line in open(gtf):
               f = line.rstrip("\n").split("\t")
               if len(f) > 8 and f[2] == "transcript":
                   a = dict(re.findall(r'(\S+) "([^"]*)"', f[8]))
                   rows.append((a["gene_id"], float(a[value])))
           return pd.DataFrame(rows, columns=["gene", value]).groupby("gene")[value].sum()
       return pd.DataFrame({s: gene_values(p) for s, p in _samples(sample_list).items()})

``calculate_accuracy`` scored GO prediction from a network: for every gene
with annotations and a neighbour, the GO terms over-represented among its
annotated neighbours (hypergeometric p below ``-p``) are its predictions, and
the accuracy is the fraction of genes with at least one prediction of which
one is among the gene's own annotations. GCEN2 has no counterpart. The
evaluation of gene-function prediction in ``benchmarks/case_gtex_blood`` is a
cross-validation instead: annotations
are held out, and each term is scored by the area under the ROC curve.

Python API
----------

.. code-block:: python

   import gcen2

   # Replaces: data_norm -m tmm
   normed = gcen2.preprocessing.normalize(expr, method='tmm')

   # Replaces: network_build -m spearman -c 0.8 -p 0.001
   net = gcen2.network.build(expr, method='spearman', cor=0.8, pval=0.001)

   # Replaces: module_identify -s 0.5
   modules = gcen2.module.identify(net, similarity=0.5)

Verified Equivalence
--------------------

GCEN2 was checked against a freshly compiled GCEN v1 0.6.3 on the sample data
(32,940 genes; ``benchmarks/compare_v1_v2.py``, results in
``benchmarks/results.json``):

- **Network building** (Spearman, ``|r| >= 0.9``, ``p <= 0.001``, with and
  without ``-s``): identical edge sets (3,892,110 edges unsigned), and every
  edge line is written exactly as v1 writes it (correlation ``%f``, p-value
  ``%e``).
- **Module identification** (``-s 0.5``): identical modules (1,657 on the full
  network).
- **GO annotation** (network and module mode, 10k-gene network): identical
  queries, study/population counts, term names and namespaces; v1's two-sided
  p-values are reproduced from GCEN2's counts (relative difference < 5e-6).
  This is also a unit test (``tests/data/v1_annotate``).
- **RWR**: v1's values divided by the number of seeds match GCEN2's
  probabilities to within ``2.2e-8`` (see *RWR seed normalisation* below), and
  the top-100 genes are the same. On the full network v1 runs out of memory
  because it builds a dense ``n x n`` matrix; GCEN2 uses a sparse one.

Behaviour Differences from GCEN v1
----------------------------------

The following differences are deliberate.

FDR mode (``build -f``)
   v1 adds an FDR column but filters edges on ``p`` and ``|r|`` only. GCEN2
   also drops edges with ``FDR > pval``. To reproduce v1 output, build without
   the FDR filter and apply Benjamini-Hochberg yourself, or filter v1 output by
   FDR. GCEN2 applies the standard Benjamini-Hochberg procedure over all
   ``n(n-1)/2`` tests. In signed mode (``-s``) v1 ranks only the positive
   correlations, so its FDR values are slightly larger. Without ``-s`` the two
   agree exactly.

RWR seed normalisation
   GCEN2 normalises the seed vector to sum to 1, so the output values are true
   probabilities. v1 does not normalise it, so with *k* seeds of weight 1 every
   v1 value is *k* times larger. As a result the default ``--prob-cutoff 0.01``
   selects far fewer genes than v1's ``-p 0.01``. To select the same genes, use
   ``--prob-cutoff`` equal to v1's ``-p`` divided by the total seed weight.

RWR output
   Seed genes are not listed in the result (v1 lists them). The cutoff test is
   ``probability >= cutoff`` (v1: ``>``).

RWR weighted networks
   GCEN2 uses ``|r|`` as the edge weight. v1 ``-w`` uses the signed correlation,
   which gives negative transition probabilities on unsigned networks.

Enrichment test
   GCEN2 uses a one-sided (over-representation) hypergeometric test and keeps
   terms with ``p <= pval``. v1 uses a two-sided Fisher exact test and keeps
   terms with ``p < pval``, reporting ``e``/``p`` (enriched/purified).

Annotation (``annotate``) and ``enrich``
   Query sets, background, counts, term names and namespaces follow v1
   exactly, and output files have v1 names (``<gene>.go`` / ``.kegg``,
   ``module_<i>.go`` plus ``module_<i>.txt``).

   - Network mode: only genes without annotation are annotated, from their
     annotated neighbours, against the annotated genes adjacent to an
     unannotated gene.
   - Module mode: each module is tested against all module genes.
   - ``enrich``: only annotated genes are counted.

   Only the test differs (see *Enrichment test*); the two-sided p-values
   computed from GCEN2's counts reproduce v1's. ``-p`` works in GCEN2; v1
   0.6.3 ``annotate`` ignores ``-p`` (only ``--pval`` is recognised).

Output file formats
   ``gcen2 build`` writes a header line ``gene_a gene_b correlation p_value``
   (v1: ``#node1 node2 correlation p-value``). Numbers are written as in v1
   (correlation with 6 decimals, p-value and FDR in ``%e`` notation), so the
   edge lines themselves are identical to v1's; only their order can differ,
   because v1 writes them in thread order.
   ``gcen2 module`` writes one ``module``/``gene`` pair per line instead of one
   module per line. ``gcen2 rwr`` writes a ``gene probability`` header.
   Annotation tables use the columns ``term_id name namespace enrichment
   study_count study_n pop_count pop_n p_value fdr``, where ``enrichment`` is
   the fold enrichment (v1: ``e``/``p``) and ``fdr`` the Benjamini-Hochberg
   adjusted p-value over every term annotated in the background (new; v1
   reports only p-values). All GCEN v1 files are still accepted as
   *input*, module files included (:func:`gcen2.read_modules`).

Breaking Changes
----------------

``build -s`` (signed network)
   Same meaning as v1: only positive correlations are kept. The default keeps
   both signs, selected by ``|r|``. Early GCEN2 development versions had this
   inverted.

``network_merge`` output format
   The direction column is now a separate ``"1"`` / ``"-1"`` column, not
   appended to the correlation value string. Update any downstream scripts
   that parse the merged network file.

Standard deviation in ``filter_genes``
   Uses sample standard deviation (``ddof=1``). GCEN v1 used population
   standard deviation (``ddof=0``). Filtering results will differ slightly.

CSV files
   GCEN2 has no ``csv_to_tsv`` / ``tsv_to_csv`` commands; its commands read
   tab-separated files. Convert in Python with ``gcen2.io.read_csv`` and
   ``gcen2.io.write_tsv`` (or ``read_tsv`` and ``write_csv``). ``read_csv``
   reads quoted fields containing commas as one field, and detects
   header-less GCEN v1 matrices as ``read_tsv`` does; v1's ``csv_to_tsv``
   split such fields at the comma.

New in GCEN2
------------

Beyond the GCEN v1 commands (full list in ``CHANGELOG.md``):

- **WGCNA modules** (``gcen2 module -e``), consensus modules
  (``consensus-modules``) and module preservation (``module-preservation``).
- **Single cell**: the hdWGCNA workflow per cell type
  (``celltype-modules``: metacells, modules, cell eigengenes, Harmony, DMEs,
  traits), module projection (``project-modules``), correlation networks per
  cell type (``celltype``); see :doc:`tutorials/single_cell`.
- **Spatial**: the hdWGCNA spatial workflow (``spatial-modules``) and
  networks with spatial information (``spatial-build``, Moran's I); see
  :doc:`tutorials/spatial`.
- **Trajectories** (``trajectory``), **multi-omics** fusion and multiplex RWR
  (``fuse``, ``rwr-multiplex``), **differential co-expression** (``diff``).
- **Figures** for every analysis (``gcen2 plot``) and **HTML / PDF reports**
  (``gcen2 report``); see :doc:`tutorials/figures`.
- A **Python API** for everything (:doc:`api`).
