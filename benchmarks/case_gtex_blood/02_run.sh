#!/bin/bash
# GCEN2 runs of the case study, from the GTEx data directory after
# 01_prepare.py:
#   * WGCNA-style modules of the 11,944 genes (PC-corrected; the settings of
#     ../wgcna_gtex/06_tom_modules.py, so the same modules as blockwiseModules)
#     with eigengenes, kME, hub genes and eigengene-trait correlations;
#   * the correlation network (r >= 0.4, FDR <= 0.001, positive edges; as
#     ../wgcna_gtex/run_all.sh) and a degree-preserving shuffle of it, the
#     networks of the gene function prediction (04_rwr_cv.py).
set -e
GCEN2=${GCEN2:-gcen2}
R=runs/case W=work/case
mkdir -p $R
$GCEN2 module -e work/expr_11944_pc5.tsv -o $R/modules.tsv --signed --min-size 30 \
  --merge-cut 0.25 -t 8 --me $R/eigengenes.tsv --kme $R/kme.tsv --hubs $R/hubs.tsv \
  --traits $W/traits.tsv --trait-cor $R/trait_cor.tsv --soft-power $R/soft_power.tsv \
  --dendrogram $R/dendrogram.tsv --tom-hubs $R/tom_hubs.tsv > $R/module.log
$GCEN2 build -i work/expr_11944_pc5.tsv -o $R/network.tsv -m pearson -c 0.4 -p 0.001 \
  -f -s -t 8 > $R/build.log
$GCEN2 shuffle -i $R/network.tsv -o $R/network_shuffled.tsv --seed 42 > $R/shuffle.log
tail -n 1 $R/module.log $R/build.log $R/shuffle.log
