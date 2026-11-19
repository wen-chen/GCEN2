#!/bin/bash
# Slide-seqV2 (no Visium grid): gcen2 spatial-modules with knn metaspots
# (k = 100 beads) and hdWGCNA with MetacellsByGroups on the spatial
# coordinates (same k, max_shared, gene fraction), each timed twice; then
# both with metacell seeds 1 and 2 (untimed) for the seed-to-seed agreement.
# Run from the spatial data directory after
#   python 01_prepare_slideseq.py raw/slideseqv2/slideseqv2.h5ad work/slideseqv2
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
GCEN2=${GCEN2:-gcen2}
S=$(dirname "$(readlink -f "$0")")
K=100; FRACTION=0.01
mkdir -p runs
rm -f runs/time_slideseq.txt
for rep in 1 2; do
  /usr/bin/time -f "gcen2 slideseqv2 $rep %e s %M KB" -a -o runs/time_slideseq.txt \
    $GCEN2 spatial-modules -i work/slideseqv2.h5ad -o runs/gcen2/slideseqv2 -t 8 \
    -k $K --fraction $FRACTION --plot > runs/gcen2_slideseqv2.log 2>&1
  /usr/bin/time -f "hdwgcna slideseqv2 $rep %e s %M KB" -a -o runs/time_slideseq.txt \
    $R "$S/02_hdwgcna_knn.R" work/slideseqv2_export runs/hdwgcna/slideseqv2 $K \
    $(( (2 * K + 3) / 7 )) $FRACTION 42 > runs/hdwgcna_slideseqv2.log 2>&1
done
for seed in 1 2; do
  $GCEN2 spatial-modules -i work/slideseqv2.h5ad -o runs/gcen2/slideseqv2_seed$seed -t 8 \
    -k $K --fraction $FRACTION --seed $seed > /dev/null 2>&1
  $R "$S/02_hdwgcna_knn.R" work/slideseqv2_export runs/hdwgcna/slideseqv2_seed$seed $K \
    $(( (2 * K + 3) / 7 )) $FRACTION $seed > /dev/null 2>&1
done
cat runs/time_slideseq.txt
${PYTHON:-python} "$S/05_compare_slideseq.py" runs
