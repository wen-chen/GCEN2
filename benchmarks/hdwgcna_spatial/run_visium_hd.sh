#!/bin/bash
# Visium HD mouse brain, 16 um bins (a square grid: knn metaspots for both
# tools; hdWGCNA's MetacellsByGroups on the spatial coordinates, within the
# k-means regions): each tool twice, then metacell seeds 1 and 2 (untimed),
# then 05_compare_slideseq.py. Run from the spatial data directory after
#   python 01_prepare_visium_hd.py raw/visium_hd_mouse_brain/binned_outputs/square_016um work/visium_hd
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
GCEN2=${GCEN2:-gcen2}
PY=${PYTHON:-python}
S=$(dirname "$(readlink -f "$0")")
K=${K:-64}; FRACTION=${FRACTION:-0.05}
SHARED=$(( (2 * K + 3) / 7 ))
mkdir -p runs
rm -f runs/time_visium_hd.txt
for rep in 1 2; do
  /usr/bin/time -f "gcen2 visium_hd $rep %e s %M KB" -a -o runs/time_visium_hd.txt \
    $GCEN2 spatial-modules -i work/visium_hd.h5ad --groupby region -o runs/gcen2/visium_hd \
    -t 8 -k $K --fraction $FRACTION --plot > runs/gcen2_visium_hd.log 2>&1
  /usr/bin/time -f "hdwgcna visium_hd $rep %e s %M KB" -a -o runs/time_visium_hd.txt \
    $R "$S/02_hdwgcna_knn.R" work/visium_hd_export runs/hdwgcna/visium_hd $K $SHARED \
    $FRACTION 42 region > runs/hdwgcna_visium_hd.log 2>&1
done
for seed in 1 2; do
  $GCEN2 spatial-modules -i work/visium_hd.h5ad --groupby region -t 8 -k $K \
    --fraction $FRACTION --seed $seed -o runs/gcen2/visium_hd_seed$seed > /dev/null 2>&1
  $R "$S/02_hdwgcna_knn.R" work/visium_hd_export runs/hdwgcna/visium_hd_seed$seed $K \
    $SHARED $FRACTION $seed region > /dev/null 2>&1
done
cat runs/time_visium_hd.txt
$PY "$S/05_compare_slideseq.py" runs visium_hd time_visium_hd.txt
