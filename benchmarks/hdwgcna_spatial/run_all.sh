#!/bin/bash
# Both tools on both samples (DLPFC 151673 by cortical layer; 10x mouse brain
# sagittal anterior by k-means region), each run twice (timing; gcen2 then
# once more, untimed, with --write-metaspots for the comparison).
# Run from the spatial data directory after 01_prepare.py:
#   python 01_prepare.py dlpfc raw/dlpfc_151673 work/dlpfc
#   python 01_prepare.py tenx raw/mouse_brain_anterior work/mouse_anterior
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
GCEN2=${GCEN2:-gcen2}
S=$(dirname "$(readlink -f "$0")")
mkdir -p runs
rm -f runs/time.txt
for name in dlpfc mouse_anterior; do
  for rep in 1 2; do
    /usr/bin/time -f "gcen2 $name $rep %e s %M KB" -a -o runs/time.txt \
      $GCEN2 spatial-modules -i work/$name.h5ad --groupby region -o runs/gcen2/$name -t 8 \
      --plot > runs/gcen2_$name.log 2>&1
    /usr/bin/time -f "hdwgcna $name $rep %e s %M KB" -a -o runs/time.txt \
      $R "$S/02_hdwgcna.R" work/${name}_export runs/hdwgcna/$name region \
      > runs/hdwgcna_$name.log 2>&1
  done
  # untimed: also the metaspot matrix, and the modules before the reassignment
  # (hdWGCNA's ConstructNetwork), for 04_compare.py
  $GCEN2 spatial-modules -i work/$name.h5ad --groupby region -o runs/gcen2/$name -t 8 \
    --write-metaspots > /dev/null 2>&1
  $GCEN2 spatial-modules -i work/$name.h5ad --groupby region -t 8 --no-reassign \
    -o runs/gcen2/${name}_noreassign > /dev/null 2>&1
  echo "$name: $(tail -1 runs/gcen2_$name.log) | $(tail -1 runs/hdwgcna_$name.log)"
done
cat runs/time.txt
${PYTHON:-python} "$S/04_compare.py" runs dlpfc mouse_anterior
