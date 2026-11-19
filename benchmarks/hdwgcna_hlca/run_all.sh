#!/bin/bash
# HLCA core: for a 50,000-cell subset and the whole atlas,
#   * gcen2 celltype-modules on the 10 largest cell types (twice) and
#     hdWGCNA's per-type workflow on the same types (once: ~25 min each);
#   * gcen2 celltype-modules on all cell types (twice);
#   * both gcen2 runs again with --backed (expression left on disk).
# Then 04_compare.py. Run from the HLCA data directory after
#   python 01_prepare.py raw/hlca_core.h5ad work 50000
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
GCEN2=${GCEN2:-gcen2}
PY=${PYTHON:-python}
S=$(dirname "$(readlink -f "$0")")
mkdir -p runs
rm -f runs/time.txt
groups=()
while read -r t; do groups+=(--group "$t"); done < work/types.txt
for size in 50k all; do
  for rep in 1 2; do
    /usr/bin/time -f "gcen2_top10 $size $rep %e s %M KB" -a -o runs/time.txt \
      $GCEN2 celltype-modules -i work/hlca_$size.h5ad --groupby cell_type \
      --use-rep X_scanvi_emb -t 8 -o runs/gcen2/${size}_top10 "${groups[@]}" \
      > runs/gcen2_${size}_top10.log 2>&1
    /usr/bin/time -f "gcen2_all $size $rep %e s %M KB" -a -o runs/time.txt \
      $GCEN2 celltype-modules -i work/hlca_$size.h5ad --groupby cell_type \
      --use-rep X_scanvi_emb -t 8 -o runs/gcen2/${size}_all > runs/gcen2_${size}_all.log 2>&1
    # the same with the expression left on disk (one cell type read at a time)
    /usr/bin/time -f "gcen2_top10_backed $size $rep %e s %M KB" -a -o runs/time.txt \
      $GCEN2 celltype-modules -i work/hlca_$size.h5ad --groupby cell_type --backed \
      --use-rep X_scanvi_emb -t 8 -o runs/gcen2/${size}_top10_backed "${groups[@]}" \
      > /dev/null 2>&1
    /usr/bin/time -f "gcen2_all_backed $size $rep %e s %M KB" -a -o runs/time.txt \
      $GCEN2 celltype-modules -i work/hlca_$size.h5ad --groupby cell_type --backed \
      --use-rep X_scanvi_emb -t 8 -o runs/gcen2/${size}_all_backed > /dev/null 2>&1
  done
  /usr/bin/time -f "hdwgcna_top10 $size 1 %e s %M KB" -a -o runs/time.txt \
    $R "$S/02_hdwgcna.R" work/hlca_$size.h5 work/types.txt runs/hdwgcna/$size \
    > runs/hdwgcna_$size.log 2>&1
  echo "$size: $(head -1 runs/gcen2_${size}_all.log)"
done
cat runs/time.txt
$PY "$S/04_compare.py" runs work/types.txt 50k all
