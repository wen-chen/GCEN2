#!/bin/bash
# GCEN2 (gcen2 trajectory) and the R implementation on the paul15 cells.
# Run from the data directory after 01_prepare.py; then 03_compare.py.
# REPS runs of each (default 2); 03_compare.py averages them.
R=${RSCRIPT:-Rscript}   # an R with WGCNA and data.table (../environment-r.yml)
GCEN2=${GCEN2:-gcen2}
S=$(dirname "$(readlink -f "$0")")
d=runs/trajectory
mkdir -p $d
rm -f $d/time.txt
T="/usr/bin/time -a -o $d/time.txt -f"
for rep in $(seq ${REPS:-2}); do
# windows of a fifth of the cells, metacells of 5, log2, |r| >= 0.5, FDR <= 0.001
$T "gcen2 %e s %M KB" $GCEN2 trajectory -i work/paul15_pt.h5ad --pseudotime dpt_pseudotime \
  -o $d/gcen2 -w 0.2 --metacell 5 -m pearson -l log2 -c 0.5 -p 0.001 -f -t 8 > $d/gcen2.log 2>&1
$T "r_reference %e s %M KB" $R "$S/02_reference.R" work/paul15_expr.tsv work/paul15_pt.tsv \
  $d/r 0.2 5 0.5 0.001 8 > $d/r.log 2>&1
cat $d/gcen2.log $d/r.log
done
echo ALL DONE
