#!/bin/bash
# GCEN2 (rwr-multiplex, fuse) against RandomWalkRestartMH and an R
# implementation of the fusion, on the three layers of 01_prepare.py.
# Run from the data directory; then 04_compare.py. REPS runs of each
# (default 2); 04_compare.py averages them.
R=${RSCRIPT:-Rscript}   # an R with RandomWalkRestartMH (see ../environment-r.yml)
GCEN2=${GCEN2:-gcen2}
S=$(dirname "$(readlink -f "$0")")
L="work/coexpression.tsv work/ppi.tsv work/regulation.tsv"
mkdir -p runs/multilayer
d=runs/multilayer
rm -f $d/time.txt
T="/usr/bin/time -a -o $d/time.txt -f"
for rep in $(seq ${REPS:-2}); do
# RandomWalkRestartMH's defaults: r = 0.7, delta = 0.5, geometric mean
$T "gcen2_rwrm %e s %M KB" $GCEN2 rwr-multiplex $(printf -- "-n %s " $L) \
  --name coexpression --name ppi --name regulation -g work/seeds.txt \
  -o $d/gcen2_rwrm.tsv --gamma 0.7 --delta 0.5 > $d/gcen2.log 2>&1
$T "rwrmh %e s %M KB" $R "$S/02_rwrmh.R" work/seeds.txt $d/rwrmh.tsv 0.7 0.5 $L \
  > $d/rwrmh.log 2>&1
$T "gcen2_fuse %e s %M KB" $GCEN2 fuse $(printf -- "-i %s " $L) \
  --name coexpression --name ppi --name regulation -o $d/gcen2_fused.tsv -t 8 >> $d/gcen2.log 2>&1
$T "r_fuse %e s %M KB" $R "$S/03_fuse.R" $d/r_fused.tsv $L > $d/r_fuse.log 2>&1
cat $d/gcen2.log $d/rwrmh.log $d/r_fuse.log
done
echo ALL DONE
