#!/bin/bash
# GCEN2 (build of each condition, every pair, then diff) and DGCA (ddcorAll),
# Pearson and Spearman, on GTEx whole blood, males (A) vs females (B). Run from the GTEx data
# directory after 01_prepare.py; then 03_compare.py. REPS runs of each
# (default 2); 03_compare.py averages them.
R=${RSCRIPT:-Rscript}   # an R with DGCA 1.0.3 (see ../environment-r.yml)
GCEN2=${GCEN2:-gcen2}
S=$(dirname "$(readlink -f "$0")")
mkdir -p runs/dgca
rm -f runs/dgca/time.txt
T="/usr/bin/time -a -o runs/dgca/time.txt -f"
for rep in $(seq ${REPS:-2}); do
for m in pearson spearman; do
for n in 2000 5000; do
  d=runs/dgca/${m}_$n; mkdir -p $d
  na=$(head -1 work/dc_${n}_male.tsv | awk -F'\t' '{print NF - 1}')
  nb=$(head -1 work/dc_${n}_female.tsv | awk -F'\t' '{print NF - 1}')
  # -c 0 -p 1: every gene pair, as DGCA reports
  $T "gcen2_build_A ${m}_$n %e s %M KB" $GCEN2 build -i work/dc_${n}_male.tsv \
    -o $d/A.network -m $m -c 0 -p 1 -t 8 > $d/gcen2.log 2>&1
  $T "gcen2_build_B ${m}_$n %e s %M KB" $GCEN2 build -i work/dc_${n}_female.tsv \
    -o $d/B.network -m $m -c 0 -p 1 -t 8 >> $d/gcen2.log 2>&1
  $T "gcen2_diff ${m}_$n %e s %M KB" $GCEN2 diff -A $d/A.network -B $d/B.network \
    -o $d/gcen2_diff.tsv --n-a $na --n-b $nb -m $m -t 8 >> $d/gcen2.log 2>&1
  $T "dgca ${m}_$n %e s %M KB" $R "$S/02_dgca.R" work/dc_${n}_male.tsv \
    work/dc_${n}_female.tsv $d/dgca.tsv 8 $m > $d/dgca.log 2>&1
  echo "run $rep, $m, $n genes: $(tail -1 $d/gcen2.log) | $(tail -1 $d/dgca.log)"
done
done
done
echo ALL DONE
