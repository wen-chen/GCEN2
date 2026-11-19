#!/bin/bash
# PyWGCNA (01_pywgcna.py) and WGCNA's step-by-step pipeline (02_wgcna_stepwise.R)
# on the six GTEx matrices of ../wgcna_gtex, REPS times each (default 2), with
# wall time and peak memory from /usr/bin/time. Run from the GTEx data
# directory, after ../wgcna_gtex/run_all.sh (runs/wgcna/ is the reference).
#   PYTHON   the Python of the PyWGCNA environment (default: python)
#   RSCRIPT  an R with WGCNA 1.74 (default: Rscript)
S=$(dirname "$(readlink -f "$0")")
PY=${PYTHON:-python}
R=${RSCRIPT:-Rscript}
export OPENBLAS_NUM_THREADS=8 OMP_NUM_THREADS=8 MKL_NUM_THREADS=8
mkdir -p runs/pywgcna runs/wgcna_step
T=runs/pywgcna/time.txt
: > $T
jobs=("pc5_5000|expr_5000_pc5.tsv|tom" "pc5_5000_A|expr_5000_pc5_A.tsv|-"
      "pc5_5000_B|expr_5000_pc5_B.tsv|-" "raw_5000|expr_5000.tsv|tom"
      "pc5_10000|expr_10000_pc5.tsv|-" "pc5_11944|expr_11944_pc5.tsv|-")
for rep in $(seq ${REPS:-2}); do
  for j in "${jobs[@]}"; do
    IFS='|' read -r name expr tom <<< "$j"
    /usr/bin/time -a -o $T -f "pywgcna $name $rep %e %M" \
      $PY -I "$S/01_pywgcna.py" work/$expr runs/pywgcna/$name $tom > runs/pywgcna/$name.log 2>&1
    /usr/bin/time -a -o $T -f "wgcna_step $name $rep %e %M" \
      $R "$S/02_wgcna_stepwise.R" work/$expr runs/wgcna_step/$name > runs/wgcna_step/$name.log 2>&1
    echo "rep $rep $name: $(tail -n 2 $T | awk '{printf "%s %ss %.1fGB  ", $1, $4, $5/1048576}')"
  done
done
echo ALL DONE
