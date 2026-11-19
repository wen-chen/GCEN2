#!/bin/bash
# Thread scaling of the main GCEN2 commands: wall time and peak memory with
# -t 1, 2, 4, 8 and 16 (REPS runs each, default 2), on an idle machine.
#
#   build     GTEx whole blood, 11,944 genes x 803 samples, r >= 0.4 with FDR
#             (../wgcna_gtex: work/expr_11944_pc5.tsv)
#   module_e  WGCNA-style modules of the same matrix (gcen2 module -e)
#   diff      differential co-expression of two 12.5M-edge networks
#             (../dgca_gtex: runs/dgca/pearson_5000/{A,B}.network)
#   celltype  hdWGCNA workflow on the 78,816 alveolar macrophages of the whole
#             Human Lung Cell Atlas, expression left on disk (--backed)
#             (../hdwgcna_hlca: work/hlca_all.h5ad)
#
# Usage: GTEX_DIR=<gtex> HLCA_DIR=<hlca> bash 01_threads.sh OUT_DIR
# Then 02_summary.py OUT_DIR.
set -e
GCEN2=${GCEN2:-gcen2}
OUT=$(readlink -f "$1")
mkdir -p "$OUT"
T="$OUT/time.txt"
: > "$T"
for rep in $(seq ${REPS:-2}); do
for t in 1 2 4 8 16; do
  cd "$GTEX_DIR"
  /usr/bin/time -a -o "$T" -f "build $t $rep %e %M" $GCEN2 build -i work/expr_11944_pc5.tsv \
    -o "$OUT/network.tsv" -m pearson -c 0.4 -p 0.001 -f -s -t $t > /dev/null
  /usr/bin/time -a -o "$T" -f "module_e $t $rep %e %M" $GCEN2 module -e work/expr_11944_pc5.tsv \
    -o "$OUT/modules.tsv" --signed --min-size 30 --merge-cut 0.25 -t $t > /dev/null
  d=runs/dgca/pearson_5000
  /usr/bin/time -a -o "$T" -f "diff $t $rep %e %M" $GCEN2 diff -A $d/A.network -B $d/B.network \
    -o "$OUT/diff.tsv" --n-a 541 --n-b 262 -m pearson -t $t > /dev/null
  rm -f "$OUT/diff.tsv"
  cd "$HLCA_DIR"
  /usr/bin/time -a -o "$T" -f "celltype $t $rep %e %M" $GCEN2 celltype-modules \
    -i work/hlca_all.h5ad --groupby cell_type --group "alveolar macrophage" --backed \
    --use-rep X_scanvi_emb -t $t -o "$OUT/celltype" > /dev/null 2>&1
  echo "rep $rep, $t threads: $(tail -n 4 "$T" | awk '{printf "%s %ss  ", $1, $4}')"
done
done
