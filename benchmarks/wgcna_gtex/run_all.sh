#!/bin/bash
# GCEN2 (gcen2 build + gcen2 module) and WGCNA on the GTEx matrices.
# Lines: "run name|matrix|gcen2 cutoff|tom"
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
GCEN2=${GCEN2:-gcen2}
S=$(dirname "$(readlink -f "$0")")
# run from the GTEx data directory (work/ from 01_prepare.py)
mkdir -p runs
rm -f runs/time.txt
jobs=("pc5_5000|work/expr_5000_pc5.tsv|0.4|tom"
      "pc5_5000_A|work/expr_5000_pc5_A.tsv|0.4|-"
      "pc5_5000_B|work/expr_5000_pc5_B.tsv|0.4|-"
      "raw_5000|work/expr_5000.tsv|0.8|tom"
      "pc5_10000|work/expr_10000_pc5.tsv|0.4|-"
      "pc5_11944|work/expr_11944_pc5.tsv|0.4|-")
for j in "${jobs[@]}"; do
  IFS='|' read -r name expr cut tom <<< "$j"
  d=runs/gcen2/$name; mkdir -p $d
  /usr/bin/time -f "gcen2_build $name %e s %M KB" -a -o runs/time.txt \
    $GCEN2 build -i $expr -o $d/network.tsv -m pearson -c $cut -p 0.001 -f -s -t 8 > $d/log 2>&1
  /usr/bin/time -f "gcen2_module $name %e s %M KB" -a -o runs/time.txt \
    $GCEN2 module -i $d/network.tsv -o $d/modules.tsv -t 8 >> $d/log 2>&1
  /usr/bin/time -f "wgcna $name %e s %M KB" -a -o runs/time.txt \
    $R "$S/02_wgcna.R" $expr runs/wgcna/$name - $tom > runs/wgcna_$name.log 2>&1
  echo "$name: $(tail -2 $d/log | tr '\n' ' ') | $(tail -1 runs/wgcna_$name.log)"
done
echo ALL DONE
