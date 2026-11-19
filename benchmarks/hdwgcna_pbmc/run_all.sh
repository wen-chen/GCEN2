#!/bin/bash
# Both tools on every run: six cell types (all cells and two split halves)
# and metacell seeds 1 and 2 for CD14 Mono, CD4 T and B. Run from the data
# directory (data/, work/; see README) after 01_preprocess.py and the split files
# (work/splits/<type>_A.txt, _B.txt: a random half of the type's cells each).
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
PY=${PYTHON:-python}
S=$(dirname "$(readlink -f "$0")")
rm -f runs/time.txt
jobs=()
for ct in "CD14 Mono" "CD4 T" "B" "CD8 T" "NK" "FCGR3A Mono"; do
  s=${ct// /_}
  jobs+=("$s|$ct|-|42" "${s}_A|$ct|work/splits/${s}_A.txt|42" "${s}_B|$ct|work/splits/${s}_B.txt|42")
done
for ct in "CD14 Mono" "CD4 T" "B"; do
  s=${ct// /_}
  for seed in 1 2; do jobs+=("${s}_seed$seed|$ct|-|$seed"); done
done
mkdir -p runs
for j in "${jobs[@]}"; do
  IFS='|' read -r name ct cells seed <<< "$j"
  /usr/bin/time -f "gcen2 $name %e s %M KB" -a -o runs/time.txt \
    $PY -I "$S/03_gcen2.py" "$ct" runs/gcen2/$name "$cells" "$seed" > runs/gcen2_$name.log 2>&1
  /usr/bin/time -f "hdwgcna $name %e s %M KB" -a -o runs/time.txt \
    $R "$S/02_hdwgcna.R" "$ct" runs/hdwgcna/$name "$cells" "$seed" > runs/hdwgcna_$name.log 2>&1
  echo "$name: $(tail -1 runs/gcen2_$name.log) | $(tail -1 runs/hdwgcna_$name.log)"
done
