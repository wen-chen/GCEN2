#!/bin/bash
# HLCA natural killer cells over 4 datasets: hdWGCNA's consensus workflow with
# harmonized eigengenes, FindAllDMEs and ModuleTraitCorrelation against gcen2
# celltype-modules --consensus-by --harmonize-by --dme-by --trait. Each tool twice (timed), then metacell seeds 1 and 2
# (untimed), then 03_compare.py. Run from the consensus data directory after
#   python 01_prepare.py ../raw/hlca_core.h5ad nk "natural killer cell" \
#     Meyer_2019 Krasnow_2020 Banovich_Kropski_2020 Lafyatis_Rojas_2019_10Xv2
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
GCEN2=${GCEN2:-gcen2}
PY=${PYTHON:-python}
S=$(dirname "$(readlink -f "$0")")
OPTS="--groupby cell_type --consensus-by dataset --harmonize-by donor_id --use-rep X_scanvi_emb --n-pcs 30 -k 25 --max-shared 10 --dme-by dataset --trait age --trait BMI --trait sex --trait-groupby dataset -t 8"
mkdir -p runs
rm -f runs/time.txt
for rep in 1 2; do
  /usr/bin/time -f "gcen2 $rep %e s %M KB" -a -o runs/time.txt \
    $GCEN2 celltype-modules -i nk.h5ad $OPTS -o runs/gcen2 > runs/gcen2.log 2>&1
  /usr/bin/time -f "hdwgcna $rep %e s %M KB" -a -o runs/time.txt \
    $R "$S/02_hdwgcna.R" nk_export runs/hdwgcna 42 > runs/hdwgcna.log 2>&1
done
for seed in 1 2; do
  $GCEN2 celltype-modules -i nk.h5ad $OPTS --seed $seed -o runs/gcen2_seed$seed > /dev/null 2>&1
  $R "$S/02_hdwgcna.R" nk_export runs/hdwgcna_seed$seed $seed > /dev/null 2>&1
done
cat runs/time.txt
$PY -I "$S/03_compare.py" nk runs runs/time.txt
