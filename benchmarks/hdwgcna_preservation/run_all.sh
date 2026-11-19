#!/bin/bash
# HLCA natural killer cells: the modules of Meyer_2019 projected onto, and
# tested for preservation in, Krasnow_2020. hdWGCNA (reference network,
# ProjectModules, query metacells, ModulePreservation with 500 permutations)
# against gcen2 celltype-modules + project-modules + module-preservation.
# Each tool twice (timed), then 03_compare.py. Run from the consensus data directory
# after ../hdwgcna_consensus/01_prepare.py (nk.h5ad, nk_export/) and
#   python 01_prepare.py nk.h5ad ref_meyer.h5ad Meyer_2019
R=${RSCRIPT:-Rscript}   # an R with the packages of environment-r.yml
GCEN2=${GCEN2:-gcen2}
PY=${PYTHON:-python}
S=$(dirname "$(readlink -f "$0")")
N=${N_PERM:-500}
NB="--use-rep X_scanvi_emb --n-pcs 30 -k 25 --max-shared 10"
mkdir -p pres
rm -f pres/time.txt
gcen2_all() {
  $GCEN2 celltype-modules -i ref_meyer.h5ad --groupby cell_type $NB --write-metacells \
    -o pres/gcen2/ref -t 8 &&
  $GCEN2 project-modules -i nk.h5ad -m pres/gcen2/ref/natural_killer_cell/modules.tsv \
    --groupby dataset --group Krasnow_2020 -o pres/gcen2/projected &&
  $GCEN2 module-preservation -r pres/gcen2/ref/natural_killer_cell/metacells.tsv \
    -m pres/gcen2/ref/natural_killer_cell/modules.tsv -q nk.h5ad --groupby dataset \
    --group Krasnow_2020 $NB --permutations $N -o pres/gcen2/preservation -t 8
}
export -f gcen2_all
export GCEN2 NB N
for rep in 1 2; do
  /usr/bin/time -f "gcen2 $rep %e s %M KB" -a -o pres/time.txt \
    bash -c gcen2_all > pres/gcen2.log 2>&1
  /usr/bin/time -f "hdwgcna $rep %e s %M KB" -a -o pres/time.txt \
    $R "$S/02_hdwgcna.R" nk_export pres/hdwgcna Meyer_2019 Krasnow_2020 $N > pres/hdwgcna.log 2>&1
done
cat pres/time.txt
$PY -I "$S/03_compare.py" nk.h5ad pres $N pres/time.txt
