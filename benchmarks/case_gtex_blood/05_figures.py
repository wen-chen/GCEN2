"""Figure panels of the case study, all drawn with gcen2.plot.

Run from the GTEx data directory after 03_modules.py and 04_rwr_cv.py.
Writes results/fig_*.{pdf,png} (next to this script):

  fig_module_trait    eigengene-donor trait correlations (gcen2 plot module-trait)
  fig_enrich_dot      top GO biological processes of the cell-type modules
                      (gcen2 plot enrich-dot, BP terms of results/module_go.tsv)
  fig_hub_ME18        hub genes of the naive T-cell module (gcen2 plot hub-network)
  fig_rwr_cv          AUROC per GO term: RWR against the controls
"""

from pathlib import Path

import pandas as pd

from gcen2.plot import enrichment, wgcna
from gcen2.plot._common import figure, finish

OUT = Path(__file__).resolve().parent / "results"
R = Path("runs/case")
FMT = ["pdf", "png"]

wgcna.module_trait(R / "trait_cor.tsv", annotate="stars", output=OUT / "fig_module_trait",
                   formats=FMT)

go = pd.read_csv(OUT / "module_go.tsv", sep="\t")
# modules of blood cell types and processes, in the order of the text
cell_modules = [2, 7, 11, 3, 8, 18, 19, 16, 6, 13, 10]
tabs = {str(m): go[go["module"] == m] for m in cell_modules}
enrichment.enrich_dot(tabs, modules=cell_modules, top=2, max_modules=len(cell_modules),
                      output=OUT / "fig_enrich_dot", formats=FMT)

wgcna.hub_network(R / "tom_hubs.tsv", R / "kme.tsv", module="18",
                  output=OUT / "fig_hub_ME18", formats=FMT)

cv = pd.read_csv(OUT / "rwr_cv.tsv", sep="\t")
labels = {"rwr": "RWR", "rwr_weighted": "RWR, |r| weights", "neighbours": "neighbour voting",
          "rwr_shuffled": "RWR, shuffled network", "degree": "degree only"}
fig = figure("single", 60)
ax = fig.add_subplot(111)
for i, m in enumerate(labels):      # median and quartiles, 5th-95th percentile whiskers
    q5, q25, q50, q75, q95 = cv[m].quantile([0.05, 0.25, 0.5, 0.75, 0.95])
    ax.hlines(i, q5, q95, color="black", lw=0.6)
    ax.hlines(i, q25, q75, color="black", lw=4)
    ax.plot(q50, i, "o", color="white", ms=3, mec="black", mew=0.6)
ax.set_yticks(range(len(labels)), list(labels.values()))
ax.axvline(0.5, color="grey", lw=0.5, ls="--")
ax.invert_yaxis()
ax.set_xlabel(f"AUROC per GO term (n = {len(cv):,}, 3-fold CV)")
finish(fig, OUT / "fig_rwr_cv", FMT)
