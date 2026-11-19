"""The hdWGCNA module workflow under metacell seeds, over all genes.

Run from the PBMC data directory after 05_tom_modules.py. For the three cell
types run with metacell seeds 42 (the default), 1 and 2 (B, CD14 Mono, CD4 T):
the adjusted Rand index over all genes (unassigned genes as one group, as in
05_tom_modules.py) of

  hdwgcna_seeds   hdWGCNA's modules between its seeds (runs/hdwgcna/<run>/modules.tsv)
  gcen2_seeds     GCEN2's (`gcen2 module -e --preset hdwgcna` on GCEN2's own
                  metacells, work/tom_modules/<run>_own.modules.tsv) between seeds
  between_tools   GCEN2 on its own metacells against hdWGCNA, all 3 x 3 seed pairs

seed_stability.tsv has the ARI over the genes in modules in both runs, which is
higher. Writes results/seed_ari_all_genes.tsv.
"""

from itertools import combinations, product
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path.cwd()
TYPES = ["B", "CD14_Mono", "CD4_T"]
SEEDS = {"42": "", "1": "_seed1", "2": "_seed2"}


def hdwgcna(run):
    m = pd.read_csv(ROOT / "runs/hdwgcna" / run / "modules.tsv", sep="\t").set_index("gene")
    return m["module"].where(m["module"] != "grey", "unassigned").astype(str)


def gcen2(run, genes):
    m = pd.read_csv(ROOT / "work/tom_modules" / f"{run}_own.modules.tsv", sep="\t")
    m = m.set_index("gene")["module"].astype(str).reindex(genes)
    return m.fillna("unassigned")


rows = []
for t in TYPES:
    h = {s: hdwgcna(t + suf) for s, suf in SEEDS.items()}
    genes = h["42"].index
    g = {s: gcen2(t + suf, genes) for s, suf in SEEDS.items()}
    for (a, b) in combinations(SEEDS, 2):
        rows.append({"cell_type": t, "comparison": "hdwgcna_seeds", "seeds": f"{a}-{b}",
                     "ARI": adjusted_rand_score(h[a].reindex(genes), h[b].reindex(genes))})
        rows.append({"cell_type": t, "comparison": "gcen2_seeds", "seeds": f"{a}-{b}",
                     "ARI": adjusted_rand_score(g[a], g[b])})
    for a, b in product(SEEDS, SEEDS):
        rows.append({"cell_type": t, "comparison": "between_tools", "seeds": f"{a}-{b}",
                     "ARI": adjusted_rand_score(g[a], h[b].reindex(genes))})
df = pd.DataFrame(rows)
df["ARI"] = df["ARI"].round(3)
out = Path(__file__).resolve().parent / "results/seed_ari_all_genes.tsv"
df.to_csv(out, sep="\t", index=False)
print(df.pivot_table(index="cell_type", columns="comparison", values="ARI",
                     aggfunc=["median", "min", "max"]).round(2).to_string())
