"""knn metaspots (no grid needed) against hdWGCNA's grid metaspots on the two
Visium samples.

For each sample, against the grid run of run_all.sh (runs/gcen2/<sample>,
the same modules as hdWGCNA):
  * knn: `gcen2 spatial-modules --metaspots knn` (k = 7, max_shared = 2,
    metaspots from obsm["spatial"], the pixel coordinates) with seeds 42, 1
    and 2: module ARI with the grid run, between the knn seeds, and the best
    |r| of each grid spot eigengene with the knn ones;
  * for scale, the grid method itself on 95% of the spots (two random
    subsets): its module ARI with the full grid run.
ARI is over all genes, unassigned genes as one group.

usage: python 06_knn_vs_grid.py RUNS_DIR SAMPLE... (from the spatial data directory)
Writes RUNS_DIR/knn_vs_grid.tsv.
"""

import subprocess
import sys
import warnings
from itertools import combinations
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from gcen2.spatial import spatial_modules

runs = Path(sys.argv[1])
GCEN2 = str(Path(sys.executable).with_name("gcen2"))


def mods(d: Path) -> pd.Series:
    return pd.read_csv(d / "modules.tsv", sep="\t").set_index("gene")["module"].astype(str)


def ari(a: pd.Series, b: pd.Series) -> float:
    g = a.index.union(b.index)
    return round(adjusted_rand_score(a.reindex(g).fillna("unassigned"),
                                     b.reindex(g).fillna("unassigned")), 3)


rows = []
for name in sys.argv[2:]:
    grid = mods(runs / "gcen2" / name)
    me_grid = pd.read_csv(runs / "gcen2" / name / "spot_eigengenes.tsv", sep="\t", index_col=0)
    knn, best_r, n_ms = {}, [], []
    for seed in (42, 1, 2):
        d = runs / "gcen2_knn" / f"{name}_seed{seed}"
        subprocess.run([GCEN2, "spatial-modules", "-i", f"work/{name}.h5ad", "--groupby",
                        "region", "-o", str(d), "-t", "8", "--metaspots", "knn",
                        "--seed", str(seed)], check=True, capture_output=True)
        knn[seed] = mods(d)
        n_ms.append(int(pd.read_csv(d / "summary.tsv", sep="\t")["n_metaspots"].iloc[0]))
        me = pd.read_csv(d / "spot_eigengenes.tsv", sep="\t", index_col=0).reindex(me_grid.index)
        r = np.corrcoef(me_grid.to_numpy().T, me.to_numpy().T)[: me_grid.shape[1],
                                                                me_grid.shape[1]:]
        best_r.append(float(np.abs(r).max(axis=1).min()))
    adata = ad.read_h5ad(f"work/{name}.h5ad")
    drop = []
    for seed in (1, 2):
        keep = np.random.default_rng(seed).random(adata.n_obs) >= 0.05
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            res = spatial_modules(adata[keep].copy(), "region", method="grid", threads=8,
                                  spot_level=False)
        lab = pd.Series(res.labels - 1, index=res.genes)
        drop.append(ari(grid, lab[lab >= 0].astype(str)))
    rows.append({
        "sample": name,
        "grid_metaspots": int(pd.read_csv(runs / "gcen2" / name / "summary.tsv",
                                          sep="\t")["n_metaspots"].iloc[0]),
        "knn_metaspots": ",".join(map(str, n_ms)),
        "ARI_knn_vs_grid": ",".join(str(ari(grid, knn[s])) for s in knn),
        "ARI_knn_seeds": ",".join(str(ari(knn[a], knn[b])) for a, b in combinations(knn, 2)),
        "ARI_grid_95pct_vs_grid": ",".join(map(str, drop)),
        "spot_ME_min_best_r": round(min(best_r), 3),
    })
    print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv(runs / "knn_vs_grid.tsv", sep="\t", index=False)
print(df.to_string(index=False))
