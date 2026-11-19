"""GCEN2 against hdWGCNA's ProjectModules and ModulePreservation.

1. Projection on hdWGCNA's reference modules: `gcen2.single_cell.project_modules`
   (grey included, as hdWGCNA) on the query cells against projected_MEs.tsv.
2. Preservation on hdWGCNA's own metacells and modules:
   `gcen2.preservation.module_preservation` (same permutations and seed)
   against preservation_obs.tsv / preservation_Z.tsv, and the time against
   ModulePreservation (t_preservation).
3. Whole workflow (each tool with its own metacells, modules and query
   metacells): every GCEN2 module matched to the hdWGCNA module with which it
   shares the most genes; their Zsummary.pres side by side.

usage: python 03_compare.py NK_H5AD RUNS_DIR [N_PERM] [TIME_FILE]
(from the consensus data directory: python 03_compare.py nk.h5ad pres 500 pres/time.txt)
Writes RUNS_DIR/compare.tsv (metric, value) and RUNS_DIR/zsummary.tsv.
"""

import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from gcen2.preservation import module_preservation
from gcen2.single_cell import project_modules

h5ad, runs = sys.argv[1], Path(sys.argv[2])
n_perm = int(sys.argv[3]) if len(sys.argv) > 3 else 500
time_file = Path(sys.argv[4]) if len(sys.argv) > 4 else None
hd = runs / "hdwgcna"
out = {}


def load(f):
    d = pd.read_csv(f, sep="\t", index_col=0)
    d.index = d.index.astype(str).str.strip()
    d = d.loc[:, [c for c in d.columns if c != "moduleSize.1"]]  # R's cbind repeats it
    return d.apply(lambda c: pd.to_numeric(c.astype(str).str.strip(), errors="coerce"))


timing = pd.read_csv(hd / "timing.tsv", sep="\t").iloc[0]
mods = pd.read_csv(hd / "ref_modules.tsv", sep="\t").set_index("gene")["module"].astype(str)
out["genes"] = int(timing["genes"])
out["modules_hdwgcna"] = int(timing["modules"])

# 1. projection
adata = ad.read_h5ad(h5ad)
query_set = "Krasnow_2020"
cells = (adata.obs["dataset"].astype(str) == query_set).to_numpy()
t0 = time.time()
me, _, table = project_modules(adata, mods, cells=cells, include_unassigned=True)
out["project_gcen2_s"] = round(time.time() - t0, 2)
out["project_hdwgcna_s"] = round(float(timing["t_project"]), 2)
pm = load(hd / "projected_MEs.tsv")
cols = [c for c in pm.columns if c in me.columns]
out["projected_modules"] = f"{len(cols)} / {pm.shape[1]}"
out["projected_ME_max_abs_diff"] = float(np.abs(me.loc[pm.index, cols].to_numpy()
                                                - pm[cols].to_numpy()).max())

# 2. preservation on hdWGCNA's metacells
ref = pd.read_csv(hd / "ref_datExpr.tsv", sep="\t", index_col=0)
qry = pd.read_csv(hd / "query_datExpr.tsv", sep="\t", index_col=0)
t0 = time.time()
obs, z = module_preservation(ref, qry, mods, n_permutations=n_perm, threads=8)
out["preservation_gcen2_s"] = round(time.time() - t0, 1)
out["preservation_hdwgcna_s"] = round(float(timing["t_preservation"]), 1)
ro, rz = load(hd / "preservation_obs.tsv"), load(hd / "preservation_Z.tsv")
for name, a, b in (("obs", obs, ro), ("Z", z, rz)):
    cols = [c for c in b.columns if c in a.columns]
    out[f"{name}_same_modules"] = sorted(a.index) == sorted(b.index)
    out[f"{name}_columns"] = f"{len(cols)} / {b.shape[1]}"
    x, y = a.loc[b.index, cols].to_numpy(float), b[cols].to_numpy(float)
    out[f"{name}_max_rel_diff"] = float(np.nanmax(np.abs(x - y) / np.maximum(np.abs(y), 1e-12)))
out["Zsummary_hdwgcna"] = ",".join(f"{m}:{v:.2f}" for m, v in rz["Zsummary.pres"].items())
out["Zsummary_gcen2_same_input"] = ",".join(f"{m}:{v:.2f}"
                                            for m, v in z["Zsummary.pres"].items())

# 3. whole workflow: GCEN2's own modules and metacells
gdir = runs / "gcen2"
gref = next(p for p in (gdir / "ref").iterdir() if p.is_dir())
gmods = pd.read_csv(gref / "modules.tsv", sep="\t")
gz = pd.read_csv(gdir / "preservation" / "preservation_Z.tsv", sep="\t", index_col=0)
rows = []
for k, g in gmods.groupby("module"):
    genes = set(g["gene"])
    overlap = {m: len(genes & set(mods.index[mods == m])) for m in set(mods) - {"grey"}}
    best = max(overlap, key=overlap.get)
    union = len(genes | set(mods.index[mods == best]))
    rows.append({"gcen2_module": f"ME{k}", "genes": len(genes), "hdwgcna_module": best,
                 "jaccard": round(overlap[best] / union, 3),
                 "Zsummary_gcen2": gz.loc[f"ME{k}", "Zsummary.pres"],
                 "Zsummary_hdwgcna": rz.loc[best, "Zsummary.pres"]})
zs = pd.DataFrame(rows)
zs.to_csv(runs / "zsummary.tsv", sep="\t", index=False, float_format="%.4g")
print(zs.to_string(index=False))
out["modules_gcen2"] = len(zs)

if time_file is not None and time_file.exists():
    t = pd.read_csv(time_file, sep=" ", header=None, usecols=[0, 2, 4],
                    names=["tool", "s", "kb"])
    for tool, g in t.groupby("tool"):
        out[f"{tool}_s"] = round(g["s"].mean(), 1)
        out[f"{tool}_GB"] = round(g["kb"].mean() / 1024 ** 2, 2)
    if {"gcen2_s", "hdwgcna_s"} <= set(out):
        out["speedup"] = round(out["hdwgcna_s"] / out["gcen2_s"], 1)
        out["memory_ratio"] = round(out["hdwgcna_GB"] / out["gcen2_GB"], 1)

df = pd.DataFrame({"metric": list(out), "value": [str(v) for v in out.values()]})
df.to_csv(runs / "compare.tsv", sep="\t", index=False)
print(df.to_string(index=False))
