"""GCEN2 against hdWGCNA's consensus workflow with harmonized eigengenes.

1. Genes: hdWGCNA's (SetupForWGCNA, then goodGenesMS over the sets) and
   GCEN2's (the same fraction rule; genes constant in a set unassigned).
2. Network step: `gcen2.consensus.consensus_tom_modules` with
   HDWGCNA_CONSENSUS on hdWGCNA's own metacell sets (multiExpr_<set>.tsv, in
   the order of sets.tsv): the powers and the modules against
   modules_network.tsv (ARI), and the time against TestSoftPowersConsensus +
   ConstructNetwork (t_network).
3. Cell level, on hdWGCNA's modules: `cell_eigengenes` (batches = the
   harmonizing column) against MEs.tsv (max |difference|) and hMEs.tsv
   (Pearson r per module), next to GCEN2 against itself with other Harmony
   seeds (hdWGCNA's own harmonized eigengenes do not vary: Seurat's RunPCA
   sets the seed to 42 before each Harmony run); kME against
   modules_final.tsv.
4. FindAllDMEs and ModuleTraitCorrelation on hdWGCNA's harmonized
   eigengenes: `find_all_dmes` and `module_trait_cor` against DMEs.tsv and
   trait_cor.tsv (hdWGCNA's p-values, 2 * (1 - pt(...)), lose digits below
   about 1e-12 and are 0 below about 1e-16; compared where above 1e-12).
5. Whole workflow (each tool with its own metacells): ARI of the final
   modules between the tools and between seeds of each tool.

usage: python 03_compare.py WORK_PREFIX RUNS_DIR [TIME_FILE]
(from the consensus data directory: python 03_compare.py nk runs runs/time.txt)
Writes RUNS_DIR/compare.tsv (metric, value).
"""

import itertools
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

from gcen2.consensus import HDWGCNA_CONSENSUS, consensus_tom_modules
from gcen2.single_cell import cell_eigengenes

prefix, runs = sys.argv[1], Path(sys.argv[2])
time_file = Path(sys.argv[3]) if len(sys.argv) > 3 else None
hd, gc = runs / "hdwgcna", runs / "gcen2"
harm_by = "donor_id"
out = {}

adata = ad.read_h5ad(prefix + ".h5ad")
ctype = str(adata.obs["cell_type"].iloc[0])
gdir = next(p for p in gc.iterdir() if p.is_dir())

# 1. genes
hd_genes = Path(hd / "genes.txt").read_text().split()
g_kme = pd.read_csv(gdir / "kme.tsv", sep="\t", index_col=0)
g_all = list(g_kme.index)
g_good = list(g_kme.index[g_kme.drop(columns="module").notna().any(axis=1)]) \
    if g_kme.shape[1] > 1 else g_all
out["genes_hdwgcna"] = len(hd_genes)
out["genes_gcen2"] = len(g_good)
out["same_genes"] = sorted(hd_genes) == sorted(g_good)

# 2. network step on hdWGCNA's metacells
sets = pd.read_csv(hd / "sets.tsv", sep="\t")
mats = []
for s in sets["set"]:
    m = pd.read_csv(hd / f"multiExpr_{''.join(c if c.isalnum() else '_' for c in s)}.tsv",
                    sep="\t", index_col=0)
    mats.append(m.T)
t0 = time.time()
res = consensus_tom_modules(mats, set_names=list(sets["set"]), threads=8, **HDWGCNA_CONSENSUS)
out["network_step_gcen2_s"] = round(time.time() - t0, 2)
out["powers_hdwgcna"] = ",".join(str(p) for p in sets["power"])
out["powers_gcen2"] = ",".join(f"{p:g}" for p in res.metadata["set_powers"].values())
net = pd.read_csv(hd / "modules_network.tsv", sep="\t").set_index("gene")["module"].astype(str)
lab = pd.Series(res.labels, index=res.genes).loc[net.index]
out["network_step_ARI"] = round(adjusted_rand_score(net, lab), 6)
out["network_modules"] = f"{len(set(net) - {'grey'})} / {int(res.labels.max())}"
timing = pd.read_csv(hd / "timing.tsv", sep="\t").iloc[0]
out["network_step_hdwgcna_s"] = round(float(timing["t_network"]), 2)

# 3. cell level on hdWGCNA's modules
me_hd = pd.read_csv(hd / "MEs.tsv", sep="\t", index_col=0)
hme_hd = pd.read_csv(hd / "hMEs.tsv", sep="\t", index_col=0)
colors = [c for c in me_hd.columns if c != "grey"]
code = {c: i + 1 for i, c in enumerate(colors)}
genes = list(net.index)
var = {g: i for i, g in enumerate(adata.var_names)}
labels = np.array([code.get(m, 0) for m in net], dtype=np.int64)
cells = adata.obs_names
batches = adata.obs[harm_by].astype(str).to_numpy()
t0 = time.time()
me, kme, hme = cell_eigengenes(adata.X, [var[g] for g in genes], labels,
                               list(range(1, len(colors) + 1)), batches=batches, seed=0)
out["cell_level_gcen2_s"] = round(time.time() - t0, 2)
me = pd.DataFrame(me, index=cells, columns=colors)
hme = pd.DataFrame(hme, index=cells, columns=colors)
out["ME_max_abs_diff"] = float(np.abs(me.to_numpy() - me_hd.loc[cells, colors].to_numpy()).max())


def rmin(a, b):
    return min(np.corrcoef(a[c], b.loc[a.index, c])[0, 1] for c in colors)


out["hME_r_gcen2_vs_hdwgcna_min"] = round(rmin(hme, hme_hd), 4)
out["hME_r_gcen2_vs_hdwgcna"] = ",".join(
    f"{np.corrcoef(hme[c], hme_hd.loc[cells, c])[0, 1]:.4f}" for c in colors)
rs = []
for s in (1, 2):  # Harmony's own variation: other seeds, the same modules
    _, _, h = cell_eigengenes(adata.X, [var[g] for g in genes], labels,
                              list(range(1, len(colors) + 1)), batches=batches, seed=s)
    h = pd.DataFrame(h, index=cells, columns=colors)
    rs.append(rmin(hme, h))
    out[f"hME_r_gcen2_seed{s}_vs_hdwgcna_min"] = round(rmin(h, hme_hd), 4)
out["hME_r_gcen2_seeds_min"] = ",".join(f"{v:.4f}" for v in rs)
out["hME_r_ME_min"] = round(rmin(hme, me), 4)  # how much harmonizing changes the eigengenes
final = pd.read_csv(hd / "modules_final.tsv", sep="\t").set_index("gene")
kme_hd = final.loc[genes, [f"kME_{c}" for c in colors]].to_numpy()
out["kME_max_abs_diff"] = round(float(np.nanmax(np.abs(kme - kme_hd))), 4)
out["kME_r"] = round(float(np.corrcoef(kme.ravel(), kme_hd.ravel())[0, 1]), 5)

# 4. DMEs and module-trait correlations on hdWGCNA's harmonized eigengenes
from gcen2.single_cell import find_all_dmes, module_trait_cor  # noqa: E402

meta = adata.obs
d = find_all_dmes(hme_hd.loc[cells], meta["dataset"].astype(str).to_numpy())
r = pd.read_csv(hd / "DMEs.tsv", sep="\t", index_col=0)
r.index = r.index.str.strip()
out["DME_rows"] = f"{len(d)} / {len(r)}"
out["DME_same_order"] = list(d.index) == list(r.index)
d = d.loc[r.index]
for c in ("p_val", "avg_log2FC", "p_val_adj"):
    a, b = d[c].to_numpy(float), r[c].to_numpy(float)
    out[f"DME_{c}_max_rel_diff"] = float(np.max(np.abs(a - b) / np.maximum(np.abs(b), 1e-300)))
out["DME_pct_same"] = bool((d[["pct.1", "pct.2"]].to_numpy(float)
                            == r[["pct.1", "pct.2"]].to_numpy(float)).all())
t = meta[["age", "BMI", "sex"]].copy()
m = module_trait_cor(hme_hd.loc[cells], t, groups=meta["dataset"].astype(str).to_numpy())
rt = pd.read_csv(hd / "trait_cor.tsv", sep="\t")
for c in ("group", "trait", "module"):
    rt[c] = rt[c].astype(str).str.strip()
for c in ("cor", "p", "fdr"):  # R's format() pads NA
    rt[c] = pd.to_numeric(rt[c].astype(str).str.strip(), errors="coerce")
x = m.merge(rt, on=["group", "trait", "module"], suffixes=("", "_r"))
out["trait_rows"] = f"{len(x)} / {len(rt)}"
out["trait_cor_max_abs_diff"] = float(np.nanmax(np.abs(x["cor"] - x["cor_r"])))
ok = x["p_r"] > 1e-12
out["trait_p_max_rel_diff"] = float(np.nanmax(np.abs(x["p"] - x["p_r"])[ok] / x["p_r"][ok]))
out["trait_p_r_zero"] = int((x["p_r"] == 0).sum())
okf = x["fdr_r"] > 1e-12
out["trait_fdr_max_rel_diff"] = float(np.nanmax(np.abs(x["fdr"] - x["fdr_r"])[okf]
                                                / x["fdr_r"][okf]))
out["stats_hdwgcna_s"] = round(float(timing["t_stats"]), 2)

# 5. whole workflow


def gcen2_final(d):
    sub = next(p for p in d.iterdir() if p.is_dir())
    m = pd.read_csv(sub / "modules.tsv", sep="\t").set_index("gene")["module"]
    return m


def hd_final(d):
    f = pd.read_csv(d / "modules_final.tsv", sep="\t").set_index("gene")["module"].astype(str)
    return f


def ari(a, b, universe):
    return adjusted_rand_score(a.reindex(universe).astype(str).fillna("none"),
                               b.reindex(universe).astype(str).fillna("none"))


universe = hd_genes
g_runs = [gc] + sorted(runs.glob("gcen2_seed*"))
h_runs = [hd] + sorted(runs.glob("hdwgcna_seed*"))
gm = [gcen2_final(d) for d in g_runs]
hm = [hd_final(d).replace("grey", np.nan) for d in h_runs]
out["ARI_gcen2_vs_hdwgcna"] = ",".join(f"{ari(a, b, universe):.3f}" for a in gm for b in hm)
out["ARI_gcen2_seeds"] = ",".join(f"{ari(a, b, universe):.3f}"
                                  for a, b in itertools.combinations(gm, 2))
out["ARI_hdwgcna_seeds"] = ",".join(f"{ari(a, b, universe):.3f}"
                                    for a, b in itertools.combinations(hm, 2))

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
