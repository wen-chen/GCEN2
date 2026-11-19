"""Compare gcen2 spatial-modules with hdWGCNA (02_hdwgcna.R) per sample:
genes, metaspots, datExpr, modules after ConstructNetwork (gcen2
--no-reassign) and after ModuleConnectivity's reassignment (adjusted Rand
index over all genes), spot-level eigengenes and kME, and time / peak
memory (runs/time.txt, median of the repeats).

usage: python 04_compare.py RUNS_DIR NAME [NAME ...]
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

runs = Path(sys.argv[1])
times = pd.read_csv(runs / "time.txt", sep=" ", header=None,
                    names=["tool", "name", "rep", "s", "_s", "kb", "_kb"])


def gcen2_labels(d: Path, genes) -> pd.Series:
    m = pd.read_csv(d / "modules.tsv", sep="\t")
    lab = pd.Series("NA", index=genes, dtype=object)
    lab[m["gene"].to_numpy()] = m["module"].astype(str).to_numpy()
    return lab


rows = []
for name in sys.argv[2:]:
    g, h = runs / "gcen2" / name, runs / "hdwgcna" / name
    hd = pd.read_csv(h / "datExpr.tsv", sep="\t", index_col=0)
    gx = pd.read_csv(g / "metaspots.tsv", sep="\t", index_col=0).T
    gk = pd.read_csv(g / "spot_kme.tsv", sep="\t", index_col=0, keep_default_na=False)
    genes = list(pd.read_csv(g / "kme.tsv", sep="\t", index_col=0).index)
    glab = gcen2_labels(g, genes)
    glab0 = gcen2_labels(runs / "gcen2" / f"{name}_noreassign", genes)
    hmod = pd.read_csv(h / "modules.tsv", sep="\t").set_index("gene")
    hm = hmod["module"].astype(str)
    hm0 = pd.read_csv(h / "modules_network.tsv", sep="\t").set_index("gene")["module"]
    gs = pd.read_csv(g / "spot_eigengenes.tsv", sep="\t", index_col=0)
    hs = pd.read_csv(h / "MEs.tsv", sep="\t", index_col=0).drop(columns="grey",
                                                                  errors="ignore")
    gs = gs.loc[hs.index]
    # each hdWGCNA module <-> the gcen2 module with the same genes (ConstructNetwork)
    match = {}
    for m in sorted(set(hm0) - {"grey"}):
        match[m] = glab0[hm0.index[hm0 == m]].value_counts().index[0]
    me_diff = max(float(np.abs(gs[f"ME{match[m]}"] - hs[m]).max()) for m in match)
    me_cor = min(float(np.corrcoef(gs[f"ME{match[m]}"], hs[m])[0, 1]) for m in match)
    kme_diff = max(float(np.abs(gk.loc[genes, f"kME{match[m]}"].astype(float)
                                - hmod.loc[genes, f"kME_{m}"]).max()) for m in match)
    t = times[times["name"] == name].groupby("tool")[["s", "kb"]].median()
    rows.append({
        "sample": name, "genes": len(genes), "same_genes": list(hd.columns) == genes,
        "metaspots": len(hd), "same_metaspots": list(hd.index) == list(gx.index),
        "max_abs_datExpr_diff": float(np.abs(gx.loc[hd.index, hd.columns].to_numpy()
                                             - hd.to_numpy()).max()),
        "modules": len(match),
        "ARI_network": adjusted_rand_score(hm0.loc[genes].astype(str), glab0.loc[genes]),
        "reassigned_genes": int((glab != glab0).sum()),
        "ARI_final": adjusted_rand_score(hm.loc[genes], glab.loc[genes]),
        "spot_ME_min_r": me_cor, "spot_ME_max_abs_diff": me_diff,
        "spot_kME_max_abs_diff": kme_diff,
        "gcen2_s": t.loc["gcen2", "s"], "gcen2_GB": t.loc["gcen2", "kb"] / 2**20,
        "hdwgcna_s": t.loc["hdwgcna", "s"], "hdwgcna_GB": t.loc["hdwgcna", "kb"] / 2**20,
    })
out = pd.DataFrame(rows)
out["speedup"] = out["hdwgcna_s"] / out["gcen2_s"]
out["memory_ratio"] = out["hdwgcna_GB"] / out["gcen2_GB"]
with pd.option_context("display.float_format", "{:.6g}".format):
    print(out.T.to_string())
out.to_csv(runs / "compare.tsv", sep="\t", index=False)
