"""PyWGCNA 2.2.1 on a GTEx matrix, with WGCNA's settings of ../wgcna_gtex.

usage: python 01_pywgcna.py EXPR_TSV OUT_DIR [tom|-] [POWER]

EXPR_TSV is genes x samples (as for WGCNA and GCEN2). PyWGCNA's own pipeline
(WGCNA.preprocess + WGCNA.findModules: pickSoftThreshold, adjacency, TOM,
average linkage, cutreeHybrid, mergeCloseModules) with a signed network and
signed TOM, the power vector of WGCNA's run (1-10, 12-20) and WGCNA's
RsquaredCut 0.85 (PyWGCNA's default is 0.9), minModuleSize 30,
MEDissThres 0.25 and deepSplit 2 (findModules' default cutreeHybrid
arguments otherwise: pamRespectsDendro False). TPMcutoff is switched off so
that PyWGCNA keeps every gene of the input (the matrices are already filtered,
log-transformed and, for most runs, have 5 PCs regressed out); no sample is
cut (PyWGCNA's default). With POWER, only that power is offered to
pickSoftThreshold, which fixes it.

Writes to OUT_DIR: modules.tsv (gene, module; grey = unassigned), info.tsv
(power, modules, genes in modules, seconds per step) and, with "tom",
tom_top.tsv (the 300,000 gene pairs with the largest TOM, as 02_wgcna.R of
../wgcna_gtex writes).
"""

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import PyWGCNA  # noqa: E402

expr_file, out = sys.argv[1], Path(sys.argv[2])
want_tom = len(sys.argv) > 3 and sys.argv[3] == "tom"
powers = [int(sys.argv[4])] if len(sys.argv) > 4 else list(range(1, 11)) + list(range(12, 21, 2))
out.mkdir(parents=True, exist_ok=True)

t0 = time.perf_counter()
x = pd.read_csv(expr_file, sep="\t", index_col=0)
w = PyWGCNA.WGCNA(name="gtex", geneExp=x.T, TPMcutoff=-np.inf,
                  powers=powers, RsquaredCut=0.85,
                  networkType="signed", TOMType="signed",
                  minModuleSize=30, MEDissThres=0.25, outputPath=str(out) + "/")
t_read = time.perf_counter() - t0

t0 = time.perf_counter()
w.preprocess(show=False)
t_pre = time.perf_counter() - t0
t0 = time.perf_counter()
w.findModules()
t_mod = time.perf_counter() - t0
plt.close("all")

colors = w.datExpr.var["moduleColors"]
pd.DataFrame({"gene": colors.index, "module": colors.to_numpy()}).to_csv(
    out / "modules.tsv", sep="\t", index=False)
if want_tom:
    tom = w.TOM.to_numpy()
    iu = np.triu_indices_from(tom, k=1)
    v = tom[iu]
    top = np.argsort(-v, kind="stable")[:300000]
    g = w.TOM.columns.to_numpy()
    pd.DataFrame({"gene_a": g[iu[0][top]], "gene_b": g[iu[1][top]], "tom": v[top]}).to_csv(
        out / "tom_top.tsv", sep="\t", index=False, float_format="%.9g")
n_mod = int((colors.unique() != "grey").sum())
info = {"samples": w.datExpr.shape[0], "genes": w.datExpr.shape[1], "power": int(w.power),
        "modules": n_mod, "module_genes": int((colors != "grey").sum()),
        "t_read": round(t_read, 2), "t_preprocess": round(t_pre, 2),
        "t_findModules": round(t_mod, 2)}
pd.DataFrame([info]).to_csv(out / "info.tsv", sep="\t", index=False)
print(info)
