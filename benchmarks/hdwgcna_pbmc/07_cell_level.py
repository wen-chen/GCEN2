"""Cell-level eigengenes and kME: GCEN2 against hdWGCNA's ModuleEigengenes and
ModuleConnectivity on the same modules.

For each cell type, hdWGCNA's run (02_hdwgcna.R in RUNS/<type>) gives its
modules after ConstructNetwork (modules.tsv), its eigengenes of the type's
cells (MEs.tsv), and its modules after ModuleConnectivity with the kME
(modules_final.tsv). GCEN2's `gcen2.single_cell.cell_eigengenes` is run on
the same cells (work/pbmc.h5ad, counts) with hdWGCNA's modules, then
hdWGCNA's ReassignModules rule (`reassign_negative_kme`). Compared: the
eigengenes, the kME and the final module of every gene.

usage: python 07_cell_level.py RUNS_DIR (from the PBMC data directory)
Writes RUNS_DIR/cell_level.tsv.
"""

import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from gcen2.single_cell import cell_eigengenes
from gcen2.spatial.modules import reassign_negative_kme

runs = Path(sys.argv[1])
adata = ad.read_h5ad("work/pbmc.h5ad")
var = {g: i for i, g in enumerate(adata.var_names)}
rows = []
for ct in ("CD14 Mono", "CD4 T", "B", "CD8 T", "NK", "FCGR3A Mono"):
    d = runs / ct.replace(" ", "_")
    net = pd.read_csv(d / "modules.tsv", sep="\t").set_index("gene")["module"].astype(str)
    final = pd.read_csv(d / "modules_final.tsv", sep="\t").set_index("gene")
    me_hd = pd.read_csv(d / "MEs.tsv", sep="\t", index_col=0)
    colors = [c for c in me_hd.columns if c != "grey"]
    code = {c: i + 1 for i, c in enumerate(colors)}
    genes = list(net.index)
    labels = np.array([code.get(m, 0) for m in net], dtype=np.int64)
    cells = np.flatnonzero((adata.obs["cell_type"] == ct).to_numpy())
    me, kme = cell_eigengenes(adata.X[cells], [var[g] for g in genes], labels,
                              list(range(1, len(colors) + 1)))
    me = pd.DataFrame(me, index=adata.obs_names[cells], columns=colors)
    me_hd = me_hd.loc[me.index, colors]
    kme_hd = final.loc[genes, [f"kME_{c}" for c in colors]].to_numpy()
    new, moved = reassign_negative_kme(labels, kme, list(range(1, len(colors) + 1)))
    names = np.array(["grey"] + colors, dtype=object)
    same_final = bool((names[new] == final.loc[genes, "module"].astype(str).to_numpy()).all())
    rows.append({
        "cell_type": ct, "cells": len(cells), "genes": len(genes), "modules": len(colors),
        "ME_max_abs_diff": float(np.abs(me.to_numpy() - me_hd.to_numpy()).max()),
        "ME_max_abs": float(np.abs(me_hd.to_numpy()).max()),
        "kME_max_abs_diff": float(np.nanmax(np.abs(kme - kme_hd))),
        "reassigned": len(moved), "same_final_modules": same_final,
    })
    print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv(runs / "cell_level.tsv", sep="\t", index=False)
print(df.to_string(index=False))
