"""gcen2 on one PBMC cell type, with the hdWGCNA settings (k = 25,
max_shared = 10, at most 1000 metacells, the shared PCA and gene list).

usage: 03_gcen2.py CELL_TYPE OUT_DIR [CELLS_FILE] [SEED] [KERNEL]

CELLS_FILE "-": all cells of the type; SEED default 42; KERNEL the correlation
kernel of gcen2.network.build (auto, the default; pairwise: GCEN v1's).

Writes to OUT_DIR: metacells.tsv (metacells x genes, log-normalised),
network.tsv (edges), modules.tsv (gene, module), members.tsv (metacell
member barcodes), timing.tsv.
"""

import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

import gcen2
from gcen2.single_cell import build_by_celltype, metacells

WORK = Path("work")  # run from the PBMC data directory
ct, out = sys.argv[1], Path(sys.argv[2])
cells_file = sys.argv[3] if len(sys.argv) > 3 and sys.argv[3] != "-" else None
seed = int(sys.argv[4]) if len(sys.argv) > 4 else 42
kernel = sys.argv[5] if len(sys.argv) > 5 else "auto"
out.mkdir(parents=True, exist_ok=True)
THREADS = 8
PARAMS = dict(k=25, max_shared=10, max_metacells=1000, min_cells=100, use_rep="X_pca",
              n_pcs=30, seed=seed)

t0 = time.perf_counter()
adata = ad.read_h5ad(WORK / "pbmc.h5ad")
genes = (WORK / f"genes_{ct.replace(' ', '_')}.txt").read_text().split()
cells = (Path(cells_file).read_text().split() if cells_file
         else adata.obs_names[adata.obs["cell_type"] == ct].tolist())
sub = adata[cells, genes].copy()
del adata
t_read = time.perf_counter() - t0

# the network as gcen2 celltype builds it (metacells, effective-sample p-values, FDR)
t0 = time.perf_counter()
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    net = build_by_celltype(sub, "cell_type", method="pearson", cor=0.3, pval=0.001,
                            signed=True, threads=THREADS, kernel=kernel, **PARAMS)[ct]
t_network = time.perf_counter() - t0
t0 = time.perf_counter()
mods = gcen2.module.identify(net, similarity=0.5, threads=THREADS)
t_modules = time.perf_counter() - t0

# the metacell matrix (same seed, same metacells) for the correlation comparisons
meta = metacells(sub, groupby="cell_type", threads=THREADS, **PARAMS)
pd.DataFrame(meta.X, index=meta.obs_names, columns=meta.var_names).to_csv(
    out / "metacells.tsv", sep="\t", float_format="%.6g")
members = np.asarray(sub.obs_names)[meta.uns["metacell_members"]]
pd.DataFrame({"metacell": meta.obs_names,
              "cells": [",".join(m) for m in members]}).to_csv(out / "members.tsv", sep="\t",
                                                               index=False)
net.edges.to_csv(out / "network.tsv", sep="\t", index=False)
rows = [(g, i) for i, m in enumerate(mods.modules) for g in sorted(m)]
pd.DataFrame(rows, columns=["gene", "module"]).to_csv(out / "modules.tsv", sep="\t",
                                                      index=False)
pd.DataFrame([{"cells": len(cells), "genes": len(genes), "metacells": meta.n_obs,
               "p_value_samples": net.metadata["p_value_samples"],
               "edges": net.edge_count, "modules": len(mods.modules),
               "t_read": t_read, "t_network": t_network, "t_modules": t_modules}]).to_csv(
    out / "timing.tsv", sep="\t", index=False)
print(f"{ct}: {len(cells)} cells, {len(genes)} genes, {meta.n_obs} metacells, "
      f"{net.edge_count} edges, {len(mods.modules)} modules; network {t_network:.1f}s, "
      f"modules {t_modules:.1f}s")
