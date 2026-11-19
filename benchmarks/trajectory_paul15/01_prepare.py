"""Paul et al. 2015 mouse myeloid differentiation (2,730 cells) with a
diffusion pseudotime, as in the scanpy DPT tutorial.

Run from a data directory of your choice (scanpy downloads its copy of the
data to data/paul15/ if it is not there); outputs in work/:
  paul15_pt.h5ad     counts per 10,000 of the 2,000 most variable genes, cells
                     x genes, with obs["dpt_pseudotime"] (GCEN2's input)
  paul15_expr.tsv    the same matrix, genes x cells, 17 significant digits
  paul15_pt.tsv      cell, pseudotime (the R implementation's input)
"""

from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sc.settings.datasetdir = Path("data")
WORK = Path("work")
WORK.mkdir(exist_ok=True)

adata = sc.datasets.paul15()
adata.X = adata.X.astype(np.float64)
sc.pp.normalize_total(adata, target_sum=1e4)
x = adata.copy()
sc.pp.log1p(x)
sc.pp.highly_variable_genes(x, n_top_genes=2000, flavor="seurat")
sc.pp.pca(x, n_comps=50, random_state=0)
sc.pp.neighbors(x, n_neighbors=20, random_state=0)
sc.tl.diffmap(x, random_state=0)
x.uns["iroot"] = int(np.flatnonzero(x.obs["paul15_clusters"] == "7MEP")[0])
sc.tl.dpt(x)

out = adata[:, x.var["highly_variable"].to_numpy()].copy()
out.obs = pd.DataFrame({"dpt_pseudotime": x.obs["dpt_pseudotime"].to_numpy(),
                        "cluster": x.obs["paul15_clusters"].astype(str).to_numpy()},
                       index=out.obs_names)
out.uns = {}
out.write_h5ad(WORK / "paul15_pt.h5ad")
pd.DataFrame(out.X.T, index=out.var_names, columns=out.obs_names).to_csv(
    WORK / "paul15_expr.tsv", sep="\t", float_format="%.17g")
out.obs[["dpt_pseudotime"]].to_csv(WORK / "paul15_pt.tsv", sep="\t", float_format="%.17g")
print(f"{out.n_obs} cells x {out.n_vars} genes; "
      f"{int(np.isfinite(out.obs['dpt_pseudotime']).sum())} with a finite pseudotime")
