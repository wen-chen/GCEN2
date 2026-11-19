"""The effective sample size of metacell correlations, observed against
predicted.

Run from the PBMC data directory. For each of the six cell types: GCEN2's
metacells (03_gcen2.py's settings) of the type's cells plus 200 null genes
(as in 09_strategies.py: a gene's fractions of the library permuted across
the cells, Poisson counts of fraction x library size). The correlation of a
null gene with a gene y over the metacells has variance s(y) / M, where
s(y) = y'Cy / y'y and C the metacells' overlap (shared cells / k), so y has
M / s(y) effective samples (`gcen2.single_cell.metacell_effective_samples`).
Observed: 1 / var(artanh r) + 3 over the 200 null genes, per gene y, as the
median in four bins of s(y). Writes results/effective_samples.tsv.
"""

import sys
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

from gcen2.single_cell import metacell_effective_samples, metacells

HERE = Path(__file__).resolve().parent
TYPES = sys.argv[1:] or ["CD14 Mono", "CD4 T", "B", "CD8 T", "NK", "FCGR3A Mono"]
META = dict(k=25, max_shared=10, max_metacells=1000, min_cells=100, use_rep="X_pca",
            n_pcs=30, seed=42)
BINS = [(0, 50), (50, 90), (90, 99), (99, 100)]

warnings.simplefilter("ignore")
adata = ad.read_h5ad("work/pbmc.h5ad")
rows = []
for ct in TYPES:
    genes = Path(f"work/genes_{ct.replace(' ', '_')}.txt").read_text().split()
    sub = adata[adata.obs["cell_type"] == ct, genes].copy()
    sub.X = sp.csr_matrix(sub.X)
    rng = np.random.default_rng(0)
    src = rng.choice(len(genes), 200, replace=False)
    lib = np.asarray(sub.X.sum(axis=1)).ravel()
    frac = sub.X[:, src].toarray() / np.maximum(lib, 1)[:, None]
    null = np.column_stack([rng.poisson(frac[rng.permutation(len(lib)), j] * lib)
                            for j in range(200)])
    full = ad.AnnData(sp.hstack([sub.X, sp.csr_matrix(null)]).tocsr(), obs=sub.obs.copy(),
                      obsm={"X_pca": sub.obsm["X_pca"]})
    full.var_names = genes + [f"null_{genes[j]}" for j in src]
    m = metacells(full, groupby="cell_type", threads=8, **META)
    x = np.asarray(m.X, dtype=np.float64).T  # genes x metacells
    pred = metacell_effective_samples(pd.DataFrame(x, index=full.var_names),
                                      m.uns["metacell_members"], 25).to_numpy()[:len(genes)]
    z = x - x.mean(axis=1, keepdims=True)
    z /= np.maximum(np.linalg.norm(z, axis=1, keepdims=True), 1e-300)
    r = z[len(genes):] @ z[:len(genes)].T  # null x real
    obs = 1 / np.arctanh(np.clip(r, -0.999999, 0.999999)).var(axis=0) + 3
    s = m.n_obs / pred
    for lo, hi in BINS:
        sel = (s >= np.percentile(s, lo)) & (s <= np.percentile(s, hi))
        rows.append({"cell_type": ct, "cells": full.n_obs, "metacells": m.n_obs,
                     "cells_per_k": full.n_obs // 25, "s_percentiles": f"{lo}-{hi}",
                     "observed": np.median(obs[sel]), "predicted": np.median(pred[sel])})
        print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv(HERE / "results/effective_samples.tsv", sep="\t", index=False, float_format="%.4g")
print(df.round(1).to_string(index=False))
