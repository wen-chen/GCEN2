"""A Visium HD sample at one bin size as AnnData, plus an export for
02_hdwgcna_knn.R.

Visium HD bins lie on a square grid, so `gcen2 spatial-modules` takes knn
metaspots for them (``--metaspots auto``: obs array_row / array_col are not a
hexagonal Visium grid), and hdWGCNA's MetaspotsByGroups does not apply.

usage: python 01_prepare_visium_hd.py BIN_DIR OUT_PREFIX [N_REGIONS]
BIN_DIR: one bin size of Space Ranger's binned outputs, e.g.
binned_outputs/square_016um, with filtered_feature_bc_matrix.h5 and
spatial/tissue_positions.parquet (or a tissue_positions.csv converted from it
by any Python with pyarrow). region = k-means (N_REGIONS, default 12, seed 0)
of 30 principal components of log-normalized highly variable genes, as
01_prepare.py does for the 10x Visium sample.
Writes OUT_PREFIX.h5ad (counts; obs array_row, array_col, region;
obsm["spatial"] full-resolution pixels x, y) and OUT_PREFIX_export/
(matrix.h5: a copy of filtered_feature_bc_matrix.h5; genes.tsv: the gene
names used here; meta.tsv: barcode, x, y, region).
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.cluster import KMeans

raw, prefix = Path(sys.argv[1]), sys.argv[2]
k = int(sys.argv[3]) if len(sys.argv) > 3 else 12
adata = sc.read_10x_h5(raw / "filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
csv = raw / "spatial" / "tissue_positions.csv"
pos = pd.read_csv(csv) if csv.exists() else pd.read_parquet(
    raw / "spatial" / "tissue_positions.parquet")
pos = pos.set_index("barcode").loc[adata.obs_names]
assert (pos["in_tissue"] == 1).all()
adata.obs["array_row"] = pos["array_row"].to_numpy()
adata.obs["array_col"] = pos["array_col"].to_numpy()
adata.obsm["spatial"] = pos[["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy(dtype=float)

x = adata.copy()
sc.pp.normalize_total(x, target_sum=1e4)
sc.pp.log1p(x)
sc.pp.highly_variable_genes(x, n_top_genes=2000)
x = x[:, x.var["highly_variable"]].copy()
sc.pp.scale(x, max_value=10)
sc.tl.pca(x, n_comps=30, random_state=0)
lab = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(x.obsm["X_pca"])
adata.obs["region"] = pd.Categorical([f"R{v:02d}" for v in lab])
adata.X = adata.X.astype(np.float32)
print(f"{adata.n_obs} bins, {adata.n_vars} genes, {adata.X.nnz:,} non-zero, median "
      f"{np.median(np.asarray(adata.X.sum(axis=1)).ravel()):.0f} counts per bin")
print(adata.obs["region"].value_counts().sort_index().to_string())
adata.write_h5ad(prefix + ".h5ad", compression="gzip")

out = Path(prefix + "_export")
out.mkdir(parents=True, exist_ok=True)
shutil.copy(raw / "filtered_feature_bc_matrix.h5", out / "matrix.h5")
# gene names as made unique here (Seurat's make.unique would differ), in the
# order of the h5 file
(out / "genes.tsv").write_text("\n".join(adata.var_names) + "\n")
pd.DataFrame({"x": adata.obsm["spatial"][:, 0], "y": adata.obsm["spatial"][:, 1],
              "region": adata.obs["region"].astype(str).to_numpy()},
             index=adata.obs_names.rename("barcode")).to_csv(out / "meta.tsv", sep="\t")
