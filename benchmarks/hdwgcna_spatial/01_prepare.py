"""Visium samples as AnnData (raw counts; obs array_row / array_col, the
Visium grid, and a ``region`` column; obsm["spatial"], full-resolution
pixels x / y), plus an export for 02_hdwgcna.R (counts.mtx, genes.tsv,
barcodes.tsv, meta.tsv with row / col / imagerow / imagecol / region, and
spatial/ with the image files Seurat's Read10X_Image reads).

usage: python 01_prepare.py dlpfc RAW_DIR OUT_PREFIX
       python 01_prepare.py tenx  RAW_DIR OUT_PREFIX [N_REGIONS]

dlpfc: spatialLIBD DLPFC sample 151673 (Maynard et al. 2021); region = the
       manual cortical layer annotation (spots without one are dropped).
tenx:  a Space Ranger output (filtered_feature_bc_matrix.h5, spatial/);
       region = k-means (N_REGIONS, default 12, seed 0) of the first 30
       principal components of log-normalized highly variable genes.
Writes OUT_PREFIX.h5ad and the directory OUT_PREFIX_export/.
"""
import shutil
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp

kind, raw, prefix = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
adata = sc.read_10x_h5(raw / "filtered_feature_bc_matrix.h5")
adata.var_names_make_unique()
cols = ["in_tissue", "array_row", "array_col", "pxl_row", "pxl_col"]
if kind == "dlpfc":
    pos_file = raw / "tissue_positions_list.txt"
    image_files = [raw / "scalefactors_json.json", raw / "tissue_lowres_image.png"]
else:
    pos_file = raw / "spatial" / "tissue_positions_list.csv"
    image_files = [raw / "spatial" / "scalefactors_json.json",
                   raw / "spatial" / "tissue_lowres_image.png"]
pos = pd.read_csv(pos_file, header=None, index_col=0, names=["barcode"] + cols)
pos = pos.loc[adata.obs_names]
assert (pos["in_tissue"] == 1).all()
adata.obs["array_row"] = pos["array_row"].to_numpy()
adata.obs["array_col"] = pos["array_col"].to_numpy()
adata.obsm["spatial"] = pos[["pxl_col", "pxl_row"]].to_numpy(dtype=float)
n0 = adata.n_obs

if kind == "dlpfc":
    lay = pd.read_csv(raw / "barcode_level_layer_map.tsv", sep="\t", header=None,
                      names=["barcode", "sample", "layer"], dtype=str)
    lay = lay[lay["sample"] == "151673"].set_index("barcode")["layer"]
    adata.obs["region"] = lay.reindex(adata.obs_names).to_numpy()
    adata = adata[adata.obs["region"].notna()].copy()
elif kind == "tenx":
    from sklearn.cluster import KMeans

    k = int(sys.argv[4]) if len(sys.argv) > 4 else 12
    x = adata.copy()
    sc.pp.normalize_total(x, target_sum=1e4)
    sc.pp.log1p(x)
    sc.pp.highly_variable_genes(x, n_top_genes=2000)
    x = x[:, x.var["highly_variable"]].copy()
    sc.pp.scale(x, max_value=10)
    sc.tl.pca(x, n_comps=30, random_state=0)
    lab = KMeans(n_clusters=k, n_init=10, random_state=0).fit_predict(x.obsm["X_pca"])
    adata.obs["region"] = [f"R{v:02d}" for v in lab]
else:
    sys.exit(f"unknown kind {kind}")
adata.obs["region"] = adata.obs["region"].astype("category")
adata.X = adata.X.astype(np.float32)
print(f"{n0} spots, {adata.n_obs} kept; {adata.n_vars} genes")
print(adata.obs["region"].value_counts().sort_index().to_string())
adata.write_h5ad(prefix + ".h5ad", compression="gzip")

out = Path(prefix + "_export")
(out / "spatial").mkdir(parents=True, exist_ok=True)
scipy.io.mmwrite(out / "counts.mtx", sp.csr_matrix(adata.X).T.astype(np.int64).tocoo())
(out / "genes.tsv").write_text("\n".join(adata.var_names) + "\n")
(out / "barcodes.tsv").write_text("\n".join(adata.obs_names) + "\n")
pd.DataFrame({"region": adata.obs["region"].astype(str), "row": adata.obs["array_row"],
              "col": adata.obs["array_col"],
              "imagerow": pos.loc[adata.obs_names, "pxl_row"],
              "imagecol": pos.loc[adata.obs_names, "pxl_col"]},
             index=adata.obs_names).to_csv(out / "meta.tsv", sep="\t")
pd.read_csv(pos_file, header=None).to_csv(out / "spatial" / "tissue_positions_list.csv",
                                          header=False, index=False)
for f in image_files:
    shutil.copy(f, out / "spatial" / f.name)
