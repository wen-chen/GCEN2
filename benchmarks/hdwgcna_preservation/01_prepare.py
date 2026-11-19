"""The reference dataset of the NK cells (../hdwgcna_consensus/01_prepare.py
output) as its own AnnData, for gcen2 celltype-modules.

usage: python 01_prepare.py NK_H5AD OUT_H5AD DATASET
"""
import sys

import anndata as ad

src, out, dataset = sys.argv[1], sys.argv[2], sys.argv[3]
a = ad.read_h5ad(src)
r = a[(a.obs["dataset"].astype(str) == dataset).to_numpy()].copy()
print(f"{r.n_obs} cells of {dataset}")
r.write_h5ad(out, compression="gzip")
