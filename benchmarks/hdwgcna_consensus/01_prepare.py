"""One cell type of the Human Lung Cell Atlas core over several datasets, as
inputs for GCEN2 and hdWGCNA's consensus workflow.

From the CELLxGENE file (see ../hdwgcna_hlca/01_prepare.py): raw counts, gene
symbols, the cells of CELL_TYPE in the given datasets (obs: cell_type,
dataset, donor_id; and the traits age (age_or_mean_of_age_range), BMI and
sex) and the scANVI embedding (the metacell neighbours).

  OUT.h5ad          counts (X, float32 CSR), obs cell_type / dataset /
                    donor_id / age / BMI / sex, obsm X_scanvi_emb
  OUT_export/       for R: matrix.h5 (10x HDF5, genes x cells), genes.tsv,
                    meta.tsv (barcode, cell_type, dataset, donor_id, age,
                    BMI, sex),
                    rep.tsv (cells x 30 scANVI dimensions)

usage: python 01_prepare.py HLCA_H5AD OUT "CELL_TYPE" DATASET [DATASET ...]
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

src, out, cell_type, datasets = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4:]
a = ad.read_h5ad(src, backed="r")
obs = a.obs
keep = np.flatnonzero((obs["cell_type"].astype(str) == cell_type).to_numpy()
                      & obs["dataset"].astype(str).isin(datasets).to_numpy())
print(f"{keep.size} cells of {cell_type} in {len(datasets)} datasets", flush=True)
raw = a.raw if a.raw is not None else a
X = sp.csr_matrix(raw.X[keep], dtype=np.float32)
var = raw.var
genes = pd.Index(var["feature_name"].astype(str).to_numpy() if "feature_name" in var
                 else np.asarray(var.index, dtype=str))
if not genes.is_unique:
    genes = pd.Index(ad.utils.make_index_unique(genes))
meta = pd.DataFrame({c: obs[c].astype(str).to_numpy()[keep]
                     for c in ("cell_type", "dataset", "donor_id")},
                    index=np.asarray(obs.index, dtype=str)[keep])
meta["age"] = obs["age_or_mean_of_age_range"].to_numpy(dtype=np.float64)[keep]
meta["BMI"] = obs["BMI"].to_numpy(dtype=np.float64)[keep]
meta["sex"] = pd.Categorical(obs["sex"].astype(str).to_numpy()[keep])
rep = np.asarray(a.obsm["X_scanvi_emb"][keep], dtype=np.float64)
print(meta.groupby("dataset").size().to_string(), flush=True)
print(f"{meta['donor_id'].nunique()} donors", flush=True)

b = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes))
b.obsm["X_scanvi_emb"] = rep
b.write_h5ad(out + ".h5ad", compression="gzip")

e = Path(out + "_export")
e.mkdir(parents=True, exist_ok=True)
csc = X.T.tocsc()  # genes x cells, one column per cell
with h5py.File(e / "matrix.h5", "w") as f:
    m = f.create_group("matrix")
    m["data"] = csc.data
    m["indices"] = csc.indices.astype(np.int64)
    m["indptr"] = csc.indptr.astype(np.int64)
    m["shape"] = np.array(csc.shape, dtype=np.int64)
    m["barcodes"] = np.array(meta.index, dtype="S")
    ft = m.create_group("features")
    ft["name"] = np.array(genes, dtype="S")
    ft["id"] = np.array(genes, dtype="S")
    ft["feature_type"] = np.array(["Gene Expression"] * len(genes), dtype="S")
    ft["genome"] = np.array(["GRCh38"] * len(genes), dtype="S")
(e / "genes.tsv").write_text("\n".join(genes) + "\n")
meta.rename_axis("barcode").to_csv(e / "meta.tsv", sep="\t")
pd.DataFrame(rep, index=meta.index, columns=[f"rep_{i + 1}" for i in range(rep.shape[1])]) \
    .rename_axis("barcode").to_csv(e / "rep.tsv", sep="\t", float_format="%.10g")
