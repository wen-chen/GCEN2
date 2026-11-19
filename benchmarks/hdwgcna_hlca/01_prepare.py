"""The Human Lung Cell Atlas core (Sikkema et al. 2023, CELLxGENE dataset
066943a2-fdac-4b29-b348-40cede398e4e) as inputs for both tools.

From the CELLxGENE file: raw counts (``raw.X``), gene symbols
(``var["feature_name"]``), ``obs["cell_type"]`` (50 cell types) and the
integrated scANVI embedding (``obsm["X_scanvi_emb"]``, used by both tools to
find metacell neighbours, as hdWGCNA's ``reduction``).

For every size (``all`` and random subsets of N cells, seed 0):
  work/hlca_<size>.h5ad   counts (X, float32 CSR), obs cell_type, obsm X_scanvi_emb
  work/hlca_<size>.h5     the same for R, cells sorted by cell type: a 10x
                          HDF5 matrix plus /obs/cell_type and /obsm/rep
and work/types.txt: the 10 largest cell types of the full atlas.

usage: python 01_prepare.py HLCA_H5AD WORK_DIR [N ...]
"""
import sys
from pathlib import Path

import anndata as ad
import h5py
import numpy as np
import pandas as pd
import scipy.sparse as sp

src, work = Path(sys.argv[1]), Path(sys.argv[2])
sizes = [int(x) for x in sys.argv[3:]]
work.mkdir(parents=True, exist_ok=True)

a = ad.read_h5ad(src)
X = a.raw.X if a.raw is not None else a.X
X = sp.csr_matrix(X, dtype=np.float32)
step = 1 << 24
for i in range(0, X.data.size, step):
    v = X.data[i:i + step]
    assert (v >= 0).all() and (v == np.round(v)).all(), "raw.X is not counts"
var = a.raw.var if a.raw is not None else a.var
symbols = var["feature_name"].astype(str).to_numpy() if "feature_name" in var else \
    np.asarray(var.index, dtype=str)
genes = pd.Index(symbols)
if not genes.is_unique:  # make unique as Seurat would
    genes = pd.Index(ad.utils.make_index_unique(genes))
rep_key = "X_scanvi_emb"
assert rep_key in a.obsm, f"no {rep_key}; obsm has {list(a.obsm)}"
rep = np.asarray(a.obsm[rep_key], dtype=np.float64)
types = a.obs["cell_type"].astype(str).to_numpy()
barcodes = np.asarray(a.obs_names, dtype=str)
print(f"{a.n_obs} cells, {len(genes)} genes, {X.nnz:,} non-zero counts, "
      f"{len(set(types))} cell types, {rep_key} {rep.shape[1]} dims", flush=True)
del a

top = pd.Series(types).value_counts()
# plain lines (names such as "CD4-positive, alpha-beta T cell" contain commas)
(work / "types.txt").write_text("\n".join(top.head(10).index) + "\n")
print(top.head(10).to_string(), flush=True)


def write(name: str, rows: np.ndarray):
    x = X[rows]
    obs = pd.DataFrame({"cell_type": types[rows]}, index=barcodes[rows])
    out = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=genes))
    out.obsm[rep_key] = rep[rows]
    out.write_h5ad(work / f"hlca_{name}.h5ad")
    # for R: cells sorted by cell type, so that a type is one contiguous slice;
    # 10x layout (genes x cells CSC = cells x genes CSR), indptr also as
    # float64 (hdf5r reads int64 as integer64)
    order = np.argsort(obs["cell_type"].to_numpy(), kind="stable")
    xs = x[order]
    with h5py.File(work / f"hlca_{name}.h5", "w") as h:
        m = h.create_group("matrix")
        m.create_dataset("data", data=xs.data.astype(np.int32))
        m.create_dataset("indices", data=xs.indices.astype(np.int32))
        m.create_dataset("indptr", data=xs.indptr.astype(np.int64))
        m.create_dataset("indptr_f", data=xs.indptr.astype(np.float64))
        m.create_dataset("shape", data=np.array([xs.shape[1], xs.shape[0]], dtype=np.int32))
        m.create_dataset("barcodes", data=obs.index.to_numpy()[order].astype("S"))
        f = m.create_group("features")
        f.create_dataset("id", data=genes.to_numpy().astype("S"))
        f.create_dataset("name", data=genes.to_numpy().astype("S"))
        f.create_dataset("feature_type", data=np.full(len(genes), b"Gene Expression"))
        f.create_dataset("genome", data=np.full(len(genes), b"GRCh38"))
        h.create_dataset("obs/cell_type", data=obs["cell_type"].to_numpy()[order].astype("S"))
        h.create_dataset("obsm/rep", data=rep[rows][order])
    print(f"hlca_{name}: {len(rows)} cells, {x.nnz:,} non-zero", flush=True)


write("all", np.arange(len(types)))
for n in sizes:
    rows = np.sort(np.random.default_rng(0).choice(len(types), size=n, replace=False))
    write(f"{n // 1000}k", rows)
