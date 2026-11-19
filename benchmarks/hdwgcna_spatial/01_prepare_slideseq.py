"""Slide-seqV2 mouse hippocampus (Stickels et al. 2021; squidpy's
``slideseqv2`` dataset, puck 200115_08) as AnnData with raw counts, plus an
export for 02_hdwgcna_knn.R.

squidpy ships log-normalized values of 4,000 highly variable genes
(``log1p(counts / total * s)``). The counts are recovered exactly per bead:
the smallest non-zero ``expm1`` value of a bead is one count, so dividing
by it gives integers (the largest deviation is 5e-4, float32 rounding). The
counts cover these 4,000 genes only.

usage: python 01_prepare_slideseq.py SLIDESEQV2_H5AD OUT_PREFIX
Writes OUT_PREFIX.h5ad (X: counts; obs: cluster, x, y; obsm["spatial"]) and
OUT_PREFIX_export/ (counts.mtx, genes.tsv, barcodes.tsv, meta.tsv).
"""
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp

src, prefix = Path(sys.argv[1]), sys.argv[2]
a = ad.read_h5ad(src)
E = sp.csr_matrix(a.X, dtype=np.float64)
E.data = np.expm1(E.data)
lengths = np.diff(E.indptr)
assert (lengths > 0).all(), "a bead without counts"
unit = np.minimum.reduceat(E.data, E.indptr[:-1])
C = sp.diags(1.0 / unit) @ E
err = np.abs(C.data - np.round(C.data)).max()
assert err < 0.01, f"not integer counts (max deviation {err})"
C.data = np.round(C.data)
C = sp.csr_matrix(C, dtype=np.float32)
C.eliminate_zeros()

obs = pd.DataFrame({"cluster": a.obs["cluster"].astype(str).to_numpy(),
                    "x": a.obsm["spatial"][:, 0], "y": a.obsm["spatial"][:, 1]},
                   index=a.obs_names.astype(str))
out = ad.AnnData(X=C, obs=obs, var=pd.DataFrame(index=a.var_names.astype(str)))
out.obsm["spatial"] = np.asarray(a.obsm["spatial"], dtype=np.float64)
out.write_h5ad(f"{prefix}.h5ad", compression="gzip")

exp = Path(f"{prefix}_export")
exp.mkdir(parents=True, exist_ok=True)
scipy.io.mmwrite(exp / "counts.mtx", sp.csr_matrix(C.T, dtype=np.int64))
(exp / "genes.tsv").write_text("\n".join(out.var_names) + "\n")
(exp / "barcodes.tsv").write_text("\n".join(out.obs_names) + "\n")
obs.rename_axis("barcode").to_csv(exp / "meta.tsv", sep="\t")
print(f"{out.n_obs} beads, {out.n_vars} genes, max deviation from integers {err:.1e}, "
      f"median {np.median(np.asarray(C.sum(axis=1)).ravel()):.0f} counts per bead; "
      f"clusters: {obs['cluster'].value_counts().to_dict()}")
