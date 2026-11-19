"""Rebuild data/example from its sources (see data/example/README.md).

    python tools/make_example_data.py --gcen-v1 <GCEN v1 0.6.3 sample dir> \
        --pbmc <pbmc.h5ad> --pbmc-umap <umap.tsv> --visium <mouse_anterior.h5ad>

pbmc.h5ad and mouse_anterior.h5ad are the prepared data of
benchmarks/hdwgcna_pbmc/01_preprocess.py and benchmarks/hdwgcna_spatial/
01_prepare.py; umap.tsv holds x, y by cell (UMAP of their X_pca, 30 PCs,
umap-learn random_state 0)."""

from __future__ import annotations

import argparse
import gzip
import shutil
from pathlib import Path

import numpy as np

OUT = Path(__file__).resolve().parents[1] / "data" / "example"
BULK = ["gene_expr.tsv", "gene_kegg.assoc", "K2ko.tsv", "rwr_seed_genes.list",
        "enrichment_gene.list", "background_gene.list", "gene_list.txt"]
GZIP = ["go-basic.obo", "gene_go.assoc"]


def portable(a):
    """Strings of obs / var (index and columns) as object arrays, so the file
    is written with the string-array encoding that older anndata (0.8-0.10,
    e.g. on Python 3.9) reads; pandas 3's str dtype would be written as
    nullable-string-array, which they cannot read."""
    import pandas as pd

    for df in (a.obs, a.var):
        df.index = pd.Index(np.asarray(df.index, dtype=object), dtype=object,
                            name=df.index.name)
        for c in df.columns:
            if pd.api.types.is_string_dtype(df[c].dtype) and \
                    not isinstance(df[c].dtype, pd.CategoricalDtype):
                df[c] = df[c].astype(object)
    return a


def bulk(src: Path) -> None:
    d = OUT / "bulk"
    d.mkdir(parents=True, exist_ok=True)
    for f in BULK:
        shutil.copyfile(src / f, d / f)
    for f in GZIP:
        # no file name or time in the header (as gzip -n): the same bytes every time
        with open(src / f, "rb") as fi, open(d / f"{f}.gz", "wb") as raw, \
                gzip.GzipFile(filename="", mode="wb", compresslevel=9, fileobj=raw,
                              mtime=0) as fo:
            shutil.copyfileobj(fi, fo)


def single_cell(h5ad: Path, umap: Path) -> None:
    import anndata as ad
    import pandas as pd
    import scipy.sparse as sp

    a = ad.read_h5ad(h5ad)
    rng = np.random.default_rng(0)
    idx = []
    for t, n in (("CD14 Mono", 800), ("B", 600), ("NK", 400)):
        cells = np.flatnonzero(a.obs["cell_type"].to_numpy() == t)
        idx += sorted(rng.choice(cells, n, replace=False))
    b = a[np.sort(idx)].copy()
    detected = np.asarray((b.X > 0).mean(axis=0)).ravel()
    c = b[:, detected >= 0.1].copy()
    c.X = sp.csr_matrix(c.X.astype(np.int32))
    xy = pd.read_csv(umap, sep="\t", index_col=0)
    c.obsm["X_umap"] = xy.loc[c.obs_names, ["x", "y"]].to_numpy(np.float32)
    c.obsm["X_pca"] = c.obsm["X_pca"].astype(np.float32)
    for k in ("cell_type", "cluster"):
        c.obs[k] = c.obs[k].astype(str).astype("category")
    (OUT / "single_cell").mkdir(parents=True, exist_ok=True)
    portable(c).write_h5ad(OUT / "single_cell" / "pbmc_subset.h5ad", compression="gzip")


def spatial(h5ad: Path) -> None:
    import anndata as ad
    import scipy.sparse as sp

    s = ad.read_h5ad(h5ad)
    x = s.X if sp.issparse(s.X) else sp.csr_matrix(s.X)
    detected = np.asarray((x > 0).mean(axis=0)).ravel()
    s = s[:, np.sort(np.argsort(-detected, kind="stable")[:2000])].copy()
    s.X = sp.csr_matrix(s.X.astype(np.int32))
    for k in list(s.uns):
        del s.uns[k]
    (OUT / "spatial").mkdir(parents=True, exist_ok=True)
    portable(s).write_h5ad(OUT / "spatial" / "mouse_brain_anterior_subset.h5ad",
                           compression="gzip")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gcen-v1", type=Path, help="GCEN v1 0.6.3 sample data directory")
    ap.add_argument("--pbmc", type=Path)
    ap.add_argument("--pbmc-umap", type=Path)
    ap.add_argument("--visium", type=Path)
    a = ap.parse_args()
    if a.gcen_v1:
        bulk(a.gcen_v1)
    if a.pbmc:
        single_cell(a.pbmc, a.pbmc_umap)
    if a.visium:
        spatial(a.visium)


if __name__ == "__main__":
    main()
