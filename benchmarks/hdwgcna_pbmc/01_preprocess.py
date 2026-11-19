"""PBMC 10k v3 (10x Genomics): QC, normalisation, PCA, clustering and
marker-based cell-type labels. Writes the shared inputs for gcen2 and
hdWGCNA so both use the same cells, labels, PCA and gene lists.

Run from the PBMC data directory (data/ holds the 10x file); outputs in work/:
  pbmc.h5ad                counts in X (int), obs cell_type, obsm X_pca,
                           var highly_variable
  counts.mtx, genes.tsv, barcodes.tsv, meta.tsv, pca.tsv   (for Seurat)
  genes_<type>.txt         genes expressed in >= 5% of the type's cells
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.io
import scipy.sparse as sp
from sklearn.cluster import KMeans

DATA = Path("data/pbmc_10k_v3_filtered_feature_bc_matrix.h5")  # relative to the data directory
WORK = Path("work")
WORK.mkdir(parents=True, exist_ok=True)
TYPES = ["CD14 Mono", "CD4 T", "CD8 T", "NK", "B", "FCGR3A Mono", "DC", "Platelet"]
MARKERS = {
    "CD4 T": ["CD3E", "IL7R", "CCR7", "LTB", "CD4"],
    "CD8 T": ["CD3E", "CD8A", "CD8B", "GZMK"],
    "NK": ["NKG7", "GNLY", "KLRD1", "PRF1", "FCGR3A"],
    "B": ["MS4A1", "CD79A", "CD79B", "BANK1"],
    "CD14 Mono": ["CD14", "LYZ", "S100A8", "S100A9", "VCAN"],
    "FCGR3A Mono": ["FCGR3A", "MS4A7", "LST1", "CDKN1C"],
    "DC": ["FCER1A", "CST3", "CLEC10A", "CD1C"],
    "Platelet": ["PPBP", "PF4", "GP9"],
}

adata = sc.read_10x_h5(DATA)
adata.var_names_make_unique()
adata.var["mt"] = adata.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], inplace=True, percent_top=None)
keep = ((adata.obs["n_genes_by_counts"] >= 200) & (adata.obs["n_genes_by_counts"] <= 6000)
        & (adata.obs["pct_counts_mt"] < 15))
adata = adata[keep.to_numpy()].copy()
sc.pp.filter_genes(adata, min_cells=3)
print(f"after QC: {adata.n_obs} cells x {adata.n_vars} genes")

counts = adata.X.copy()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000)
lognorm = adata.X.copy()
sc.pp.scale(adata, max_value=10)
sc.tl.pca(adata, n_comps=30, mask_var="highly_variable", random_state=0)
adata.X = lognorm

km = KMeans(n_clusters=20, n_init=10, random_state=0).fit(adata.obsm["X_pca"][:, :20])
adata.obs["cluster"] = km.labels_.astype(str)

# cluster -> type by the mean z-scored marker expression
z = {}
for t, genes in MARKERS.items():
    g = [x for x in genes if x in adata.var_names]
    sc.tl.score_genes(adata, g, score_name=f"score_{t}", random_state=0)
    z[t] = adata.obs.groupby("cluster", observed=True)[f"score_{t}"].mean()
scores = pd.DataFrame(z)
scores = (scores - scores.mean()) / scores.std()
cluster_type = scores.idxmax(axis=1)
print(pd.concat([scores.round(2), cluster_type.rename("type")], axis=1).to_string())
adata.obs["cell_type"] = adata.obs["cluster"].map(cluster_type).astype(str)
print(adata.obs["cell_type"].value_counts().to_string())

# CD4/CD8 marker sanity check
for t in ("CD4 T", "CD8 T", "NK", "B", "CD14 Mono"):
    sub = adata[adata.obs["cell_type"] == t]
    means = {g: float(np.asarray(sub[:, g].X.mean())) for g in
             ("CD3E", "CD8A", "IL7R", "NKG7", "MS4A1", "CD14", "LYZ") if g in adata.var_names}
    print(t, {k: round(v, 2) for k, v in means.items()})

out = sc.AnnData(X=counts.astype(np.float32), obs=adata.obs[["cell_type", "cluster",
                                                                   "total_counts"]].copy(),
                 var=adata.var[["highly_variable"]].copy())
out.obsm["X_pca"] = adata.obsm["X_pca"]
out.write_h5ad(WORK / "pbmc.h5ad")

# for Seurat: genes x cells counts, identical labels and PCA
scipy.io.mmwrite(WORK / "counts.mtx", sp.csc_matrix(counts.T).astype(np.int32))
pd.Series(out.var_names).to_csv(WORK / "genes.tsv", index=False, header=False)
pd.Series(out.obs_names).to_csv(WORK / "barcodes.tsv", index=False, header=False)
out.obs[["cell_type"]].to_csv(WORK / "meta.tsv", sep="\t")
pd.DataFrame(out.obsm["X_pca"], index=out.obs_names,
             columns=[f"PC_{i + 1}" for i in range(30)]).to_csv(WORK / "pca.tsv", sep="\t")

# per-type gene lists (hdWGCNA SetupForWGCNA gene_select="fraction", 0.05)
X = sp.csr_matrix(counts)
for t in sorted(set(out.obs["cell_type"])):
    m = (out.obs["cell_type"] == t).to_numpy()
    frac = np.asarray((X[m] > 0).mean(axis=0)).ravel()
    genes = out.var_names[frac >= 0.05]
    (WORK / f"genes_{t.replace(' ', '_')}.txt").write_text("\n".join(genes) + "\n")
    print(f"{t}: {m.sum()} cells, {len(genes)} genes in >= 5%")
sys.exit(0)
