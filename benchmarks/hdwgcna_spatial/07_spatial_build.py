"""`gcen2 spatial-build` on DLPFC 151673: what the three modes add to a
co-expression network of the spots.

Run from the spatial data directory (after 01_prepare.py). The spots of
work/dlpfc.h5ad (3,611, counts), genes with counts in >= 5% of the spots,
log-normalised, Spearman, |r| >= CUT and p <= 0.05 (the command's defaults
but the cutoff), 6 nearest neighbours, 8 threads:

  plain      gcen2.network.build on the spots (no spatial information)
  filtered   `-m filtered`: the plain edges between genes with significant
             positive spatial autocorrelation (Moran's I test, FDR <= 0.05)
  i>0        `-m filtered --moran-fdr 1`: the plain edges between genes with
             Moran's I > 0 (the rule before the test was added)
  weighted   `-m weighted`: every spot replaced by the Gaussian-weighted mean
             of itself and its neighbours, then the network
  regions    `-m regions --region-key region`: one network per cortical layer
             (L1-L6, WM; manual annotation); edges pooled over the layers

Measured: edges and genes; precision of the edges for STRING interactions
(experimental or database >= 700, the PBMC benchmark's gold standard, path
in STRING_GOLD) against all gene pairs; false edges: 200 null genes, each a
gene's fractions of the spots' counts permuted across the spots, then Poisson
counts of fraction x library size (no spatial pattern of its own, no
co-expression; as the PBMC benchmark's null genes), are added first, and the fraction of the pairs with a null
gene that become edges is counted; Moran's I of the real and the null genes.

usage: 07_spatial_build.py [CUT ...]   (default 0.6 0.3)
Writes results/spatial_build.tsv and results/spatial_build_moran.tsv.
"""

import os
import sys
import time
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

import gcen2
from gcen2 import spatial
from gcen2.spatial.core import _get_expression

HERE = Path(__file__).resolve().parent
CUTS = [float(c) for c in sys.argv[1:]] or [0.6, 0.3]
N_NULL = 200
KW = dict(method="spearman", threads=8)

adata = ad.read_h5ad("work/dlpfc.h5ad")
adata.X = sp.csr_matrix(adata.X)
rng = np.random.default_rng(0)
detected = np.asarray((adata.X > 0).sum(axis=0)).ravel()
cand = np.flatnonzero(detected >= round(0.05 * adata.n_obs))
src = rng.choice(cand, N_NULL, replace=False)
xs = adata.X[:, src].toarray()
library = np.asarray(adata.X.sum(axis=1)).ravel()
frac = xs / np.maximum(library, 1)[:, None]
null = np.empty_like(xs)
for j in range(N_NULL):  # fractions permuted across spots, Poisson(fraction x library)
    null[:, j] = rng.poisson(frac[rng.permutation(len(library)), j] * library)
null_names = [f"null_{adata.var_names[j]}" for j in src]
full = ad.AnnData(sp.hstack([adata.X, sp.csr_matrix(null)]).tocsr(), obs=adata.obs.copy(),
                  obsm={"spatial": adata.obsm["spatial"]})
full.var_names = list(adata.var_names) + null_names
is_null = set(null_names)

X, genes = _get_expression(full, None, 0.05, None)
real = [g for g in genes if g not in is_null]
n_real = len(real)
n_null = len(genes) - n_real
null_pairs = n_null * n_real + n_null * (n_null - 1) // 2
gold = pd.read_csv(os.environ.get("STRING_GOLD", "../gcen-pbmc/work/string_gold.tsv"),
                   sep="\t")
real_set = set(real)
gold = gold[gold["a"].isin(real_set) & gold["b"].isin(real_set)]
gold_pairs = set(zip(gold["a"], gold["b"])) | set(zip(gold["b"], gold["a"]))
baseline = len(gold) / (n_real * (n_real - 1) / 2)


def summary(mode, cut, edges, t, spots=len(full)):
    e_null = edges["gene_a"].isin(is_null) | edges["gene_b"].isin(is_null)
    re = edges[~e_null]
    hits = sum((a, b) in gold_pairs for a, b in zip(re["gene_a"], re["gene_b"]))
    prec = hits / len(re) if len(re) else np.nan
    return {"mode": mode, "cor": cut, "spots": spots, "edges": len(re),
            "genes": len(set(re["gene_a"]) | set(re["gene_b"])),
            "positive": float((re["correlation"] > 0).mean()) if len(re) else np.nan,
            "string_precision": prec, "string_enrichment": prec / baseline,
            "null_edges": int(e_null.sum()), "null_edge_rate": e_null.sum() / null_pairs,
            "t": t}


rows = []
for cut in CUTS:
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        t0 = time.perf_counter()
        plain = gcen2.network.build(pd.DataFrame(X.T, index=genes), cor=cut, pval=0.05, **KW)
        rows.append(summary("plain", cut, plain.edges, time.perf_counter() - t0))
        for mode, fn, extra in (
                ("filtered", spatial.build_spatial_filtered, {}),
                ("i>0", spatial.build_spatial_filtered, {"moran_fdr": 1.0}),
                ("weighted", spatial.build_spatial_weighted, {})):
            t0 = time.perf_counter()
            net = fn(full, n_neighbors=6, cor_cutoff=cut, p_cutoff=0.05, **extra, **KW)
            rows.append(summary(mode, cut, net.edges, time.perf_counter() - t0))
        t0 = time.perf_counter()
        nets = spatial.build_spatial_regions(full, region_key="region", cor=cut, pval=0.05, **KW)
        t = time.perf_counter() - t0
    pooled = pd.concat([n.edges for n in nets.values()]).drop_duplicates(["gene_a", "gene_b"])
    rows.append(summary("regions", cut, pooled, t))
    for name, n in sorted(nets.items()):
        rows.append(summary(f"regions:{name}", cut, n.edges, np.nan,
                            spots=int((full.obs["region"].astype(str) == name).sum())))
    for r in rows[-len(nets) - 5:]:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()},
              flush=True)

df = pd.DataFrame(rows)
df["string_baseline"] = baseline
df.to_csv(HERE / "results/spatial_build.tsv", sep="\t", index=False, float_format="%.6g")
print(df.round(4).to_string(index=False))

mi = spatial.spatial_autocorrelation(full, n_neighbors=6, smoothed=True)
mi["null"] = mi["gene"].isin(is_null)
mi.to_csv(HERE / "results/spatial_build_moran.tsv", sep="\t", index=False,
          float_format="%.6g")
for flag, g in mi.groupby("null"):
    print("null" if flag else "real", len(g), "genes; I > 0:", round((g.morans_i > 0).mean(), 3),
          "FDR <= 0.05:", round(((g.morans_i > 0) & (g.fdr <= 0.05)).mean(), 3),
          "median I", round(g.morans_i.median(), 4), "median I smoothed",
          round(g.morans_i_smoothed.median(), 4))
print(mi.head(15).to_string(index=False))
