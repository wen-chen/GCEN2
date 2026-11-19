"""How single cells become correlation samples: metacells against pseudo-replicates,
single cells and MAGIC imputation.

Run from the PBMC data directory (after 01_preprocess.py). For each of the six
cell types, on its gene list (work/genes_<type>.txt), with GCEN2's own code
paths and the cutoffs of 03_gcen2.py (Pearson, r >= 0.3, p <= 0.001, positive
edges only, 8 threads) but without FDR, so that the null edges measure the
p-values themselves:

  metacell     `gcen2 celltype` default: kNN metacells (k = 25, max_shared = 10,
               at most 1,000, the shared PCA), p-values on every pair's
               effective samples from the metacells' overlap
  metacell_nk  the same with `--effective-samples cells`: p-values on cells / k
  pseudo10     `--strategy pseudo_replicates --replicates 10 -l log2 --layer norm`:
  pseudo50     the cells split at random into 10 (50) groups, their
               library-size-normalised counts (counts / total x 1e4) summed,
               log2(x + 1). (On raw counts the sums vary with the groups'
               library sizes, and the top pairs are ribosomal protein genes.)
  cells        every cell a sample (log-normalised: counts / total x 1e4, log1p)
  magic        `build_dropout_aware(zero_handling="impute")`: MAGIC 3.0.0
               (defaults) on the same log-normalised cells, then Pearson

Measured:
  * false edges: 200 null genes are added. Each takes a gene of the list,
    permutes its fraction of every cell's counts (count / library size)
    across the type's cells and draws Poisson counts with mean fraction x
    library size: no co-expression with anything, and the expected count
    proportional to the library size, as for a gene whose expression does not
    vary. (Permuting the counts within bins of library size instead keeps
    E[count | library size], which is not proportional for genes such as the
    ribosomal ones, whose share differs between cell states with different
    library sizes; after normalisation such null genes follow the states.) They are
    added before the network is built; the fraction of the pairs with a null
    gene that become edges. Valid p-values give at most 0.0005 (p <= 0.001,
    positive side only).
  * ranking: STRING enrichment of the top 1,000 / 10,000 gene pairs by
    correlation (real genes; work/string_gold.tsv, as 04_compare.py).
  * reproducibility: correlation of the gene-gene correlations, and overlap
    of the top 10,000 pairs, of the two disjoint halves of the cells
    (work/splits, as 04_compare.py).

MAGIC is not installed with GCEN2 (optional extra); point MAGIC_LIB at a
directory with magic-impute 3.0.0 and its dependencies (`pip install --target`)
when it is not in the environment.

usage: 09_strategies.py [CELL_TYPE ...]   (default: all six)
Writes results/strategies.tsv.
"""

import importlib.util
import os
import sys
import time
import warnings
from pathlib import Path

if os.environ.get("MAGIC_LIB"):
    sys.path.append(os.environ["MAGIC_LIB"])

import anndata as ad
import magic
import numpy as np
import pandas as pd
import scipy.sparse as sp

import gcen2
from gcen2.single_cell import build_by_celltype, build_dropout_aware, metacells
from gcen2.single_cell.core import _partition_cells

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location("compare", HERE / "04_compare.py")
compare = importlib.util.module_from_spec(spec)
spec.loader.exec_module(compare)

WORK = Path("work")
TYPES = sys.argv[1:] or ["CD14 Mono", "CD4 T", "B", "CD8 T", "NK", "FCGR3A Mono"]
THREADS = 8
CUT = dict(method="pearson", cor=0.3, pval=0.001, signed=True, fdr=False, threads=THREADS)
META = dict(k=25, max_shared=10, max_metacells=1000, min_cells=100, use_rep="X_pca",
            n_pcs=30, seed=42)
N_NULL = 200
STRATEGIES = ["metacell", "metacell_nk", "pseudo10", "pseudo50", "cells", "magic"]


class CaptureMAGIC(magic.MAGIC):
    """magic.MAGIC that keeps its output, so the ranking uses the very matrix
    build_dropout_aware correlated (MAGIC's PCA is randomised)."""
    last = None

    def fit_transform(self, X, *args, **kwargs):
        out = super().fit_transform(X, *args, **kwargs)
        CaptureMAGIC.last = out
        return out


magic.MAGIC = CaptureMAGIC  # gcen2 calls magic.MAGIC()


def dense(x):
    return np.asarray(x.toarray() if sp.issparse(x) else x, dtype=np.float64)


def null_counts(x, library, rng):
    """Null genes from the columns of x (cells x genes counts): each column's
    fractions of the library permuted across the cells, then Poisson counts
    of fraction x library size (see the module docstring)."""
    frac = x / np.maximum(library, 1)[:, None]
    out = np.empty_like(x)
    for j in range(x.shape[1]):
        out[:, j] = rng.poisson(frac[rng.permutation(len(library)), j] * library)
    return out


def lognorm(counts):
    """cells x genes counts -> log1p(counts / total x 1e4)."""
    tot = counts.sum(axis=1, keepdims=True)
    tot[tot == 0] = 1
    return np.log1p(counts / tot * 1e4)


def samples(strategy, sub, genes):
    """samples x genes matrix of the strategy (MAGIC: None, see magic_matrix)."""
    if strategy.startswith("metacell"):
        return dense(metacells(sub, groupby="cell_type", threads=THREADS, **META)[:, genes].X)
    counts = dense(sub[:, genes].X)
    if strategy.startswith("pseudo"):
        norm = dense(sub[:, genes].layers["norm"])
        parts = _partition_cells(len(counts), int(strategy[6:]), np.random.default_rng(42))
        return np.log2(np.vstack([norm[p].sum(axis=0) for p in parts]) + 1)
    if strategy == "cells":
        return lognorm(counts)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        return np.asarray(magic.MAGIC(verbose=0).fit_transform(
            pd.DataFrame(lognorm(counts), columns=genes)), dtype=np.float64)


def upper_cor(x):
    """Upper triangle of the gene-gene Pearson correlations of samples x genes
    (genes constant in x: 0)."""
    x = x - x.mean(axis=0)
    norm = np.sqrt((x * x).sum(axis=0))
    z = np.divide(x, norm, out=np.zeros_like(x), where=norm > 0)
    c = z.T @ z
    return c[np.triu_indices(len(c), 1)].astype(np.float32)


def network(strategy, full, genes_all):
    """GCEN2's network on real + null genes, its p-value sample count and time."""
    t0 = time.perf_counter()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        if strategy.startswith("metacell"):
            eff = "cells" if strategy == "metacell_nk" else "overlap"
            net = build_by_celltype(full, "cell_type", effective=eff, **META, **CUT)[CT]
        elif strategy.startswith("pseudo"):
            net = build_by_celltype(full, "cell_type", strategy="pseudo_replicates",
                                    pseudo_replicates=int(strategy[6:]), min_cells=100,
                                    seed=42, log="log2", layer="norm", **CUT)[CT]
        else:
            x = pd.DataFrame(lognorm(dense(full.X)).T, index=genes_all)
            if strategy == "cells":
                net = gcen2.network.build(x, **CUT)
            else:
                net = build_dropout_aware(x, zero_handling="impute", **CUT)
    t = time.perf_counter() - t0
    n_p = net.metadata.get("p_value_samples", net.metadata.get("n_samples", full.n_obs))
    return net.edges, n_p, t


adata = ad.read_h5ad(WORK / "pbmc.h5ad")
gold = compare.string_gold()
rows = []
for CT in TYPES:
    tag = CT.replace(" ", "_")
    genes = (WORK / f"genes_{tag}.txt").read_text().split()
    cells = adata.obs_names[adata.obs["cell_type"] == CT]
    sub = adata[cells, genes].copy()
    sub.X = sp.csr_matrix(sub.X)

    # null genes: real genes' counts permuted across cells of similar library size
    rng = np.random.default_rng(0)
    src = rng.choice(len(genes), N_NULL, replace=False)
    null = null_counts(dense(sub.X[:, src]), np.asarray(sub.X.sum(axis=1)).ravel(), rng)
    null_names = [f"null_{genes[j]}" for j in src]
    full = ad.AnnData(sp.hstack([sub.X, sp.csr_matrix(null)]).tocsr(), obs=sub.obs.copy(),
                      obsm={"X_pca": sub.obsm["X_pca"]})
    full.var_names = genes + null_names
    tot = np.asarray(full.X.sum(axis=1)).ravel()  # real + null genes
    full.layers["norm"] = sp.diags(1e4 / np.maximum(tot, 1)) @ full.X
    sub.layers["norm"] = full[:, genes].layers["norm"]
    genes_all = list(full.var_names)
    is_null = set(null_names)
    n_real = len(genes)
    null_pairs = N_NULL * n_real + N_NULL * (N_NULL - 1) // 2

    iu = np.triu_indices(n_real, 1)
    is_gold = compare.gold_mask(np.asarray(genes), iu, gold)
    halves = {h: (WORK / "splits" / f"{tag}_{h}.txt").read_text().split() for h in "AB"}

    for s in STRATEGIES:
        CaptureMAGIC.last = None
        edges, n_p, t_net = network(s, full, genes_all)
        e_null = edges["gene_a"].isin(is_null) | edges["gene_b"].isin(is_null)
        if s == "magic":  # the matrix build_dropout_aware correlated, real genes
            x = np.asarray(CaptureMAGIC.last, dtype=np.float64)[:, :n_real]
        else:
            x = samples(s, sub, genes)
        v = upper_cor(x)
        prec = compare.precision_at(v, is_gold)
        del v
        hv = [upper_cor(samples(s, sub[halves[h]], genes)) for h in "AB"]
        top = [set(np.argsort(-h)[:10000].tolist()) for h in hv]
        row = {"cell_type": CT, "strategy": s, "cells": len(cells), "genes": n_real,
               "samples": len(x), "p_value_samples": n_p,
               "real_edges": int((~e_null).sum()), "null_edges": int(e_null.sum()),
               "null_edge_rate": e_null.sum() / null_pairs,
               "P@1000": prec["P@1000"], "enr@1000": prec["enr@1000"],
               "P@10000": prec["P@10000"], "enr@10000": prec["enr@10000"],
               "baseline": prec["baseline"],
               "split_cor_of_cor": float(np.corrcoef(hv[0], hv[1])[0, 1]),
               "split_top10k": len(top[0] & top[1]) / 10000, "t_network": t_net}
        rows.append(row)
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in row.items()},
              flush=True)
        del hv, top

df = pd.DataFrame(rows)
out = HERE / "results/strategies.tsv"
df.to_csv(out, sep="\t", index=False, float_format="%.6g")
print(df.drop(columns=["baseline", "genes"]).round(4).to_string(index=False))
