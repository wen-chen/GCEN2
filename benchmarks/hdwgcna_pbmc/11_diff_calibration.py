"""Calibration of the differential co-expression test between cell types
(`gcen2 celltype --diff`): the two random halves of a cell type as two "cell
types", which differ only by chance.

Run from the PBMC data directory (after 01_preprocess.py and the split files).
For each of the six cell types: 2,000 genes drawn at random from its gene
list (work/genes_<type>.txt); the cells of each half (work/splits) give
metacells and a network with every gene pair (as 03_gcen2.py: k = 25,
max_shared = 10, at most 1,000 metacells, the shared PCA, Pearson; no cutoff:
r >= 0, p <= 1, no FDR, both signs); `write_celltype_networks(diff=True)`
tests every pair, Fisher's z, on the sample sizes of the networks' own
p-values:

  overlap   `gcen2 celltype` default: every pair on the smaller of its two
            genes' effective samples (metacell_effective_samples), in each half
  cells     `--effective-samples cells`: n_cells // k for every pair

The halves have the same correlations, so the fraction of pairs with p <= 0.05
(0.01) should be 0.05 (0.01). Pairs are also grouped by their effective
samples (overlap, smaller of the two halves) into five bins of equal size, to
see whether either rule is off for pairs of smooth genes (few effective
samples) or noisy ones.

usage: 11_diff_calibration.py [CELL_TYPE ...]   (default: all six)
Writes results/diff_calibration.tsv.
"""

import shutil
import sys
import tempfile
import warnings
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd

from gcen2.single_cell import build_by_celltype, write_celltype_networks

HERE = Path(__file__).resolve().parent
WORK = Path("work")
TYPES = sys.argv[1:] or ["CD14 Mono", "CD4 T", "B", "CD8 T", "NK", "FCGR3A Mono"]
N_GENES = 2000
META = dict(k=25, max_shared=10, max_metacells=1000, min_cells=20, use_rep="X_pca",
            n_pcs=30, seed=42)
NET = dict(method="pearson", cor=0.0, pval=1.0, fdr=False, signed=False, threads=8)
N_BINS = 5

adata = ad.read_h5ad(WORK / "pbmc.h5ad")
rows = []
for ct in TYPES:
    tag = ct.replace(" ", "_")
    genes = (WORK / f"genes_{tag}.txt").read_text().split()
    genes = sorted(np.random.default_rng(0).choice(genes, N_GENES, replace=False).tolist())
    half = {}
    for h in "AB":
        for c in (WORK / "splits" / f"{tag}_{h}.txt").read_text().split():
            half[c] = h
    sub = adata[list(half), genes].copy()
    sub.obs["half"] = [half[c] for c in sub.obs_names]
    res = {}
    for eff in ("overlap", "cells"):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            nets = build_by_celltype(sub, "half", effective=eff, **META, **NET)
        tmp = Path(tempfile.mkdtemp(prefix="diffcal_", dir=WORK))
        try:
            write_celltype_networks({h: nets[h] for h in "AB"}, tmp, diff=True, pval=0.05,
                                    threads=8)
            summary = pd.read_csv(tmp / "diff" / "summary.tsv", sep="\t")
            d = pd.read_csv(tmp / "diff" / summary["file"][0], sep="\t",
                            usecols=["gene_a", "gene_b", "p_value", "category"])
        finally:
            shutil.rmtree(tmp)
        d = d[d["category"].str.startswith("shared")].reset_index(drop=True)
        res[eff] = (d, summary.iloc[0], nets)
    d_o, s_o, nets_o = res["overlap"]
    d_c, s_c, _ = res["cells"]
    assert (d_o["gene_a"].to_numpy() == d_c["gene_a"].to_numpy()).all()
    # every pair's effective samples: the smaller of the two halves' (overlap)
    n_pair = np.minimum.reduce([
        np.minimum(g[d_o["gene_a"]].to_numpy(), g[d_o["gene_b"]].to_numpy())
        for g in (nets_o[h].metadata["gene_effective_samples"] for h in "AB")])
    edges = np.unique(np.quantile(n_pair, np.linspace(0, 1, N_BINS + 1)))
    bins = np.clip(np.searchsorted(edges, n_pair, side="right") - 1, 0, len(edges) - 2)
    base = {"cell_type": ct, "metacells_A": nets_o["A"].metadata["n_samples"],
            "metacells_B": nets_o["B"].metadata["n_samples"]}
    for eff, (d, s, _) in res.items():
        p = d["p_value"].to_numpy()
        ok = np.isfinite(p)
        groups = [("all", np.ones(len(p), bool), "", "")] + [
            (f"bin{b + 1}", bins == b, edges[b], edges[b + 1]) for b in range(len(edges) - 1)]
        for name, m, lo, hi in groups:
            rows.append({**base, "effective": eff, "n_a": int(s["n_a"]), "n_b": int(s["n_b"]),
                         "pairs": name, "pair_n_eff_from": lo, "pair_n_eff_to": hi,
                         "tested": int((m & ok).sum()), "untested": int((m & ~ok).sum()),
                         "p<=0.05": float((p[m & ok] <= 0.05).mean()),
                         "p<=0.01": float((p[m & ok] <= 0.01).mean())})
            print({k: (round(v, 4) if isinstance(v, float) else v)
                   for k, v in rows[-1].items()}, flush=True)

df = pd.DataFrame(rows)
df.to_csv(HERE / "results/diff_calibration.tsv", sep="\t", index=False, float_format="%.6g")
print(df.round(4).to_string(index=False))
