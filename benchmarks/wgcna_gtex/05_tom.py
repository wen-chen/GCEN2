"""GCEN2's TOM edge score (gcen2 build --tom) against WGCNA's.

For the 5,000-gene runs (corrected and uncorrected): the soft power picked
like pickSoftThreshold, STRING precision of the top pairs by TOM over all
gene pairs, and the TOM of WGCNA's top 300,000 pairs (runs/wgcna/*/tom_top.tsv,
from blockwiseModules) next to GCEN2's. Signed network, as in 02_wgcna.R.
Run after run_all.sh; prints one line per data set.
"""

import os
import time

import numpy as np
import pandas as pd

from gcen2._matrix_io import read_matrix
from gcen2.core import _gcen_core

ROOT = "."  # run from the GTEx data directory
# STRING gold standard written by ../hdwgcna_pbmc/04_compare.py (work/string_gold.tsv
# of the PBMC data directory): GCEN2_STRING_GOLD, else work/string_gold.tsv here
GOLD = pd.read_csv(os.environ.get("GCEN2_STRING_GOLD", "work/string_gold.tsv"), sep="\t")

for data, expr in (("pc5", "expr_5000_pc5.tsv"), ("raw", "expr_5000.tsv")):
    m, genes, _ = read_matrix(f"{ROOT}/work/{expr}")
    n = len(genes)
    idx = {g: i for i, g in enumerate(genes)}
    t0 = time.time()
    ij, _, _, _, tom, info = _gcen_core.correlation.network_edges(
        m, "pearson", 8, 0.0, 1.0, True, False, True, 0.0)  # all pairs, signed, power auto
    seconds = time.time() - t0
    code = ij[:, 0].astype(np.int64) * n + ij[:, 1]
    a, b = GOLD["a"].map(idx), GOLD["b"].map(idx)
    ok = a.notna() & b.notna()
    gold_codes = (np.minimum(a[ok], b[ok]).astype(np.int64) * n
                  + np.maximum(a[ok], b[ok]).astype(np.int64)).to_numpy()
    is_gold = np.isin(code, gold_codes)
    order = np.argsort(-tom, kind="stable")
    prec = [is_gold[order[:k]].mean() for k in (1000, 10000, 100000)]

    w = pd.read_csv(f"{ROOT}/runs/wgcna/{data}_5000/tom_top.tsv", sep="\t")
    wi, wj = w["gene_a"].map(idx).to_numpy(), w["gene_b"].map(idx).to_numpy()
    wcode = np.minimum(wi, wj).astype(np.int64) * n + np.maximum(wi, wj)
    mine = pd.Series(tom, index=code).loc[wcode].to_numpy()
    overlap = len(np.intersect1d(code[order[:10000]], wcode[:10000])) / 10000
    print(f"{data}: power {info['power']:g} ({seconds:.1f} s); STRING P@1k/10k/100k by TOM "
          + " / ".join(f"{p:.2%}" for p in prec)
          + f"; vs WGCNA TOM: max relative difference {np.max(np.abs(mine / w['tom'] - 1)):.1e}, "
            f"top-10k overlap {overlap:.3f}")
