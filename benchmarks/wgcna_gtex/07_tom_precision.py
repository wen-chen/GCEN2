"""Topological overlap in single against double precision: the modules of
`gcen2 module -e` on GTEx whole blood (5,000 genes, 5 PCs regressed out).

Run from the GTEx data directory (after 01_prepare.py and run_all.sh).
gcen2.tom_modules.tom_modules with the settings of 06_tom_modules.py (signed,
power auto, minModuleSize 30, mergeCutHeight 0.25, deepSplit 2), with its
1 - TOM matrix replaced; everything after it (average-linkage tree, dynamic
hybrid tree cut with PAM, kME reassignment, merging) is the same code:

  double          as `gcen2 module -e` (checked against
                  runs/gcen2/pc5_5000/tom_modules.tsv)
  TOM float       the TOM computed in single precision, as `gcen2 build --tom`
                  scores edges (adjacency and its products in float;
                  network.build(all_edges=True, tom=True) on all pairs, same
                  power), 1 - TOM then in double
  TOM to float    the double TOM rounded to single precision, 1 - TOM in double
  1 - TOM float   the double 1 - TOM rounded to single precision, i.e. the
                  dissimilarity matrix stored in float32

Measured per variant: the largest absolute difference of 1 - TOM from double;
the adjusted Rand index of the modules after the tree cut and of the final
modules with double's; genes in another module (modules matched by their
largest overlap; unassigned genes as one group).

Writes results/tom_precision.tsv (one row per variant),
results/tom_precision_modules.tsv (gene and the final module of each variant)
and results/tom_precision_error.tsv (histogram of the absolute differences of
1 - TOM from double, per variant).
"""

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

import gcen2
import gcen2.tom_modules as tm
from gcen2.core import require_extension

ROOT = Path.cwd()
HERE = Path(__file__).resolve().parent
THREADS = 8
OPTS = dict(method="pearson", signed=True, min_size=30, merge_cut=0.25, deep_split=2,
            threads=THREADS)

expr = pd.read_csv(ROOT / "work/expr_5000_pc5.tsv", sep="\t", index_col=0)
genes = list(expr.index)
core = require_extension()
kept = {}


class Ext:
    """The extension, with 1 - TOM replaced by transform(1 - TOM) and kept."""

    def __init__(self, transform=None):
        self.transform = transform
        self.module = self

    def __getattr__(self, name):
        return getattr(core.module, name)

    def tom_dissimilarity(self, *args):
        diss, info = core.module.tom_dissimilarity(*args)
        if self.transform is not None:
            diss[...] = self.transform(diss)
        kept["diss"] = np.array(diss, copy=True)
        return diss, info


orig_cut = tm.cutree_hybrid


def cut_and_keep(*a, **k):
    lab = orig_cut(*a, **k)
    kept["cut"] = lab.copy()
    return lab


def single_tom(power):
    """1 - TOM of all pairs from the single-precision edge TOM (build --tom)."""
    net = gcen2.network.build(expr, method="pearson", signed=True, all_edges=True, tom=True,
                              power=power, threads=THREADS)
    pos = pd.Index(genes)
    a = pos.get_indexer(net.edges["gene_a"])
    b = pos.get_indexer(net.edges["gene_b"])
    assert (a >= 0).all() and (b >= 0).all()
    assert len(net.edges) == len(genes) * (len(genes) - 1) // 2
    tom = np.zeros((len(genes), len(genes)))
    tom[a, b] = net.edges["tom"].to_numpy()
    tom[b, a] = tom[a, b]
    diss = 1.0 - tom
    np.fill_diagonal(diss, 0.0)
    return diss


def f32(x):
    return x.astype(np.float32).astype(np.float64)


def moved(x, y):
    """Genes whose module in y is not the y-module overlapping x's module most."""
    best = pd.crosstab(x, y).idxmax(axis=1)
    return int((pd.Series(y) != best.reindex(x).to_numpy()).sum())


tm.cutree_hybrid = cut_and_keep
runs = {}
tm.require_extension = lambda: Ext()
r = tm.tom_modules(expr, **OPTS)
power = r.power
runs["double"] = (r.labels.copy(), kept["cut"].copy(), kept.pop("diss"))
diss_float = single_tom(power)
for name, transform in (("TOM float", lambda d: diss_float),
                        ("TOM to float", lambda d: 1.0 - f32(1.0 - d)),
                        ("1 - TOM float", f32)):
    tm.require_extension = lambda t=transform: Ext(t)
    r = tm.tom_modules(expr, **OPTS)
    runs[name] = (r.labels.copy(), kept["cut"].copy(), kept.pop("diss"))
tm.cutree_hybrid = orig_cut

ref = pd.read_csv(ROOT / "runs/gcen2/pc5_5000/tom_modules.tsv", sep="\t").set_index("gene")
d_lab, d_cut, d_diss = runs["double"]
assert adjusted_rand_score(ref["module"].reindex(genes).fillna(-1), d_lab) == 1.0
iu = np.triu_indices(len(genes), 1)
bins = np.arange(-12, -5.75, 0.25)
rows, hist = [], {"log10_abs_from": bins[:-1], "log10_abs_to": bins[1:]}
for name, (lab, cut, diss) in runs.items():
    err = np.abs(diss[iu] - d_diss[iu])
    if name != "double":
        hist[name] = np.histogram(np.log10(np.maximum(err, 1e-12)), bins=bins)[0]
    rows.append({"variant": name, "power": power, "diss_abs_diff_max": float(err.max()),
                 "diss_abs_diff_median": float(np.median(err)),
                 "cut_modules": len(set(cut) - {0}), "cut_ARI": adjusted_rand_score(d_cut, cut),
                 "cut_genes_moved": moved(d_cut, cut), "modules": len(set(lab) - {0}),
                 "final_ARI": adjusted_rand_score(d_lab, lab),
                 "final_genes_moved": moved(d_lab, lab), "unassigned": int((lab == 0).sum())})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(HERE / "results/tom_precision.tsv", sep="\t", index=False,
                          float_format="%.6g")
pd.DataFrame(hist).to_csv(HERE / "results/tom_precision_error.tsv", sep="\t", index=False)
pd.DataFrame({"gene": genes, **{n: v[0] for n, v in runs.items()}}).to_csv(
    HERE / "results/tom_precision_modules.tsv", sep="\t", index=False)
