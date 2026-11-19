"""Cross-validated gene function prediction with RWR on the GTEx blood network.

Run from the GTEx data directory after 02_run.sh. For every GO
biological-process term with 20-200 annotated genes in the network
(annotation propagated along is_a; IEP, HEP, ND and NOT left out by
01_prepare.py), the annotated genes are split into 3 folds; each fold in turn
is hidden, the other two are the seeds, and the genes are ranked by score. The
AUROC is that of the hidden genes against the genes annotated to other BP
terms but not this one (unannotated genes are not counted as negatives); a
term's AUROC is the mean over its folds.

Scores:
  rwr            gcen2 RWR (restart 0.5, unweighted, as `gcen2 rwr`)
  rwr_weighted   RWR with |r| as edge weights (gcen2.rwr.predict(weighted_network=True))
  neighbours     neighbour voting: the fraction of a gene's neighbours that are seeds
  rwr_shuffled   RWR on a degree-preserving shuffle of the network (gcen2 shuffle):
                 same degrees, random partners
  degree         the gene's degree alone (no seeds)

The degree-only and shuffled-network controls measure how much of the
performance comes from which genes are hubs rather than from who is connected
to whom (Gillis & Pavlidis, PLoS One 2011). The RWR here calls the C++ kernel
behind gcen2.rwr.predict directly with the transition matrix built once; the
script checks on the first term that the probabilities equal predict()'s.

Writes results/rwr_cv.tsv (one row per term) and results/rwr_cv_summary.tsv.
"""

import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.sparse import diags
from scipy.stats import rankdata, wilcoxon

from gcen2.annotation import load_association, load_go
from gcen2.core import extension
from gcen2.io import read_network
from gcen2.rwr import predict

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
R, W = Path("runs/case"), Path("work/case")
GAMMA, FOLDS, MIN_GENES, MAX_GENES = 0.5, 3, 20, 200
rng = np.random.default_rng(42)
core = extension()

go = load_go("data/go/go-basic.obo")
assoc = load_association(W / "go_assoc.tsv", go)
bp = {t.id for t in go.terms.values() if t.namespace == "biological_process"}


def transition(net, genes, weighted):
    """Column-stochastic matrix as in gcen2.rwr.predict."""
    a = net.sparse_adjacency(genes, weighted=weighted)
    col = np.asarray(a.sum(axis=0)).ravel()
    col[col == 0] = 1.0
    return (a @ diags(1.0 / col)).tocsr(), a


def rwr(w, p0):
    pt, _ = core.rwr.random_walk_restart(w.indptr.astype(np.int64), w.indices.astype(np.int32),
                                         w.data, p0, GAMMA, 1_000_000, 1e-10)
    return np.asarray(pt)


net = read_network(R / "network.tsv")
genes = sorted(set(net.edges["gene_a"]) | set(net.edges["gene_b"]))
idx = {g: i for i, g in enumerate(genes)}
n = len(genes)
w_unw, adj = transition(net, genes, weighted=False)
w_wt, _ = transition(net, genes, weighted=True)
shuf = read_network(R / "network_shuffled.tsv")
assert sorted(set(shuf.edges["gene_a"]) | set(shuf.edges["gene_b"])) == genes
w_shuf, _ = transition(shuf, genes, weighted=False)
degree = np.asarray(adj.sum(axis=0)).ravel()

# genes with any BP annotation form the universe; terms by their genes in it
gene_bp = {g: assoc[g] & bp for g in genes if assoc.get(g, set()) & bp}
universe = np.array(sorted(idx[g] for g in gene_bp))
term_genes: dict[str, list[int]] = {}
for g, ts in gene_bp.items():
    for t in ts:
        term_genes.setdefault(t, []).append(idx[g])
terms = sorted(t for t, gs in term_genes.items() if MIN_GENES <= len(gs) <= MAX_GENES)
print(f"{n:,} genes in the network, {len(universe):,} with a BP annotation, "
      f"{len(terms):,} terms with {MIN_GENES}-{MAX_GENES} genes")


def auroc(score, pos, neg):
    s = np.concatenate([score[pos], score[neg]])
    r = rankdata(s)
    npos = len(pos)
    return (r[:npos].sum() - npos * (npos + 1) / 2) / (npos * len(neg))


checked = False
rows = []
t0 = time.time()
for k, term in enumerate(terms):
    pos_all = np.array(sorted(term_genes[term]))
    folds = np.array_split(rng.permutation(pos_all), FOLDS)
    res = {m: [] for m in ("rwr", "rwr_weighted", "neighbours", "rwr_shuffled", "degree")}
    in_term = np.zeros(n, bool)
    in_term[pos_all] = True
    for f in range(FOLDS):
        hidden = folds[f]
        seeds = np.concatenate([folds[j] for j in range(FOLDS) if j != f])
        neg = universe[~in_term[universe]]
        p0 = np.zeros(n)
        p0[seeds] = 1.0 / len(seeds)
        s_rwr = rwr(w_unw, p0)
        if not checked:
            ref = predict(net, [genes[i] for i in seeds], gamma=GAMMA, prob_cutoff=0.0)
            ref = ref.ranked_genes.set_index("gene")["probability"]
            mine = pd.Series(s_rwr, index=genes).drop([genes[i] for i in seeds])
            assert np.array_equal(mine.loc[ref.index].to_numpy(), ref.to_numpy())
            checked = True
        seedvec = np.zeros(n)
        seedvec[seeds] = 1.0
        scores = {"rwr": s_rwr, "rwr_weighted": rwr(w_wt, p0),
                  "neighbours": (adj @ seedvec) / np.maximum(degree, 1),
                  "rwr_shuffled": rwr(w_shuf, p0), "degree": degree.astype(float)}
        for m, s in scores.items():
            res[m].append(auroc(s, hidden, neg))
    rows.append({"term": term, "name": go.terms[term].name, "genes": len(pos_all),
                 **{m: round(float(np.mean(v)), 4) for m, v in res.items()}})
    if (k + 1) % 200 == 0:
        print(f"{k + 1:,} terms, {time.time() - t0:.0f} s", flush=True)

df = pd.DataFrame(rows)
df.to_csv(OUT / "rwr_cv.tsv", sep="\t", index=False)
methods = ["rwr", "rwr_weighted", "neighbours", "rwr_shuffled", "degree"]
summary = pd.DataFrame({
    "method": methods,
    "median_AUROC": [round(df[m].median(), 3) for m in methods],
    "mean_AUROC": [round(df[m].mean(), 3) for m in methods],
    "terms_AUROC_ge_0.7": [int((df[m] >= 0.7).sum()) for m in methods],
    "terms_better_than_shuffled": [int((df[m] > df["rwr_shuffled"]).sum()) for m in methods],
    "terms_better_than_degree": [int((df[m] > df["degree"]).sum()) for m in methods],
    "median_gain_over_shuffled": [round(float((df[m] - df["rwr_shuffled"]).median()), 3)
                                  for m in methods],
    # two-sided Wilcoxon signed-rank test of RWR against each other score
    "wilcoxon_p_vs_rwr": [float("nan") if m == "rwr" else
                          float(f"{wilcoxon(df['rwr'], df[m]).pvalue:.1e}") for m in methods],
})
summary.insert(1, "terms", len(df))
summary.to_csv(OUT / "rwr_cv_summary.tsv", sep="\t", index=False)
print(summary.to_string(index=False))
print(f"{time.time() - t0:.0f} s")
