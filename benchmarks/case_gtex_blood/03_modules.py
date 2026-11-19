"""What the GTEx whole-blood modules are: GO enrichment, hub genes, donor traits.

Run from the GTEx data directory after 02_run.sh. For every module of
runs/case/modules.tsv: GO biological-process enrichment against the 11,944
genes (gcen2.annotation.enrich: one-sided hypergeometric test, GO annotation
propagated along is_a; its fdr column: Benjamini-Hochberg over every GO term
annotated in the background), the five hub genes (highest kME) and the eigengene correlations with
the donor traits, in all 803 samples and in each split half
(../wgcna_gtex: work/half_A.txt, half_B.txt), which checks that an
association holds in two independent sets of donors.

Writes to results/ (next to this script):
  module_go.tsv          every module x term with FDR <= 0.05
  module_summary.tsv     module, genes, hub genes, top 3 terms, trait correlations
  trait_halves.tsv       eigengene-trait correlations in each half
  pc_traits.tsv          the top 6 expression PCs of the uncorrected data
                         (work/expr_11944.tsv) against the traits: what the
                         5 regressed PCs removed
"""

from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import pearsonr

from gcen2.annotation import enrich, load_association, load_go

OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
R, W = Path("runs/case"), Path("work/case")

go = load_go("data/go/go-basic.obo")
assoc = load_association(W / "go_assoc.tsv", go)
background = (W / "genes.txt").read_text().split()
modules = pd.read_csv(R / "modules.tsv", sep="\t")
hubs = pd.read_csv(R / "hubs.tsv", sep="\t")
tc = pd.read_csv(R / "trait_cor.tsv", sep="\t")
traits = pd.read_csv(W / "traits.tsv", sep="\t", index_col=0)
me = pd.read_csv(R / "eigengenes.tsv", sep="\t", index_col=0)
if me.shape[0] != traits.shape[0]:      # modules as rows
    me = me.T
me = me.loc[traits.index]

go_rows, summary = [], []
for mod, d in modules.groupby("module"):
    r = enrich(d["gene"].tolist(), background, go, assoc, pval=1.0).results
    r = r[r["namespace"] == "biological_process"]
    r = r.sort_values(["p_value", "term_id"]).loc[lambda x: x["fdr"] <= 0.05]
    r.insert(0, "module", mod)
    go_rows.append(r)
    row = {"module": f"ME{mod}", "genes": len(d),
           "hub_genes": ",".join(hubs.loc[hubs["module"] == mod, "gene"].head(5)),
           "top_terms": "; ".join(f"{n} ({e:.1f}x, FDR {f:.1e})" for n, e, f in
                                  zip(r["name"].head(3), r["enrichment"].head(3),
                                      r["fdr"].head(3)))}
    for t in traits.columns:
        c = tc[(tc["module"] == f"ME{mod}") & (tc["trait"] == t)].iloc[0]
        row[f"r_{t}"], row[f"p_{t}"] = round(c["cor"], 3), float(f"{c['p']:.2e}")
    summary.append(row)
pd.concat(go_rows).to_csv(OUT / "module_go.tsv", sep="\t", index=False)
pd.DataFrame(summary).to_csv(OUT / "module_summary.tsv", sep="\t", index=False)

halves = []
for h in ("A", "B"):
    samples = [s for s in Path(f"work/half_{h}.txt").read_text().split() if s in traits.index]
    for mod in me.columns:
        for t in traits.columns:
            x, y = me.loc[samples, mod].to_numpy(), traits.loc[samples, t].to_numpy()
            ok = ~np.isnan(y)
            c, p = pearsonr(x[ok], y[ok])
            halves.append({"half": h, "module": mod, "trait": t, "n": int(ok.sum()),
                           "cor": round(c, 3), "p": float(f"{p:.2e}")})
pd.DataFrame(halves).to_csv(OUT / "trait_halves.tsv", sep="\t", index=False)

# what the 5 regressed PCs carried: the top PCs of the uncorrected data vs the traits
x = pd.read_csv("work/expr_11944.tsv", sep="\t", index_col=0)[traits.index]
z = (x.to_numpy() - x.to_numpy().mean(axis=1, keepdims=True)).T
u, s, _ = np.linalg.svd(z, full_matrices=False)
pcs = []
for i in range(6):
    row = {"pc": i + 1, "variance": round(float(s[i] ** 2 / (s ** 2).sum()), 3)}
    for t in traits.columns:
        ok = ~traits[t].isna().to_numpy()
        row[t] = round(float(np.corrcoef(u[ok, i], traits[t].to_numpy()[ok])[0, 1]), 2)
    pcs.append(row)
pd.DataFrame(pcs).to_csv(OUT / "pc_traits.tsv", sep="\t", index=False)
print(pd.DataFrame(summary)[["module", "genes", "hub_genes", "r_age", "p_age",
                             "r_female", "p_female"]].to_string(index=False))
