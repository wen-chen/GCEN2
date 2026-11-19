"""Compare GCEN2 and WGCNA on GTEx whole blood.

  * correlations: GCEN2's Pearson (all pairs of the first 300 genes, Python
    API) against WGCNA::cor;
  * STRING (experimental/database >= 700): precision of the top-N pairs by
    Pearson r (the same for both tools), by GCEN2's Spearman r and by
    WGCNA's TOM; and of the pairs within modules;
  * modules: number, genes, split-half ARI (samples halves, each corrected
    on its own), ARI between the tools;
  * time and memory (runs/time.txt).
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import rankdata
from sklearn.metrics import adjusted_rand_score

import gcen2

ROOT = Path.cwd()  # run from the GTEx data directory
RES = ROOT / "results"
RES.mkdir(exist_ok=True)
# STRING gold standard written by ../hdwgcna_pbmc/04_compare.py (work/string_gold.tsv
# of the PBMC data directory): GCEN2_STRING_GOLD, else work/string_gold.tsv here
GOLD = pd.read_csv(os.environ.get("GCEN2_STRING_GOLD", "work/string_gold.tsv"), sep="\t")
TOP_N = [1000, 10000, 100000]


def load(expr):
    x = pd.read_csv(ROOT / "work" / expr, sep="\t", index_col=0)
    return x.index.to_numpy(), x.to_numpy(np.float64)


def corr(X):
    Z = X - X.mean(axis=1, keepdims=True)
    Z /= np.linalg.norm(Z, axis=1, keepdims=True)
    return Z @ Z.T


def gold_matrix(genes):
    idx = {g: i for i, g in enumerate(genes)}
    a = GOLD["a"].map(idx)
    b = GOLD["b"].map(idx)
    ok = a.notna() & b.notna()
    G = np.zeros((len(genes), len(genes)), dtype=bool)
    i = np.minimum(a[ok], b[ok]).astype(int)
    j = np.maximum(a[ok], b[ok]).astype(int)
    G[i, j] = True
    return G


def precision(scores, g):
    o = np.argsort(-scores, kind="stable")
    return {f"P@{k}": g[o[:k]].mean() for k in TOP_N}


def gcen2_modules(path):
    m = pd.read_csv(path, sep="\t")
    m = m.drop_duplicates("gene")
    return dict(zip(m["gene"], m["module"].astype(str)))


def wgcna_modules(path):
    m = pd.read_csv(path, sep="\t")
    m = m[m["module"] != "grey"]
    return dict(zip(m["gene"], m["module"].astype(str)))


def module_stats(lab, genes, G):
    la = np.array([lab.get(x, "") for x in genes])
    iu = np.triu_indices(len(genes), 1)
    a, b = la[iu[0]], la[iu[1]]
    same = (a == b) & (a != "")
    g = G[iu]
    sizes = pd.Series([v for v in la if v]).value_counts()
    return {"modules": len(sizes), "module_genes": int((la != "").sum()),
            "largest": int(sizes.iloc[0]) if len(sizes) else 0,
            "within_P": g[same].mean(), "within_enr": g[same].mean() / g.mean()}


def ari(l1, l2):
    common = sorted(set(l1) & set(l2))
    return adjusted_rand_score([l1[x] for x in common], [l2[x] for x in common]), len(common)


def main():
    out = []
    # 1. correlations against WGCNA::cor
    genes, X = load("expr_5000_pc5.tsv")
    chk = pd.read_csv(ROOT / "runs/wgcna/pc5_5000/cor_check.tsv", sep="\t")
    df = pd.DataFrame(X[:300], index=genes[:300])
    net = gcen2.network.build(df, method="pearson", all_edges=True)
    mine = net.edges.set_index(["gene_a", "gene_b"])["correlation"]
    key = list(zip(chk["gene_a"], chk["gene_b"]))
    diff = np.abs(mine.reindex(key).to_numpy() - chk["cor"].to_numpy())
    print(f"GCEN2 vs WGCNA::cor, {len(key)} pairs: max |diff| {diff.max():.2e}")

    # 2. STRING precision and modules
    rows = []
    for data, expr in (("pc5", "expr_5000_pc5.tsv"), ("raw", "expr_5000.tsv")):
        genes, X = load(expr)
        G = gold_matrix(genes)
        iu = np.triu_indices(len(genes), 1)
        g = G[iu]
        base = g.mean()
        pear = corr(X)[iu]
        spear = corr(rankdata(X, axis=1))[iu]
        rows.append({"data": data, "ranking": "Pearson r (both)", **precision(pear, g)})
        rows.append({"data": data, "ranking": "GCEN2 Spearman r", **precision(spear, g)})
        tom = pd.read_csv(ROOT / f"runs/wgcna/{data}_5000/tom_top.tsv", sep="\t")
        idx = {x: i for i, x in enumerate(genes)}
        a = tom["gene_a"].map(idx).to_numpy()
        b = tom["gene_b"].map(idx).to_numpy()
        hit = G[np.minimum(a, b), np.maximum(a, b)]
        rows.append({"data": data, "ranking": "WGCNA TOM",
                     **{f"P@{k}": hit[:k].mean() for k in TOP_N}})
        for r in rows[-3:]:
            r["baseline"] = base
        run = f"{data}_5000"
        lg = gcen2_modules(ROOT / f"runs/gcen2/{run}/modules.tsv")
        lw = wgcna_modules(ROOT / f"runs/wgcna/{run}/modules.tsv")
        for tool, lab in (("GCEN2", lg), ("WGCNA", lw)):
            row = {"data": data, "tool": tool, **module_stats(lab, genes, G)}
            if data == "pc5":
                la = (gcen2_modules if tool == "GCEN2" else wgcna_modules)(
                    ROOT / f"runs/{tool.lower()}/pc5_5000_A/modules.tsv")
                lb = (gcen2_modules if tool == "GCEN2" else wgcna_modules)(
                    ROOT / f"runs/{tool.lower()}/pc5_5000_B/modules.tsv")
                row["split_ARI"], row["split_genes"] = ari(la, lb)
            out.append(row)
        a_ri, n_c = ari(lg, lw)
        out.append({"data": data, "tool": "GCEN2 vs WGCNA", "between_ARI": a_ri,
                    "between_genes": n_c})
    # split-half agreement of the correlations themselves
    _, XA = load("expr_5000_pc5_A.tsv")
    _, XB = load("expr_5000_pc5_B.tsv")
    iu = np.triu_indices(len(XA), 1)
    print(f"split-half correlation of the Pearson r (pc5, 5,000 genes): "
          f"{np.corrcoef(corr(XA)[iu], corr(XB)[iu])[0, 1]:.3f}")

    pd.set_option("display.width", 250)
    s = pd.DataFrame(rows)
    m = pd.DataFrame(out)
    s.to_csv(RES / "string.tsv", sep="\t", index=False)
    m.to_csv(RES / "modules.tsv", sep="\t", index=False)
    print("\n== STRING precision\n" + s.round(4).to_string(index=False))
    print("\n== modules\n" + m.round(4).to_string(index=False))

    # 3. time and memory
    t = pd.read_csv(ROOT / "runs/time.txt", sep=" ", header=None,
                    names=["step", "run", "sec", "s", "kb", "KB"], comment="C")
    t = t[["step", "run", "sec", "kb"]]
    g2 = t[t.step.str.startswith("gcen2")].groupby("run").agg(
        gcen2_s=("sec", "sum"), gcen2_GB=("kb", "max"))
    w = t[t.step == "wgcna"].set_index("run")[["sec", "kb"]].rename(
        columns={"sec": "wgcna_s", "kb": "wgcna_GB"})
    tt = g2.join(w)
    tt["gcen2_GB"] /= 1048576
    tt["wgcna_GB"] /= 1048576
    tt["speedup"] = tt["wgcna_s"] / tt["gcen2_s"]
    tt.to_csv(RES / "time.tsv", sep="\t")
    print("\n== time and memory\n" + tt.round(2).to_string())


if __name__ == "__main__":
    main()
