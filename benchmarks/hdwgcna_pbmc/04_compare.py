"""Compare gcen2 and hdWGCNA runs on PBMC cell types.

usage: 04_compare.py [TOOL ...]   (default: gcen2 hdwgcna, whichever exist)

Per cell type and tool:
  * STRING gold standard (experimental or database channel >= 700, so the
    co-expression channel plays no part): precision of the top-N gene pairs
    and its enrichment over all pairs, ranking by correlation on the
    metacells (both tools) and, for hdWGCNA, by TOM; same for the pairs
    inside modules.
  * split-half reproducibility: correlation of the gene-gene correlations
    of two disjoint halves of the cells, overlap of their top 10,000 pairs,
    adjusted Rand index of their modules.
  * seed stability: the same comparisons between the default run (seed
    42) and runs with other metacell seeds.
Between tools: correlation of their gene-gene correlations, module ARI.
Writes results/*.tsv and prints them.
"""

import gzip
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path.cwd()  # run from the PBMC data directory
RUNS = ROOT / "runs"
RES = ROOT / "results"
RES.mkdir(exist_ok=True)
TYPES = ["CD14 Mono", "CD4 T", "B", "CD8 T", "NK", "FCGR3A Mono"]
SEEDS = [1, 2]  # runs <type>_seed<n>, compared with the default seed (42)
TOP_N = [1000, 10000, 100000]


def string_gold(min_score=700):
    cache = ROOT / "work" / "string_gold.tsv"
    if cache.exists():
        return pd.read_csv(cache, sep="\t")
    names = {}
    with gzip.open(ROOT / "data/string/9606.protein.info.v12.0.txt.gz", "rt") as fh:
        next(fh)
        for line in fh:
            f = line.split("\t", 2)
            names[f[0]] = f[1]
    rows = []
    with gzip.open(ROOT / "data/string/9606.protein.links.detailed.v12.0.txt.gz", "rt") as fh:
        next(fh)
        for line in fh:
            p1, p2, _, _, _, _, exp, db, _, _ = line.split()
            if int(exp) >= min_score or int(db) >= min_score:
                a, b = names.get(p1), names.get(p2)
                if a and b and a < b:
                    rows.append((a, b))
    gold = pd.DataFrame(rows, columns=["a", "b"]).drop_duplicates()
    gold.to_csv(cache, sep="\t", index=False)
    return gold


def corr_matrix(path):
    """Gene-gene Pearson correlations of a metacells x genes table."""
    df = pd.read_csv(path, sep="\t", index_col=0)
    x = df.to_numpy(np.float64)
    x = x - x.mean(axis=0)
    norm = np.sqrt((x * x).sum(axis=0))
    keep = norm > 0
    z = x[:, keep] / norm[keep]
    return z.T @ z, df.columns[keep].to_numpy()


def upper(c):
    iu = np.triu_indices(len(c), 1)
    return c[iu], iu


def gold_mask(genes, iu, gold):
    idx = {g: i for i, g in enumerate(genes)}
    n = len(genes)
    ga = gold["a"].map(idx)
    gb = gold["b"].map(idx)
    ok = ga.notna() & gb.notna()
    i = np.minimum(ga[ok], gb[ok]).astype(np.int64)
    j = np.maximum(ga[ok], gb[ok]).astype(np.int64)
    codes = set((i * n + j).tolist())
    return np.fromiter(((a * n + b) in codes for a, b in zip(*iu)), bool, len(iu[0]))


def precision_at(scores, is_gold):
    order = np.argsort(-scores, kind="stable")
    base = is_gold.mean()
    return {f"P@{n}": is_gold[order[:n]].mean() for n in TOP_N} | {
        f"enr@{n}": is_gold[order[:n]].mean() / base for n in TOP_N} | {"baseline": base}


def module_labels(path, genes):
    m = pd.read_csv(path, sep="\t")
    col = "module"
    m = m[m[col].astype(str) != "grey"]
    first = m.drop_duplicates("gene")  # gcen2 modules may overlap: first one
    lab = dict(zip(first["gene"], first[col].astype(str)))
    return lab


def module_string(lab, genes, iu, is_gold):
    """STRING precision of the gene pairs within the same module."""
    g = np.asarray(genes)
    la = np.array([lab.get(x, "") for x in g])
    a, b = la[iu[0]], la[iu[1]]
    same = (a == b) & (a != "")
    return {"module_genes": int((la != "").sum()), "n_modules": len(set(lab.values())),
            "within_pairs": int(same.sum()),
            "within_P": is_gold[same].mean() if same.any() else np.nan,
            "within_enr": (is_gold[same].mean() / is_gold.mean()) if same.any() else np.nan}


def ari(lab1, lab2):
    common = sorted(set(lab1) & set(lab2))
    if len(common) < 2:
        return np.nan, len(common)
    return adjusted_rand_score([lab1[g] for g in common], [lab2[g] for g in common]), len(common)


def compare_runs(tool, run_a, run_b):
    """cor of the gene-gene correlations, top-10k overlap and module ARI of
    two runs of one tool, or None when one is missing."""
    ha, ea = tool_files(tool, run_a)
    hb, eb = tool_files(tool, run_b)
    if not (ea.exists() and eb.exists()):
        return None
    ca, ga = corr_matrix(ea)
    cb, gb = corr_matrix(eb)
    common = np.intersect1d(ga, gb)
    pa = pd.Index(ga).get_indexer(common)
    pb = pd.Index(gb).get_indexer(common)
    va, _ = upper(ca[np.ix_(pa, pa)])
    vb, _ = upper(cb[np.ix_(pb, pb)])
    top_a = set(np.argsort(-va)[:10000])
    top_b = set(np.argsort(-vb)[:10000])
    a_ri, n_common = ari(module_labels(ha / "modules.tsv", ga),
                         module_labels(hb / "modules.tsv", gb))
    return {"cor_of_cor": np.corrcoef(va, vb)[0, 1],
            "top10k_overlap": len(top_a & top_b) / 10000,
            "module_ARI": a_ri, "module_genes_both": n_common}


def tool_files(tool, run):
    d = RUNS / tool / run
    expr = d / ("metacells.tsv" if tool == "gcen2" else "datExpr.tsv")
    return d, expr


def main(tools):
    gold = string_gold()
    print(f"STRING gold pairs (experimental/database >= 700): {len(gold)}")
    rows, rep_rows, cross_rows, seed_rows = [], [], [], []
    for ct in TYPES:
        s = ct.replace(" ", "_")
        cache = {}
        for tool in tools:
            d, expr = tool_files(tool, s)
            if not expr.exists():
                continue
            c, genes = corr_matrix(expr)
            v, iu = upper(c)
            is_gold = gold_mask(genes, iu, gold)
            cache[tool] = (c, genes, v)
            timing = pd.read_csv(d / "timing.tsv", sep="\t").iloc[0].to_dict()
            row = {"cell_type": ct, "tool": tool, "ranking": "correlation",
                   **precision_at(v, is_gold)}
            lab = module_labels(d / "modules.tsv", genes)
            row |= module_string(lab, genes, iu, is_gold)
            row |= {k: timing[k] for k in timing if k.startswith("t_") or k in (
                "metacells", "power", "edges")}
            rows.append(row)
            if tool == "hdwgcna" and (d / "tom_top.tsv").exists():
                tom = pd.read_csv(d / "tom_top.tsv", sep="\t")
                idx = {g: i for i, g in enumerate(genes)}
                n = len(genes)
                i = tom["gene_a"].map(idx).to_numpy()
                j = tom["gene_b"].map(idx).to_numpy()
                lo, hi = np.minimum(i, j), np.maximum(i, j)
                gold_codes = set((iu[0][is_gold].astype(np.int64) * n
                                  + iu[1][is_gold]).tolist())
                hit = np.fromiter(((a * n + b) in gold_codes for a, b in zip(lo, hi)), bool)
                rows.append({"cell_type": ct, "tool": tool, "ranking": "TOM",
                             **{f"P@{k}": hit[:k].mean() for k in TOP_N},
                             **{f"enr@{k}": hit[:k].mean() / is_gold.mean() for k in TOP_N},
                             "baseline": is_gold.mean()})
            # split halves
            r = compare_runs(tool, f"{s}_A", f"{s}_B")
            if r:
                rep_rows.append({"cell_type": ct, "tool": tool, **r})
            for seed in SEEDS:
                r = compare_runs(tool, s, f"{s}_seed{seed}")
                if r:
                    seed_rows.append({"cell_type": ct, "tool": tool, "seed": seed, **r})
        if "gcen2" in cache and "hdwgcna" in cache:
            (c1, g1, _), (c2, g2, _) = cache["gcen2"], cache["hdwgcna"]
            common = np.intersect1d(g1, g2)
            p1 = pd.Index(g1).get_indexer(common)
            p2 = pd.Index(g2).get_indexer(common)
            v1, _ = upper(c1[np.ix_(p1, p1)])
            v2, _ = upper(c2[np.ix_(p2, p2)])
            l1 = module_labels(RUNS / "gcen2" / s / "modules.tsv", g1)
            l2 = module_labels(RUNS / "hdwgcna" / s / "modules.tsv", g2)
            a_ri, n_common = ari(l1, l2)
            cross_rows.append({"cell_type": ct, "genes": len(common),
                               "cor_of_cor": np.corrcoef(v1, v2)[0, 1],
                               "top10k_overlap": len(set(np.argsort(-v1)[:10000])
                                                     & set(np.argsort(-v2)[:10000])) / 10000,
                               "module_ARI": a_ri, "module_genes_both": n_common})
    pd.set_option("display.width", 250)
    pd.set_option("display.max_columns", 40)
    for name, rr in (("string", rows), ("split_half", rep_rows), ("seed_stability", seed_rows),
                     ("between_tools", cross_rows)):
        df = pd.DataFrame(rr)
        df.to_csv(RES / f"{name}.tsv", sep="\t", index=False)
        print(f"\n== {name}\n{df.round(4).to_string(index=False)}")


if __name__ == "__main__":
    main(sys.argv[1:] or ["gcen2", "hdwgcna"])
