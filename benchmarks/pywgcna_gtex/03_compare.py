"""PyWGCNA against WGCNA (step by step and blockwiseModules) and GCEN2.

Run from the GTEx data directory after run_all.sh, with the STRING gold
standard in GCEN2_STRING_GOLD (as ../wgcna_gtex/04_compare.py). For each of
the six runs:

  power         chosen by PyWGCNA, by WGCNA (step by step; blockwiseModules
                and GCEN2 choose the same, ../wgcna_gtex/results/tom_modules.tsv)
  modules       number, genes in modules (grey = unassigned), smallest module
  ARI           adjusted Rand index over all genes, unassigned genes as one
                group: PyWGCNA vs WGCNA step by step (the pipeline PyWGCNA
                reimplements), PyWGCNA vs blockwiseModules, WGCNA step by step
                vs blockwiseModules (how much the two WGCNA pipelines differ);
                GCEN2 vs blockwiseModules is in ../wgcna_gtex (1.0 in all six)
  STRING        within-module STRING precision over all pairs of genes in the
                same module, as enrichment over the baseline of all pairs
  time, memory  mean wall time and peak memory over the repetitions

and for pc5_5000 and raw_5000 PyWGCNA's TOM against blockwiseModules' over
WGCNA's 300,000 top pairs (largest relative difference; overlap of the top
10,000). Writes results/compare.tsv and results/tom.tsv (next to this script).
"""

import os
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
REF = HERE.parent / "wgcna_gtex" / "results"
RUNS = ["pc5_5000", "pc5_5000_A", "pc5_5000_B", "raw_5000", "pc5_10000", "pc5_11944"]
GOLD = pd.read_csv(os.environ.get("GCEN2_STRING_GOLD", "work/string_gold.tsv"), sep="\t")


def modules(path):
    m = pd.read_csv(path, sep="\t").drop_duplicates("gene").set_index("gene")["module"]
    return m.astype(str)


def string_enrichment(m):
    """Within-module STRING precision over the baseline of all pairs."""
    genes = m.index
    idx = {g: i for i, g in enumerate(genes)}
    a, b = GOLD["a"].map(idx), GOLD["b"].map(idx)
    ok = a.notna() & b.notna()
    a, b = a[ok].astype(int).to_numpy(), b[ok].astype(int).to_numpy()
    n = len(genes)
    base = len(set(zip(np.minimum(a, b), np.maximum(a, b)))) / (n * (n - 1) / 2)
    lab = m.to_numpy()
    assigned = lab != "grey"
    sizes = pd.Series(lab[assigned]).value_counts()
    within = int((sizes * (sizes - 1) // 2).sum())
    pairs = set(zip(np.minimum(a, b), np.maximum(a, b)))
    hits = sum(1 for p, q in pairs if assigned[p] and lab[p] == lab[q])
    return hits / within / base


def summary(m):
    sizes = m[m != "grey"].value_counts()
    return len(sizes), int(sizes.sum()), int(sizes.min())


times = pd.read_csv("runs/pywgcna/time.txt", sep=" ",
                    names=["tool", "run", "rep", "seconds", "peak_kb"])
times = times.groupby(["tool", "run"]).agg(seconds=("seconds", "mean"),
                                          peak_kb=("peak_kb", "max"), reps=("rep", "size"))
ref = pd.read_csv(REF / "tom_modules.tsv", sep="\t").set_index("run")

rows = []
for run in RUNS:
    py = modules(f"runs/pywgcna/{run}/modules.tsv")
    st = modules(f"runs/wgcna_step/{run}/modules.tsv").reindex(py.index)
    bw = modules(f"runs/wgcna/{run}/modules.tsv").reindex(py.index)
    pinfo = pd.read_csv(f"runs/pywgcna/{run}/info.tsv", sep="\t").iloc[0]
    sinfo = pd.read_csv(f"runs/wgcna_step/{run}/timing.tsv", sep="\t").iloc[0]
    row = {"run": run, "power_pywgcna": int(pinfo["power"]), "power_wgcna": int(sinfo["power"])}
    for name, m in (("pywgcna", py), ("wgcna_step", st), ("blockwise", bw)):
        k, g, smallest = summary(m)
        row |= {f"{name}_modules": k, f"{name}_genes": g, f"{name}_smallest": smallest,
                f"{name}_STRING_enrichment": round(string_enrichment(m), 1)}
    row |= {"ARI_pywgcna_vs_step": round(adjusted_rand_score(st, py), 3),
            "ARI_pywgcna_vs_blockwise": round(adjusted_rand_score(bw, py), 3),
            "ARI_step_vs_blockwise": round(adjusted_rand_score(bw, st), 3)}
    for tool in ("pywgcna", "wgcna_step"):
        t = times.loc[(tool, run)]
        row |= {f"{tool}_s": round(t["seconds"], 1), f"{tool}_GB": round(t["peak_kb"] / 1048576, 2)}
    row |= {"blockwise_s": round(ref.loc[run, "wgcna_s"], 1),
            "blockwise_GB": round(ref.loc[run, "wgcna_GB"], 2),
            "gcen2_s": round(ref.loc[run, "gcen2_s"], 1),
            "gcen2_GB": round(ref.loc[run, "gcen2_GB"], 2)}
    rows.append(row)
    print(row, flush=True)
df = pd.DataFrame(rows)
df.to_csv(OUT / "compare.tsv", sep="\t", index=False)

tom_rows = []
for run in ("pc5_5000", "raw_5000"):
    w = pd.read_csv(f"runs/wgcna/{run}/tom_top.tsv", sep="\t")
    p = pd.read_csv(f"runs/pywgcna/{run}/tom_top.tsv", sep="\t")
    for d in (w, p):
        d.index = pd.MultiIndex.from_arrays([np.minimum(d.gene_a, d.gene_b),
                                             np.maximum(d.gene_a, d.gene_b)])
    j = w["tom"].to_frame("wgcna").join(p["tom"].rename("pywgcna"), how="left")
    covered = j["pywgcna"].notna()
    rel = ((j["pywgcna"] - j["wgcna"]).abs() / j["wgcna"])[covered]
    tom_rows.append({"run": run, "pairs": len(j), "in_pywgcna_top": int(covered.sum()),
                     "max_relative_difference": float(f"{rel.max():.2e}"),
                     "top10000_overlap": len(set(w.index[:10000]) & set(p.index[:10000]))})
    print(tom_rows[-1], flush=True)
pd.DataFrame(tom_rows).to_csv(OUT / "tom.tsv", sep="\t", index=False)
