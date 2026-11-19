"""Parameter sensitivity on GTEx whole blood (5,000 genes, 5 PCs regressed out).

Run from the GTEx data directory (../wgcna_gtex: work/expr_5000_pc5.tsv) with
the STRING gold standard of ../hdwgcna_pbmc in GCEN2_STRING_GOLD (as
../wgcna_gtex/04_compare.py). One parameter at a time from the defaults of the
other benchmarks:

  correlation network   gcen2 build -m pearson -p 0.001 -f -s -c C, C = 0.3 ... 0.8,
                        then gcen2 module -s 0.5
  GCEN modules          gcen2 module -s S on the r >= 0.4 network, S = 0.3 ... 0.7
  WGCNA-style modules   gcen2 module -e --signed (defaults: power auto, deepSplit 2,
                        min size 30, merge cut 0.25), varying the power, deepSplit,
                        the minimum size and the merge cut

Per setting: edges, genes, STRING precision of the edges; modules, genes in
modules, largest module, STRING precision within modules and its enrichment
over all pairs; for WGCNA-style modules the adjusted Rand index with the
default, over all 5,000 genes (unassigned genes as one group) and over the
genes in modules in both.

Writes results/sensitivity_network.tsv, sensitivity_gcen_modules.tsv,
sensitivity_wgcna_modules.tsv (next to this script).
"""

import os
import shutil
import subprocess
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

GCEN2 = os.environ.get("GCEN2", "gcen2")
OUT = Path(__file__).resolve().parent / "results"
OUT.mkdir(exist_ok=True)
EXPR = "work/expr_5000_pc5.tsv"
GOLD = pd.read_csv(os.environ.get("GCEN2_STRING_GOLD", "work/string_gold.tsv"), sep="\t")

genes = pd.read_csv(EXPR, sep="\t", index_col=0, usecols=[0]).index.to_numpy()
idx = {g: i for i, g in enumerate(genes)}
n = len(genes)
a, b = GOLD["a"].map(idx), GOLD["b"].map(idx)
ok = a.notna() & b.notna()
gold = set(zip(np.minimum(a[ok], b[ok]).astype(int), np.maximum(a[ok], b[ok]).astype(int)))
base_rate = len(gold) / (n * (n - 1) / 2)
tmp = Path(tempfile.mkdtemp(prefix="gcen2_sens_"))


def run(*args):
    subprocess.run([GCEN2, *map(str, args)], check=True, capture_output=True, text=True)


def edge_stats(net):
    e = pd.read_csv(net, sep="\t", usecols=[0, 1])
    i, j = e.iloc[:, 0].map(idx).to_numpy(), e.iloc[:, 1].map(idx).to_numpy()
    pairs = zip(np.minimum(i, j), np.maximum(i, j))
    hits = sum(p in gold for p in pairs)
    return {"edges": len(e), "genes": len(set(i) | set(j)),
            "edge_STRING_precision": round(hits / max(len(e), 1), 4)}


def labels(path):
    m = pd.read_csv(path, sep="\t").drop_duplicates("gene")
    lab = np.full(n, "", dtype=object)
    lab[m["gene"].map(idx).to_numpy()] = m["module"].astype(str).to_numpy()
    return lab


def module_stats(lab):
    sizes = pd.Series(lab[lab != ""]).value_counts()
    within = hits = 0
    for mod in sizes.index:
        members = np.flatnonzero(lab == mod)
        k = len(members)
        within += k * (k - 1) // 2
        s = set(members.tolist())
        hits += sum(1 for p, q in gold if p in s and q in s and lab[q] == mod)
    prec = hits / within if within else float("nan")
    return {"modules": len(sizes), "module_genes": int((lab != "").sum()),
            "largest": int(sizes.iloc[0]) if len(sizes) else 0,
            "within_STRING_precision": round(prec, 4),
            "within_STRING_enrichment": round(prec / base_rate, 1)}


# 1. correlation cutoff
rows = []
for c in (0.3, 0.4, 0.5, 0.6, 0.7, 0.8):
    net, mod = tmp / f"net_{c}.tsv", tmp / f"mod_{c}.tsv"
    run("build", "-i", EXPR, "-o", net, "-m", "pearson", "-c", c, "-p", 0.001, "-f", "-s",
        "-t", 8)
    run("module", "-i", net, "-o", mod, "-s", 0.5, "-t", 8)
    rows.append({"cutoff_r": c, **edge_stats(net), **module_stats(labels(mod))})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "sensitivity_network.tsv", sep="\t", index=False)

# 2. GCEN module similarity (r >= 0.4 network)
rows = []
for s in (0.3, 0.4, 0.5, 0.6, 0.7):
    mod = tmp / f"gmod_{s}.tsv"
    run("module", "-i", tmp / "net_0.4.tsv", "-o", mod, "-s", s, "-t", 8)
    rows.append({"similarity": s, **module_stats(labels(mod))})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "sensitivity_gcen_modules.tsv", sep="\t", index=False)

# 3. WGCNA-style modules, one parameter at a time
default = {"power": "auto", "deep-split": 2, "min-size": 30, "merge-cut": 0.25}
grid = [("power", v) for v in (6, 8, 10, 12, 14, 16)] + \
       [("deep-split", v) for v in (0, 1, 3, 4)] + \
       [("min-size", v) for v in (20, 50, 100)] + \
       [("merge-cut", v) for v in (0.15, 0.35)]


def wgcna(params, name):
    out = tmp / f"w_{name}.tsv"
    args = ["module", "-e", EXPR, "-o", out, "--signed", "-t", 8]
    for k, v in params.items():
        args += [f"--{k}", v]
    run(*args)
    return labels(out)


ref = wgcna(default, "default")
rows = [{"parameter": "default", "value": "", **module_stats(ref), "ARI_vs_default": 1.0,
         "ARI_vs_default_assigned": 1.0, "genes_assigned_both": int((ref != "").sum())}]
for k, v in grid:
    lab = wgcna({**default, k: v}, f"{k}_{v}")
    ari = adjusted_rand_score(np.where(ref == "", "-", ref), np.where(lab == "", "-", lab))
    both = (ref != "") & (lab != "")
    rows.append({"parameter": k, "value": v, **module_stats(lab),
                 "ARI_vs_default": round(ari, 3),
                 "ARI_vs_default_assigned": round(adjusted_rand_score(ref[both], lab[both]), 3),
                 "genes_assigned_both": int(both.sum())})
    print(rows[-1], flush=True)
pd.DataFrame(rows).to_csv(OUT / "sensitivity_wgcna_modules.tsv", sep="\t", index=False)
shutil.rmtree(tmp)
