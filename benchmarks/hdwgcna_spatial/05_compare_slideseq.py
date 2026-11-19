"""Slide-seqV2: GCEN2 (knn metaspots) against hdWGCNA (MetacellsByGroups on
the spatial coordinates), from run_slideseq.sh.

* genes: the same gene list?
* network step: `gcen2 module -e --preset hdwgcna` on hdWGCNA's own metacell
  matrix (datExpr.tsv) against hdWGCNA's ConstructNetwork modules (ARI);
* whole workflow: GCEN2's modules (its own metaspots) against hdWGCNA's,
  next to each tool's agreement with itself under another metacell seed
  (both sample centres at random), and the best |r| of each GCEN2 bead
  eigengene with hdWGCNA's;
* time and peak memory (median of the two timed runs).
ARI is over all genes, unassigned genes as one group.

usage: python 05_compare_slideseq.py RUNS_DIR [SAMPLE] [TIME_FILE]
(defaults slideseqv2, time_slideseq.txt; also used for Visium HD)
Writes RUNS_DIR/<SAMPLE>_compare.tsv.
"""

import subprocess
import sys
from itertools import combinations
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import adjusted_rand_score

runs = Path(sys.argv[1])
name = sys.argv[2] if len(sys.argv) > 2 else "slideseqv2"
time_file = sys.argv[3] if len(sys.argv) > 3 else "time_slideseq.txt"
GCEN2 = str(Path(sys.executable).with_name("gcen2"))
genes = pd.Index((runs / f"hdwgcna/{name}/genes.txt").read_text().split())


def gcen2_mods(d: Path) -> pd.Series:
    m = pd.read_csv(d / "modules.tsv", sep="\t").set_index("gene")["module"].astype(str)
    return m.reindex(genes).fillna("unassigned")


def hd_mods(d: Path, name="modules.tsv") -> pd.Series:
    m = pd.read_csv(d / name, sep="\t").set_index("gene")["module"].astype(str)
    return m.replace("grey", "unassigned").reindex(genes).fillna("unassigned")


def ari(a, b) -> float:
    return round(adjusted_rand_score(a, b), 3)


row = {}
g_genes = set(pd.read_csv(runs / f"gcen2/{name}/kme.tsv", sep="\t")["gene"])
row["genes"] = len(genes)
row["same_genes"] = g_genes == set(genes)

# network step on hdWGCNA's datExpr
work = runs / f"{name}_tom"
work.mkdir(exist_ok=True)
expr = work / "datExpr_genes_x_metacells.tsv"
pd.read_csv(runs / f"hdwgcna/{name}/datExpr.tsv", sep="\t", index_col=0).T.to_csv(
    expr, sep="\t", float_format="%.17g")
p = subprocess.run(["/usr/bin/time", "-f", "%e %M", GCEN2, "module", "-e", str(expr), "-o",
                    str(work / "modules.tsv"), "--preset", "hdwgcna", "-t", "8"],
                   capture_output=True, text=True, check=True)
sec, kb = p.stderr.strip().splitlines()[-1].split()
hd_net = hd_mods(runs / f"hdwgcna/{name}", "modules_network.tsv")
row["network_step_ARI"] = ari(gcen2_mods(work), hd_net)
row["network_step_gcen2_s"] = float(sec)
row["network_step_hdwgcna_s"] = float(pd.read_csv(runs / f"hdwgcna/{name}/timing.tsv",
                                                  sep="\t")["t_network"].iloc[0])

# whole workflow
g = {s: gcen2_mods(runs / f"gcen2/{name}{s}") for s in ("", "_seed1", "_seed2")}
h = {s: hd_mods(runs / f"hdwgcna/{name}{s}") for s in ("", "_seed1", "_seed2")}
row["gcen2_modules"] = int((g[""].unique() != "unassigned").sum())
row["hdwgcna_modules"] = int((h[""].unique() != "unassigned").sum())
row["ARI_gcen2_vs_hdwgcna"] = ",".join(str(ari(g[a], h[b])) for a in g for b in h)
row["ARI_gcen2_seeds"] = ",".join(str(ari(g[a], g[b])) for a, b in combinations(g, 2))
row["ARI_hdwgcna_seeds"] = ",".join(str(ari(h[a], h[b])) for a, b in combinations(h, 2))
me_g = pd.read_csv(runs / f"gcen2/{name}/spot_eigengenes.tsv", sep="\t", index_col=0)
me_h = pd.read_csv(runs / f"hdwgcna/{name}/MEs.tsv", sep="\t", index_col=0)
me_h = me_h.drop(columns=[c for c in me_h.columns if c == "grey"]).reindex(me_g.index)
r = np.corrcoef(me_g.to_numpy().T, me_h.to_numpy().T)[: me_g.shape[1], me_g.shape[1]:]
row["bead_ME_best_abs_r"] = ",".join(f"{x:.2f}" for x in np.abs(r).max(axis=1))

t = pd.read_csv(runs / time_file, sep=" ", header=None,
                names=["tool", "sample", "rep", "sec", "s", "kb", "KB"])
for tool in ("gcen2", "hdwgcna"):
    x = t[t.tool == tool]
    row[f"{tool}_s"] = round(float(x.sec.median()), 2)
    row[f"{tool}_GB"] = round(float(x.kb.median()) / 1024**2, 3)
row["speedup"] = round(row["hdwgcna_s"] / row["gcen2_s"], 1)
row["memory_ratio"] = round(row["hdwgcna_GB"] / row["gcen2_GB"], 1)
out = pd.Series(row)
out.to_frame(name).rename_axis("metric").to_csv(runs / f"{name}_compare.tsv",
                                                        sep="\t")
print(out.to_string())
