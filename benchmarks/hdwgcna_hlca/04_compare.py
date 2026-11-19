"""HLCA: GCEN2 celltype-modules against hdWGCNA's per-type workflow, from
run_all.sh.

Per size and cell type (the 10 largest):
  * genes: does GCEN2 select the same genes as hdWGCNA's
    SetupForWGCNA(gene_select = "fraction")?
  * network step: `gcen2 module -e --preset hdwgcna` on hdWGCNA's own
    metacell matrix (datExpr.tsv) against hdWGCNA's modules (ARI);
  * whole workflow: GCEN2's modules (its own metacells) against hdWGCNA's
    (ARI; both sample metacell centres at random), and, for scale, GCEN2's
    modules between metacell seeds 42, 1 and 2 (runs/gcen2/<size>_top10_seed<s>,
    when present) and those seeds against hdWGCNA; with hdWGCNA's seed runs
    (runs/hdwgcna/<size>_seed<s>, 02_hdwgcna.R ... <s>), hdWGCNA's modules
    between its seeds 42, 1 and 2 and GCEN2's three seeds against its seeds
    1 and 2 (with the above, all 9 pairs of seeds);
ARI over all genes, unassigned genes as one group. Time and peak memory of
the whole runs come from runs/time.txt (GCEN2: mean of the two runs).

usage: python 04_compare.py RUNS_DIR TYPES_FILE SIZE...
Writes RUNS_DIR/compare.tsv and RUNS_DIR/time.tsv.
"""

import re
import subprocess
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

runs = Path(sys.argv[1])
types = Path(sys.argv[2]).read_text().splitlines()
sizes = sys.argv[3:]
GCEN2 = str(Path(sys.executable).with_name("gcen2"))


def labels(path: Path, genes, grey="grey") -> pd.Series:
    m = pd.read_csv(path, sep="\t").set_index("gene")["module"].astype(str)
    return m.replace(grey, "unassigned").reindex(genes).fillna("unassigned")


def hd_final_of(d: Path, genes) -> pd.Series:
    """hdWGCNA's modules after ModuleConnectivity's reassignment when run."""
    f = d / "modules_final.tsv"
    return labels(f if f.exists() else d / "modules.tsv", genes)


rows = []
for size in sizes:
    ct = pd.read_csv(runs / "gcen2" / f"{size}_top10" / "celltypes.tsv", sep="\t")
    gdir = dict(zip(ct["cell_type"], ct["dir"]))
    hd_time = pd.read_csv(runs / "hdwgcna" / size / "timing.tsv", sep="\t").set_index(
        "cell_type")
    for t in types:
        hd = runs / "hdwgcna" / size / re.sub(r"[^A-Za-z0-9]+", "_", t)
        g = runs / "gcen2" / f"{size}_top10" / gdir[t]
        genes = (hd / "genes.txt").read_text().split()
        g_genes = set(pd.read_csv(g / "kme.tsv", sep="\t")["gene"])
        expr = hd / "datExpr_genes_x_metacells.tsv"
        pd.read_csv(hd / "datExpr.tsv", sep="\t", index_col=0).T.to_csv(
            expr, sep="\t", float_format="%.17g")
        p = subprocess.run(["/usr/bin/time", "-f", "%e %M", GCEN2, "module", "-e", str(expr),
                            "-o", str(hd / "gcen2_modules.tsv"), "--preset", "hdwgcna",
                            "-t", "8"], capture_output=True, text=True, check=True)
        sec, kb = p.stderr.strip().splitlines()[-1].split()
        expr.unlink()
        hd_mod = labels(hd / "modules.tsv", genes)  # ConstructNetwork
        # whole workflow: after ModuleConnectivity's reassignment when run
        # (GCEN2's modules.tsv is after its reassignment by default)
        final = hd / "modules_final.tsv"
        hd_final = labels(final, genes) if final.exists() else hd_mod
        row = ct.set_index("cell_type").loc[t]
        rows.append({
            "size": size, "cell_type": t, "cells": int(row["n_cells"]),
            "genes": len(genes), "same_genes": g_genes == set(genes),
            "metacells_gcen2": int(row["n_metacells"]),
            "metacells_hdwgcna": int(hd_time.loc[t, "metacells"]),
            "network_step_ARI": round(adjusted_rand_score(
                labels(hd / "gcen2_modules.tsv", genes), hd_mod), 3),
            "network_step_gcen2_s": float(sec),
            "network_step_hdwgcna_s": float(hd_time.loc[t, "t_network"]),
            "modules_gcen2": int(row["modules"]),
            "modules_hdwgcna": int(hd_time.loc[t, "modules"]),
            "ARI_gcen2_vs_hdwgcna": round(adjusted_rand_score(
                labels(g / "modules.tsv", genes), hd_final), 3),
        })
        seeds = [g] + [runs / "gcen2" / f"{size}_top10_seed{s}" / gdir[t] for s in (1, 2)]
        if all((d / "modules.tsv").exists() for d in seeds):
            mods = [labels(d / "modules.tsv", genes) for d in seeds]
            rows[-1]["ARI_gcen2_seeds"] = ",".join(
                str(round(adjusted_rand_score(a, b), 3))
                for i, a in enumerate(mods) for b in mods[i + 1:])
            rows[-1]["ARI_gcen2_seeds_vs_hdwgcna"] = ",".join(
                str(round(adjusted_rand_score(m, hd_final), 3)) for m in mods[1:])
            hd_seeds = [runs / "hdwgcna" / f"{size}_seed{s}" / hd.name for s in (1, 2)]
            if all((d / "modules.tsv").exists() for d in hd_seeds):
                hmods = [hd_final] + [hd_final_of(d, genes) for d in hd_seeds]
                rows[-1]["ARI_hdwgcna_seeds"] = ",".join(
                    str(round(adjusted_rand_score(a, b), 3))
                    for i, a in enumerate(hmods) for b in hmods[i + 1:])
                rows[-1]["ARI_gcen2_seeds_vs_hdwgcna_seeds"] = ",".join(
                    str(round(adjusted_rand_score(m, h), 3)) for m in mods for h in hmods[1:])
        print(rows[-1], flush=True)
df = pd.DataFrame(rows)
df.to_csv(runs / "compare.tsv", sep="\t", index=False)

t = pd.read_csv(runs / "time.txt", sep=" ", header=None,
                names=["run", "size", "rep", "sec", "s", "kb", "KB"])
tt = (t.groupby(["run", "size"])
      .agg(seconds=("sec", "mean"), peak_GB=("kb", lambda x: x.max() / 1024**2), n=("sec", "size"))
      .reset_index())
tt.to_csv(runs / "time.tsv", sep="\t", index=False, float_format="%.3f")
print(df.to_string(index=False))
print(tt.to_string(index=False))
