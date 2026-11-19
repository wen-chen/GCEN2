"""GCEN2's WGCNA-style modules (gcen2 module -e) against blockwiseModules.

For the six runs of run_all.sh: `gcen2 module -e` with WGCNA's settings
there (signed, power auto, minModuleSize 30, mergeCutHeight 0.25, deepSplit
2), timed with /usr/bin/time, and the adjusted Rand index of its modules
with WGCNA's (runs/wgcna/<run>/modules.tsv; unassigned genes as one group,
over all genes). Writes results/tom_modules.tsv.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path.cwd()  # run from the GTEx data directory
GCEN2 = sys.argv[1] if len(sys.argv) > 1 else "gcen2"
RUNS = {"pc5_5000": "expr_5000_pc5.tsv", "pc5_5000_A": "expr_5000_pc5_A.tsv",
        "pc5_5000_B": "expr_5000_pc5_B.tsv", "raw_5000": "expr_5000.tsv",
        "pc5_10000": "expr_10000_pc5.tsv", "pc5_11944": "expr_11944_pc5.tsv"}

times = pd.read_csv(ROOT / "runs/time.txt", sep=" ", header=None,
                    names=["step", "run", "sec", "s", "kb", "KB"], comment="C")
rows = []
for run, expr in RUNS.items():
    out = ROOT / "runs/gcen2" / run / "tom_modules.tsv"
    p = subprocess.run(["/usr/bin/time", "-f", "%e %M", GCEN2, "module", "-e",
                        str(ROOT / "work" / expr), "-o", str(out), "--signed",
                        "--min-size", "30", "--merge-cut", "0.25", "-t", "8"],
                       capture_output=True, text=True, check=True)
    sec, kb = p.stderr.split()[-2:]
    power = next(line for line in p.stdout.splitlines() if line.startswith("TOM power"))
    mine = pd.read_csv(out, sep="\t").set_index("gene")["module"]
    w = pd.read_csv(ROOT / f"runs/wgcna/{run}/modules.tsv", sep="\t").set_index("gene")["module"]
    a = mine.reindex(w.index).fillna(-1).astype(int).astype(str)
    wt = times[(times.step == "wgcna") & (times.run == run)].iloc[0]
    rows.append({"run": run, "power": power.split()[2],
                 "gcen2_modules": mine.nunique(), "wgcna_modules": w[w != "grey"].nunique(),
                 "gcen2_genes": len(mine), "wgcna_genes": int((w != "grey").sum()),
                 "ARI": adjusted_rand_score(a, w), "gcen2_s": float(sec),
                 "gcen2_GB": int(kb) / 1048576, "wgcna_s": wt.sec, "wgcna_GB": wt.kb / 1048576})
    print(rows[-1], flush=True)
df = pd.DataFrame(rows)
(ROOT / "results").mkdir(exist_ok=True)
df.to_csv(ROOT / "results/tom_modules.tsv", sep="\t", index=False)
print(df.round(3).to_string(index=False))
