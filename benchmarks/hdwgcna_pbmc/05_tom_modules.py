"""GCEN2's WGCNA-style modules (gcen2 module -e --preset hdwgcna) against
hdWGCNA's ConstructNetwork modules.

For every run of run_all.sh:
  * on hdWGCNA's own metacell matrix (runs/hdwgcna/<run>/datExpr.tsv):
    `gcen2 module -e --preset hdwgcna -t 8`, timed with /usr/bin/time, its
    soft power (TestSoftPowers' rule: bicor, R^2 >= 0.8, power > 3) and the
    adjusted Rand index of its modules with hdWGCNA's (over all genes,
    unassigned genes as one group);
  * on GCEN2's metacells (runs/gcen2/<run>/metacells.tsv, same settings,
    sampled independently): the same command, i.e. the whole hdWGCNA module
    workflow without R, and its ARI with hdWGCNA's modules; for scale, the
    ARI of hdWGCNA's modules between metacell seeds is in seed_stability.tsv.
hdWGCNA's time is t_network of timing.tsv (SetDatExpr, TestSoftPowers,
ConstructNetwork). Writes results/tom_modules.tsv.
"""

import subprocess
import sys
from pathlib import Path

import pandas as pd
from sklearn.metrics import adjusted_rand_score

ROOT = Path.cwd()  # run from the PBMC data directory
RUNS = ROOT / "runs"
WORK = ROOT / "work" / "tom_modules"
WORK.mkdir(parents=True, exist_ok=True)
GCEN2 = str(Path(sys.executable).with_name("gcen2"))


def ari(a: pd.Series, b: pd.Series) -> float:
    a = a.astype(str)
    b = b.reindex(a.index).fillna("unassigned").astype(str)
    return adjusted_rand_score(a, b)


def gcen2_modules(src: Path, name: str):
    """(modules by gene, power, seconds, GB) of gcen2 module -e on a
    metacells x genes table."""
    expr = WORK / f"{name}.tsv"
    if not expr.exists():
        pd.read_csv(src, sep="\t", index_col=0).T.to_csv(expr, sep="\t", float_format="%.17g")
    out = WORK / f"{name}.modules.tsv"
    p = subprocess.run(["/usr/bin/time", "-f", "%e %M", GCEN2, "module", "-e", str(expr), "-o",
                        str(out), "--preset", "hdwgcna", "-t", "8"],
                       capture_output=True, text=True, check=True)
    sec, kb = p.stderr.strip().splitlines()[-1].split()
    power = next(x for x in p.stdout.splitlines() if x.startswith("TOM power")).split()[2]
    mods = pd.read_csv(out, sep="\t").set_index("gene")["module"]
    return mods, float(power), float(sec), int(kb) / 2**20


rows = []
for run in sorted(p.name for p in (RUNS / "hdwgcna").iterdir() if p.is_dir()):
    h = RUNS / "hdwgcna" / run
    hdw = pd.read_csv(h / "modules.tsv", sep="\t").set_index("gene")["module"]
    hdw = hdw.where(hdw != "grey", "unassigned")
    t = pd.read_csv(h / "timing.tsv", sep="\t").iloc[0]
    mods, power, sec, gb = gcen2_modules(h / "datExpr.tsv", f"{run}_hdw")
    row = {"run": run, "metacells": int(t.metacells), "genes": int(t.genes),
           "hdwgcna_power": t.power, "gcen2_power": power,
           "hdwgcna_modules": int(t.modules), "gcen2_modules": mods.nunique(),
           "ARI": round(ari(hdw, mods), 6), "gcen2_s": sec, "gcen2_GB": round(gb, 2),
           "hdwgcna_network_s": t.t_network}
    own = RUNS / "gcen2" / run / "metacells.tsv"
    if own.exists():
        m2, p2, s2, _ = gcen2_modules(own, f"{run}_own")
        row.update({"own_power": p2, "own_modules": m2.nunique(),
                    "own_ARI": round(ari(hdw, m2), 3), "own_s": s2})
    rows.append(row)
    print(row, flush=True)

df = pd.DataFrame(rows)
df.to_csv(ROOT / "results/tom_modules.tsv", sep="\t", index=False)
print(df.to_string(index=False))
