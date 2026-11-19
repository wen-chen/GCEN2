"""gcen2 trajectory vs the R implementation, window by window and edge by
edge; writes results/compare.tsv. Run from the data directory after
run_all.sh."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUNS = Path("runs/trajectory")


def pairs(df):
    a, b = df["gene_a"].astype(str), df["gene_b"].astype(str)
    swap = a > b
    return pd.Index(a.where(~swap, b) + "\t" + b.where(~swap, a))


def main() -> int:
    res = {}
    k = 0
    same_edges, max_dr = True, 0.0
    while (RUNS / f"gcen2/window_{k}.network").exists():
        g = pd.read_csv(RUNS / f"gcen2/window_{k}.network", sep="\t")
        r = pd.read_csv(RUNS / f"r/window_{k}.network", sep="\t")
        g.index, r.index = pairs(g), pairs(r)
        same_edges &= g.index.sort_values().equals(r.index.sort_values())
        common = g.index.intersection(r.index)
        if len(common):
            max_dr = max(max_dr, float(np.abs(g.loc[common, "correlation"]
                                              - r.loc[common, "correlation"]).max()))
        res[f"edges_window_{k}"] = len(g)
        k += 1
    res["windows"] = k
    res["same_edges_every_window"] = bool(same_edges)
    res["window_max_abs_dr"] = max_dr

    g = pd.read_csv(RUNS / "gcen2/edge_dynamics.tsv", sep="\t")
    r = pd.read_csv(RUNS / "r/edge_dynamics.tsv", sep="\t")
    g.index, r.index = pairs(g), pairs(r)
    common = g.index.intersection(r.index)
    g, r = g.loc[common], r.loc[common]
    cor = [f"cor_{i}" for i in range(k)]
    res.update({
        "dynamics_edges_gcen2": len(g.index), "dynamics_edges_r": len(r.index),
        "dynamics_edges_shared": len(common),
        "dynamics_max_abs_dr": float(np.abs(g[cor].to_numpy() - r[cor].to_numpy()).max()),
        "same_pattern_first_last": bool((g[["pattern", "n_windows", "first", "last"]]
                                         == r[["pattern", "n_windows", "first", "last"]])
                                        .all().all()),
        "max_abs_dtrend": float(np.abs(g["trend"] - r["trend"]).max()),
    })
    for p, n in g["pattern"].value_counts().items():
        res[f"pattern_{p}"] = int(n)
    runs = {}  # mean over the runs
    for line in (RUNS / "time.txt").read_text().splitlines():
        name, s, _, kb, _ = line.split()
        runs.setdefault(name, []).append((float(s), float(kb) / 1024**2))
    for name, v in runs.items():
        res[f"{name}_s"] = round(sum(x[0] for x in v) / len(v), 2)
        res[f"{name}_GB"] = round(sum(x[1] for x in v) / len(v), 3)
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    lines = ["metric\tvalue"] + [f"{m}\t{v:.4g}" if isinstance(v, float) else f"{m}\t{v}"
                                 for m, v in res.items()]
    (out / "compare.tsv").write_text("\n".join(lines) + "\n")
    (out / "time.txt").write_text((RUNS / "time.txt").read_text())
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
