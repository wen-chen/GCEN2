"""GCEN2 vs RandomWalkRestartMH (RWR-M) and vs the R fusion; writes
results/compare.tsv. Run from the data directory after run_all.sh.

RandomWalkRestartMH strips "_1" from node names when it reports them
(``gsub("_1", "", ...)``), e.g. hCG_1984214 becomes hCG984214; those genes
are matched after the same substitution.
"""

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUNS = Path("runs/multilayer")


def times() -> dict:
    """step: (seconds, peak GB), the mean over the runs."""
    runs = {}
    for line in (RUNS / "time.txt").read_text().splitlines():
        name, s, _, kb, _ = line.split()
        runs.setdefault(name, []).append((float(s), float(kb) / 1024**2))
    return {k: tuple(sum(x) / len(v) for x in zip(*v)) for k, v in runs.items()}


def rwr() -> dict:
    g = pd.read_csv(RUNS / "gcen2_rwrm.tsv", sep="\t")
    r = pd.read_csv(RUNS / "rwrmh.tsv", sep="\t")
    g["key"] = [re.sub("_1", "", x, count=1) for x in g["gene"]]
    seeds = set(Path("work/seeds.txt").read_text().split())
    r = r[~r["NodeNames"].isin(seeds)].set_index("NodeNames")
    m = g.join(r, on="key", how="inner")
    top = r.sort_values("Score", ascending=False, kind="stable").index
    return {
        "rwr_genes_gcen2": len(g), "rwr_genes_rwrmh": len(r), "rwr_genes_matched": len(m),
        "rwr_genes_renamed_by_rwrmh": int((g["key"] != g["gene"]).sum()),
        "rwr_max_abs_diff": float(np.abs(m["probability"] - m["Score"]).max()),
        "rwr_top100_shared": len(set(g["key"][:100]) & set(top[:100])),
    }


def fusion() -> dict:
    g = pd.read_csv(RUNS / "gcen2_fused.tsv", sep="\t", usecols=[0, 1, 2, 3])
    r = pd.read_csv(RUNS / "r_fused.tsv", sep="\t")
    # R's pmin/pmax order the pair by the locale's collation, GCEN2 by bytes
    swap = r["gene_a"] > r["gene_b"]
    r.loc[swap, ["gene_a", "gene_b"]] = r.loc[swap, ["gene_b", "gene_a"]].to_numpy()
    m = g.merge(r, on=["gene_a", "gene_b"], how="outer", indicator=True)
    both = m["_merge"] == "both"
    return {
        "fuse_edges_gcen2": len(g), "fuse_edges_r": len(r), "fuse_edges_shared": int(both.sum()),
        "fuse_max_abs_score_diff": float(np.abs(m.loc[both, "score_x"]
                                                - m.loc[both, "score_y"]).max()),
        "fuse_same_n_layers": bool((m.loc[both, "n_layers_x"]
                                    == m.loc[both, "n_layers_y"]).all()),
    }


def main() -> int:
    t = times()
    res = {**rwr(), **fusion()}
    for key, name in (("rwr_gcen2", "gcen2_rwrm"), ("rwr_rwrmh", "rwrmh"),
                      ("fuse_gcen2", "gcen2_fuse"), ("fuse_r", "r_fuse")):
        res[f"{key}_s"], res[f"{key}_GB"] = round(t[name][0], 2), round(t[name][1], 3)
    genes = set()
    for name in ("coexpression", "ppi", "regulation"):
        layer = pd.read_csv(f"work/{name}.tsv", sep="\t", usecols=[0, 1])
        res[f"edges_{name}"] = len(layer)
        genes |= set(layer["gene_a"]) | set(layer["gene_b"])
    res["genes"] = len(genes)
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    lines = ["metric\tvalue"] + [f"{k}\t{v:.4g}" if isinstance(v, float) else f"{k}\t{v}"
                                 for k, v in res.items()]
    (out / "compare.tsv").write_text("\n".join(lines) + "\n")
    (out / "time.txt").write_text((RUNS / "time.txt").read_text())
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
