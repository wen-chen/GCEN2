"""GCEN2 diff vs DGCA ddcorAll, pair by pair; writes results/compare.tsv.

Run from the GTEx data directory after run_all.sh. GCEN2 writes correlations
with 6 decimals (the GCEN v1 network format), so its z and p are those of the
rounded r; ``max_abs_dz_same_r`` applies DGCA's formula to the rounded r to
show that the formulas agree and the rest is file precision. DGCA's z is
B - A, GCEN2's A - B: compared with the sign flipped. Pearson and Spearman
(variance 1.06 / (n - 3) for both tools); GCEN2's fdr against DGCA's
Benjamini-Hochberg pValDiff_adj.
"""

import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
RUNS = Path("runs/dgca")


def time_table() -> dict:
    """(step, genes): (seconds, peak GB), the mean over the runs."""
    runs = {}
    for line in (RUNS / "time.txt").read_text().splitlines():
        f = line.split()
        if len(f) == 6 and f[3] == "s":
            runs.setdefault((f[0], f[1]), []).append((float(f[2]), float(f[4]) / 1024**2))
    return {k: tuple(sum(x) / len(v) for x in zip(*v)) for k, v in runs.items()}


def pair_key(df, a, b):
    x, y = df[a].astype(str), df[b].astype(str)
    swap = x > y
    return pd.Index(x.where(~swap, y) + "\t" + y.where(~swap, x))


def compare(m: str, n: int, times: dict) -> dict:
    d = RUNS / f"{m}_{n}"
    g = pd.read_csv(d / "gcen2_diff.tsv", sep="\t")
    r = pd.read_csv(d / "dgca.tsv", sep="\t")
    g.index = pair_key(g, "gene_a", "gene_b")
    r.index = pair_key(r, "Gene1", "Gene2")
    common = g.index.intersection(r.index)
    g, r = g.loc[common], r.loc[common]
    na, nb = (pd.read_csv(f"work/dc_{n}_{s}.tsv", sep="\t", nrows=0).shape[1] - 1
              for s in ("male", "female"))

    v = 1.06 if m == "spearman" else 1.0

    def z_dgca(ra, rb):  # DGCA's dCorrs, as A - B
        ra, rb = np.minimum(ra, 0.9999), np.minimum(rb, 0.9999)
        return (np.arctanh(ra) - np.arctanh(rb)) / np.sqrt(v / (na - 3) + v / (nb - 3))

    zg, zr = g["z_score"].to_numpy(), -r["zScoreDiff"].to_numpy()
    pg, pr = g["p_value"].to_numpy(), r["pValDiff"].to_numpy()
    ok = (pg > 1e-300) & (pr > 1e-300)
    qg, qr = g["fdr"].to_numpy(), r["pValDiff_adj"].to_numpy()
    okq = (qg > 1e-300) & (qr > 1e-300)
    sig_g, sig_r = pg <= 0.05, pr <= 0.05
    gb = [times[(f"gcen2_{s}", f"{m}_{n}")] for s in ("build_A", "build_B", "diff")]
    dg = times[("dgca", f"{m}_{n}")]
    return {
        "method": m, "genes": n, "samples_A": na, "samples_B": nb,
        "pairs_gcen2": len(g.index), "pairs_dgca": len(r.index),
        "pairs_expected": n * (n - 1) // 2, "pairs_common": len(common),
        "max_abs_dr_A": np.abs(g["cor_a"] - r["A_cor"]).max(),
        "max_abs_dr_B": np.abs(g["cor_b"] - r["B_cor"]).max(),
        "max_abs_dz": np.abs(zg - zr).max(),
        "max_abs_dz_same_r": np.abs(zg - z_dgca(g["cor_a"].to_numpy(),
                                                 g["cor_b"].to_numpy())).max(),
        "max_rel_dp": (np.abs(pg[ok] - pr[ok]) / pr[ok]).max(),
        "significant_gcen2": int(sig_g.sum()), "significant_dgca": int(sig_r.sum()),
        "significant_differ": int((sig_g != sig_r).sum()),
        "max_rel_dfdr": (np.abs(qg[okq] - qr[okq]) / qr[okq]).max(),
        "fdr_significant_gcen2": int((qg <= 0.05).sum()),
        "fdr_significant_dgca": int((qr <= 0.05).sum()),
        "fdr_significant_differ": int(((qg <= 0.05) != (qr <= 0.05)).sum()),
        "gcen2_s": round(sum(t[0] for t in gb), 2),
        "gcen2_GB": round(max(t[1] for t in gb), 3),
        "dgca_s": round(dg[0], 2), "dgca_GB": round(dg[1], 3),
    }


def main() -> int:
    times = time_table()
    cols = {f"{m}_{n}": compare(m, n, times) for m in ("pearson", "spearman")
            for n in (2000, 5000)}
    lines = ["metric\t" + "\t".join(cols)]
    for m in cols["pearson_2000"]:
        lines.append("\t".join([m] + [f"{c[m]:.4g}" if isinstance(c[m], float) else str(c[m])
                                      for c in cols.values()]))
    out = HERE / "results"
    out.mkdir(exist_ok=True)
    (out / "compare.tsv").write_text("\n".join(lines) + "\n")
    (out / "time.txt").write_text((RUNS / "time.txt").read_text())
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
