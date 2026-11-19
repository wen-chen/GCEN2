"""GTEx v10 whole blood (803 samples): log2(TPM + 1) of the expressed genes,
as genes x samples TSVs shared by GCEN2 and WGCNA.

Run from the GTEx data directory (data/ holds the GTEx file); outputs in work/:
  expr_<N>.tsv     the N most variable genes (N = 5000, 10000 and all 11,944), gene
                   symbols as row names, samples as columns
  half_A.txt, half_B.txt   a random half of the samples each
"""

from pathlib import Path

import numpy as np
import pandas as pd

DATA = Path("data/gene_tpm_v10_whole_blood.gct.gz")  # relative to the data directory
WORK = Path("work")
WORK.mkdir(parents=True, exist_ok=True)

tpm = pd.read_csv(DATA, sep="\t", skiprows=2)
genes = tpm["Description"].astype(str)
x = tpm.drop(columns=["Name", "Description"]).to_numpy(np.float64)
print(f"{x.shape[0]} genes x {x.shape[1]} samples")

expressed = np.median(x, axis=1) >= 1.0
x, genes = x[expressed], genes[expressed].to_numpy()
# one row per symbol: the most expressed
order = np.argsort(-x.mean(axis=1), kind="stable")
_, first = np.unique(genes[order], return_index=True)
keep = np.sort(order[first])
x, genes = x[keep], genes[keep]
logx = np.log2(x + 1.0)
print(f"expressed (median TPM >= 1), unique symbols: {len(genes)}")

var_order = np.argsort(-logx.var(axis=1), kind="stable")
samples = tpm.columns[2:]
for n in (5000, 10000, len(genes)):  # len(genes): all expressed genes
    idx = np.sort(var_order[:n])
    pd.DataFrame(logx[idx], index=genes[idx], columns=samples).to_csv(
        WORK / f"expr_{n}.tsv", sep="\t", float_format="%.6g")
    print(f"expr_{n}.tsv written")

rng = np.random.default_rng(0)
perm = rng.permutation(len(samples))
half = len(samples) // 2
(WORK / "half_A.txt").write_text("\n".join(samples[np.sort(perm[:half])]) + "\n")
(WORK / "half_B.txt").write_text("\n".join(samples[np.sort(perm[half:])]) + "\n")


# Confounders (Parsana et al., Genome Biol 2019): the top expression PCs of
# GTEx carry technical variation that correlates most genes. Regress out the
# top 5 PCs (5 gave the best STRING precision among 1, 3, 5, 10, 20 PCs on
# the 5,000 genes, so absolute precisions on corrected data are optimistic;
# both tools get the same input). Each split half is corrected on its own.
def remove_pcs(df, k=5):
    x = df.to_numpy(np.float64)
    mean = x.mean(axis=1, keepdims=True)
    xc = x - mean
    u, s, vt = np.linalg.svd(xc, full_matrices=False)
    return pd.DataFrame(xc - (u[:, :k] * s[:k]) @ vt[:k] + mean, index=df.index,
                        columns=df.columns)


halves = {h: (WORK / f"half_{h}.txt").read_text().split() for h in "AB"}
for n in (5000, 10000, len(genes)):
    df = pd.read_csv(WORK / f"expr_{n}.tsv", sep="\t", index_col=0)
    remove_pcs(df).to_csv(WORK / f"expr_{n}_pc5.tsv", sep="\t", float_format="%.6g")
    if n == 5000:
        for h, cols in halves.items():
            df[cols].to_csv(WORK / f"expr_{n}_{h}.tsv", sep="\t", float_format="%.6g")
            remove_pcs(df[cols]).to_csv(WORK / f"expr_{n}_pc5_{h}.tsv", sep="\t",
                                        float_format="%.6g")
    print(f"expr_{n}_pc5.tsv written")
