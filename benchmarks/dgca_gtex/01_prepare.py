"""GTEx v10 whole blood split by sex, for differential co-expression.

Run from the GTEx data directory after ``wgcna_gtex/01_prepare.py`` (which
writes work/expr_5000.tsv, log2(TPM + 1) of the 5,000 most variable genes).
data/ must also hold the GTEx subject phenotypes
(GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt). Outputs in work/:

  dc_<N>_male.tsv, dc_<N>_female.tsv   the N most variable genes (N = 2000 and
                                       5000), genes x samples, one file per sex
"""

from pathlib import Path

import numpy as np
import pandas as pd

WORK = Path("work")
pheno = pd.read_csv("data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", sep="\t",
                    index_col="SUBJID")
expr = pd.read_csv(WORK / "expr_5000.tsv", sep="\t", index_col=0)
subject = expr.columns.str.split("-").str[:2].str.join("-")
sex = pheno["SEX"].reindex(subject).map({1: "male", 2: "female"}).to_numpy()
print({s: int((sex == s).sum()) for s in ("male", "female")})

order = np.argsort(-expr.var(axis=1).to_numpy(), kind="stable")
for n in (2000, 5000):
    sub = expr.iloc[np.sort(order[:n])]
    for s in ("male", "female"):
        sub.loc[:, sex == s].to_csv(WORK / f"dc_{n}_{s}.tsv", sep="\t")
    print(f"dc_{n}_*.tsv written")
