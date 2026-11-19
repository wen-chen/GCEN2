"""GO annotation and donor traits for the GTEx whole-blood case study.

Run from the GTEx data directory after ``../wgcna_gtex/01_prepare.py`` (which
writes work/expr_11944_pc5.tsv). Inputs in data/:

  GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt   donor phenotypes
  go/go-basic.obo, go/HUMAN-uniprot.gaf.gz                GO release 2026-08-05

Outputs in work/case/:

  go_assoc.tsv   gene symbol -> GO terms (GCEN v1 association format). NOT
                 qualifiers and annotations inferred from expression (IEP, HEP)
                 or without data (ND) are left out, so that co-expression is
                 not scored against annotations derived from co-expression.
  traits.tsv     samples x numeric donor traits: female (0/1), age (midpoint
                 of the 10-year bracket), hardy (Hardy death scale 0-4) and
                 ventilator (Hardy 0: died on a ventilator, 0/1)
  genes.txt      the 11,944 genes (background of the enrichment)
"""

import gzip
from pathlib import Path

import pandas as pd

WORK = Path("work")
OUT = WORK / "case"
OUT.mkdir(parents=True, exist_ok=True)
SKIP_EVIDENCE = {"IEP", "HEP", "ND"}

terms: dict[str, set[str]] = {}
with gzip.open("data/go/HUMAN-uniprot.gaf.gz", "rt") as fh:
    for line in fh:
        if line.startswith("!"):
            continue
        f = line.rstrip("\n").split("\t")
        if "NOT" in f[3].split("|") or f[6] in SKIP_EVIDENCE:
            continue
        terms.setdefault(f[2], set()).add(f[4])
with open(OUT / "go_assoc.tsv", "w") as fh:
    for gene in sorted(terms):
        fh.write(gene + "\t" + "\t".join(sorted(terms[gene])) + "\n")

expr_genes = pd.read_csv(WORK / "expr_11944_pc5.tsv", sep="\t", index_col=0,
                         usecols=[0]).index
(OUT / "genes.txt").write_text("\n".join(expr_genes) + "\n")
print(f"{len(terms):,} annotated genes; {sum(g in terms for g in expr_genes):,} "
      f"of the {len(expr_genes):,} expressed genes")

samples = pd.read_csv(WORK / "expr_11944_pc5.tsv", sep="\t", index_col=0, nrows=0).columns
pheno = pd.read_csv("data/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt", sep="\t",
                    index_col="SUBJID")
p = pheno.reindex(samples.str.split("-").str[:2].str.join("-"))
age = p["AGE"].str.split("-", expand=True).astype(float).mean(axis=1) + 0.5
traits = pd.DataFrame({"female": (p["SEX"] == 2).astype(int).to_numpy(),
                       "age": age.to_numpy(),
                       "hardy": p["DTHHRDY"].to_numpy(),
                       "ventilator": (p["DTHHRDY"] == 0).astype(float).to_numpy()},
                      index=samples)
traits.loc[p["DTHHRDY"].isna().to_numpy(), "ventilator"] = float("nan")
traits.index.name = "sample"
traits.to_csv(OUT / "traits.tsv", sep="\t")
print(traits.describe().loc[["count", "mean"]].round(2).to_string())
