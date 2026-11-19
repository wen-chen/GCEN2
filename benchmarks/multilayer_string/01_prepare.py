"""Three gene network layers for the multi-layer benchmark.

Usage (from a data directory of your choice; outputs in work/):
  python 01_prepare.py <GTEx expr_11944_pc5.tsv> <STRING directory> <TRRUST file>

  coexpression.tsv  GTEx whole blood co-expression (``gcen2 build -m pearson
                    -c 0.5 -p 0.001 -f -s`` on the 11,944 genes with 5 PCs
                    removed, from wgcna_gtex/01_prepare.py); weight = r
  ppi.tsv           STRING v12 human, combined score >= 700, gene symbols;
                    weight = score / 1000
  regulation.tsv    TRRUST v2 human transcription factor -> target pairs;
                    weight = 1
  seeds.txt         interferon-response genes, the seeds of the random walks

Every layer is undirected, without self-loops or repeated pairs, as
``gene_a, gene_b, weight`` with a header.
"""

import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

expr, string_dir, trrust = sys.argv[1], Path(sys.argv[2]), sys.argv[3]
WORK = Path("work")
WORK.mkdir(exist_ok=True)


def write(name, a, b, w):
    a, b = np.asarray(a, dtype=object), np.asarray(b, dtype=object)
    lo, hi = np.where(a < b, a, b), np.where(a < b, b, a)
    df = pd.DataFrame({"gene_a": lo, "gene_b": hi, "weight": w})
    df = df[df.gene_a != df.gene_b].sort_values("weight", ascending=False)
    df = df.drop_duplicates(["gene_a", "gene_b"]).sort_values(["gene_a", "gene_b"])
    df.to_csv(WORK / name, sep="\t", index=False)
    genes = len(set(df.gene_a) | set(df.gene_b))
    print(f"{name}: {len(df):,} edges, {genes:,} genes")


net = WORK / "coexpression.network"
subprocess.run(["gcen2", "build", "-i", expr, "-o", str(net), "-m", "pearson", "-c", "0.5",
                "-p", "0.001", "-f", "-s", "-t", "8"], check=True)
co = pd.read_csv(net, sep="\t")
write("coexpression.tsv", co.gene_a, co.gene_b, co.correlation.to_numpy())

info = pd.read_csv(string_dir / "9606.protein.info.v12.0.txt.gz", sep="\t")
symbol = dict(zip(info["#string_protein_id"], info["preferred_name"]))
links = pd.read_csv(string_dir / "9606.protein.links.v12.0.txt.gz", sep=" ")
links = links[links.combined_score >= 700]
write("ppi.tsv", links.protein1.map(symbol), links.protein2.map(symbol),
      links.combined_score.to_numpy() / 1000)

reg = pd.read_csv(trrust, sep="\t", header=None, names=["tf", "target", "mode", "pmid"])
write("regulation.tsv", reg.tf, reg.target, np.ones(len(reg)))

(WORK / "seeds.txt").write_text("\n".join(
    ["IFIT1", "IFIT3", "ISG15", "MX1", "OAS1", "RSAD2", "IFI44L", "STAT1", "IRF7"]) + "\n")
