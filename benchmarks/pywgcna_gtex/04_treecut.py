"""PyWGCNA's gene tree and dynamic tree cut against WGCNA's (04_treecut.R).

Run with the PyWGCNA environment from the GTEx data directory, after
04_treecut.R. The same steps as there with PyWGCNA's functions (adjacency,
TOMsimilarity, scipy average linkage as in findModules, cutreeHybrid with and
without the PAM stage). cutreeHybrid returns two label columns: "Name", the
branch labels (0 = unassigned), and "Value", the labels renumbered by module
size, which findModules passes on to labels2colors.

Writes results/treecut.tsv (next to this script): tree heights; per stage the
unassigned genes and clusters of each, the ARI over all genes (unassigned as
one group) and over the genes both assign; and what the renumbering does to
the unassigned genes.
"""

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import PyWGCNA  # noqa: E402
from scipy.cluster.hierarchy import linkage  # noqa: E402
from scipy.spatial.distance import squareform  # noqa: E402
from sklearn.metrics import adjusted_rand_score as ari  # noqa: E402

W = PyWGCNA.WGCNA
x = pd.read_csv("work/expr_5000_pc5.tsv", sep="\t", index_col=0).T
adj = W.adjacency(x, power=10, adjacencyType="signed")
d = 1 - W.TOMsimilarity(np.array(adj), TOMType="signed")
d.index = d.columns = x.columns
tree = linkage(squareform(d.round(8).values, checks=False), method="average")
r_h = np.loadtxt("runs/treecut/r_heights.txt")
r = pd.read_csv("runs/treecut/r_labels.tsv", sep="\t", index_col=0).loc[x.columns]

rows = [{"stage": "tree", "measure": "largest difference of the merge heights",
         "value": float(f"{np.abs(np.sort(tree[:, 2]) - r_h).max():.2e}")}]
for stage, col, kw in (("tree cut (pamStage FALSE)", "nopam", {"pamStage": False}),
                       ("with PAM stage (pamRespectsDendro FALSE)", "pam",
                        {"pamRespectsDendro": False})):
    lab = W.cutreeHybrid(dendro=tree, distM=d.round(8), minClusterSize=30, deepSplit=2, **kw)
    ref = r[col].to_numpy()
    name, value = lab["Name"].to_numpy(), lab["Value"].to_numpy()
    both = (ref != 0) & (name != 0)
    zero_value = value[name == 0]
    rows += [
        {"stage": stage, "measure": "unassigned genes, WGCNA", "value": int((ref == 0).sum())},
        {"stage": stage, "measure": "unassigned genes, PyWGCNA", "value": int((name == 0).sum())},
        {"stage": stage, "measure": "clusters, WGCNA", "value": int(len(set(ref) - {0}))},
        {"stage": stage, "measure": "clusters, PyWGCNA", "value": int(len(set(name) - {0}))},
        {"stage": stage, "measure": "ARI, all genes", "value": round(ari(ref, name), 3)},
        {"stage": stage, "measure": "genes assigned by both", "value": int(both.sum())},
        {"stage": stage, "measure": "ARI, genes assigned by both",
         "value": round(ari(ref[both], name[both]), 3)},
        {"stage": stage, "measure": "renumbered: label of the unassigned genes",
         "value": ",".join(map(str, sorted(set(zero_value)))) if len(zero_value) else "-"},
        {"stage": stage, "measure": "renumbered: genes labelled 0 (unassigned)",
         "value": int((value == 0).sum())},
        {"stage": stage, "measure": "renumbered: branch labels of the genes labelled 0",
         "value": ",".join(map(str, sorted(set(name[value == 0]))))},
    ]
out = Path(__file__).resolve().parent / "results"
out.mkdir(exist_ok=True)
df = pd.DataFrame(rows)
df.to_csv(out / "treecut.tsv", sep="\t", index=False)
print(df.to_string(index=False))
