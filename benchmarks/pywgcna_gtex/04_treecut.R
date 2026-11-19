# WGCNA's gene tree and dynamic tree cut on the corrected 5,000 genes, for
# 04_treecut.py (which runs PyWGCNA's on the same data and compares).
# Signed adjacency (power 10, the power both choose), signed TOM, average
# linkage; cutreeDynamic (hybrid, deepSplit 2, minClusterSize 30) with the
# PAM stage as PyWGCNA's findModules (pamRespectsDendro FALSE) and without
# it (pamStage FALSE: the tree cut proper).
#
# usage: Rscript 04_treecut.R   (from the GTEx data directory)
# Writes runs/treecut/r_heights.txt and r_labels.tsv.

suppressPackageStartupMessages(library(WGCNA))
enableWGCNAThreads(nThreads = 8)
dir.create("runs/treecut", recursive = TRUE, showWarnings = FALSE)
x <- read.delim("work/expr_5000_pc5.tsv", row.names = 1, check.names = FALSE)
datExpr <- t(as.matrix(x))
adj <- adjacency(datExpr, power = 10, type = "signed")
d <- 1 - TOMsimilarity(adj, TOMType = "signed", verbose = 0)
tree <- hclust(as.dist(d), method = "average")
writeLines(sprintf("%.17g", sort(tree$height)), "runs/treecut/r_heights.txt")
pam <- cutreeDynamic(dendro = tree, distM = d, deepSplit = 2, pamRespectsDendro = FALSE,
                     minClusterSize = 30, verbose = 0)
nopam <- cutreeDynamic(dendro = tree, distM = d, deepSplit = 2, pamStage = FALSE,
                       minClusterSize = 30, verbose = 0)
write.table(data.frame(gene = colnames(datExpr), pam = pam, nopam = nopam),
            "runs/treecut/r_labels.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
