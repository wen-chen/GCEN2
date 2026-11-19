# WGCNA 1.74 / dynamicTreeCut 1.63-1 / fastcluster 1.3.0 reference for gcen2's
# TOM modules (tests/test_tom_modules.py).
# usage: Rscript make_modules_reference.R   (run in this directory)
#
# modexpr.tsv: 200 genes x 30 samples, modules of 50, 40, 30, 20 and 12 genes
# (two sub-groups in the first) plus 48 noise genes (seeded). For power 6,
# signed and unsigned (TOMType = networkType, Pearson):
#   diss_<type>.tsv.gz   1 - TOMsimilarityFromExpr (full matrix)
#   hclust_<type>.tsv    fastcluster::hclust(as.dist(diss), "average"): merge, height
#   cutree_<type>.tsv    cutreeDynamic(hybrid, deepSplit 2, minClusterSize 10,
#                        cutHeight 0.995) labels
#   modules_<type>.tsv   blockwiseModules(power 6, minModuleSize 10,
#                        mergeCutHeight 0.25, numericLabels): unmerged, merged
suppressPackageStartupMessages(library(WGCNA))
set.seed(11)
n_s <- 30
sizes <- c(50, 40, 30, 20, 12)
base <- matrix(rnorm(length(sizes) * n_s), length(sizes))
sub <- matrix(rnorm(2 * n_s), 2)
rows <- list()
for (m in seq_along(sizes)) for (g in 1:sizes[m]) {
  v <- base[m, ] * runif(1, 0.5, 1.5)
  if (m == 1) v <- v + 0.8 * sub[1 + (g > 25), ]
  if (m == 3 && g %% 2 == 0) v <- -v  # anti-correlated half
  rows[[length(rows) + 1]] <- v + rnorm(n_s, sd = 0.8)
}
for (g in 1:48) rows[[length(rows) + 1]] <- rnorm(n_s)
x <- do.call(rbind, rows)
x <- x[sample(nrow(x)), ]
rownames(x) <- sprintf("m%03d", 1:nrow(x)); colnames(x) <- sprintf("s%02d", 1:n_s)
write.table(format(x, digits = 17), "modexpr.tsv", sep = "\t", quote = FALSE, col.names = NA)
x <- as.matrix(read.delim("modexpr.tsv", row.names = 1))
d <- t(x)
for (type in c("signed", "unsigned")) {
  diss <- 1 - TOMsimilarityFromExpr(d, power = 6, networkType = type, TOMType = type,
                                    verbose = 0)
  gz <- gzfile(sprintf("diss_%s.tsv.gz", type), "w")
  write.table(format(diss, digits = 17), gz, sep = "\t", quote = FALSE,
              row.names = FALSE, col.names = FALSE)
  close(gz)
  tree <- fastcluster::hclust(as.dist(diss), method = "average")
  write.table(data.frame(a = tree$merge[, 1], b = tree$merge[, 2],
                         height = sprintf("%.17g", tree$height)),
              sprintf("hclust_%s.tsv", type), sep = "\t", quote = FALSE, row.names = FALSE)
  lab <- cutreeDynamic(tree, distM = diss, deepSplit = 2, minClusterSize = 10,
                       cutHeight = 0.995, verbose = 0)
  write.table(data.frame(gene = rownames(x), label = lab), sprintf("cutree_%s.tsv", type),
              sep = "\t", quote = FALSE, row.names = FALSE)
  net <- blockwiseModules(d, power = 6, networkType = type, TOMType = type,
                          maxBlockSize = 10000, minModuleSize = 10, mergeCutHeight = 0.25,
                          numericLabels = TRUE, verbose = 0)
  write.table(data.frame(gene = rownames(x), unmerged = net$unmergedColors,
                         merged = net$colors),
              sprintf("modules_%s.tsv", type), sep = "\t", quote = FALSE, row.names = FALSE)
  cat(type, ": cutree", max(lab), "clusters,", sum(lab == 0), "unassigned; modules",
      max(net$unmergedColors), "->", max(net$colors), "\n")
}
