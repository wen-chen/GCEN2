# WGCNA 1.74 reference for gcen2's TOM edge scores (tests/test_tom.py).
# usage: Rscript make_reference.R   (run in this directory; rewrites the .tsv files)
#
# expr.tsv: 80 genes x 36 samples, four co-expressed groups plus noise (seeded).
# sft_<type>.tsv: pickSoftThreshold table (default powers); tom_<type>.tsv:
# TOMsimilarityFromExpr (Pearson, TOMType = networkType, TOMDenom = "min") at
# power 6, upper triangle.
suppressPackageStartupMessages(library(WGCNA))
set.seed(7)
n_g <- 80; n_s <- 36
base <- matrix(rnorm(4 * n_s), 4)
x <- base[rep(1:4, length.out = n_g), ] * runif(n_g, 0.2, 1.5) + matrix(rnorm(n_g * n_s), n_g)
x[5, ] <- -x[1, ]  # a negatively correlated gene
rownames(x) <- sprintf("g%02d", 1:n_g); colnames(x) <- sprintf("s%02d", 1:n_s)
write.table(format(x, digits = 17), "expr.tsv", sep = "\t", quote = FALSE, col.names = NA)
x <- as.matrix(read.delim("expr.tsv", row.names = 1))
d <- t(x)
for (type in c("signed", "unsigned")) {
  sft <- pickSoftThreshold(d, networkType = type, verbose = 0)
  fi <- sft$fitIndices
  write.table(data.frame(power = fi$Power, r2 = sprintf("%.17g", fi$SFT.R.sq),
                         slope = sprintf("%.17g", fi$slope), mean_k = sprintf("%.17g", fi$mean.k.),
                         estimate = ifelse(is.na(sft$powerEstimate), NA, sft$powerEstimate)),
              sprintf("sft_%s.tsv", type), sep = "\t", quote = FALSE, row.names = FALSE)
  tom <- TOMsimilarityFromExpr(d, power = 6, networkType = type, TOMType = type, verbose = 0)
  ut <- which(upper.tri(tom), arr.ind = TRUE)
  write.table(data.frame(gene_a = rownames(x)[ut[, 1]], gene_b = rownames(x)[ut[, 2]],
                         tom = sprintf("%.9g", tom[ut])),
              sprintf("tom_%s.tsv", type), sep = "\t", quote = FALSE, row.names = FALSE)
}
