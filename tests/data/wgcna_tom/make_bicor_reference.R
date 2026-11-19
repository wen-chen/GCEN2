# WGCNA 1.74 reference for gcen2's bicor soft-threshold table (the power
# selection of hdWGCNA's TestSoftPowers: pickSoftThreshold(corFnc = bicor)).
# usage: Rscript make_bicor_reference.R   (run in this directory)
#
# bicor_expr.tsv: the first 60 genes of modexpr.tsv, plus 4 genes with a zero
# MAD (mostly zeros: bicor falls back to Pearson) and 2 with outliers.
#   bicor.tsv           bicor(datExpr) (genes x genes)
#   sft_bicor_signed.tsv pickSoftThreshold(networkType = "signed",
#                        corFnc = bicor, powers 1..10, 12..30 by 2)$fitIndices
suppressPackageStartupMessages(library(WGCNA))
x <- as.matrix(read.delim("modexpr.tsv", row.names = 1))[1:60, ]
set.seed(9)
n <- ncol(x)
z <- matrix(0, 4, n)
for (i in 1:4) z[i, sample(n, 6 + i)] <- rexp(6 + i)
o <- x[1:2, ]
o[1, 3] <- 25; o[2, c(5, 9)] <- c(-30, 18)
x <- rbind(x, z, o)
rownames(x)[61:66] <- c(sprintf("zero%d", 1:4), "out1", "out2")
write.table(format(x, digits = 17), "bicor_expr.tsv", sep = "\t", quote = FALSE,
            col.names = NA)
x <- as.matrix(read.delim("bicor_expr.tsv", row.names = 1))
d <- t(x)
b <- bicor(d)
write.table(format(b, digits = 17), "bicor.tsv", sep = "\t", quote = FALSE, col.names = NA)
sft <- pickSoftThreshold(d, powerVector = c(1:10, seq(12, 30, 2)), networkType = "signed",
                         corFnc = bicor, verbose = 0)$fitIndices
write.table(format(sft, digits = 17), "sft_bicor_signed.tsv", sep = "\t", quote = FALSE,
            row.names = FALSE)
