# DGCA 1.0.3: differential correlation of every gene pair between two
# conditions (ddcorAll, no permutations, Benjamini-Hochberg), written with
# fwrite.
# Usage: Rscript 02_dgca.R <A.tsv> <B.tsv> <out.tsv> [threads] [pearson|spearman]
# A and B are genes x samples TSVs with the same genes (01_prepare.py).
suppressPackageStartupMessages({
  library(DGCA)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
threads <- if (length(args) >= 4) as.integer(args[4]) else 8
method <- if (length(args) >= 5) args[5] else "pearson"
WGCNA::allowWGCNAThreads(threads)

t0 <- Sys.time()
read <- function(p) {
  d <- fread(p, sep = "\t", data.table = FALSE)
  m <- as.matrix(d[, -1])
  rownames(m) <- d[[1]]
  m
}
a <- read(args[1])
b <- read(args[2])
stopifnot(identical(rownames(a), rownames(b)))
design <- makeDesign(c(rep("A", ncol(a)), rep("B", ncol(b))))
# corr_cutoff 0.9999: the clipping GCEN2 applies before Fisher's z
# (DGCA's default is 0.99).
res <- ddcorAll(inputMat = cbind(a, b), design = design, compare = c("A", "B"),
                corrType = method, adjust = "BH", nPerms = 0,
                corr_cutoff = 0.9999, classify = TRUE)
fwrite(res, args[3], sep = "\t")
cat(sprintf("DGCA %s (%s): %d genes, %d pairs, %.1f s\n", packageVersion("DGCA"), method,
            nrow(a), nrow(res), as.numeric(Sys.time() - t0, units = "secs")))
