# WGCNA 1.74 reference for the module eigengenes, kME and module-trait
# correlations of gcen2's TOM modules (tests/test_tom_modules.py).
# usage: Rscript make_downstream_reference.R   (run in this directory, after
# make_modules_reference.R)
#
# traits.tsv: 30 samples x 3 traits (seeded; "batch" has two missing values).
# For the merged blockwiseModules labels in modules_<type>.tsv:
#   me_<type>.tsv     moduleEigengenes(eigengenes) incl. ME0, and varExplained
#   kme_<type>.tsv    signedKME(datExpr, MEs) and corPvalueStudent p-values
#   trait_<type>.tsv  corAndPvalue(MEs, traits, use = "p"): cor, p, nObs
suppressPackageStartupMessages(library(WGCNA))
x <- as.matrix(read.delim("modexpr.tsv", row.names = 1))
d <- t(x)
set.seed(5)
traits <- data.frame(row.names = rownames(d),
                     weight = d[, 1] + rnorm(nrow(d), sd = 0.7),
                     age = round(runif(nrow(d), 20, 70)),
                     batch = rep(c(0, 1), length.out = nrow(d)))
traits$batch[c(4, 17)] <- NA
write.table(traits, "traits.tsv", sep = "\t", quote = FALSE, col.names = NA, na = "NA")
f17 <- function(m) format(m, digits = 17)
for (type in c("signed", "unsigned")) {
  lab <- read.delim(sprintf("modules_%s.tsv", type))$merged
  me <- moduleEigengenes(d, lab)
  MEs <- me$eigengenes
  write.table(f17(as.matrix(MEs)), sprintf("me_%s.tsv", type), sep = "\t", quote = FALSE,
              col.names = NA)
  write.table(data.frame(module = colnames(MEs), var = sprintf("%.17g", unlist(me$varExplained))),
              sprintf("varexp_%s.tsv", type), sep = "\t", quote = FALSE, row.names = FALSE)
  kme <- signedKME(d, MEs)
  p <- corPvalueStudent(as.matrix(kme), nrow(d))
  colnames(p) <- sub("^kME", "p.", colnames(kme))
  write.table(f17(cbind(as.matrix(kme), p)), sprintf("kme_%s.tsv", type), sep = "\t",
              quote = FALSE, col.names = NA)
  cp <- corAndPvalue(MEs, traits, use = "p")
  out <- data.frame(module = rep(colnames(MEs), ncol(traits)),
                    trait = rep(colnames(traits), each = ncol(MEs)),
                    cor = sprintf("%.17g", as.vector(cp$cor)),
                    p = sprintf("%.17g", as.vector(cp$p)),
                    n = as.vector(cp$nObs))
  write.table(out, sprintf("trait_%s.tsv", type), sep = "\t", quote = FALSE, row.names = FALSE)
}
