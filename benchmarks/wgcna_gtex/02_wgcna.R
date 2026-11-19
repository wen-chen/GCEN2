# WGCNA 1.74 on a genes x samples log2(TPM + 1) matrix: the standard
# pipeline (pickSoftThreshold, then blockwiseModules, signed network).
#
# usage: Rscript 02_wgcna.R EXPR_TSV OUT_DIR [SAMPLES_FILE|-] [tom]
#
# Writes to OUT_DIR: modules.tsv (gene, module; grey = unassigned),
# timing.tsv and, with "tom", tom_top.tsv (the gene pairs with the largest TOM).

suppressPackageStartupMessages(library(WGCNA))
args <- commandArgs(trailingOnly = TRUE)
expr_file <- args[1]; out <- args[2]
samples_file <- if (length(args) > 2 && args[3] != "-") args[3] else NA
want_tom <- length(args) > 3 && args[4] == "tom"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
threads <- 8
enableWGCNAThreads(nThreads = threads)
options(stringsAsFactors = FALSE)
n_top <- 300000

t0 <- proc.time()[["elapsed"]]
x <- read.delim(expr_file, row.names = 1, check.names = FALSE)
datExpr <- t(as.matrix(x))
if (!is.na(samples_file)) datExpr <- datExpr[readLines(samples_file), ]
t_read <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)),
                         networkType = "signed", verbose = 0)
power <- sft$powerEstimate
# WGCNA FAQ: no power reaches the scale-free fit -> 12 for a signed network
# with more than 40 samples
if (is.na(power)) power <- 12
t_power <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
net <- blockwiseModules(datExpr, power = power, networkType = "signed",
                        TOMType = "signed", maxBlockSize = 30000, minModuleSize = 30,
                        mergeCutHeight = 0.25, numericLabels = FALSE,
                        saveTOMs = want_tom, saveTOMFileBase = file.path(out, "TOM"),
                        nThreads = threads, verbose = 0)
t_network <- proc.time()[["elapsed"]] - t0

write.table(data.frame(gene = colnames(datExpr), module = net$colors),
            file.path(out, "modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
if (want_tom) {
  load(net$TOMFiles[1])  # TOM (dist)
  TOM <- as.matrix(TOM)
  g <- colnames(datExpr)[net$blockGenes[[1]]]
  ut <- which(upper.tri(TOM))
  top <- ut[order(-TOM[ut])[seq_len(min(n_top, length(ut)))]]
  rr <- (top - 1) %% nrow(TOM) + 1; cc <- (top - 1) %/% nrow(TOM) + 1
  write.table(data.frame(gene_a = g[rr], gene_b = g[cc], tom = TOM[top]),
              file.path(out, "tom_top.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  unlink(net$TOMFiles)
}
if (want_tom) {
  # WGCNA's correlations of the first 300 genes, to check GCEN2's against
  cm <- WGCNA::cor(datExpr[, 1:300])
  ut <- which(upper.tri(cm), arr.ind = TRUE)
  write.table(data.frame(gene_a = colnames(cm)[ut[, 1]], gene_b = colnames(cm)[ut[, 2]],
                         cor = sprintf("%.17g", cm[ut])),
              file.path(out, "cor_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
}
n_mod <- length(setdiff(unique(net$colors), "grey"))
write.table(data.frame(samples = nrow(datExpr), genes = ncol(datExpr), power = power,
                       modules = n_mod, module_genes = sum(net$colors != "grey"),
                       t_read = t_read, t_power = t_power, t_network = t_network),
            file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d samples x %d genes: power %s, %d modules (%d genes); power %.1fs, network %.1fs\n",
            nrow(datExpr), ncol(datExpr), power, n_mod, sum(net$colors != "grey"),
            t_power, t_network))
