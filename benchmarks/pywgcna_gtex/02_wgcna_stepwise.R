# WGCNA 1.74, step by step (the WGCNA tutorial's pipeline, which PyWGCNA
# reimplements): pickSoftThreshold, adjacency, TOMsimilarity, average-linkage
# hclust, cutreeDynamic (hybrid, deepSplit 2, pamRespectsDendro FALSE as
# PyWGCNA's findModules), labels2colors, mergeCloseModules. Signed network
# and TOM, minimum module size 30, merge cut 0.25; the power as 02_wgcna.R of
# ../wgcna_gtex (pickSoftThreshold, 12 when no power reaches the fit).
#
# usage: Rscript 02_wgcna_stepwise.R EXPR_TSV OUT_DIR
# Writes to OUT_DIR: modules.tsv (gene, module; grey = unassigned), timing.tsv.

suppressPackageStartupMessages(library(WGCNA))
args <- commandArgs(trailingOnly = TRUE)
expr_file <- args[1]; out <- args[2]
dir.create(out, recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 8)
options(stringsAsFactors = FALSE)

t0 <- proc.time()[["elapsed"]]
x <- read.delim(expr_file, row.names = 1, check.names = FALSE)
datExpr <- t(as.matrix(x))
t_read <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
sft <- pickSoftThreshold(datExpr, powerVector = c(1:10, seq(12, 20, 2)),
                         networkType = "signed", verbose = 0)
power <- sft$powerEstimate
if (is.na(power)) power <- 12
t_power <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
adj <- adjacency(datExpr, power = power, type = "signed")
TOM <- TOMsimilarity(adj, TOMType = "signed", verbose = 0)
rm(adj)
dissTOM <- 1 - TOM
rm(TOM)
geneTree <- hclust(as.dist(dissTOM), method = "average")
dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM, deepSplit = 2,
                             pamRespectsDendro = FALSE, minClusterSize = 30, verbose = 0)
rm(dissTOM)
dynamicColors <- labels2colors(dynamicMods)
merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = 0.25, verbose = 0)
colors <- merge$colors
t_network <- proc.time()[["elapsed"]] - t0

write.table(data.frame(gene = colnames(datExpr), module = colors),
            file.path(out, "modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
n_mod <- length(setdiff(unique(colors), "grey"))
write.table(data.frame(samples = nrow(datExpr), genes = ncol(datExpr), power = power,
                       modules = n_mod, module_genes = sum(colors != "grey"),
                       t_read = t_read, t_power = t_power, t_network = t_network),
            file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d samples x %d genes: power %s, %d modules (%d genes); power %.1fs, network %.1fs\n",
            nrow(datExpr), ncol(datExpr), power, n_mod, sum(colors != "grey"),
            t_power, t_network))
