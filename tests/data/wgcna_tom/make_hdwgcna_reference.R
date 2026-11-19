# WGCNA 1.74 reference for gcen2's hdwgcna preset: what hdWGCNA 0.4.12
# ConstructNetwork runs on a datExpr (blockwiseConsensusModules on one set
# with its defaults), with the power of TestSoftPowers (pickSoftThreshold,
# corFnc = bicor, powers 1..10, 12..30; lowest power > 3 with SFT.R.sq >= 0.8).
# usage: Rscript make_hdwgcna_reference.R   (run in this directory)
#   hdwgcna_modules.tsv  gene, merged labels (minModuleSize 10 instead of 50
#                        for the 200 genes of modexpr.tsv); power in hdwgcna_power.txt
suppressPackageStartupMessages(library(WGCNA))
x <- as.matrix(read.delim("modexpr.tsv", row.names = 1))
d <- t(x)
sft <- pickSoftThreshold(d, powerVector = c(1:10, seq(12, 30, 2)), networkType = "signed",
                         corFnc = bicor, verbose = 0)$fitIndices
power <- min(subset(sft, SFT.R.sq >= 0.8 & Power > 3)$Power)
writeLines(as.character(power), "hdwgcna_power.txt")
tmp <- tempfile(); dir.create(tmp)
net <- blockwiseConsensusModules(
  list(set = list(data = d)), power = power, maxBlockSize = 30000, randomSeed = 12345,
  corType = "pearson", consensusQuantile = 0.3, networkType = "signed", TOMType = "signed",
  TOMDenom = "min", calibrationQuantile = 0.95, sampleForCalibration = TRUE,
  sampleForCalibrationFactor = 1000, useDiskCache = TRUE, cacheDir = tmp, deepSplit = 4,
  pamStage = FALSE, detectCutHeight = 0.995, minModuleSize = 10, mergeCutHeight = 0.2,
  saveConsensusTOMs = FALSE, numericLabels = TRUE, verbose = 0)
write.table(data.frame(gene = rownames(x), module = net$colors), "hdwgcna_modules.tsv",
            sep = "\t", quote = FALSE, row.names = FALSE)
cat("power", power, "modules", max(net$colors), "unassigned", sum(net$colors == 0), "\n")
