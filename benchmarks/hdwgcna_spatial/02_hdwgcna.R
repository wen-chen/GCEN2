# hdWGCNA (v0.4.12) spatial workflow (the "hdWGCNA in spatial transcriptomics"
# tutorial) on a Visium sample exported by 01_prepare.py:
# SetupForWGCNA(gene_select = "fraction") -> MetaspotsByGroups -> NormalizeMetacells
# -> SetDatExpr (all metaspots) -> TestSoftPowers -> ConstructNetwork
# -> ModuleEigengenes (spots) -> ModuleConnectivity.
#
# usage: Rscript 02_hdwgcna.R EXPORT_DIR OUT_DIR [GROUP_BY]
# Writes genes.txt, datExpr.tsv (metaspots x genes), modules_network.tsv
# (ConstructNetwork), modules.tsv (after ModuleConnectivity's reassignment, with
# the kME_<module> columns), MEs.tsv (spots x modules), timing.tsv to OUT_DIR.
for (i in 1:10) {
  if (tryCatch({suppressPackageStartupMessages(library(enrichR)); TRUE},
               error = function(e) FALSE)) break
  Sys.sleep(5)
}
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(WGCNA); library(hdWGCNA)
})
args <- commandArgs(trailingOnly = TRUE)
W <- args[1]; out <- args[2]
group_by <- if (length(args) > 2) args[3] else "region"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(42)
# hdWGCNA 0.4.12 MetaspotsByGroups checks its input with GetAssayData(slot = ),
# which SeuratObject 5 made defunct; patch that check only (the aggregation,
# ConstructMetaspots, already uses LayerData)
f <- hdWGCNA::MetaspotsByGroups
src <- paste(deparse(f), collapse = "\n")
old <- "GetAssayData(seurat_obj, assay = assay, \n            slot = slot)"
stopifnot(grepl(old, src, fixed = TRUE))
src <- sub(old, "SeuratObject::LayerData(seurat_obj, assay = assay, layer = slot)", src,
           fixed = TRUE)
f <- eval(parse(text = src)); environment(f) <- asNamespace("hdWGCNA")
MetaspotsByGroups <- f
enableWGCNAThreads(nThreads = 8)

t0 <- proc.time()[["elapsed"]]
counts <- readMM(file.path(W, "counts.mtx"))
rownames(counts) <- readLines(file.path(W, "genes.tsv"))
colnames(counts) <- readLines(file.path(W, "barcodes.tsv"))
meta <- read.delim(file.path(W, "meta.tsv"), row.names = 1, check.names = FALSE)
so <- CreateSeuratObject(as(counts, "CsparseMatrix"), assay = "Spatial",
                         meta.data = meta)
img <- Read10X_Image(file.path(W, "spatial"), assay = "Spatial")
so[["slice1"]] <- img[colnames(so)]
so <- NormalizeData(so, verbose = FALSE)
t_read <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
so <- SetupForWGCNA(so, gene_select = "fraction", fraction = 0.05, wgcna_name = "vis")
so <- MetaspotsByGroups(so, group.by = group_by, ident.group = group_by, assay = "Spatial")
so <- NormalizeMetacells(so)
t_metaspots <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
so <- SetDatExpr(so, group.by = NULL, group_name = NULL, assay = "Spatial")
so <- TestSoftPowers(so)
so <- ConstructNetwork(so, tom_name = "vis", tom_outdir = file.path(out, "TOM"),
                       overwrite_tom = TRUE)
t_network <- proc.time()[["elapsed"]] - t0
m0 <- GetModules(so)
write.table(data.frame(gene = m0$gene_name, module = m0$module),
            file.path(out, "modules_network.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

t0 <- proc.time()[["elapsed"]]
so <- ModuleEigengenes(so, verbose = FALSE)
so <- ModuleConnectivity(so)
t_me <- proc.time()[["elapsed"]] - t0

power <- tryCatch(GetWGCNAParams(so)$power, error = function(e) NA)
modules <- GetModules(so)
writeLines(GetWGCNAGenes(so), file.path(out, "genes.txt"))
write.table(data.frame(gene = modules$gene_name, module = modules$module,
                       modules[, grep("^kME_", colnames(modules))], check.names = FALSE),
            file.path(out, "modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
datExpr <- GetDatExpr(so)
write.table(datExpr, file.path(out, "datExpr.tsv"), sep = "\t", quote = FALSE, col.names = NA)
write.table(GetMEs(so, harmonized = FALSE), file.path(out, "MEs.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)
write.table(data.frame(spots = ncol(so), genes = ncol(datExpr), metaspots = nrow(datExpr),
                       power = power,
                       modules = length(setdiff(unique(as.character(modules$module)), "grey")),
                       t_read = t_read, t_metaspots = t_metaspots, t_network = t_network,
                       t_me = t_me),
            file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d metaspots, %d genes, power %s, %d modules; metaspots %.1fs, network %.1fs, MEs %.1fs\n",
            nrow(datExpr), ncol(datExpr), power,
            length(setdiff(unique(as.character(modules$module)), "grey")),
            t_metaspots, t_network, t_me))
