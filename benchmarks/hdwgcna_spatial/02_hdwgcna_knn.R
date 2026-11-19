# hdWGCNA (v0.4.12) on spatial data without a Visium grid (Slide-seqV2,
# exported by 01_prepare_slideseq.py). hdWGCNA's MetaspotsByGroups needs the
# Visium grid, so this is the hdWGCNA workflow with metacells built in the
# spatial coordinates instead (MetacellsByGroups on a "spatial" reduction):
# SetupForWGCNA(gene_select = "fraction") -> MetacellsByGroups (k nearest
# beads, summed) -> NormalizeMetacells -> SetDatExpr -> TestSoftPowers
# -> ConstructNetwork -> ModuleEigengenes (beads) -> ModuleConnectivity.
#
# Also used for Visium HD bins (01_prepare_visium_hd.py: EXPORT_DIR/matrix.h5,
# genes.tsv, meta.tsv with a region column, GROUP_BY = region).
#
# usage: Rscript 02_hdwgcna_knn.R EXPORT_DIR OUT_DIR [K] [MAX_SHARED] [FRACTION] [SEED] [GROUP_BY]
# Writes genes.txt, datExpr.tsv (metacells x genes), modules_network.tsv,
# modules.tsv (after ModuleConnectivity, with kME_<module>), MEs.tsv and
# timing.tsv to OUT_DIR.
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
k <- if (length(args) > 2) as.integer(args[3]) else 100L
max_shared <- if (length(args) > 3) as.integer(args[4]) else round(2 * k / 7)
fraction <- if (length(args) > 4) as.numeric(args[5]) else 0.01
seed <- if (length(args) > 5) as.integer(args[6]) else 42L
group_by <- if (length(args) > 6) args[7] else "all"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 8)

t0 <- proc.time()[["elapsed"]]
if (file.exists(file.path(W, "matrix.h5"))) {
  counts <- Read10X_h5(file.path(W, "matrix.h5"), use.names = TRUE)
  rownames(counts) <- readLines(file.path(W, "genes.tsv"))
} else {
  counts <- readMM(file.path(W, "counts.mtx"))
  rownames(counts) <- readLines(file.path(W, "genes.tsv"))
  colnames(counts) <- readLines(file.path(W, "barcodes.tsv"))
}
meta <- read.delim(file.path(W, "meta.tsv"), row.names = 1, check.names = FALSE)
counts <- counts[, rownames(meta)]
meta$all <- "all"
so <- CreateSeuratObject(as(counts, "CsparseMatrix"), meta.data = meta)
so <- NormalizeData(so, verbose = FALSE)
xy <- as.matrix(meta[colnames(so), c("x", "y")])
colnames(xy) <- c("spatial_1", "spatial_2")
so[["spatial"]] <- CreateDimReducObject(embeddings = xy, key = "spatial_", assay = "RNA")
t_read <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
so <- SetupForWGCNA(so, gene_select = "fraction", fraction = fraction, wgcna_name = "ss")
# MetacellsByGroups subsets the object per group, which resets R's RNG, so
# seed inside ConstructMetacells (see ../hdwgcna_pbmc/02_hdwgcna.R)
trace("ConstructMetacells", bquote(set.seed(.(seed))), print = FALSE,
      where = asNamespace("hdWGCNA"))
so <- MetacellsByGroups(so, group.by = group_by, ident.group = group_by, reduction = "spatial",
                        dims = 1:2, k = k, max_shared = max_shared, min_cells = 50,
                        target_metacells = 100000, max_iter = 100000, mode = "sum")
so <- NormalizeMetacells(so)
t_metaspots <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
so <- SetDatExpr(so, group.by = NULL, group_name = NULL, assay = "RNA")
so <- TestSoftPowers(so)
so <- ConstructNetwork(so, tom_name = "ss", tom_outdir = file.path(out, "TOM"),
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
nmod <- length(setdiff(unique(as.character(modules$module)), "grey"))
write.table(data.frame(beads = ncol(so), genes = ncol(datExpr), metacells = nrow(datExpr),
                       power = power, modules = nmod, t_read = t_read,
                       t_metaspots = t_metaspots, t_network = t_network, t_me = t_me),
            file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d metacells, %d genes, power %s, %d modules; metacells %.1fs, network %.1fs, MEs %.1fs\n",
            nrow(datExpr), ncol(datExpr), power, nmod, t_metaspots, t_network, t_me))
