# hdWGCNA (v0.4.12) on one PBMC cell type, with the inputs gcen2 uses: the
# same cells, labels, PCA (as the Seurat "pca" reduction) and gene list.
#
# usage: Rscript 02_hdwgcna.R CELL_TYPE OUT_DIR [CELLS_FILE|-] [SEED]
#
# Writes to OUT_DIR: datExpr.tsv (metacells x genes), modules.tsv (gene,
# module after ConstructNetwork; grey = unassigned), modules_final.tsv (after
# ModuleConnectivity's reassignment, with kME_<module>), MEs.tsv (cells x
# modules, ModuleEigengenes), tom_top.tsv (gene pairs with the largest TOM),
# members.tsv (metacell member barcodes), timing.tsv.

# enrichR (a dependency of hdWGCNA) contacts maayanlab.cloud when attached,
# which times out now and then: attach it first, with retries
for (i in 1:10) {
  if (tryCatch({suppressPackageStartupMessages(library(enrichR)); TRUE},
               error = function(e) FALSE)) break
  Sys.sleep(5)
}
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(WGCNA); library(hdWGCNA)
})
args <- commandArgs(trailingOnly = TRUE)
ct <- args[1]; out <- args[2]
cells_file <- if (length(args) > 2 && args[3] != "-") args[3] else NA
seed <- if (length(args) > 3) as.integer(args[4]) else 42L
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(seed)
threads <- 8
enableWGCNAThreads(nThreads = threads)
W <- "work"  # run from the PBMC data directory
n_top <- 300000

t0 <- proc.time()[["elapsed"]]
counts <- readMM(file.path(W, "counts.mtx"))
rownames(counts) <- readLines(file.path(W, "genes.tsv"))
colnames(counts) <- readLines(file.path(W, "barcodes.tsv"))
meta <- read.delim(file.path(W, "meta.tsv"), row.names = 1, check.names = FALSE)
pca <- as.matrix(read.delim(file.path(W, "pca.tsv"), row.names = 1, check.names = FALSE))
cells <- if (is.na(cells_file)) rownames(meta)[meta$cell_type == ct] else readLines(cells_file)
genes <- readLines(file.path(W, paste0("genes_", gsub(" ", "_", ct), ".txt")))

so <- CreateSeuratObject(as(counts[, cells], "CsparseMatrix"),
                         meta.data = meta[cells, , drop = FALSE])
so <- NormalizeData(so, verbose = FALSE)
so[["pca"]] <- CreateDimReducObject(embeddings = pca[cells, ], key = "PC_", assay = "RNA")
t_read <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
so <- SetupForWGCNA(so, features = genes, wgcna_name = ct)
# Metacell centres are drawn with sample(), but MetacellsByGroups subsets the
# Seurat object per group just before, and subsetting a Seurat object resets
# R's RNG (set.seed(1) and set.seed(2) give the same runif() afterwards), so
# a seed set by the caller has no effect. Seed inside ConstructMetacells.
trace("ConstructMetacells", bquote(set.seed(.(seed))), print = FALSE,
      where = asNamespace("hdWGCNA"))
so <- MetacellsByGroups(so, group.by = "cell_type", ident.group = "cell_type",
                        reduction = "pca", dims = 1:30, k = 25, max_shared = 10,
                        min_cells = 100, target_metacells = 1000, mode = "average")
so <- NormalizeMetacells(so)
t_metacells <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
so <- SetDatExpr(so, group_name = ct, group.by = "cell_type", assay = "RNA",
                 layer = "data")
so <- TestSoftPowers(so, networkType = "signed")
so <- ConstructNetwork(so, tom_name = gsub(" ", "_", ct),
                       tom_outdir = file.path(out, "TOM"), overwrite_tom = TRUE)
t_network <- proc.time()[["elapsed"]] - t0

power <- tryCatch(GetWGCNAParams(so)$power, error = function(e) NULL)
if (is.null(power)) power <- NA
modules <- GetModules(so)
write.table(data.frame(gene = modules$gene_name, module = modules$module),
            file.path(out, "modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
datExpr <- GetDatExpr(so)
write.table(datExpr, file.path(out, "datExpr.tsv"), sep = "\t", quote = FALSE,
            col.names = NA)
mc <- GetMetacellObject(so)
write.table(data.frame(metacell = colnames(mc), cells = mc$cells_merged),
            file.path(out, "members.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

# eigengenes and kME on the single cells (the object holds this type's cells
# only); ModuleConnectivity also reassigns genes with a negative kME
t0 <- proc.time()[["elapsed"]]
so <- ModuleEigengenes(so, verbose = FALSE)
so <- ModuleConnectivity(so)
t_me <- proc.time()[["elapsed"]] - t0
mf <- GetModules(so)
write.table(data.frame(gene = mf$gene_name, module = mf$module,
                       mf[, grep("^kME_", colnames(mf))], check.names = FALSE),
            file.path(out, "modules_final.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(GetMEs(so, harmonized = FALSE), file.path(out, "MEs.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)

TOM <- GetTOM(so)
g <- colnames(TOM)
ut <- which(upper.tri(TOM))
vals <- TOM[ut]
top <- ut[order(-vals)[seq_len(min(n_top, length(ut)))]]
rr <- (top - 1) %% nrow(TOM) + 1; cc <- (top - 1) %/% nrow(TOM) + 1
write.table(data.frame(gene_a = g[rr], gene_b = g[cc], tom = TOM[top]),
            file.path(out, "tom_top.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

write.table(data.frame(cells = length(cells), genes = ncol(datExpr),
                       metacells = nrow(datExpr), power = power,
                       modules = length(setdiff(unique(as.character(modules$module)), "grey")),
                       t_read = t_read, t_metacells = t_metacells, t_network = t_network,
                       t_me = t_me),
            file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%s: %d metacells, %d genes, power %s, network %.1fs (metacells %.1fs)\n",
            ct, nrow(datExpr), ncol(datExpr), power, t_network, t_metacells))
