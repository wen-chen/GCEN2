# hdWGCNA (v0.4.12) per cell type of the HLCA core, the single-cell workflow
# with the settings of ../hdwgcna_pbmc: for each cell type, its cells only
# -> NormalizeData -> SetupForWGCNA(gene_select = "fraction", 0.05) ->
# MetacellsByGroups (the scANVI embedding as the reduction, 30 dims, k 25,
# max_shared 10, up to 1000 metacells) -> NormalizeMetacells -> SetDatExpr
# -> TestSoftPowers(signed) -> ConstructNetwork -> ModuleEigengenes and
# ModuleConnectivity on the type's cells.
#
# The atlas has more non-zero counts than an R sparse matrix holds (2^31 - 1),
# so the cells of each type are read on their own: work/hlca_<size>.h5 is
# sorted by cell type (01_prepare.py), and a type is one slice.
#
# usage: Rscript 02_hdwgcna.R H5 TYPES_FILE OUT_DIR [SEED]
# Writes OUT_DIR/<type>/{genes.txt, datExpr.tsv, modules.tsv (ConstructNetwork),
# modules_final.tsv (after ModuleConnectivity, with kME)} and
# OUT_DIR/timing.tsv (one row per cell type).
for (i in 1:10) {
  if (tryCatch({suppressPackageStartupMessages(library(enrichR)); TRUE},
               error = function(e) FALSE)) break
  Sys.sleep(5)
}
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(WGCNA); library(hdWGCNA); library(hdf5r)
})
args <- commandArgs(trailingOnly = TRUE)
h5 <- args[1]; types <- readLines(args[2]); out <- args[3]
seed <- if (length(args) > 3) as.integer(args[4]) else 42L
dir.create(out, recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 8)
# MetacellsByGroups subsets the object per group, which resets R's RNG: seed
# inside ConstructMetacells (see ../hdwgcna_pbmc/02_hdwgcna.R)
trace("ConstructMetacells", bquote(set.seed(.(seed))), print = FALSE,
      where = asNamespace("hdWGCNA"))

h <- H5File$new(h5, mode = "r")
ip <- h[["matrix/indptr_f"]][]
cell_type <- h[["obs/cell_type"]][]
barcodes <- h[["matrix/barcodes"]][]
genes <- h[["matrix/features/name"]][]

read_type <- function(t) {
  cols <- which(cell_type == t)
  i0 <- min(cols); i1 <- max(cols)
  stopifnot(length(cols) == i1 - i0 + 1)  # contiguous
  a <- ip[i0] + 1; b <- ip[i1 + 1]
  x <- sparseMatrix(i = h[["matrix/indices"]][a:b] + 1L, p = ip[i0:(i1 + 1)] - ip[i0],
                    x = as.numeric(h[["matrix/data"]][a:b]),
                    dims = c(length(genes), length(cols)),
                    dimnames = list(genes, barcodes[cols]))
  rep <- t(h[["obsm/rep"]][, i0:i1])
  rownames(rep) <- barcodes[cols]
  colnames(rep) <- paste0("rep_", seq_len(ncol(rep)))
  list(counts = x, rep = rep)
}

rows <- list()
for (t in types) {
  name <- gsub("[^A-Za-z0-9]+", "_", t)
  d <- file.path(out, name)
  dir.create(d, recursive = TRUE, showWarnings = FALSE)
  t0 <- proc.time()[["elapsed"]]
  r <- read_type(t)
  so <- CreateSeuratObject(r$counts, meta.data = data.frame(
    cell_type = rep(t, ncol(r$counts)), row.names = colnames(r$counts)))
  so <- NormalizeData(so, verbose = FALSE)
  so[["rep"]] <- CreateDimReducObject(embeddings = r$rep, key = "rep_", assay = "RNA")
  rm(r); gc(verbose = FALSE)
  t_read <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  so <- SetupForWGCNA(so, gene_select = "fraction", fraction = 0.05, wgcna_name = name)
  so <- MetacellsByGroups(so, group.by = "cell_type", ident.group = "cell_type",
                          reduction = "rep", dims = 1:30, k = 25, max_shared = 10,
                          min_cells = 100, target_metacells = 1000, mode = "average")
  so <- NormalizeMetacells(so)
  t_metacells <- proc.time()[["elapsed"]] - t0

  t0 <- proc.time()[["elapsed"]]
  so <- SetDatExpr(so, group_name = t, group.by = "cell_type", assay = "RNA", layer = "data")
  so <- TestSoftPowers(so, networkType = "signed")
  so <- ConstructNetwork(so, tom_name = name, tom_outdir = file.path(d, "TOM"),
                         overwrite_tom = TRUE)
  t_network <- proc.time()[["elapsed"]] - t0

  power <- tryCatch(GetWGCNAParams(so)$power, error = function(e) NA)
  modules <- GetModules(so)  # after ConstructNetwork

  # eigengenes and kME on the type's single cells; ModuleConnectivity also
  # reassigns genes with a negative kME
  t0 <- proc.time()[["elapsed"]]
  so <- ModuleEigengenes(so, verbose = FALSE)
  so <- ModuleConnectivity(so)
  t_me <- proc.time()[["elapsed"]] - t0
  mf <- GetModules(so)
  write.table(data.frame(gene = mf$gene_name, module = mf$module,
                         mf[, grep("^kME_", colnames(mf))], check.names = FALSE),
              file.path(d, "modules_final.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  datExpr <- GetDatExpr(so)
  writeLines(GetWGCNAGenes(so), file.path(d, "genes.txt"))
  write.table(data.frame(gene = modules$gene_name, module = modules$module),
              file.path(d, "modules.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(datExpr, file.path(d, "datExpr.tsv"), sep = "\t", quote = FALSE, col.names = NA)
  unlink(file.path(d, "TOM"), recursive = TRUE)
  nmod <- length(setdiff(unique(as.character(modules$module)), "grey"))
  rows[[t]] <- data.frame(cell_type = t, cells = ncol(so), genes = ncol(datExpr),
                          metacells = nrow(datExpr), power = power, modules = nmod,
                          t_read = t_read, t_metacells = t_metacells, t_network = t_network,
                          t_me = t_me)
  cat(sprintf("%s: %d cells, %d metacells, %d genes, power %s, %d modules; read %.1fs, metacells %.1fs, network %.1fs, MEs + kME %.1fs\n",
              t, ncol(so), nrow(datExpr), ncol(datExpr), power, nmod, t_read, t_metacells,
              t_network, t_me))
  rm(so); gc(verbose = FALSE)
}
write.table(do.call(rbind, rows), file.path(out, "timing.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)
