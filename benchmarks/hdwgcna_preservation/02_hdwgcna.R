# hdWGCNA (v0.4.12) module projection and preservation (the "Projecting
# modules to new datasets" and "Module preservation and reproducibility"
# tutorials) on HLCA natural killer cells: the modules of one dataset (REF)
# projected onto, and tested for preservation in, another (QUERY).
#   reference: NormalizeData -> SetupForWGCNA(fraction 0.05) ->
#     MetacellsByGroups (scANVI embedding, 30 dims, k 25, max_shared 10,
#     up to 1000) -> NormalizeMetacells -> SetDatExpr -> TestSoftPowers ->
#     ConstructNetwork -> ModuleEigengenes -> ModuleConnectivity
#   query: NormalizeData -> ProjectModules -> MetacellsByGroups (as above) ->
#     NormalizeMetacells -> SetDatExpr -> ModulePreservation (N_PERM
#     permutations, seed 12345)
#
# usage: Rscript 02_hdwgcna.R EXPORT_DIR OUT_DIR [REF] [QUERY] [N_PERM]
# EXPORT_DIR: the export of ../hdwgcna_consensus/01_prepare.py.
# Writes to OUT_DIR: ref_modules.tsv (gene, module after ModuleConnectivity),
# ref_datExpr.tsv and query_datExpr.tsv (metacells x genes), projected_MEs.tsv
# (query cells), preservation_obs.tsv, preservation_Z.tsv and timing.tsv.
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
ref_set <- if (length(args) > 2) args[3] else "Meyer_2019"
query_set <- if (length(args) > 3) args[4] else "Krasnow_2020"
n_perm <- if (length(args) > 4) as.integer(args[5]) else 500L
dir.create(out, recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 8)
trace("ConstructMetacells", quote(set.seed(42)), print = FALSE, where = asNamespace("hdWGCNA"))

t0 <- proc.time()[["elapsed"]]
counts <- Read10X_h5(file.path(W, "matrix.h5"), use.names = TRUE)
rownames(counts) <- readLines(file.path(W, "genes.tsv"))
meta <- read.delim(file.path(W, "meta.tsv"), row.names = 1, check.names = FALSE)
rep <- as.matrix(read.delim(file.path(W, "rep.tsv"), row.names = 1, check.names = FALSE))
colnames(rep) <- paste0("rep_", seq_len(ncol(rep)))
make_obj <- function(set) {
  cells <- rownames(meta)[meta$dataset == set]
  so <- CreateSeuratObject(counts[, cells], meta.data = meta[cells, ])
  so <- NormalizeData(so, verbose = FALSE)
  so[["rep"]] <- CreateDimReducObject(embeddings = rep[cells, ], key = "rep_", assay = "RNA")
  so
}
ref <- make_obj(ref_set)
query <- make_obj(query_set)
type <- unique(meta$cell_type)
t_read <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
ref <- SetupForWGCNA(ref, gene_select = "fraction", fraction = 0.05, wgcna_name = "ref")
ref <- MetacellsByGroups(ref, group.by = "cell_type", ident.group = "cell_type", reduction = "rep",
                         dims = 1:30, k = 25, max_shared = 10, min_cells = 100,
                         target_metacells = 1000, mode = "average")
ref <- NormalizeMetacells(ref)
ref <- SetDatExpr(ref, group_name = type, group.by = "cell_type", assay = "RNA", layer = "data")
ref <- TestSoftPowers(ref, networkType = "signed")
ref <- ConstructNetwork(ref, tom_name = "ref", tom_outdir = file.path(out, "TOM"),
                        overwrite_tom = TRUE)
ref <- ModuleEigengenes(ref, verbose = FALSE)
ref <- ModuleConnectivity(ref)
t_ref <- proc.time()[["elapsed"]] - t0
unlink(file.path(out, "TOM"), recursive = TRUE)
mods <- GetModules(ref)
write.table(data.frame(gene = mods$gene_name, module = mods$module), file.path(out, "ref_modules.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(GetDatExpr(ref), file.path(out, "ref_datExpr.tsv"), sep = "\t", quote = FALSE,
            col.names = NA)

t0 <- proc.time()[["elapsed"]]
query <- ProjectModules(seurat_obj = query, seurat_ref = ref, wgcna_name = "ref",
                        wgcna_name_proj = "projected", assay = "RNA")
t_project <- proc.time()[["elapsed"]] - t0
write.table(format(GetMEs(query, harmonized = FALSE, wgcna_name = "projected"), digits = 17),
            file.path(out, "projected_MEs.tsv"), sep = "\t", quote = FALSE, col.names = NA)

t0 <- proc.time()[["elapsed"]]
query <- MetacellsByGroups(query, group.by = "cell_type", ident.group = "cell_type",
                           reduction = "rep", dims = 1:30, k = 25, max_shared = 10,
                           min_cells = 100, target_metacells = 1000, mode = "average",
                           wgcna_name = "projected")
query <- NormalizeMetacells(query, wgcna_name = "projected")
query <- SetDatExpr(query, group_name = type, group.by = "cell_type", assay = "RNA",
                    layer = "data", wgcna_name = "projected")
t_query_metacells <- proc.time()[["elapsed"]] - t0
write.table(GetDatExpr(query, wgcna_name = "projected"), file.path(out, "query_datExpr.tsv"),
            sep = "\t", quote = FALSE, col.names = NA)

t0 <- proc.time()[["elapsed"]]
query <- ModulePreservation(query, seurat_ref = ref, name = "pres", n_permutations = n_perm,
                            wgcna_name = "projected", wgcna_name_ref = "ref", verbose = 0)
t_preservation <- proc.time()[["elapsed"]] - t0
mp <- GetModulePreservation(query, "pres", wgcna_name = "projected")
write.table(format(mp$obs, digits = 17), file.path(out, "preservation_obs.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)
write.table(format(mp$Z, digits = 17), file.path(out, "preservation_Z.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)
nmod <- length(setdiff(unique(as.character(mods$module)), "grey"))
write.table(data.frame(ref_cells = ncol(ref), query_cells = ncol(query),
                       genes = length(GetWGCNAGenes(ref, "ref")), modules = nmod,
                       permutations = n_perm, t_read = t_read, t_ref = t_ref,
                       t_project = t_project, t_query_metacells = t_query_metacells,
                       t_preservation = t_preservation),
            file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d modules; reference %.1fs, projection %.1fs, query metacells %.1fs, preservation %.1fs\n",
            nmod, t_ref, t_project, t_query_metacells, t_preservation))
