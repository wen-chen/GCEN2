# hdWGCNA (v0.4.12) consensus network of one cell type over several datasets,
# with harmonized module eigengenes (the "Consensus network analysis" and
# ModuleEigengenes(group.by.vars) parts of the hdWGCNA tutorials), on the
# export of 01_prepare.py:
#   NormalizeData -> SetupForWGCNA(gene_select = "fraction", 0.05)
#   -> MetacellsByGroups(group.by = c(cell_type, CONSENSUS_BY); the scANVI
#      embedding, 30 dims, k 25, max_shared 10, up to 1000 metacells per group)
#   -> NormalizeMetacells -> SetMultiExpr(multi.group.by = CONSENSUS_BY)
#   -> TestSoftPowersConsensus -> ConstructNetwork(consensus = TRUE)
#   -> ScaleData -> ModuleEigengenes(group.by.vars = HARMONIZE_BY)
#   -> ModuleConnectivity (on the harmonized eigengenes, as hdWGCNA does)
#   -> FindAllDMEs(group.by = CONSENSUS_BY)
#   -> ModuleTraitCorrelation(traits = age, BMI, sex; group.by = CONSENSUS_BY).
#
# The harmonized eigengenes do not depend on SEED: Seurat's RunPCA (called
# per module by ModuleEigengenes) sets the seed to 42 just before Harmony.
#
# usage: Rscript 02_hdwgcna.R EXPORT_DIR OUT_DIR [SEED] [CONSENSUS_BY] [HARMONIZE_BY]
# Writes to OUT_DIR: genes.txt, sets.tsv (set, metacells, power),
# multiExpr_<set>.tsv (metacells x genes of each set), modules_network.tsv
# (after ConstructNetwork), MEs.tsv, hMEs.tsv, modules_final.tsv (after
# ModuleConnectivity, with kME), DMEs.tsv, trait_cor.tsv (group, trait,
# module, cor, p, fdr) and timing.tsv.
for (i in 1:10) {
  if (tryCatch({suppressPackageStartupMessages(library(enrichR)); TRUE},
               error = function(e) FALSE)) break
  Sys.sleep(5)
}
suppressPackageStartupMessages({
  library(Matrix); library(Seurat); library(WGCNA); library(hdWGCNA); library(harmony)
})
args <- commandArgs(trailingOnly = TRUE)
W <- args[1]; out <- args[2]
seed <- if (length(args) > 2) as.integer(args[3]) else 42L
cons_by <- if (length(args) > 3) args[4] else "dataset"
harm_by <- if (length(args) > 4) args[5] else "donor_id"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
enableWGCNAThreads(nThreads = 8)

t0 <- proc.time()[["elapsed"]]
counts <- Read10X_h5(file.path(W, "matrix.h5"), use.names = TRUE)
rownames(counts) <- readLines(file.path(W, "genes.tsv"))
meta <- read.delim(file.path(W, "meta.tsv"), row.names = 1, check.names = FALSE)
meta$sex <- factor(meta$sex)  # ModuleTraitCorrelation takes numeric or factor traits
rep <- as.matrix(read.delim(file.path(W, "rep.tsv"), row.names = 1, check.names = FALSE))
counts <- counts[, rownames(meta)]
so <- CreateSeuratObject(counts, meta.data = meta)
so <- NormalizeData(so, verbose = FALSE)
colnames(rep) <- paste0("rep_", seq_len(ncol(rep)))
so[["rep"]] <- CreateDimReducObject(embeddings = rep[colnames(so), ], key = "rep_", assay = "RNA")
type <- unique(meta$cell_type)
stopifnot(length(type) == 1)
t_read <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
so <- SetupForWGCNA(so, gene_select = "fraction", fraction = 0.05, wgcna_name = "cons")
# MetacellsByGroups subsets the object per group, which resets R's RNG: seed
# inside ConstructMetacells (see ../hdwgcna_pbmc/02_hdwgcna.R)
trace("ConstructMetacells", bquote(set.seed(.(seed))), print = FALSE,
      where = asNamespace("hdWGCNA"))
so <- MetacellsByGroups(so, group.by = c("cell_type", cons_by), ident.group = "cell_type",
                        reduction = "rep", dims = 1:30, k = 25, max_shared = 10,
                        min_cells = 100, target_metacells = 1000, mode = "average")
so <- NormalizeMetacells(so)
t_metacells <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
# the sets in sorted (C locale) order: ConstructNetwork takes the powers of
# TestSoftPowersConsensus in that order (dplyr::group_split) and pairs them
# with the sets by position, so with SetMultiExpr's default order (first
# appearance) a set can get another set's power
so <- SetMultiExpr(so, group_name = type, group.by = "cell_type", multi.group.by = cons_by,
                   multi_groups = sort(unique(as.character(meta[[cons_by]])), method = "radix"),
                   assay = "RNA")
so <- TestSoftPowersConsensus(so)
so <- ConstructNetwork(so, consensus = TRUE, tom_name = "cons",
                       tom_outdir = file.path(out, "TOM"), overwrite_tom = TRUE)
t_network <- proc.time()[["elapsed"]] - t0
unlink(file.path(out, "TOM"), recursive = TRUE)

multi <- GetMultiExpr(so)
powers <- GetWGCNAParams(so)$power
writeLines(GetWGCNAGenes(so), file.path(out, "genes.txt"))
write.table(data.frame(set = names(multi), metacells = sapply(multi, function(m) nrow(m$data)),
                       power = powers),
            file.path(out, "sets.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
for (s in names(multi)) {
  write.table(multi[[s]]$data, file.path(out, paste0("multiExpr_", gsub("[^A-Za-z0-9]+", "_", s),
                                                    ".tsv")),
              sep = "\t", quote = FALSE, col.names = NA)
}
m0 <- GetModules(so)
write.table(data.frame(gene = m0$gene_name, module = m0$module),
            file.path(out, "modules_network.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

t0 <- proc.time()[["elapsed"]]
so <- ScaleData(so, features = GetWGCNAGenes(so), verbose = FALSE)
set.seed(seed)
so <- ModuleEigengenes(so, group.by.vars = harm_by, verbose = FALSE)
so <- ModuleConnectivity(so, group.by = "cell_type", group_name = type)
t_me <- proc.time()[["elapsed"]] - t0

t0 <- proc.time()[["elapsed"]]
dmes <- FindAllDMEs(so, group.by = cons_by)
so <- ModuleTraitCorrelation(so, traits = c("age", "BMI", "sex"), group.by = cons_by)
t_stats <- proc.time()[["elapsed"]] - t0
write.table(format(dmes, digits = 17), file.path(out, "DMEs.tsv"), sep = "\t", quote = FALSE,
            col.names = NA)
mt <- GetModuleTraitCorrelation(so)
tc <- do.call(rbind, lapply(names(mt$cor), function(g) {
  r <- as.matrix(mt$cor[[g]]); p <- as.matrix(mt$pval[[g]]); f <- as.matrix(mt$fdr[[g]])
  f <- f[rownames(r), colnames(r), drop = FALSE]
  data.frame(group = g, trait = rep(rownames(r), ncol(r)),
             module = rep(colnames(r), each = nrow(r)), cor = as.vector(r),
             p = as.vector(p), fdr = as.vector(f))
}))
write.table(format(tc, digits = 17), file.path(out, "trait_cor.tsv"), sep = "\t", quote = FALSE,
            row.names = FALSE)

write.table(GetMEs(so, harmonized = FALSE), file.path(out, "MEs.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)
write.table(GetMEs(so, harmonized = TRUE), file.path(out, "hMEs.tsv"), sep = "\t",
            quote = FALSE, col.names = NA)
mf <- GetModules(so)
write.table(data.frame(gene = mf$gene_name, module = mf$module,
                       mf[, grep("^kME_", colnames(mf))], check.names = FALSE),
            file.path(out, "modules_final.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
nmod <- length(setdiff(unique(as.character(m0$module)), "grey"))
write.table(data.frame(cells = ncol(so), genes = length(GetWGCNAGenes(so)),
                       metacells = sum(sapply(multi, function(m) nrow(m$data))),
                       powers = paste(powers, collapse = ","), modules = nmod,
                       t_read = t_read, t_metacells = t_metacells, t_network = t_network,
                       t_me = t_me, t_stats = t_stats),
            file.path(out, "timing.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("%d cells, %d sets, powers %s, %d modules; metacells %.1fs, network %.1fs, MEs + kME %.1fs\n",
            ncol(so), length(multi), paste(powers, collapse = ","), nmod, t_metacells,
            t_network, t_me))

