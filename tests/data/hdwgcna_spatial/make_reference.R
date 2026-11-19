# hdWGCNA 0.4.12 reference for gcen2's metaspots and spot-level module
# eigengenes (tests/test_spatial_modules.py).
# usage: Rscript make_reference.R   (run in this directory)
#
# A synthetic Visium-like sample: hexagonal grid (row + col even), rows 0..19,
# cols 0..41, with holes; groups A (left), B (upper right), C (lower right) and
# D (a small block, < 50 spots: skipped); 30 genes, Poisson counts.
#   counts.tsv      spots x genes counts, with row, col, group
#   members.tsv     MetaspotsByGroups(group.by = "group"): metaspot (centre
#                   barcode), its neighbours in hdWGCNA's order
#   metaspots.tsv   NormalizeMetacells data (metaspots x genes)
#   spot_me.tsv     ComputeModuleEigengene's computation (ScaleData, RunPCA,
#                   sign by average expression) on the spots for modules
#                   m1 = g1..g10, m2 = g11..g20
suppressPackageStartupMessages({library(Matrix); library(Seurat); library(hdWGCNA)})
set.seed(7)
grid <- expand.grid(row = 0:19, col = 0:41)
grid <- grid[(grid$row + grid$col) %% 2 == 0, ]
grid <- grid[runif(nrow(grid)) > 0.08, ]
grid$group <- ifelse(grid$col < 18, "A", ifelse(grid$row < 11, "B", "C"))
grid$group[grid$col >= 36 & grid$row >= 15] <- "D"
n <- nrow(grid)
bc <- sprintf("S%04d-1", seq_len(n))
rownames(grid) <- bc
f1 <- sin(grid$col / 6) + rnorm(n, sd = 0.3)
f2 <- cos(grid$row / 4) + rnorm(n, sd = 0.3)
mu <- matrix(log(2), n, 30)
mu[, 1:10] <- mu[, 1:10] + 0.8 * f1
mu[, 11:20] <- mu[, 11:20] + 0.8 * f2
counts <- matrix(rpois(n * 30, exp(mu)), n, 30, dimnames = list(bc, sprintf("g%d", 1:30)))
write.table(cbind(grid, counts), "counts.tsv", sep = "\t", quote = FALSE, col.names = NA)

# Seurat object with a Visium image (MetaspotsByGroups needs one): a blank
# low-resolution image and the grid as tissue positions
d <- file.path(tempdir(), "spatial"); dir.create(d, showWarnings = FALSE)
png(file.path(d, "tissue_lowres_image.png"), 60, 60); grid::grid.newpage(); invisible(dev.off())
writeLines('{"tissue_hires_scalef": 0.1, "tissue_lowres_scalef": 0.03, "fiducial_diameter_fullres": 100, "spot_diameter_fullres": 60}',
           file.path(d, "scalefactors_json.json"))
write.table(data.frame(bc, 1, grid$row, grid$col, grid$row * 100 + 100, grid$col * 100 + 100),
            file.path(d, "tissue_positions_list.csv"), sep = ",", quote = FALSE,
            row.names = FALSE, col.names = FALSE)
so <- CreateSeuratObject(as(t(counts), "CsparseMatrix"), assay = "Spatial",
                         meta.data = data.frame(group = grid$group, row = grid$row,
                                                col = grid$col, imagerow = grid$row * 100 + 100,
                                                imagecol = grid$col * 100 + 100, row.names = bc))
img <- Read10X_Image(d, assay = "Spatial")
so[["slice1"]] <- img[colnames(so)]
so <- NormalizeData(so, verbose = FALSE)

# SeuratObject 5 made GetAssayData(slot = ) defunct; MetaspotsByGroups uses it
# for an input check only (see benchmarks/hdwgcna_spatial/02_hdwgcna.R)
src <- paste(deparse(hdWGCNA::MetaspotsByGroups), collapse = "\n")
old <- "GetAssayData(seurat_obj, assay = assay, \n            slot = slot)"
stopifnot(grepl(old, src, fixed = TRUE))
f <- eval(parse(text = sub(old, "SeuratObject::LayerData(seurat_obj, assay = assay, layer = slot)",
                           src, fixed = TRUE)))
environment(f) <- asNamespace("hdWGCNA")

so <- SetupForWGCNA(so, gene_select = "all", wgcna_name = "t")
so <- suppressWarnings(f(so, group.by = "group", ident.group = "group", assay = "Spatial"))
so <- NormalizeMetacells(so)
mo <- GetMetacellObject(so)
nb <- mo@misc$spot_neighbors
write.table(data.frame(metaspot = names(nb),
                       neighbours = vapply(nb, paste, "", collapse = ",")),
            "members.tsv", sep = "\t", quote = FALSE, row.names = FALSE)
x <- as.matrix(LayerData(mo, layer = "data"))
write.table(format(t(x), digits = 17), "metaspots.tsv", sep = "\t", quote = FALSE,
            col.names = NA)

# ComputeModuleEigengene's steps on the spots
me <- sapply(list(m1 = sprintf("g%d", 1:10), m2 = sprintf("g%d", 11:20)), function(g) {
  x <- LayerData(so, layer = "data")[g, ]
  cs <- CreateSeuratObject(x, assay = "Spatial")
  cs <- SetAssayData(cs, layer = "data", new.data = x)
  cs <- ScaleData(cs, features = g, verbose = FALSE)
  pc <- suppressWarnings(RunPCA(cs, features = g, verbose = FALSE))@reductions$pca@cell.embeddings[, 1]
  if (cor(Matrix::colMeans(x), pc) < 0) pc <- -pc
  pc
})
write.table(format(me, digits = 17), "spot_me.tsv", sep = "\t", quote = FALSE, col.names = NA)
cat(nrow(grid), "spots,", length(nb), "metaspots\n")
