# Weighted network fusion in R (data.table), following the definition of
# gcen2 fuse with --normalize max and equal layer weights: every layer's
# weights divided by its largest weight; per undirected edge, the mean over
# the layers (0 where absent) and the number of layers holding it. There is
# no published R package for this definition; this is the straightforward
# implementation used to check GCEN2 and to time it.
# Usage: Rscript 03_fuse.R <out.tsv> <layer.tsv>...
suppressPackageStartupMessages(library(data.table))
args <- commandArgs(trailingOnly = TRUE)
files <- args[-1]
setDTthreads(8)

t0 <- Sys.time()
edges <- rbindlist(lapply(seq_along(files), function(k) {
  e <- fread(files[k], sep = "\t", select = 1:3, col.names = c("a", "b", "w"))
  e[, w := abs(w)]
  e[, `:=`(gene_a = pmin(a, b), gene_b = pmax(a, b))]
  e <- e[gene_a != gene_b, .(w = max(w)), by = .(gene_a, gene_b)]
  e[, `:=`(w = w / max(w), layer = k)]
}))
fused <- edges[, .(score = sum(w) / length(files), n_layers = .N), by = .(gene_a, gene_b)]
setorder(fused, -score, gene_a, gene_b)
fwrite(fused, args[1], sep = "\t")
cat(sprintf("R data.table fusion: %d layers, %d edges, %.1f s\n", length(files),
            nrow(fused), as.numeric(Sys.time() - t0, units = "secs")))
