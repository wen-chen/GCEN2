# RandomWalkRestartMH: random walk with restart on the multiplex network
# (RWR-M) of the layers, unweighted, global geometric-mean score per gene.
# Usage: Rscript 02_rwrmh.R <seeds.txt> <out.tsv> <r> <delta> <layer.tsv>...
suppressPackageStartupMessages({
  library(RandomWalkRestartMH)
  library(igraph)
  library(data.table)
})
args <- commandArgs(trailingOnly = TRUE)
seeds <- readLines(args[1])
r <- as.numeric(args[3])
delta <- as.numeric(args[4])
files <- args[-(1:4)]

t0 <- Sys.time()
layers <- lapply(files, function(f) {
  e <- fread(f, sep = "\t", select = 1:2)
  graph_from_data_frame(e, directed = FALSE)  # no weight attribute: unweighted
})
names(layers) <- sub("\\.tsv$", "", basename(files))
mp <- create.multiplex(layers)
adj <- normalize.multiplex.adjacency(compute.adjacency.matrix(mp, delta = delta))
seeds <- seeds[seeds %in% mp$Pool_of_Nodes]
res <- Random.Walk.Restart.Multiplex(adj, mp, seeds, r = r, MeanType = "Geometric",
                                     DispResults = "Alphabetic")
fwrite(res$RWRM_Results, args[2], sep = "\t")
cat(sprintf("RandomWalkRestartMH %s: %d genes, %d layers, %d seeds, %.1f s\n",
            packageVersion("RandomWalkRestartMH"), mp$Number_of_Nodes_Multiplex,
            length(layers), length(seeds), as.numeric(Sys.time() - t0, units = "secs")))
