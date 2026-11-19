# RandomWalkRestartMH reference scores for tests/test_multilayer.py.
# Rscript make_reference.R (from this directory; RandomWalkRestartMH, git
# commit 3d61974 of github.com/alberto-valdeolivas/RandomWalkRestartMH)
suppressPackageStartupMessages({library(RandomWalkRestartMH); library(igraph)})
seeds <- readLines("seeds.txt")
for (weighted in c(FALSE, TRUE)) {
  layers <- lapply(c(a = "a.tsv", b = "b.tsv", c = "c.tsv"), function(f) {
    e <- read.delim(f)
    if (!weighted) e$weight <- NULL
    graph_from_data_frame(e, directed = FALSE)
  })
  for (delta in c(0.5, 0.2)) {
    mp <- create.multiplex(layers)
    adj <- normalize.multiplex.adjacency(compute.adjacency.matrix(mp, delta = delta))
    res <- Random.Walk.Restart.Multiplex(adj, mp, seeds, r = 0.7, DispResults = "Alphabetic")
    out <- sprintf("%s_delta%s.tsv", if (weighted) "weighted" else "unweighted", delta)
    write.table(format(res$RWRM_Results, digits = 17), out, sep = "\t", quote = FALSE,
                row.names = FALSE)
  }
}
