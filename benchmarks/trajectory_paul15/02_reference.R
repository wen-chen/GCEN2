# Co-expression networks in sliding pseudotime windows, in R with WGCNA's
# correlation functions, following the definition of gcen2 trajectory:
# cells ordered by pseudotime; windows of w cells every s cells (a last
# window ends at the last cell); consecutive cells averaged into metacells;
# log2(x + 1); Pearson r, Student p, Benjamini-Hochberg over all gene pairs;
# an edge has |r| >= cutoff and p, FDR <= cutoff. Then every edge of any
# window: its r in each window, the windows holding it, its pattern (stable,
# gained, lost, transient, intermittent), last - first r and the correlation
# of its r with the window centres (median pseudotime). There is no published
# package for this analysis; this is the straightforward R implementation,
# used to check GCEN2 and to time it.
# Usage: Rscript 02_reference.R <expr.tsv> <pt.tsv> <outdir> <window fraction>
#        <metacell> <cor> <pval> [threads]
suppressPackageStartupMessages({
  library(WGCNA)
  library(data.table)
})
a <- commandArgs(trailingOnly = TRUE)
out <- a[3]
frac <- as.numeric(a[4]); meta <- as.integer(a[5])
cut_r <- as.numeric(a[6]); cut_p <- as.numeric(a[7])
threads <- if (length(a) >= 8) as.integer(a[8]) else 8
allowWGCNAThreads(threads)
setDTthreads(threads)
dir.create(out, showWarnings = FALSE, recursive = TRUE)

t0 <- Sys.time()
d <- fread(a[1], sep = "\t", data.table = FALSE)
genes <- d[[1]]
x <- t(as.matrix(d[, -1]))               # cells x genes
pt <- fread(a[2], sep = "\t", data.table = FALSE)
time <- pt[[2]][match(rownames(x) <- colnames(d)[-1], pt[[1]])]
keep <- is.finite(time)
x <- x[keep, , drop = FALSE]; time <- time[keep]
o <- order(time, method = "radix")      # stable, as numpy's argsort(kind="stable")
x <- x[o, , drop = FALSE]; time <- time[o]
n <- length(time)
w <- max(1L, as.integer(round(frac * n)))
s <- max(1L, w %/% 2L)
starts <- seq(0L, n - w, by = s)
if (tail(starts, 1) + w < n) starts <- c(starts, n - w)

G <- length(genes)
ut <- which(upper.tri(diag(G)), arr.ind = TRUE)    # every pair, gene i < gene j
ut <- ut[order(ut[, 1], ut[, 2]), ]
n_pairs <- nrow(ut)
cors <- matrix(0, n_pairs, length(starts))
present <- matrix(FALSE, n_pairs, length(starts))
centers <- numeric(length(starts))
for (k in seq_along(starts)) {
  cells <- starts[k] + seq_len(w)
  block <- x[cells, , drop = FALSE]
  grp <- (seq_len(w) - 1L) %/% meta
  pooled <- rowsum(block, grp) / as.vector(table(grp))  # metacells x genes
  r <- WGCNA::cor(log2(pooled + 1))
  r[is.na(r)] <- 0
  rv <- r[ut]
  p <- corPvalueStudent(rv, nrow(pooled))
  q <- p.adjust(p, "BH", n = n_pairs)
  edge <- abs(rv) >= cut_r & p <= cut_p & q <= cut_p
  edge[is.na(edge)] <- FALSE
  cors[, k] <- pmin(pmax(rv, -1), 1)
  present[, k] <- edge
  centers[k] <- median(time[cells])
  fwrite(data.table(gene_a = genes[ut[edge, 1]], gene_b = genes[ut[edge, 2]],
                    correlation = rv[edge], p_value = p[edge], fdr = q[edge]),
         file.path(out, sprintf("window_%d.network", k - 1L)), sep = "\t")
}
u <- which(rowSums(present) > 0)
pr <- present[u, , drop = FALSE]
cr <- cors[u, , drop = FALSE]
L <- ncol(pr)
first <- max.col(pr, ties.method = "first") - 1L
last <- L - max.col(pr[, L:1, drop = FALSE], ties.method = "first")
pattern <- ifelse(rowSums(pr) == L, "stable",
           ifelse(!pr[, 1] & pr[, L], "gained",
           ifelse(pr[, 1] & !pr[, L], "lost",
           ifelse(!pr[, 1] & !pr[, L], "transient", "intermittent"))))
trend <- if (L >= 3) as.vector(cor(t(cr), centers)) else NA_real_
dyn <- data.table(gene_a = genes[ut[u, 1]], gene_b = genes[ut[u, 2]], cr,
                  n_windows = rowSums(pr), first = first, last = last, pattern = pattern,
                  cor_change = cr[, L] - cr[, 1], trend = trend)
setnames(dyn, 3:(2 + L), paste0("cor_", seq_len(L) - 1L))
fwrite(dyn, file.path(out, "edge_dynamics.tsv"), sep = "\t")
cat(sprintf("R reference: %d cells, %d genes, %d windows of %d cells, %d edges, %.1f s\n",
            n, G, L, w, nrow(dyn), as.numeric(Sys.time() - t0, units = "secs")))
