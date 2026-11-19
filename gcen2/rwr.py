"""Random Walk with Restart (RWR) for gene function prediction."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from gcen2.core import extension
from gcen2.exceptions import GCENInputError

# pandas and scipy are imported inside predict(); predict_from_file (used by
# ``gcen2 rwr``) needs neither.
if TYPE_CHECKING:
    from gcen2.models import Network, RWRResult


def predict(
    network: Network,
    seed_genes: list[str],
    gamma: float = 0.5,
    prob_cutoff: float = 0.01,
    directed: bool = False,
    weighted_network: bool = False,
    weighted_genes: bool = False,
    gene_weights: dict[str, float] | None = None,
    max_iter: int = 1_000_000,
    tol: float = 1e-10,
) -> RWRResult:
    """Predict gene function using Random Walk with Restart.

    Parameters
    ----------
    network:
        Input co-expression network.
    seed_genes:
        List of known gene identifiers to use as seeds.
    gamma:
        Restart probability (0 < gamma < 1).
    prob_cutoff:
        Minimum probability to include a gene in results.
    directed:
        If True, treat the network as directed.
    weighted_network:
        If True, use correlation values as edge weights.
    weighted_genes:
        If True, use ``gene_weights`` as seed weights; otherwise uniform.
    gene_weights:
        Optional mapping of gene → weight when weighted_genes=True.
    max_iter:
        Maximum iteration count.
    tol:
        Convergence tolerance (L1 norm).
    """
    import pandas as pd
    from scipy.sparse import diags

    from gcen2.models import RWRResult

    _check_inputs(gamma, seed_genes)

    all_genes = sorted(set(network.edges["gene_a"]) | set(network.edges["gene_b"]))
    gene_idx = {g: i for i, g in enumerate(all_genes)}
    n = len(all_genes)

    # Sparse column-stochastic transition matrix. A dense n x n matrix needs
    # 8*n^2 bytes (6 GB for 27k genes); the sparse one scales with the edges.
    W = network.sparse_adjacency(all_genes, weighted=weighted_network, directed=directed)
    col_sums = np.asarray(W.sum(axis=0)).ravel()
    col_sums[col_sums == 0] = 1.0
    W = (W @ diags(1.0 / col_sums)).tocsr()

    P0 = _seed_vector(gene_idx, n, seed_genes,
                      gene_weights if weighted_genes and gene_weights else None)

    core = extension()
    if core is not None:
        PT, steps = core.rwr.random_walk_restart(
            W.indptr.astype(np.int64), W.indices.astype(np.int32), W.data, P0,
            gamma, max_iter, tol,
        )
    else:
        PT, steps = _rwr_python(W, P0, gamma, max_iter, tol)

    result_df = pd.DataFrame({"gene": all_genes, "probability": PT})
    result_df = result_df[result_df["gene"].isin(set(all_genes) - set(seed_genes))]
    result_df = result_df[result_df["probability"] >= prob_cutoff]
    # ties by gene name: the same order with every numpy / pandas sort algorithm
    result_df = result_df.sort_values(["probability", "gene"], ascending=[False, True],
                                      kind="mergesort").reset_index(drop=True)

    return RWRResult(ranked_genes=result_df, seed_genes=seed_genes, gamma=gamma, steps=steps)


def _check_inputs(gamma: float, seed_genes) -> None:
    if not 0.0 < gamma < 1.0:
        raise GCENInputError(f"gamma must be in (0, 1), got {gamma}.")
    if not seed_genes:
        raise GCENInputError("seed_genes must not be empty.")


def _seed_vector(gene_idx: dict, n: int, seed_genes, gene_weights=None) -> np.ndarray:
    """Restart vector: weight 1 (or ``gene_weights[g]``) per seed, normalised to sum 1."""
    P0 = np.zeros(n, dtype=np.float64)
    for g in seed_genes:
        if g in gene_idx:
            P0[gene_idx[g]] = gene_weights.get(g, 1.0) if gene_weights else 1.0
    if P0.sum() == 0:
        raise GCENInputError("None of the seed genes are present in the network.")
    return P0 / P0.sum()


def _rank_desc(prob: np.ndarray, genes) -> np.ndarray:
    """Indices sorting ``prob`` descending and ties by gene name (as
    :func:`predict`): the same order with every numpy version (a quicksort
    orders ties as its implementation happens to)."""
    return np.lexsort((np.asarray(genes, dtype=str), -np.asarray(prob)))


def predict_from_file(
    path,
    seed_genes: list[str],
    gamma: float = 0.5,
    prob_cutoff: float = 0.01,
    max_iter: int = 1_000_000,
    tol: float = 1e-10,
) -> tuple[list[str], np.ndarray]:
    """RWR on an edge-list file without pandas or scipy (used by ``gcen2 rwr``).

    Same result as ``predict(read_network(path), seed_genes, gamma,
    prob_cutoff).ranked_genes`` (undirected, unweighted), returned as
    ``(genes, probabilities)`` in the same order. The C++ extension reads
    the file and builds the adjacency (parallel edges counted, as in
    :meth:`Network.sparse_adjacency`); without it this falls back to
    :func:`predict`.
    """
    _check_inputs(gamma, seed_genes)
    core = extension()
    if core is None:
        from gcen2.io import read_network

        df = predict(read_network(path), seed_genes, gamma, prob_cutoff,
                     max_iter=max_iter, tol=tol).ranked_genes
        return df["gene"].tolist(), df["probability"].to_numpy()

    from gcen2.module import read_file_adjacency

    # Rows in descending column order: scipy's ``W @ diags(...)`` in predict()
    # stores them that way, and the same order keeps the floating-point sums
    # (and so the results) bit-identical.
    indptr, indices, names, counts = read_file_adjacency(path, with_counts=True,
                                                         descending=True)
    n = len(names)
    P0 = _seed_vector({g: i for i, g in enumerate(names)}, n, seed_genes)
    if hasattr(core.rwr, "random_walk_restart_counts"):
        # W = A / column sums is formed inside the iteration (bit-identical,
        # without a float64 copy of the matrix)
        PT, _ = core.rwr.random_walk_restart_counts(indptr, indices, counts, P0, gamma,
                                                    max_iter, tol)
    else:
        # column-stochastic W, computed as scipy's ``W @ diags(1 / col_sums)``
        col_sums = np.bincount(indices, weights=counts, minlength=n)
        col_sums[col_sums == 0] = 1.0
        data = (1.0 / col_sums)[indices]
        data *= counts
        del counts
        PT, _ = core.rwr.random_walk_restart(indptr, indices, data, P0, gamma, max_iter, tol)
    PT = np.asarray(PT)
    seeds = set(seed_genes)
    keep = np.flatnonzero(np.fromiter((g not in seeds for g in names), dtype=bool, count=n)
                          & (prob_cutoff <= PT))
    order = keep[_rank_desc(PT[keep], [names[i] for i in keep])]
    return [names[i] for i in order], PT[order]


def _rwr_python(
    W: np.ndarray,
    P0: np.ndarray,
    gamma: float,
    max_iter: int,
    tol: float,
) -> tuple[np.ndarray, int]:
    """RWR iteration ``PT <- (1-gamma) W PT + gamma P0`` until the L1 change < tol.

    ``W`` may be a dense array or a scipy sparse matrix (column-stochastic).
    Same update rule and convergence test as the C++ ``random_walk_restart``;
    used as the fallback when the extension is not built.
    """
    PT = P0.copy()
    for step in range(max_iter):
        PT_new = (1 - gamma) * W @ PT + gamma * P0
        if np.sum(np.abs(PT_new - PT)) < tol:
            return PT_new, step + 1
        PT = PT_new
    return PT, max_iter
