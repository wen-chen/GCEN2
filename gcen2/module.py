"""Module identification via Jaccard similarity + DFS clustering."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from gcen2._matrix_io import require_file
from gcen2.exceptions import GCENInputError

# pandas, scipy and gcen2.models are imported where needed, so that
# ``gcen2 module`` (modules_from_file) runs without them.
if TYPE_CHECKING:
    from gcen2.models import ModuleSet, Network


def identify(network: Network, similarity: float = 0.5, threads: int = 2) -> ModuleSet:
    """Identify gene modules from a co-expression network.

    Uses Jaccard similarity-based module detection with DFS clustering,
    identical to the algorithm in GCEN v1 ``module_identify.cpp``.

    Parameters
    ----------
    network:
        Input co-expression network.
    similarity:
        Jaccard similarity cutoff. Gene pairs above this threshold are
        considered to belong to the same module.
    threads:
        Worker threads for the C++ similarity kernel.

    Returns
    -------
    ModuleSet with singleton modules (< 2 genes) excluded.
    """
    from gcen2.models import ModuleSet

    if network.edges.empty:
        return ModuleSet(modules=[], source_network=network, similarity_cutoff=similarity)

    _check_similarity(similarity)

    from gcen2.core import extension

    core = extension()
    if core is not None:
        raw_modules = _jaccard_modules_cpp(network, similarity, threads, core)
    else:
        raw_modules = _jaccard_modules(network, similarity)

    modules = [m for m in raw_modules if len(m) >= 2]
    return ModuleSet(modules=modules, source_network=network, similarity_cutoff=similarity)


def _check_similarity(similarity: float) -> None:
    if not 0.0 < similarity < 1.0:
        raise GCENInputError(f"similarity must be in (0, 1), got {similarity}.")


def read_file_adjacency(path, with_counts: bool = False, descending: bool = False):
    """Symmetric CSR adjacency of an edge-list file, read by the C++ extension.

    Returns ``(indptr, indices, names)`` (plus ``counts``, the entry
    multiplicities, with ``with_counts``): nodes sorted by name, duplicate
    and reversed edges merged, a self-loop on the diagonal; columns ascending
    within each row (descending with ``descending``). The header is
    detected as in :func:`gcen2.io.read_network`. Needs the compiled
    extension; imports neither pandas nor scipy.
    """
    from gcen2._matrix_io import network_has_header
    from gcen2.core import require_extension
    from gcen2.exceptions import GCENIOError

    path = require_file(path)
    core = require_extension()
    try:
        return core.module.read_edge_adjacency(
            str(path), network_has_header(path), with_counts, descending)
    except (RuntimeError, UnicodeDecodeError) as exc:
        raise GCENIOError(f"Failed to read network file {path}: {exc}") from exc


def modules_from_file(path, similarity: float = 0.5, threads: int = 2) -> list[set[str]]:
    """Modules of the network in an edge-list file, read without pandas.

    Same modules, in the same order, as
    ``identify(read_network(path), similarity, threads).modules``, but the
    C++ extension reads the two gene columns and builds the adjacency
    directly (``_gcen_core.module.read_edge_adjacency``), so neither pandas
    nor scipy is imported. Used by ``gcen2 module``. Without the compiled
    extension it falls back to :func:`identify`.
    """
    _check_similarity(similarity)
    from gcen2.core import extension

    core = extension()
    if core is None:
        from gcen2.io import read_network

        return identify(read_network(path, genes_only=True), similarity, threads).modules
    indptr, indices, names = read_file_adjacency(path)
    if not names:
        return []
    labels = core.module.jaccard_modules(indptr, indices, similarity, max(1, int(threads)))
    modules = _labels_to_modules(np.asarray(labels), np.array(names, dtype=object))
    return [m for m in modules if len(m) >= 2]


def _unique_adjacency(network: Network):
    """Symmetric 0/1 CSR adjacency (duplicates and reverse edges collapsed).

    Returns ``(A, genes)``; a self-loop is kept as a diagonal entry. ``A``
    has int32 indices and bool data. Gene columns may be object or
    categorical (``read_network(..., genes_only=True)``); categorical columns
    avoid hashing every gene name again.
    """
    import pandas as pd
    from scipy.sparse import coo_matrix

    edges = network.edges
    vocab = pd.Index(pd.Categorical(edges["gene_a"]).categories).union(
        pd.Index(pd.Categorical(edges["gene_b"]).categories))
    ca = pd.Categorical(edges["gene_a"], categories=vocab).codes.astype(np.int32)
    cb = pd.Categorical(edges["gene_b"], categories=vocab).codes.astype(np.int32)
    n = len(vocab)
    rows = np.concatenate([ca, cb])
    cols = np.concatenate([cb, ca])
    del ca, cb
    # tocsr() merges duplicate (i, j) entries (logical or for bool data)
    A = coo_matrix((np.ones(len(rows), dtype=bool), (rows, cols)), shape=(n, n)).tocsr()
    del rows, cols
    A.sum_duplicates()
    A.data[:] = 1
    A.sort_indices()
    return A, vocab.to_numpy(dtype=object)


def _labels_to_modules(labels: np.ndarray, genes) -> list[set[str]]:
    """Group nodes by label (label < 0 = not in any module), sorted deterministically."""
    nodes = np.flatnonzero(labels >= 0)
    if nodes.size == 0:
        return []
    order = np.argsort(labels[nodes], kind="stable")
    nodes = nodes[order]
    bounds = np.flatnonzero(np.diff(labels[nodes])) + 1
    modules = [set(genes[g].tolist()) for g in np.split(nodes, bounds)]
    # deterministic order: largest first, ties by smallest gene name
    modules.sort(key=lambda mod: (-len(mod), min(mod)))
    return modules


def _jaccard_modules_cpp(network: Network, cutoff: float, threads: int, core) -> list[set[str]]:
    """C++ kernel (``_gcen_core.module.jaccard_modules``) on the CSR adjacency."""
    A, genes = _unique_adjacency(network)
    labels = core.module.jaccard_modules(
        A.indptr.astype(np.int64, copy=False), A.indices.astype(np.int32, copy=False),
        cutoff, max(1, int(threads))
    )
    return _labels_to_modules(np.asarray(labels), genes)


def _jaccard_modules(network: Network, cutoff: float, chunk_size: int = 500) -> list[set[str]]:
    """Connected components of the Jaccard similarity graph (GCEN v1 algorithm).

    For genes a != b with neighbour sets N(a), N(b) (a self-loop puts a gene in
    its own set) and a_b = 1 if a and b are adjacent:

        J(a, b) = (|N(a) & N(b)| + a_b) / (|N(a)| + |N(b)| - |N(a) & N(b)| - a_b)

    Pairs with J > cutoff are linked and every connected component of those
    links is a module. Since cutoff > 0, only pairs that share a neighbour or
    are adjacent can qualify, i.e. the non-zeros of A @ A + A. Those are
    enumerated in row chunks of the sparse adjacency matrix instead of testing
    all n^2/2 gene pairs (the C++ all-pairs loop needed > 7 min for 30k genes).
    """
    from scipy.sparse import coo_matrix
    from scipy.sparse.csgraph import connected_components

    A, genes = _unique_adjacency(network)
    A = A.astype(np.int64)  # values are packed below
    n = len(genes)
    deg = np.asarray(A.sum(axis=1)).ravel()

    big = np.int64(1) << 32  # pack (shared, adjacent) into one sparse value
    src, dst = [], []
    for start in range(0, n, chunk_size):
        R = A[start:start + chunk_size]
        C = (R @ A + R * big).tocoo()
        i = C.row.astype(np.int64) + start
        j = C.col.astype(np.int64)
        upper = j > i
        i, j, v = i[upper], j[upper], C.data[upper]
        shared = (v % big).astype(np.float64)
        ab = (v // big).astype(np.float64)
        denom = deg[i] + deg[j] - shared - ab
        with np.errstate(divide="ignore", invalid="ignore"):
            sim = (shared + ab) / denom
        keep = (denom > 0) & (sim > cutoff)
        src.append(i[keep])
        dst.append(j[keep])

    src = np.concatenate(src) if src else np.empty(0, dtype=np.int64)
    dst = np.concatenate(dst) if dst else np.empty(0, dtype=np.int64)
    if src.size == 0:
        return []
    S = coo_matrix((np.ones(src.size), (src, dst)), shape=(n, n))
    _, labels = connected_components(S, directed=False)
    in_sim = np.zeros(n, dtype=bool)
    in_sim[src] = True
    in_sim[dst] = True
    return _labels_to_modules(np.where(in_sim, labels, -1), genes)
