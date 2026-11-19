"""Metacells: groups of transcriptionally similar cells aggregated into one
profile, the usual input of single-cell co-expression networks (hdWGCNA
``MetacellsByGroups``, Morabito et al., Cell Rep Methods 2023).

Single cells are too sparse for stable gene-gene correlations; averaging
each cell with its nearest neighbours removes most dropout noise while
keeping the biological variation between cell states. Within each group
(cell type, sample, ...):

1. every cell's ``k`` nearest neighbours (itself included) are found in a
   reduced space (``adata.obsm[use_rep]``, ``X_pca`` when present, else a
   PCA computed here);
2. cells are visited in random order and a cell's neighbourhood becomes a
   metacell unless it shares more than ``max_shared`` cells with a metacell
   already chosen, until ``max_metacells`` are chosen;
3. each metacell's expression is the mean (or sum) of its ``k`` cells, then
   (for counts) scaled to 10,000 and log1p transformed.
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np

from gcen2._anndata import is_anndata
from gcen2.exceptions import GCENInputError

AGGREGATE = ("mean", "sum")


def metacells(
    adata: Any,
    groupby: str | list[str] | None = None,
    k: int = 25,
    max_shared: int | None = None,
    max_metacells: int = 1000,
    min_cells: int = 100,
    use_rep: str | None = None,
    n_pcs: int = 30,
    layer: str | None = None,
    aggregate: str = "mean",
    normalize: bool | None = None,
    seed: int = 42,
    threads: int = 1,
    groups: Any = None,
) -> Any:
    """Aggregate cells into metacells, separately within each group.

    Parameters
    ----------
    adata:
        AnnData (cells x genes, dense or sparse).
    groupby:
        ``adata.obs`` column(s) whose combinations define the groups (e.g.
        cell type, or ``["cell_type", "sample"]``); no grouping when None.
        Metacells never mix cells of different groups.
    k:
        Cells per metacell (the cell and its k - 1 nearest neighbours).
    max_shared:
        Maximum number of cells a metacell may share with another one
        (default ``min(10, k - 1)``; hdWGCNA uses 10 with ``k = 25``).
    max_metacells:
        Maximum metacells per group.
    min_cells:
        Groups with fewer cells are skipped (with a warning).
    use_rep:
        ``adata.obsm`` key of the space for the neighbours. Default:
        ``X_pca`` when present, else a PCA of the expression (counts are
        first scaled to 10,000 per cell and log1p transformed; values that
        are not integer counts are used as they are), on the highly
        variable genes when ``adata.var["highly_variable"]`` exists.
    n_pcs:
        Components used (of ``use_rep`` or of the PCA computed here).
    layer:
        Expression aggregated (and used for the PCA): ``adata.layers[layer]``
        instead of ``X``.
    aggregate:
        ``mean`` or ``sum`` of the member cells (the same correlations).
    normalize:
        Scale each metacell to 10,000 in total and log1p transform it (as
        hdWGCNA ``NormalizeMetacells``), so that sequencing depth does not
        correlate all genes. Default (None): when the expression holds
        integer counts.
    seed:
        Random order in which cells are tried as metacell centres.
    threads:
        Threads for the neighbour search.
    groups:
        Only these groups (values of a single ``groupby`` column; the PCA,
        if computed, then uses only their cells). The same result as
        subsetting ``adata`` to them first, without copying the expression.

    Note
    ----
    Metacells overlap and each averages ``k`` cells: correlations computed
    on them are fine, but p-values that count every metacell as an
    independent sample are far too small. :func:`build_by_celltype` computes
    them on every gene pair's effective samples from the overlap of the
    metacells (:func:`gcen2.single_cell.metacell_effective_samples`, which
    takes ``uns["metacell_members"]``); when building from metacells
    directly, prefer correlation cutoffs (or ``pval`` with that in mind).

    Returns
    -------
    AnnData
        One row per metacell (dense ``X``; log-normalised when
        ``normalize``), the genes of ``adata``. ``obs``
        holds the group column(s), ``n_cells`` and the centre cell;
        ``uns["metacell_members"]`` the member cells' positions in ``adata``
        (one row of ``k`` per metacell).
    """
    import anndata as ad
    import pandas as pd

    plan = _metacell_plan(adata, groupby, k, max_shared, max_metacells, min_cells, use_rep,
                          n_pcs, layer, seed, threads, groups, aggregate)
    X, keys, member_idx = plan["X"], plan["keys"], plan["members"]
    if normalize is None:
        normalize = _is_counts(X)
    profiles = _metacell_profiles(X, member_idx, aggregate, normalize)
    n = len(member_idx)
    obs = pd.DataFrame(plan["groups"], columns=keys) if keys else pd.DataFrame(index=range(n))
    obs["n_cells"] = k
    obs["center"] = np.asarray(adata.obs_names)[member_idx[:, 0]]
    obs.index = [f"metacell_{i}" for i in range(n)]
    out = ad.AnnData(X=profiles, obs=obs, var=adata.var.copy())
    out.uns["metacell_members"] = member_idx
    out.uns["metacell_params"] = {"k": k, "max_shared": plan["max_shared"],
                                  "aggregate": aggregate, "normalize": bool(normalize),
                                  "groupby": keys, "seed": seed}
    return out


def _metacell_plan(adata, groupby, k, max_shared, max_metacells, min_cells, use_rep, n_pcs,
                   layer, seed, threads, groups, aggregate="mean") -> dict:
    """The members of every metacell (no aggregation yet): ``X`` (the
    expression), ``keys`` (groupby columns), ``members`` (metacells x k
    positions in ``adata``, centre first), ``groups`` (the group of each
    metacell, a tuple) and ``max_shared``. Groups in order of first
    appearance; one random generator over them. ``groups``: only the cells
    whose first groupby column is one of these."""
    import pandas as pd

    from gcen2._knn import knn_indices

    if not is_anndata(adata):
        raise GCENInputError("metacells needs an AnnData object.")
    if aggregate not in AGGREGATE:
        raise GCENInputError(f"aggregate must be one of {AGGREGATE}, got '{aggregate}'.")
    if k < 2:
        raise GCENInputError(f"k must be >= 2, got {k}.")
    if max_shared is None:
        max_shared = min(10, k - 1)
    if not 0 <= max_shared < k:
        raise GCENInputError(f"max_shared must be in [0, k), got {max_shared}.")
    if max_metacells < 1:
        raise GCENInputError(f"max_metacells must be >= 1, got {max_metacells}.")
    X = _expression(adata, layer)

    keys = [] if groupby is None else [groupby] if isinstance(groupby, str) else list(groupby)
    for key in keys:
        if key not in adata.obs.columns:
            raise GCENInputError(f"groupby key '{key}' not found in adata.obs.")
    pool = None  # positions of the cells of `groups` (None: all cells)
    if groups is not None:
        if not keys:
            raise GCENInputError("groups needs a groupby column.")
        wanted = {str(g) for g in ([groups] if isinstance(groups, str) else groups)}
        pool = np.flatnonzero(adata.obs[keys[0]].astype(str).isin(wanted).to_numpy())
    if keys:
        labels = adata.obs[keys].astype(str)
        if pool is not None:
            labels = labels.iloc[pool]
        group_codes, group_values = pd.MultiIndex.from_frame(labels).factorize()
    else:
        group_codes = np.zeros(adata.n_obs, dtype=np.int64)
        group_values = [()]

    rep = _representation(adata, X, use_rep, n_pcs, seed, rows=pool)
    rng = np.random.default_rng(seed)
    members, names = [], []
    for g, value in enumerate(group_values):
        cells = np.flatnonzero(group_codes == g)
        name = value if isinstance(value, tuple) else (value,)
        if len(cells) < max(min_cells, k):
            warnings.warn(f"Group {dict(zip(keys, name)) or 'all'} has {len(cells)} cells "
                          f"(< {max(min_cells, k)}); skipped.", stacklevel=3)
            continue
        # the cell itself first, equidistant cells by index (the same
        # neighbours on every platform)
        nn = knn_indices(rep[cells], k, workers=threads)
        chosen = _select(nn, max_shared, max_metacells, rng)
        pos = cells if pool is None else pool[cells]  # positions in adata
        members.extend(pos[nn[c]] for c in chosen)
        names.extend([name] * len(chosen))
    if not members:
        raise GCENInputError("No group has enough cells for a metacell.")
    return {"X": X, "keys": keys, "members": np.vstack(members), "groups": names,
            "max_shared": max_shared}


def _metacell_profiles(X, member_idx: np.ndarray, aggregate: str, normalize: bool) -> np.ndarray:
    """Dense profiles (metacells x genes) of the members; row by row, so any
    subset of the metacells gives the same rows."""
    if not normalize:
        return _aggregate(X, member_idx, aggregate)
    # Normalised to 10,000 per metacell, which is the same for the mean and the
    # sum of the members. From the sums, x * 1e4 / total rounds once (x * 1e4 is
    # exact for counts), so equal ratios give equal values: Spearman ranks of
    # the result then do not depend on the platform's log1p (CPUs differ in its
    # last bit, which reordered near-equal values).
    sums = _aggregate(X, member_idx, "sum")
    totals = sums.sum(axis=1, keepdims=True)
    safe = np.where(totals > 0, totals, 1.0)
    return np.log1p(np.where(totals > 0, sums * 1e4 / safe, 0.0))


def _expression(adata, layer):
    if layer is not None:
        if layer not in adata.layers:
            raise GCENInputError(f"Layer '{layer}' not found in adata.layers.")
        return adata.layers[layer]
    if adata.X is None:
        raise GCENInputError("adata.X is empty; give a layer.")
    return adata.X


def _is_backed(X) -> bool:
    """An on-disk matrix of an AnnData opened with ``backed="r"``."""
    import scipy.sparse as sp

    return not (sp.issparse(X) or isinstance(X, np.ndarray)) and hasattr(X, "format") \
        and hasattr(X, "shape")


def _rows(X, idx=None):
    """Rows ``idx`` (sorted positions; all if None) of X in memory: a copy
    of those rows only, also from a backed (on-disk) matrix."""
    if _is_backed(X):
        if X.format != "csr":
            raise GCENInputError("A backed AnnData needs its X stored as CSR (rows "
                                 "are read per cell type); load it into memory instead.")
        return X[:] if idx is None else X[np.asarray(idx)]
    return X if idx is None else X[idx]


def _is_counts(X) -> bool:
    """Non-negative integers (raw counts). Checked in chunks, so that a large
    matrix needs no full-size temporaries (a backed matrix is read in
    blocks of rows)."""
    import scipy.sparse as sp

    if _is_backed(X):
        step = 20000
        return X.shape[0] > 0 and all(_is_counts(X[i:i + step])
                                      for i in range(0, X.shape[0], step))
    values = X.data if sp.issparse(X) else np.asarray(X)
    values = values.reshape(-1)
    if not values.size:
        return False
    step = 1 << 22
    for i in range(0, values.size, step):
        v = values[i:i + step]
        if not (np.all(v >= 0) and np.all(v == np.round(v))):
            return False
    return True


def _representation(adata, X, use_rep, n_pcs, seed, rows=None) -> np.ndarray:
    """The neighbour space of the cells ``rows`` (all if None)."""
    if use_rep is None and "X_pca" in adata.obsm:
        use_rep = "X_pca"
    if use_rep is not None:
        if use_rep not in adata.obsm:
            raise GCENInputError(f"'{use_rep}' not found in adata.obsm.")
        rep = np.asarray(adata.obsm[use_rep])
        if rows is not None:
            rep = rep[rows]
        return np.ascontiguousarray(rep[:, :n_pcs], dtype=np.float64)
    return _pca(adata, _rows(X, rows), n_pcs, seed)


def _pca(adata, X, n_pcs, seed) -> np.ndarray:
    """PCA scores of the log-normalised expression (see :func:`metacells`)."""
    import scipy.sparse as sp
    from scipy.sparse.linalg import LinearOperator, svds

    if "highly_variable" in adata.var.columns:
        hv = adata.var["highly_variable"].to_numpy(dtype=bool)
        if hv.any():
            X = X[:, hv]
    # one code path for dense and sparse input (same neighbours either way)
    X = sp.csr_matrix(X, dtype=np.float64)
    X.eliminate_zeros()
    if _is_counts(X):
        totals = np.asarray(X.sum(axis=1)).ravel()
        scale = np.where(totals > 0, 1e4 / np.where(totals > 0, totals, 1.0), 0.0)
        X = sp.csr_matrix(sp.diags(scale) @ X)
        X.data = np.log1p(X.data)
    n_cells, n_genes = X.shape
    n_comp = max(1, min(int(n_pcs), n_cells - 1, n_genes - 1))
    mean = np.asarray(X.mean(axis=0)).ravel()
    ones = np.ones(n_cells)
    # centred implicitly, X stays sparse
    op = LinearOperator(
        (n_cells, n_genes), dtype=np.float64,
        matvec=lambda v: X @ v - ones * (mean @ v),
        rmatvec=lambda v: X.T @ v - mean * v.sum(),
        matmat=lambda V: X @ V - np.outer(ones, mean @ V),
        rmatmat=lambda V: X.T @ V - np.outer(mean, V.sum(axis=0)))
    v0 = np.random.default_rng(seed).uniform(-1, 1, min(n_cells, n_genes))
    u, s, _ = svds(op, k=n_comp, v0=v0)
    order = np.argsort(-s)
    scores = u[:, order] * s[order]
    # fix the sign of each component (largest |loading| positive) for
    # reproducible neighbours across platforms
    signs = np.sign(scores[np.abs(scores).argmax(axis=0), np.arange(scores.shape[1])])
    signs[signs == 0] = 1
    return np.ascontiguousarray(scores * signs)


def _select(nn: np.ndarray, max_shared: int, max_metacells: int,
            rng: np.random.Generator) -> list[int]:
    """Greedy choice of neighbourhoods sharing <= max_shared cells with
    every neighbourhood already chosen, centres tried in random order."""
    n = len(nn)
    owners: list[list[int]] = [[] for _ in range(n)]  # cell -> chosen metacells
    chosen: list[int] = []
    shared = np.zeros(max_metacells, dtype=np.int64)
    for c in rng.permutation(n):
        touched = [m for cell in nn[c] for m in owners[cell]]
        ok = True
        if touched:
            np.add.at(shared, touched, 1)
            ok = shared[touched].max() <= max_shared
            shared[touched] = 0
        if not ok:
            continue
        m = len(chosen)
        chosen.append(int(c))
        for cell in nn[c]:
            owners[cell].append(m)
        if len(chosen) == max_metacells:
            break
    return chosen


def _aggregate(X, member_idx: np.ndarray, how: str) -> np.ndarray:
    """Mean or sum of the member rows of X, one row per metacell."""
    import scipy.sparse as sp

    n_meta, k = member_idx.shape
    rows = np.repeat(np.arange(n_meta), k)
    weight = 1.0 / k if how == "mean" else 1.0
    # only the member cells, in their original order (the same sums as over
    # all rows of X, without converting all of X to float64)
    used, cols = np.unique(member_idx.ravel(), return_inverse=True)
    P = sp.csr_matrix((np.full(rows.size, weight), (rows, cols.ravel())),
                      shape=(n_meta, used.size))
    Xu = _rows(X, used)
    out = P @ (sp.csr_matrix(Xu, dtype=np.float64) if sp.issparse(Xu)
               else np.asarray(Xu, dtype=np.float64))
    return np.asarray(out.toarray() if sp.issparse(out) else out, dtype=np.float64)
