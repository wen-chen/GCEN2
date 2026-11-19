"""AnnData helpers without a pandas import (used by gcen2.network on its
numpy path; re-exported by gcen2.single_cell)."""

from __future__ import annotations

from typing import Any

import numpy as np

from gcen2.exceptions import GCENInputError


def require_anndata(what: str = "single-cell functionality") -> None:
    """Fail with an install hint when anndata is missing."""
    try:
        import anndata  # noqa: F401
    except ImportError as exc:
        raise GCENInputError(f"anndata is required for {what}. "
                             "Install with: pip install gcen2[single-cell]") from exc


def is_anndata(obj: Any) -> bool:
    """True for an AnnData-like object (duck-typed, anndata not imported)."""
    return hasattr(obj, "obs") and hasattr(obj, "var_names") and hasattr(obj, "X")


def _cell_mask(adata: Any, groupby: str | None, groups: Any):
    if groupby is None:
        return None
    if groupby not in adata.obs.columns:
        raise GCENInputError(f"groupby key '{groupby}' not found in adata.obs.")
    if groups is None:
        raise GCENInputError("groups is required with groupby.")
    values = [groups] if isinstance(groups, (str, int, float)) else list(groups)
    col = adata.obs[groupby]
    mask = col.isin(values).to_numpy() | col.astype(str).isin(
        [str(v) for v in values]).to_numpy()
    if not mask.any():
        raise GCENInputError(f"No cells with {groupby} in {values}.")
    return mask


def adata_cell_names(adata: Any, groupby: str | None = None, groups: Any = None) -> list[str]:
    """Names of the cells :func:`adata_matrix` keeps (its columns)."""
    names = np.asarray(adata.obs_names, dtype=str)
    mask = _cell_mask(adata, groupby, groups)
    return (names if mask is None else names[mask]).tolist()


def adata_matrix(
    adata: Any,
    layer: str | None = None,
    use_highly_variable: bool = False,
    groupby: str | None = None,
    groups: Any = None,
) -> tuple[np.ndarray, list[str]]:
    """Genes x cells float64 matrix (C order) and gene names of an AnnData.

    Parameters
    ----------
    adata:
        AnnData with cells as rows (``X`` dense or scipy sparse).
    layer:
        ``adata.layers[layer]`` instead of ``adata.X``.
    use_highly_variable:
        Keep the genes flagged in ``adata.var["highly_variable"]``.
    groupby, groups:
        Keep the cells whose ``adata.obs[groupby]`` is ``groups`` (one value
        or a list of values).

    A sparse matrix is densified only after the cell and gene subsets, and
    directly in genes x cells order.
    """
    if not is_anndata(adata):
        raise GCENInputError("Expected an AnnData object.")
    if use_highly_variable:
        if "highly_variable" not in adata.var.columns:
            raise GCENInputError(
                "use_highly_variable=True but 'highly_variable' not in adata.var. "
                "Run sc.pp.highly_variable_genes(adata) first.")
        gene_mask = adata.var["highly_variable"].to_numpy(dtype=bool)
    else:
        gene_mask = None
    cell_mask = _cell_mask(adata, groupby, groups)
    if layer is not None:
        if layer not in adata.layers:
            raise GCENInputError(f"Layer '{layer}' not found in adata.layers.")
        X = adata.layers[layer]
    else:
        X = adata.X
    if X is None:
        raise GCENInputError("adata.X is empty; give a layer.")
    if cell_mask is not None:
        X = X[cell_mask]
    if gene_mask is not None:
        X = X[:, gene_mask]
    genes = np.asarray(adata.var_names, dtype=str)
    if gene_mask is not None:
        genes = genes[gene_mask]
    if hasattr(X, "toarray"):  # scipy sparse: transpose first (CSR -> CSC view)
        matrix = X.T.toarray(order="C").astype(np.float64, copy=False)
    else:
        matrix = np.asarray(X, dtype=np.float64).T
    return np.ascontiguousarray(matrix), genes.tolist()

