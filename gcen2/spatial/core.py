"""Spatial co-expression network analysis: three build methods."""
from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd

from gcen2._anndata import require_anndata
from gcen2.exceptions import GCENInputError
from gcen2.models import Network
from gcen2.network import build as _network_build
from gcen2.spatial.neighbors import build_neighbor_graph


def _get_expression(adata: Any, layer: str | None, fraction: float = 0.0,
                    normalize: bool | None = None) -> tuple[np.ndarray, list[str]]:
    """Dense (n_spots, n_genes) expression of the genes with a value > 0 in
    at least ``round(fraction * n_spots)`` spots (all genes if 0), and their
    names. ``normalize``: Seurat LogNormalize (log1p of counts per 10,000);
    None: when the data are raw counts."""
    from gcen2.single_cell.metacell import _expression, _is_counts
    from gcen2.single_cell.modules import fraction_genes
    from gcen2.spatial.modules import log_normalize

    if not 0 <= fraction <= 1:
        raise GCENInputError("fraction must be between 0 and 1.")
    X = _expression(adata, layer)
    counts = _is_counts(X)
    if normalize is None:
        normalize = counts
    elif normalize and not counts:
        raise GCENInputError("normalize=True needs raw counts (adata.X or a layer).")
    idx = (np.flatnonzero(fraction_genes(X, fraction)) if fraction > 0
           else np.arange(X.shape[1]))
    if len(idx) == 0:
        raise GCENInputError("no gene passes the fraction filter.")
    if normalize:
        M = log_normalize(X, idx)
    else:
        import scipy.sparse as sp

        M = X[:, idx]
        M = M.toarray() if sp.issparse(M) else np.asarray(M)
        M = np.asarray(M, dtype=float)
    var = np.asarray(adata.var_names, dtype=str)
    return M, [str(var[i]) for i in idx]


def _get_coords(adata: Any) -> np.ndarray:
    if "spatial" not in adata.obsm:
        raise GCENInputError(
            "No spatial coordinates found in adata.obsm['spatial']. "
            "Make sure you have loaded a spatial AnnData object."
        )
    return np.asarray(adata.obsm["spatial"], dtype=float)


def _morans_i(X: np.ndarray, W) -> np.ndarray:
    """Compute Moran's I for each gene (column of X).

    Parameters
    ----------
    X : (n_spots, n_genes) expression matrix
    W : (n_spots, n_spots) spatial weight matrix (scipy sparse)

    Returns
    -------
    Array of shape (n_genes,) with Moran's I values.
    """
    n = X.shape[0]
    w_sum = W.sum()
    if w_sum == 0:
        return np.zeros(X.shape[1])

    z = X - X.mean(axis=0)          # mean-center each gene
    z2_sum = (z ** 2).sum(axis=0)
    z2_sum = np.where(z2_sum == 0, 1e-12, z2_sum)

    Wz = W @ z                       # spatial lag
    numerator = (z * Wz).sum(axis=0)
    I = (n / w_sum) * (numerator / z2_sum)
    return np.asarray(I).ravel()


def _morans_test(X: np.ndarray, W) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Moran's I of every gene (column of X) with its one-sided test for
    positive spatial autocorrelation: I, z and p, from the moments of I under
    randomisation (Cliff and Ord 1981), which allow any distribution of the
    values and an asymmetric W (the kNN graph). Constant genes: z = 0, p = 1."""
    from scipy.stats import norm

    n = X.shape[0]
    I = _morans_i(X, W)
    W = W.tocsr()
    S0 = W.sum()
    if S0 == 0 or n < 4:
        return I, np.zeros_like(I), np.ones_like(I)
    Wt = W.T.tocsr()
    S1 = 0.5 * float(((W + Wt).power(2)).sum())
    S2 = float(((np.asarray(W.sum(axis=1)).ravel() + np.asarray(W.sum(axis=0)).ravel()) ** 2)
               .sum())
    z = X - X.mean(axis=0)
    m2 = (z ** 2).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        b2 = n * (z ** 4).sum(axis=0) / m2 ** 2
    e = -1.0 / (n - 1)
    var = ((n * ((n * n - 3 * n + 3) * S1 - n * S2 + 3 * S0 * S0)
            - b2 * ((n * n - n) * S1 - 2 * n * S2 + 6 * S0 * S0))
           / ((n - 1) * (n - 2) * (n - 3) * S0 * S0) - e * e)
    ok = (m2 > 0) & (var > 0)
    zs = np.zeros_like(I)
    zs[ok] = (I[ok] - e) / np.sqrt(var[ok])
    p = np.ones_like(I)
    p[ok] = norm.sf(zs[ok])
    return I, zs, p


def _bh(p: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values (R ``p.adjust(p, "BH")``)."""
    p = np.asarray(p, dtype=float)
    m = len(p)
    if m == 0:
        return p
    order = np.argsort(p, kind="stable")
    q = p[order] * m / np.arange(1, m + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(m)
    out[order] = np.minimum(q, 1.0)
    return out


def _smooth(X: np.ndarray, coords: np.ndarray, W, bandwidth: float | None = None) -> np.ndarray:
    """Each spot replaced by the Gaussian-weighted mean of itself and its
    neighbours (``W``; with ``bandwidth``, the weights recomputed from the
    distances); the expression of :func:`build_spatial_weighted`."""
    import scipy.sparse as sp

    if bandwidth is not None:
        if bandwidth <= 0:
            raise GCENInputError(f"bandwidth must be positive, got {bandwidth}")
        W = W.copy()
        W.eliminate_zeros()
        W_coo = W.tocoo()
        if len(W_coo.data) > 0:
            dists = np.sqrt(((coords[W_coo.row] - coords[W_coo.col]) ** 2).sum(axis=1))
            W_coo.data[:] = np.exp(-(dists ** 2) / (2 * bandwidth ** 2))
            W = W_coo.tocsr()
    # self-loops: each spot includes itself in the weighted average
    W_self = W + sp.eye(X.shape[0], format="csr")
    row_sums = np.asarray(W_self.sum(axis=1)).ravel()
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    W_norm = sp.diags(1.0 / row_sums) @ W_self  # row-normalised
    return np.asarray(W_norm @ X, dtype=float)


def spatial_autocorrelation(
    adata: Any,
    *,
    layer: str | None = None,
    fraction: float = 0.05,
    normalize: bool | None = None,
    n_neighbors: int | None = None,
    radius: float | None = None,
    smoothed: bool = False,
    bandwidth: float | None = None,
) -> pd.DataFrame:
    """Moran's I of every gene (the genes and neighbour graph of
    :func:`build_spatial_filtered`): a table ``gene``, ``morans_i``,
    ``pval`` (one-sided test for positive autocorrelation under
    randomisation), ``fdr`` (Benjamini-Hochberg over the genes), ``mean`` and
    ``detected`` (fraction of spots with a value > 0); with ``smoothed``, also
    ``morans_i_smoothed``, the I of the expression smoothed as in
    :func:`build_spatial_weighted`. Rows in decreasing ``morans_i``."""
    require_anndata("spatial analysis")
    coords = _get_coords(adata)
    X, gene_names = _get_expression(adata, layer, fraction, normalize)
    X = np.asarray(X, dtype=float)
    W = build_neighbor_graph(coords, n_neighbors=n_neighbors, radius=radius)
    mi, _, p = _morans_test(X, W)
    out = pd.DataFrame({"gene": list(gene_names), "morans_i": mi, "pval": p, "fdr": _bh(p),
                        "mean": X.mean(axis=0), "detected": (X > 0).mean(axis=0)})
    if smoothed:
        out["morans_i_smoothed"] = _morans_i(_smooth(X, coords, W, bandwidth), W)
    return out.sort_values("morans_i", ascending=False, kind="stable").reset_index(drop=True)


def build_spatial_filtered(
    adata: Any,
    *,
    layer: str | None = None,
    fraction: float = 0.05,
    normalize: bool | None = None,
    n_neighbors: int | None = None,
    radius: float | None = None,
    cor_cutoff: float = 0.6,
    p_cutoff: float = 0.05,
    moran_fdr: float = 0.05,
    **build_kwargs,
) -> Network:
    """Build a co-expression network filtered by spatial autocorrelation.

    Constructs a standard co-expression network, then keeps the edges
    between genes with significant positive spatial autocorrelation: Moran's
    I > 0 and a Benjamini-Hochberg FDR <= ``moran_fdr`` of its one-sided test
    (moments under randomisation, Cliff and Ord 1981). ``moran_fdr=1``: only
    I > 0, which is not a test (on a DLPFC section 82% of the genes, and 69%
    of genes permuted across spots of similar depth, have I > 0).

    Parameters
    ----------
    adata : AnnData with ``obsm['spatial']`` coordinates.
    layer : Expression layer to use (``adata.X`` if None).
    fraction : Only genes with a value > 0 in at least
        ``round(fraction * n_spots)`` spots (hdWGCNA's gene selection; 0: all genes).
    normalize : Log-normalize (Seurat LogNormalize); None: when the data are
        raw counts.
    n_neighbors : KNN neighbours. Mutually exclusive with *radius*.
    radius : Radius cutoff. Mutually exclusive with *n_neighbors*.
    cor_cutoff : Minimum absolute correlation for an edge.
    p_cutoff : Maximum p-value for correlation significance (passed to network.build).
    moran_fdr : Maximum FDR of a gene's Moran's I test (see above).
    **build_kwargs : Additional keyword arguments forwarded to
        ``gcen2.network.build`` (e.g. ``method``, ``signed``).

    Returns
    -------
    Network with spatially unsupported edges removed.
    """
    require_anndata("spatial analysis")
    if not 0 < moran_fdr <= 1:
        raise GCENInputError(f"moran_fdr must be in (0, 1], got {moran_fdr}")

    coords = _get_coords(adata)
    X, gene_names = _get_expression(adata, layer, fraction, normalize)

    expr_df = pd.DataFrame(X.T, index=gene_names, columns=list(adata.obs_names))
    net = _network_build(
        expr_df,
        cor=cor_cutoff,
        pval=p_cutoff,
        **build_kwargs,
    )

    if net.edge_count == 0:
        return net

    W = build_neighbor_graph(coords, n_neighbors=n_neighbors, radius=radius)
    mi, _, p = _morans_test(X, W)
    spatial = set(np.asarray(gene_names, dtype=object)[(mi > 0) & (_bh(p) <= moran_fdr)])

    edges = net.edges.copy()
    mask = edges["gene_a"].isin(spatial) & edges["gene_b"].isin(spatial)
    return Network(
        edges=edges[mask].reset_index(drop=True),
        is_signed=net.is_signed,
        method=net.method,
    )


def build_spatial_weighted(
    adata: Any,
    *,
    layer: str | None = None,
    fraction: float = 0.05,
    normalize: bool | None = None,
    n_neighbors: int | None = None,
    radius: float | None = None,
    bandwidth: float | None = None,
    cor_cutoff: float = 0.6,
    p_cutoff: float = 0.05,
    **build_kwargs,
) -> Network:
    """Build a co-expression network using spatially smoothed expression.

    Smooths the expression matrix by replacing each spot's expression with a
    Gaussian-weighted average of its spatial neighbours (including itself).
    The smoothed matrix is then used for standard co-expression network
    construction, which reduces dropout noise while preserving spatial signal.

    Parameters
    ----------
    adata : AnnData with ``obsm['spatial']`` coordinates.
    layer : Expression layer to use (``adata.X`` if None).
    fraction : Only genes with a value > 0 in at least
        ``round(fraction * n_spots)`` spots (hdWGCNA's gene selection; 0: all genes).
    normalize : Log-normalize (Seurat LogNormalize); None: when the data are
        raw counts.
    n_neighbors : KNN neighbours. Mutually exclusive with *radius*.
    radius : Radius cutoff. Mutually exclusive with *n_neighbors*.
    bandwidth : Gaussian kernel bandwidth (in coordinate units). Auto-estimated
        from the median neighbour distance when None.
    cor_cutoff : Minimum absolute correlation for an edge.
    p_cutoff : Maximum p-value for an edge.

    Returns
    -------
    Network built on spatially smoothed expression values.
    """
    require_anndata("spatial analysis")
    coords = _get_coords(adata)
    X, gene_names = _get_expression(adata, layer, fraction, normalize)

    W = build_neighbor_graph(coords, n_neighbors=n_neighbors, radius=radius)
    X_smooth = _smooth(X, coords, W, bandwidth)

    expr_df = pd.DataFrame(
        np.asarray(X_smooth, dtype=float).T,
        index=gene_names,
        columns=list(adata.obs_names),
    )
    net = _network_build(expr_df, cor=cor_cutoff, pval=p_cutoff, **build_kwargs)
    return net


def build_spatial_regions(
    adata: Any,
    *,
    region_key: str | None = None,
    n_regions: int = 5,
    layer: str | None = None,
    fraction: float = 0.05,
    normalize: bool | None = None,
    n_neighbors: int | None = None,
    radius: float | None = None,
    **build_kwargs,
) -> dict[str, Network]:
    """Build one co-expression network per spatial region.

    Partitions spots into regions (manual annotation or automatic KMeans
    clustering on spatial coordinates), then builds a separate co-expression
    network for each region.

    Parameters
    ----------
    adata : AnnData with ``obsm['spatial']`` coordinates.
    region_key : Column in ``adata.obs`` to use as region labels. If None,
        performs automatic KMeans clustering into *n_regions* clusters.
    n_regions : Number of clusters for automatic partitioning (used only
        when *region_key* is None).
    layer : Expression layer to use (``adata.X`` if None).
    fraction : Only genes with a value > 0 in at least
        ``round(fraction * n_spots)`` spots (hdWGCNA's gene selection; 0: all genes).
    normalize : Log-normalize (Seurat LogNormalize); None: when the data are
        raw counts.
    n_neighbors, radius : Not used: accepted so that the three spatial build
        functions take the same arguments (the regions are not a neighbour graph).
    **build_kwargs : Forwarded to ``gcen2.network.build`` (e.g. ``cor``, ``pval``).

    Returns
    -------
    dict mapping region name → Network. Regions with fewer than 3 spots
    are skipped with a UserWarning.
    """
    require_anndata("spatial analysis")
    del n_neighbors, radius  # see the docstring
    coords = _get_coords(adata)  # validates obsm['spatial'] exists; stored for KMeans

    # Determine region labels
    if region_key is not None:
        if region_key not in adata.obs.columns:
            raise GCENInputError(
                f"region_key '{region_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        labels = adata.obs[region_key].astype(str).tolist()
    else:
        if n_regions < 1:
            raise GCENInputError(f"n_regions must be >= 1, got {n_regions}")
        try:
            from sklearn.cluster import KMeans
        except ImportError as exc:
            raise GCENInputError(
                "scikit-learn is required for automatic region clustering: "
                "pip install scikit-learn"
            ) from exc
        km = KMeans(n_clusters=n_regions, random_state=42, n_init=10)
        cluster_ids = km.fit_predict(coords)
        labels = [f"region_{c}" for c in cluster_ids]

    X, gene_names = _get_expression(adata, layer, fraction, normalize)
    spot_names = list(adata.obs_names)

    result: dict[str, Network] = {}
    for region in sorted(set(labels)):
        mask = [lbl == region for lbl in labels]
        n_spots = sum(mask)
        if n_spots < 3:
            warnings.warn(
                f"Region '{region}' has only {n_spots} spots (< 3), skipping.",
                UserWarning,
                stacklevel=2,
            )
            continue

        X_region = X[mask, :]
        spot_subset = [s for s, m in zip(spot_names, mask) if m]
        expr_df = pd.DataFrame(X_region.T, index=gene_names, columns=spot_subset)
        net = _network_build(expr_df, **build_kwargs)
        result[region] = net

    return result
