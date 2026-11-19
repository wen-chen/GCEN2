"""Spatial neighbor graph construction (KNN or radius-based)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from gcen2.exceptions import GCENInputError


def build_neighbor_graph(
    coords: np.ndarray,
    n_neighbors: int | None = None,
    radius: float | None = None,
) -> sp.csr_matrix:
    """Build a Gaussian-weighted spatial neighbor graph.

    Parameters
    ----------
    coords:
        Array of shape ``(n_spots, 2)`` with 2-D spatial coordinates.
    n_neighbors:
        Number of nearest neighbors (KNN mode). Mutually exclusive with *radius*.
    radius:
        Distance cutoff (radius mode). Mutually exclusive with *n_neighbors*.

    Returns
    -------
    scipy.sparse.csr_matrix of shape ``(n_spots, n_spots)`` with Gaussian
    edge weights ``exp(-d² / (2σ²))``, where σ = median neighbor distance.
    """
    if n_neighbors is not None and radius is not None:
        raise GCENInputError("Specify either n_neighbors or radius, not both")
    if n_neighbors is None and radius is None:
        raise GCENInputError("Must specify either n_neighbors or radius")

    from scipy.spatial import KDTree

    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    tree = KDTree(coords)

    if n_neighbors is not None:
        if n_neighbors < 1:
            raise GCENInputError(f"n_neighbors must be >= 1, got {n_neighbors}")
        if n_neighbors >= n:
            raise GCENInputError(
                f"n_neighbors={n_neighbors} must be less than n_spots={n}"
            )
        # k + 1 with the point itself first, equidistant points by index (the
        # same neighbours on every platform; grids have many ties)
        from gcen2._knn import knn_indices

        indices = knn_indices(coords, n_neighbors + 1)[:, 1:]  # drop self
        diff = coords[indices] - coords[:, None, :]
        dists = np.sqrt((diff * diff).sum(axis=2))

        sigma = float(np.median(dists[dists > 0]))
        if not np.isfinite(sigma) or sigma == 0.0:
            sigma = 1.0
        weights = np.exp(-(dists.ravel() ** 2) / (2 * sigma ** 2))
        rows = np.repeat(np.arange(n), n_neighbors)
        cols = indices.ravel()

    else:  # radius mode
        if radius < 0:
            raise GCENInputError(f"radius must be non-negative, got {radius}")
        pairs = tree.query_pairs(radius, output_type="ndarray")
        if len(pairs) == 0:
            return sp.csr_matrix((n, n), dtype=float)

        i_idx, j_idx = pairs[:, 0], pairs[:, 1]
        dists_pair = np.linalg.norm(coords[i_idx] - coords[j_idx], axis=1)
        sigma = float(np.median(dists_pair))
        if not np.isfinite(sigma) or sigma == 0.0:
            sigma = 1.0
        w = np.exp(-(dists_pair ** 2) / (2 * sigma ** 2))

        # Symmetric: add both directions
        rows = np.concatenate([i_idx, j_idx])
        cols = np.concatenate([j_idx, i_idx])
        weights = np.concatenate([w, w])

    W = sp.csr_matrix((weights, (rows, cols)), shape=(n, n), dtype=float)
    W = W.maximum(W.T)
    return W
