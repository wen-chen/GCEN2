"""k nearest neighbours whose order does not depend on the platform."""
from __future__ import annotations

import numpy as np


def knn_indices(points, k: int, workers: int = 1, block: int = 4096) -> np.ndarray:
    """Indices (n, k) of the k nearest points of every point: the point itself
    first, then by increasing squared Euclidean distance and, for equal
    distances, by index.

    ``scipy.spatial.cKDTree`` returns equidistant points in an order that
    differs between platforms, and grids of spots or bins have many of them,
    so the neighbours themselves (not only their order) could differ. Here the
    k + 1 nearest are queried; where the k-th and (k + 1)-th are equally far,
    all points within that distance are gathered and the tie is broken by
    index. On a grid the squared distances are exact, so the result is the
    same everywhere.
    """
    from scipy.spatial import cKDTree

    pts = np.asarray(points, dtype=np.float64)
    if pts.ndim == 1:
        pts = pts[:, None]
    pts = np.ascontiguousarray(pts)
    n = len(pts)
    k = int(k)
    if not 1 <= k <= n:
        raise ValueError(f"k must be between 1 and the number of points ({n}), got {k}")
    tree = cKDTree(pts)
    kq = min(k + 1, n)
    _, nn = tree.query(pts, k=kq, workers=max(1, int(workers)))
    nn = np.asarray(nn, dtype=np.int64).reshape(n, kq)
    out = np.empty((n, k), dtype=np.int64)
    ambiguous = []
    for r0 in range(0, n, block):
        r1 = min(r0 + block, n)
        cand = nn[r0:r1]
        rows = np.arange(r0, r1)
        diff = pts[cand] - pts[rows, None, :]
        sq = (diff * diff).sum(axis=2)
        not_self = cand != rows[:, None]
        order = np.lexsort((cand, sq, not_self), axis=-1)
        cand = np.take_along_axis(cand, order, axis=1)
        sq = np.take_along_axis(sq, order, axis=1)
        out[r0:r1] = cand[:, :k]
        amb = cand[:, 0] != rows  # itself not among the candidates (duplicates)
        if kq > k:
            amb |= sq[:, k - 1] == sq[:, k]
        ambiguous.extend((rows[amb]).tolist())
    for i in ambiguous:
        far = float(np.sqrt(((pts[out[i, -1]] - pts[i]) ** 2).sum()))
        cand = np.asarray(tree.query_ball_point(pts[i], far * (1 + 1e-9) + 1e-300),
                          dtype=np.int64)
        if not (cand == i).any():
            cand = np.append(cand, i)
        diff = pts[cand] - pts[i]
        sq = (diff * diff).sum(axis=1)
        order = np.lexsort((cand, sq, cand != i))
        out[i] = cand[order][:k]
    return out
