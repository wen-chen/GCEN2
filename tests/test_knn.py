"""gcen2._knn.knn_indices: the k nearest points with a platform-independent
order of ties (cKDTree's order of equidistant points differs between
platforms, and grids of spots or bins have many)."""

import numpy as np
import pytest

from gcen2._knn import knn_indices


def _brute(pts, k):
    d = ((pts[:, None, :] - pts[None, :, :]) ** 2).sum(axis=2)
    n = len(pts)
    out = []
    for i in range(n):
        key = np.lexsort((np.arange(n), d[i], np.arange(n) != i))
        out.append(key[:k])
    return np.array(out)


@pytest.mark.parametrize("k", [1, 2, 5, 6, 9, 13])
def test_grid_ties_by_index(k):
    """On a grid, equidistant points are taken by index: the same as sorting
    every point by (not itself, squared distance, index)."""
    xy = np.array([(r, c) for r in range(10) for c in range(10)], dtype=float)
    got = knn_indices(xy, k, block=7)
    assert np.array_equal(got, _brute(xy, k))
    assert np.array_equal(got[:, 0], np.arange(len(xy)))


def test_duplicates_and_random_points():
    rng = np.random.default_rng(0)
    pts = rng.normal(size=(300, 3))
    pts[50:60] = pts[10]  # duplicates: distance 0 ties
    for k in (1, 4, 12, 25):
        assert np.array_equal(knn_indices(pts, k, workers=2), _brute(pts, k))


def test_same_sets_as_kdtree_without_ties():
    from scipy.spatial import cKDTree

    pts = np.random.default_rng(1).normal(size=(500, 30))
    _, nn = cKDTree(pts).query(pts, k=8)
    assert np.array_equal(knn_indices(pts, 8), nn)


def test_k_out_of_range():
    with pytest.raises(ValueError):
        knn_indices(np.zeros((3, 2)), 4)
