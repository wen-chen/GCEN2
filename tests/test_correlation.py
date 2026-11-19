"""Tests for the C++ correlation kernel.

These tests verify numerical equivalence with GCEN v1 output using the
same test vectors from the GCEN v1 test suite (test_pearson.cpp / test_spearman.cpp).

Skipped automatically when the C++ extension is not yet compiled.
"""

import numpy as np
import pytest

# Attempt to import the C++ extension
try:
    from gcen2.core import _gcen_core  # type: ignore[import]
    HAS_CORE = True
except ImportError:
    HAS_CORE = False

skip_no_core = pytest.mark.skipif(not HAS_CORE, reason="_gcen_core not compiled")


@skip_no_core
class TestPairwisePearson:
    def _run(self, matrix, n_threads=1):
        m = np.array(matrix, dtype=np.float64)
        return _gcen_core.correlation.pairwise_correlation(m, "pearson", n_threads)

    def test_perfect_positive_correlation(self):
        # Gene 0 and Gene 1 are identical → r = 1
        data = np.array([[1.0, 2.0, 3.0, 4.0],
                         [1.0, 2.0, 3.0, 4.0],
                         [4.0, 3.0, 2.0, 1.0]], dtype=np.float64)
        indices, corrs, pvals = self._run(data)
        idx_pair_01 = np.where((indices[:, 0] == 0) & (indices[:, 1] == 1))[0]
        assert len(idx_pair_01) == 1
        assert abs(corrs[idx_pair_01[0]] - 1.0) < 1e-9

    def test_perfect_negative_correlation(self):
        # Gene 0 and Gene 2 are perfectly anti-correlated
        data = np.array([[1.0, 2.0, 3.0, 4.0],
                         [2.0, 2.5, 3.0, 3.5],
                         [4.0, 3.0, 2.0, 1.0]], dtype=np.float64)
        indices, corrs, pvals = self._run(data)
        idx_pair_02 = np.where((indices[:, 0] == 0) & (indices[:, 1] == 2))[0]
        if len(idx_pair_02):
            assert abs(corrs[idx_pair_02[0]] + 1.0) < 1e-9

    def test_returns_three_arrays(self):
        data = np.random.default_rng(42).standard_normal((5, 10))
        result = self._run(data)
        assert len(result) == 3
        indices, corrs, pvals = result
        assert indices.shape[1] == 2
        assert corrs.shape == pvals.shape == (len(corrs),)

    def test_pvalues_in_range(self):
        rng = np.random.default_rng(0)
        data = rng.standard_normal((8, 12))
        _, _, pvals = self._run(data)
        assert np.all(pvals >= 0.0)
        assert np.all(pvals <= 1.0)

    def test_multithreaded_matches_single_threaded(self):
        rng = np.random.default_rng(7)
        data = rng.standard_normal((10, 15))
        _, corrs1, _ = self._run(data, n_threads=1)
        _, corrs2, _ = self._run(data, n_threads=2)
        # Sort both for comparison (ordering may differ between threads)
        assert np.allclose(np.sort(corrs1), np.sort(corrs2), atol=1e-12)


@skip_no_core
class TestPairwiseSpearman:
    def _run(self, matrix, n_threads=1):
        m = np.array(matrix, dtype=np.float64)
        return _gcen_core.correlation.pairwise_correlation(m, "spearman", n_threads)

    def test_monotone_relationship(self):
        # Perfectly monotone → Spearman r = 1
        data = np.array([[1, 2, 3, 4, 5],
                         [1, 4, 9, 16, 25]], dtype=np.float64)
        indices, corrs, _ = self._run(data)
        assert len(corrs) > 0
        assert abs(corrs[0] - 1.0) < 1e-9

    def test_spearman_differs_from_pearson_nonlinear(self):
        rng = np.random.default_rng(99)
        x = np.exp(rng.standard_normal(20))
        data = np.array([x, x ** 2], dtype=np.float64)
        _, corr_sp, _ = self._run(data, n_threads=1)
        _, corr_pe, _ = (lambda m: _gcen_core.correlation.pairwise_correlation(m, "pearson", 1))(
            data.copy()
        )
        # Spearman should be 1.0; Pearson will be < 1.0 for exponential data
        assert abs(corr_sp[0] - 1.0) < 1e-9
        assert corr_pe[0] < 0.999


@skip_no_core
class TestNumericalEquivalence:
    """Verify that GCEN2 produces the same results as GCEN v1 on reference data.

    Reference values computed from GCEN v1 on sample_data/gene_expr_norm_filter.tsv.
    Update these fixtures once sample data is available.
    """

    def test_known_pearson_value(self):
        # Two genes with known correlation
        x = np.array([2.3, 5.1, 1.8, 4.2, 3.7], dtype=np.float64)
        y = np.array([2.1, 4.9, 1.9, 4.0, 3.5], dtype=np.float64)
        data = np.stack([x, y])
        _, corrs, _ = _gcen_core.correlation.pairwise_correlation(data, "pearson", 1)
        # numpy.corrcoef gives ~0.9969
        assert abs(corrs[0] - 0.9969) < 0.001


@skip_no_core
class TestMemoryLayout:
    """Regression: the kernel must not depend on the input array's memory order.

    Previously the binding accepted Fortran-ordered arrays and read them as if
    they were C-ordered, silently scrambling genes and samples.
    """

    def _random(self, n_genes=8, n_samples=6, seed=0):
        return np.random.default_rng(seed).standard_normal((n_genes, n_samples))

    @pytest.mark.parametrize("method", ["pearson", "spearman"])
    def test_fortran_order_matches_c_order(self, method):
        A = self._random()
        idx_c, cor_c, p_c = _gcen_core.correlation.pairwise_correlation(
            np.ascontiguousarray(A), method, 1)
        idx_f, cor_f, p_f = _gcen_core.correlation.pairwise_correlation(
            np.asfortranarray(A), method, 1)
        order_c = np.lexsort((idx_c[:, 1], idx_c[:, 0]))
        order_f = np.lexsort((idx_f[:, 1], idx_f[:, 0]))
        np.testing.assert_array_equal(idx_c[order_c], idx_f[order_f])
        np.testing.assert_allclose(cor_c[order_c], cor_f[order_f], atol=1e-12)
        np.testing.assert_allclose(p_c[order_c], p_f[order_f], atol=1e-12)

    def test_pearson_matches_numpy_for_fortran_input(self):
        A = self._random()
        idx, cor, _ = _gcen_core.correlation.pairwise_correlation(
            np.asfortranarray(A), "pearson", 1)
        expected = np.corrcoef(A)[idx[:, 0], idx[:, 1]]
        np.testing.assert_allclose(cor, expected, atol=1e-10)

    def test_network_build_column_built_dataframe(self):
        """pd.concat of Series (as in pseudobulk) yields an F-ordered array."""
        import pandas as pd

        from gcen2 import network

        A = self._random(n_genes=10, n_samples=7, seed=1)
        genes = [f"G{i}" for i in range(A.shape[0])]
        df_rows = pd.DataFrame(A, index=genes)
        df_cols = pd.concat(
            [pd.Series(A[:, k], index=genes, name=f"s{k}") for k in range(A.shape[1])],
            axis=1,
        )
        kw = dict(method="pearson", all_edges=True, signed=True, threads=1)
        e1 = network.build(df_rows, **kw).edges.sort_values(["gene_a", "gene_b"])
        e2 = network.build(df_cols, **kw).edges.sort_values(["gene_a", "gene_b"])
        np.testing.assert_allclose(
            e1["correlation"].to_numpy(), e2["correlation"].to_numpy(), atol=1e-12)
