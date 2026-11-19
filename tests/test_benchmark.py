"""Benchmark tests for GCEN2 computation kernels.

Measures wall-clock time for the main performance-critical operations.
NOT run in regular CI — requires pytest-benchmark and are marked @slow.

Run with:
    pytest tests/test_benchmark.py -v -m benchmark --benchmark-sort=mean

Reference targets (SRS §5.1):
    Network build (10k genes, 100 samples) < 5 min on 8 cores.
"""

import numpy as np
import pytest

try:
    from gcen2.core import _gcen_core
    HAS_CORE = True
except ImportError:
    HAS_CORE = False

try:
    import pytest_benchmark  # noqa: F401
    HAS_BENCHMARK = True
except ImportError:
    HAS_BENCHMARK = False

skip_no_core = pytest.mark.skipif(not HAS_CORE, reason="_gcen_core not compiled")
skip_no_benchmark = pytest.mark.skipif(
    not HAS_BENCHMARK, reason="pytest-benchmark not installed"
)


def _make_expr(n_genes, n_samples, seed=42):
    rng = np.random.default_rng(seed)
    return rng.standard_normal((n_genes, n_samples)).astype(np.float64)


@skip_no_core
@skip_no_benchmark
@pytest.mark.benchmark
@pytest.mark.slow
class TestCorrelationBenchmarks:
    """Benchmark the C++ correlation kernel at various gene counts."""

    def test_bench_500_genes(self, benchmark):
        data = _make_expr(500, 100)
        indices, corrs, pvals = benchmark(
            _gcen_core.correlation.pairwise_correlation, data, "spearman", 2
        )
        assert len(corrs) <= 500 * 499 // 2

    def test_bench_1k_genes(self, benchmark):
        data = _make_expr(1000, 100)
        indices, corrs, pvals = benchmark(
            _gcen_core.correlation.pairwise_correlation, data, "spearman", 2
        )
        assert len(corrs) <= 1000 * 999 // 2

    def test_bench_2k_genes_4threads(self, benchmark):
        data = _make_expr(2000, 100)
        indices, corrs, pvals = benchmark(
            _gcen_core.correlation.pairwise_correlation, data, "spearman", 4
        )
        assert len(corrs) <= 2000 * 1999 // 2

    def test_thread_scaling(self):
        """4 threads should be faster than 1 thread (logged, not asserted)."""
        import time
        data = _make_expr(2000, 100)

        t0 = time.perf_counter()
        _gcen_core.correlation.pairwise_correlation(data, "spearman", 1)
        t1 = time.perf_counter()

        t2 = time.perf_counter()
        _gcen_core.correlation.pairwise_correlation(data, "spearman", 4)
        t3 = time.perf_counter()

        speedup = (t1 - t0) / max(t3 - t2, 1e-9)
        print(f"\n1-thread: {t1-t0:.2f}s, 4-thread: {t3-t2:.2f}s, speedup: {speedup:.2f}x")
        # Log only — don't assert speedup (CI machines vary)


@skip_no_core
@skip_no_benchmark
@pytest.mark.benchmark
@pytest.mark.slow
class TestModuleBenchmarks:

    def _make_csr(self, n, density=0.05, seed=0):
        from scipy.sparse import random as sparse_random
        A = sparse_random(n, n, density=density, random_state=seed, format="csr")
        A = ((A + A.T) > 0).astype(np.float64).tocsr()
        A.setdiag(0)
        A.eliminate_zeros()
        A.sort_indices()
        return A.indptr.astype(np.int64), A.indices.astype(np.int32)

    def test_bench_500_nodes(self, benchmark):
        indptr, indices = self._make_csr(500, density=0.05)
        labels = benchmark(_gcen_core.module.jaccard_modules, indptr, indices, 0.3, 2)
        assert len(labels) == 500

    def test_bench_1k_nodes(self, benchmark):
        indptr, indices = self._make_csr(1000, density=0.03)
        labels = benchmark(_gcen_core.module.jaccard_modules, indptr, indices, 0.3, 2)
        assert len(labels) == 1000


@skip_no_core
@skip_no_benchmark
@pytest.mark.benchmark
@pytest.mark.slow
class TestRWRBenchmarks:

    def test_bench_rwr_1k(self, benchmark):
        n = 1000
        rng = np.random.default_rng(0)
        W = np.abs(rng.standard_normal((n, n)))
        col_sums = W.sum(axis=0)
        col_sums[col_sums == 0] = 1.0
        W /= col_sums
        P0 = np.zeros(n)
        P0[:10] = 0.1

        from scipy.sparse import csr_matrix
        Ws = csr_matrix(W)
        PT, steps = benchmark(
            _gcen_core.rwr.random_walk_restart, Ws.indptr.astype(np.int64),
            Ws.indices.astype(np.int32), Ws.data, P0, 0.5, 1_000_000, 1e-10
        )
        assert steps < 1_000_000  # must converge
        assert abs(PT.sum() - 1.0) < 0.01
