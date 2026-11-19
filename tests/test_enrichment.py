"""Tests for the C++ enrichment (hypergeometric p-value) kernel.

Reference values imported from GCEN v1 test_hypergeometric.cpp.
Skipped when the C++ extension is not compiled.
"""

import numpy as np
import pytest

try:
    from gcen2.core import _gcen_core  # type: ignore[import]
    HAS_CORE = True
except ImportError:
    HAS_CORE = False

skip_no_core = pytest.mark.skipif(not HAS_CORE, reason="_gcen_core not compiled")

# Also test the Python fallback in annotation.py (always available)
from tests.reference import hypergeometric_pvalue


class TestPythonHypergeometric:
    """Test the scipy-based Python fallback."""

    def test_pvalue_range(self):
        p = hypergeometric_pvalue(100, 20, 10, 5)
        assert 0.0 <= p <= 1.0

    def test_high_enrichment_low_pvalue(self):
        # All 20 study genes are in the term → should be highly significant
        p = hypergeometric_pvalue(1000, 20, 20, 20)
        assert p < 0.05

    def test_no_enrichment_high_pvalue(self):
        # Study proportion matches background
        p = hypergeometric_pvalue(100, 10, 10, 1)
        assert p > 0.05

    def test_zero_study_count(self):
        p = hypergeometric_pvalue(100, 10, 10, 0)
        # P(X >= 0) should be 1.0
        assert abs(p - 1.0) < 1e-9


@skip_no_core
class TestCppHypergeometric:
    """Test the C++ implementation — must match scipy values."""

    def _py(self, pop_n, pop_count, study_n, study_count):
        return hypergeometric_pvalue(pop_n, pop_count, study_n, study_count)

    def _cpp(self, pop_n, pop_count, study_n, study_count):
        return _gcen_core.enrichment.hypergeometric_pvalue(
            pop_n, pop_count, study_n, study_count
        )

    @pytest.mark.parametrize("args", [
        (100, 20, 10, 5),
        (1000, 50, 30, 10),
        (500, 100, 20, 15),
        (200, 10, 10, 1),
    ])
    def test_matches_python_fallback(self, args):
        py_p = self._py(*args)
        cpp_p = self._cpp(*args)
        assert abs(py_p - cpp_p) < 1e-6, (
            f"Mismatch for args={args}: Python={py_p:.8f}, C++={cpp_p:.8f}"
        )


@skip_no_core
class TestBatchHypergeometric:
    def test_batch_matches_single(self):
        pop_n = 500
        study_n = 30
        pop_counts = np.array([20, 50, 10, 80, 5], dtype=np.int32)
        study_counts = np.array([5, 8, 2, 15, 1], dtype=np.int32)

        batch_p = _gcen_core.enrichment.batch_hypergeometric(
            pop_n, pop_counts, study_n, study_counts
        )
        for i in range(len(pop_counts)):
            single_p = _gcen_core.enrichment.hypergeometric_pvalue(
                pop_n, int(pop_counts[i]), study_n, int(study_counts[i])
            )
            assert abs(batch_p[i] - single_p) < 1e-9

    def test_batch_returns_correct_length(self):
        n = 50
        pop_counts = np.ones(n, dtype=np.int32) * 10
        study_counts = np.ones(n, dtype=np.int32) * 2
        result = _gcen_core.enrichment.batch_hypergeometric(
            100, pop_counts, 20, study_counts
        )
        assert len(result) == n

    def test_batch_pvalues_in_range(self):
        rng = np.random.default_rng(5)
        n = 30
        pop_counts = rng.integers(1, 50, size=n).astype(np.int32)
        study_counts = (pop_counts * rng.uniform(0, 0.5, size=n)).astype(np.int32)
        result = _gcen_core.enrichment.batch_hypergeometric(
            200, pop_counts, 40, study_counts
        )
        assert np.all(result >= 0.0)
        assert np.all(result <= 1.0)


def _random_annotation(seed=0, n_genes=400, n_terms=120):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    assoc = {}
    for g in genes:
        k = int(rng.integers(0, 6))
        assoc[g] = {f"T{int(t):03d}" for t in rng.choice(n_terms, size=k, replace=False)}
    # a few terms concentrated in the first genes, so that queries hit them
    for g in genes[:40]:
        assoc[g] |= {"T900", "T901"}
    return genes, assoc


def test_fdr_is_bh_over_all_background_terms():
    """fdr = R's p.adjust(p, "BH") over every term annotated in the background,
    the terms without a query gene counted with p = 1."""
    from gcen2.annotation import enrich

    try:
        from scipy.stats import false_discovery_control  # SciPy >= 1.11
    except ImportError:  # the oldest supported SciPy: R's p.adjust(p, "BH") by hand
        def false_discovery_control(p):
            p = np.asarray(p, dtype=float)
            order = np.argsort(p)[::-1]
            adj = np.minimum.accumulate(p[order] * len(p) / np.arange(len(p), 0, -1))
            out = np.empty_like(p)
            out[order] = np.minimum(adj, 1.0)
            return out

    genes, assoc = _random_annotation()
    query, background = genes[:60], genes
    full = enrich(query, background, None, assoc, pval=1.0).results.set_index("term_id")
    tested = sorted(set().union(*(assoc[g] for g in background)))
    p = full["p_value"].reindex(tested).fillna(1.0).to_numpy()
    ref = dict(zip(tested, false_discovery_control(p)))
    np.testing.assert_allclose(full["fdr"].to_numpy(), [ref[t] for t in full.index], rtol=1e-12)
    # the cutoff keeps the same adjusted values
    cut = enrich(query, background, None, assoc, pval=0.01).results.set_index("term_id")
    assert len(cut) < len(full)
    np.testing.assert_array_equal(cut["fdr"].to_numpy(), full.loc[cut.index, "fdr"].to_numpy())


def test_fdr_same_with_and_without_the_extension(monkeypatch, tmp_path):
    """The C++ kernel (fused counts, test and BH) and the numpy path write the
    same rows; p and fdr agree to rounding (the two hypergeometric codes
    differ in the last digits)."""
    import pandas as pd

    from gcen2 import annotation
    from gcen2.annotation import KEGGMapping
    from gcen2.models import ModuleSet

    genes, assoc = _random_annotation(seed=1)
    modules = ModuleSet(modules=[set(genes[:50]), set(genes[100:180]), set(genes[300:330])])
    kegg = KEGGMapping(k_to_pathways={}, pathway_names={})
    a, b = tmp_path / "a", tmp_path / "b"
    annotation.annotate_to_dir(modules, a, assoc, kegg=kegg)
    monkeypatch.setattr(annotation, "extension", lambda: None)
    annotation.annotate_to_dir(modules, b, assoc, kegg=kegg)
    files = sorted(p.name for p in a.glob("*.kegg"))
    assert files and files == sorted(p.name for p in b.glob("*.kegg"))
    for f in files:
        x, y = (pd.read_csv(d / f, sep="\t") for d in (a, b))
        assert list(x.columns)[-2:] == ["p_value", "fdr"]
        assert x["term_id"].tolist() == y["term_id"].tolist()
        np.testing.assert_allclose(x[["p_value", "fdr"]], y[["p_value", "fdr"]], rtol=1e-11)
