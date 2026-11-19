"""Smoke tests for gcen2.plotting.

These tests verify that plot functions execute without errors and return
the expected types. They do not check visual correctness.

Requires matplotlib. Skipped when matplotlib is not installed.
Uses non-interactive backend to avoid display requirements.
"""

import numpy as np
import pandas as pd
import pytest

try:
    import matplotlib
    matplotlib.use("Agg")  # non-interactive backend — MUST be set before pyplot import
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

requires_mpl = pytest.mark.skipif(not HAS_MPL, reason="matplotlib not installed")


@pytest.fixture(autouse=True)
def close_figures():
    yield
    if HAS_MPL:
        plt.close("all")


@requires_mpl
class TestNetworkStatic:
    def test_returns_figure(self, small_network):
        from gcen2.plotting import network_static
        fig = network_static(small_network)
        assert fig is not None

    def test_saves_file(self, small_network, tmp_path):
        from gcen2.plotting import network_static
        out = tmp_path / "network.png"
        network_static(small_network, output=out)
        assert out.exists()


@requires_mpl
class TestDotplot:
    def _make_enrichment(self):
        from gcen2.models import EnrichmentResult
        df = pd.DataFrame({
            "term_id": [f"GO:{i:04d}" for i in range(10)],
            "name":    [f"term_{i}" for i in range(10)],
            "namespace": ["biological_process"] * 10,
            "enrichment": np.random.default_rng(0).uniform(1, 5, 10),
            "study_count": np.arange(1, 11),
            "study_n": [50] * 10,
            "pop_count": np.arange(10, 110, 10),
            "pop_n": [1000] * 10,
            "p_value": np.logspace(-10, -1, 10),
        })
        return EnrichmentResult(results=df, ontology_type="GO",
                                gene_or_module="test_gene", p_cutoff=0.05)

    def test_dotplot_returns_figure(self):
        from gcen2.plotting import dotplot
        fig = dotplot(self._make_enrichment())
        assert fig is not None

    def test_barplot_returns_figure(self):
        from gcen2.plotting import barplot
        fig = barplot(self._make_enrichment())
        assert fig is not None

    def test_dotplot_empty_raises(self):
        from gcen2.exceptions import GCENInputError
        from gcen2.models import EnrichmentResult
        from gcen2.plotting import dotplot
        empty = EnrichmentResult(
            results=pd.DataFrame(), ontology_type="GO",
            gene_or_module="x", p_cutoff=0.05
        )
        with pytest.raises(GCENInputError):
            dotplot(empty)


@requires_mpl
class TestExpressionPlots:
    def test_boxplot_returns_figure(self, small_expr):
        from gcen2.plotting import expression_boxplot
        fig = expression_boxplot(small_expr)
        assert fig is not None

    def test_boxplot_saves_file(self, small_expr, tmp_path):
        from gcen2.plotting import expression_boxplot
        out = tmp_path / "boxplot.png"
        expression_boxplot(small_expr, output=out)
        assert out.exists()


@requires_mpl
class TestDegreeDistribution:
    def test_returns_figure(self, small_network):
        from gcen2.plotting import degree_distribution
        fig = degree_distribution(small_network)
        assert fig is not None
