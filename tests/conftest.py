"""Shared pytest fixtures for GCEN2 tests."""

import os

# figures are only saved: no GUI backend (over SSH X forwarding its display
# probe can block); an explicit MPLBACKEND wins
os.environ.setdefault("MPLBACKEND", "Agg")
# subprocesses of the tests print UTF-8 (Windows pipes default to the ANSI code
# page, which cannot encode the arrows and Greek letters gcen2 prints)
os.environ.setdefault("PYTHONIOENCODING", "utf-8")

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import pytest  # noqa: E402

from gcen2.models import ExpressionMatrix, Network


@pytest.fixture
def small_expr_df() -> pd.DataFrame:
    """Small deterministic expression matrix (10 genes × 8 samples)."""
    rng = np.random.default_rng(42)
    data = rng.exponential(scale=10.0, size=(10, 8))
    genes = [f"Gene{i}" for i in range(10)]
    samples = [f"Sample{j}" for j in range(8)]
    return pd.DataFrame(data, index=genes, columns=samples)


@pytest.fixture
def small_expr(small_expr_df) -> ExpressionMatrix:
    return ExpressionMatrix(
        data=small_expr_df,
        gene_names=np.array(small_expr_df.index),
        sample_names=np.array(small_expr_df.columns),
    )


@pytest.fixture
def small_network() -> Network:
    """Small synthetic network with 5 edges."""
    edges = pd.DataFrame(
        {
            "gene_a": ["GeneA", "GeneA", "GeneB", "GeneC", "GeneD"],
            "gene_b": ["GeneB", "GeneC", "GeneC", "GeneD", "GeneE"],
            "correlation": [0.9, 0.85, 0.8, 0.75, 0.7],
            "p_value": [0.001, 0.002, 0.003, 0.004, 0.005],
        }
    )
    return Network(edges=edges)
