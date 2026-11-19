"""Tests for gcen2.preprocessing — normalization and filtering."""

import numpy as np
import pandas as pd
import pytest

from gcen2 import preprocessing
from gcen2.exceptions import GCENInputError
from gcen2.models import ExpressionMatrix


class TestNormalization:
    def test_upqt_preserves_shape(self, small_expr_df):
        result = preprocessing.normalize_upqt(small_expr_df)
        assert result.shape == small_expr_df.shape

    def test_median_preserves_shape(self, small_expr_df):
        result = preprocessing.normalize_median(small_expr_df)
        assert result.shape == small_expr_df.shape

    def test_deseq_preserves_shape(self, small_expr_df):
        result = preprocessing.normalize_deseq(small_expr_df)
        assert result.shape == small_expr_df.shape

    def test_tmm_preserves_shape(self, small_expr_df):
        result = preprocessing.normalize_tmm(small_expr_df)
        assert result.shape == small_expr_df.shape

    def test_hkg_normalization(self, small_expr_df):
        hkg = [small_expr_df.index[0], small_expr_df.index[1]]
        result = preprocessing.normalize_hkg(small_expr_df, hkg)
        assert result.shape == small_expr_df.shape

    def test_hkg_missing_genes_raises(self, small_expr_df):
        with pytest.raises(GCENInputError):
            preprocessing.normalize_hkg(small_expr_df, ["nonexistent_gene"])

    def test_normalize_dispatcher(self, small_expr_df):
        for method in ["upqt", "median", "deseq", "tmm"]:
            result = preprocessing.normalize(small_expr_df, method=method)
            assert isinstance(result, pd.DataFrame)
            assert result.shape == small_expr_df.shape

    def test_normalize_unknown_method_raises(self, small_expr_df):
        with pytest.raises(GCENInputError):
            preprocessing.normalize(small_expr_df, method="invalid")

    def test_normalize_returns_expression_matrix(self, small_expr):
        result = preprocessing.normalize(small_expr, method="tmm")
        assert isinstance(result, ExpressionMatrix)
        assert result.is_normalized
        assert result.normalization_method == "tmm"

    def test_no_inplace_mutation(self, small_expr_df):
        original = small_expr_df.copy()
        preprocessing.normalize_tmm(small_expr_df)
        pd.testing.assert_frame_equal(small_expr_df, original)


class TestFiltering:
    def test_filter_by_percent_mean(self, small_expr_df):
        result = preprocessing.filter_genes(small_expr_df, percent_mean=0.5)
        assert len(result) <= len(small_expr_df)
        assert len(result) >= 1

    def test_filter_by_cutoff_mean(self, small_expr_df):
        threshold = small_expr_df.mean(axis=1).max() + 1
        result = preprocessing.filter_genes(small_expr_df, cutoff_mean=threshold)
        assert len(result) == 0

    def test_filter_no_op(self, small_expr_df):
        result = preprocessing.filter_genes(small_expr_df)
        assert result.shape == small_expr_df.shape

    def test_std_uses_ddof1(self, small_expr_df):
        std_ddof1 = small_expr_df.std(axis=1, ddof=1)
        std_ddof0 = small_expr_df.std(axis=1, ddof=0)
        assert not np.allclose(std_ddof1.values, std_ddof0.values)


class TestLogTransform:
    def test_log2_transform(self, small_expr_df):
        result = preprocessing.log_transform(small_expr_df, base="log2")
        expected = np.log2(small_expr_df + 1.0)
        np.testing.assert_allclose(result.values, expected.values)

    def test_invalid_base_raises(self, small_expr_df):
        with pytest.raises(GCENInputError):
            preprocessing.log_transform(small_expr_df, base="log3")
