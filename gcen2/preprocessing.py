"""Data normalization and filtering.

Python reimplementation of GCEN v1 data_norm / data_filter, ported from
src/util/norm.hpp and src/data_filter.cpp.

Bug fixes vs GCEN v1:
  - deviation() now uses sample std dev (ddof=1) instead of population std dev (ddof=0).
    Pass ddof=0 explicitly if you need the old population-std behaviour.
"""

from __future__ import annotations

from typing import Union

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.models import ExpressionMatrix


def normalize(
    expr: Union[ExpressionMatrix, pd.DataFrame],
    method: str = "tmm",
    housekeeping_genes: list[str] | None = None,
) -> Union[ExpressionMatrix, pd.DataFrame]:
    """Normalize gene expression data.

    Parameters
    ----------
    expr:
        Expression matrix (genes × samples).
    method:
        Normalization method: 'upqt', 'median', 'deseq', 'tmm', or 'hkg'.
    housekeeping_genes:
        Required when method='hkg'. List of reference gene identifiers.

    Returns
    -------
    Normalized expression matrix of the same type as the input.
    """
    df, is_wrapped = _unwrap(expr)

    dispatch = {
        "upqt":   normalize_upqt,
        "median": normalize_median,
        "deseq":  normalize_deseq,
        "tmm":    normalize_tmm,
        "hkg":    lambda d: normalize_hkg(d, housekeeping_genes or []),
    }
    if method not in dispatch:
        raise GCENInputError(
            f"Unknown normalization method '{method}'. "
            f"Choose from: {', '.join(dispatch.keys())}"
        )

    result_df = dispatch[method](df)

    if is_wrapped:
        return ExpressionMatrix(
            data=result_df,
            gene_names=np.array(result_df.index, dtype=str),
            sample_names=np.array(result_df.columns, dtype=str),
            is_normalized=True,
            normalization_method=method,
        )
    return result_df


def normalize_upqt(expr: pd.DataFrame) -> pd.DataFrame:
    """Upper quartile normalization.

    Mirrors GCEN v1 norm_upqt(): divide each sample by the 75th percentile of
    its nonzero values, then rescale so that the global mean of upper quartiles
    is preserved.
    """
    # Per-sample upper quartile of nonzero expression
    upqt = expr.apply(
        lambda col: _upper_quartile(col[col > 0].to_numpy()) if (col > 0).any() else 1.0
    )
    upqt_mean = upqt.mean()
    # factor[j] = upqt[j] / upqt_mean  →  result = expr / factor (columnwise)
    factor = upqt / upqt_mean
    factor = factor.replace(0.0, np.nan).fillna(1.0)
    return expr.divide(factor, axis=1)


def normalize_median(expr: pd.DataFrame) -> pd.DataFrame:
    """Median normalization.

    Mirrors GCEN v1 norm_median(): same idea as upqt but uses the median of
    nonzero values.
    """
    med = expr.apply(
        lambda col: np.median(col[col > 0].to_numpy()) if (col > 0).any() else 1.0
    )
    med_mean = med.mean()
    factor = med / med_mean
    factor = factor.replace(0.0, np.nan).fillna(1.0)
    return expr.divide(factor, axis=1)


def normalize_deseq(expr: pd.DataFrame) -> pd.DataFrame:
    """DESeq size-factor normalization.

    Mirrors GCEN v1 norm_deseq():
      1. For each gene, compute geometric mean across samples.
      2. For each sample, compute the ratio of its expression to the gene geometric mean.
      3. Size factor = median of these ratios across all genes with nonzero geometric mean.
      4. Divide each sample by its size factor, then rescale so the mean size factor is 1.
    """
    log_expr = np.log(expr.replace(0, np.nan).values)  # NaN for zeros
    # Geometric mean per gene = exp(mean of log values, ignoring NaN)
    log_geo_mean = np.nanmean(log_expr, axis=1)  # shape (n_genes,)

    # Ratio matrix: each gene-sample divided by gene geometric mean
    with np.errstate(invalid="ignore"):
        ratio = log_expr - log_geo_mean[:, np.newaxis]  # log-space subtraction

    # Median ratio per sample (only genes with nonzero geometric mean)
    valid_rows = np.isfinite(log_geo_mean)
    median_ratios = np.nanmedian(ratio[valid_rows, :], axis=0)  # shape (n_samples,)
    size_factors = np.exp(median_ratios)

    # Rescale so mean of size factors = 1
    size_factors = size_factors / size_factors.mean()
    result = expr.values / size_factors[np.newaxis, :]
    return pd.DataFrame(result, index=expr.index, columns=expr.columns)


def normalize_tmm(expr: pd.DataFrame) -> pd.DataFrame:
    """TMM (Trimmed Mean of M-values) normalization.

    Mirrors GCEN v1 norm_tmm() / calc_weighted_factor() exactly:
      1. Select reference sample = sample whose upper quartile is closest to the
         global mean of upper quartiles.
      2. For each sample, compute weighted trimmed mean of log₂ fold-change
         (M values) against the reference, trimming 30% of M and 5% of A.
      3. Normalize: divide each sample by 2^TMM factor.
      4. Rescale so geometric mean of factors = 1.
    """
    data = expr.values.astype(float)  # (n_genes, n_samples)
    n_samples = data.shape[1]

    # Step 1: select reference sample
    upqt_vals = np.array([
        _upper_quartile(data[:, j][data[:, j] > 0]) if (data[:, j] > 0).any() else 0.0
        for j in range(n_samples)
    ])
    ref_col = int(np.argmin(np.abs(upqt_vals - upqt_vals.mean())))

    ref = data[:, ref_col]
    ref_sum = ref.sum()

    # Step 2: compute TMM factor per sample
    factors = np.ones(n_samples)
    for j in range(n_samples):
        obs = data[:, j]
        obs_sum = obs.sum()
        if obs_sum == 0:
            continue

        # Only use genes where both obs and ref are > 0
        mask = (obs > 0) & (ref > 0)
        if mask.sum() < 3:
            continue

        obs_m = obs[mask]
        ref_m = ref[mask]
        obs_sum_m = obs_sum
        ref_sum_m = ref_sum

        M = np.log2((obs_m / obs_sum_m) / (ref_m / ref_sum_m))
        A = 0.5 * (np.log2(obs_m / obs_sum_m) + np.log2(ref_m / ref_sum_m))
        W = ((obs_sum_m - obs_m) / obs_sum_m / obs_m
             + (ref_sum_m - ref_m) / ref_sum_m / ref_m)

        n = float(len(M))
        lo_m = np.floor(n * 0.3) + 1.0
        hi_m = n + 1.0 - lo_m
        lo_a = np.floor(n * 0.05) + 1.0
        hi_a = n + 1.0 - lo_a

        M_rank = _get_rank(M)
        A_rank = _get_rank(A)
        keep = (M_rank >= lo_m) & (M_rank <= hi_m) & (A_rank >= lo_a) & (A_rank <= hi_a)

        if keep.sum() == 0:
            continue

        W_k = W[keep]
        M_k = M[keep]
        w_sum = (1.0 / W_k).sum()
        if w_sum == 0:
            continue

        factors[j] = 2.0 ** (np.sum(M_k / W_k) / w_sum)

    # Step 3: rescale so geometric mean of factors = 1
    log_factors = np.log(factors)
    factors = factors / np.exp(log_factors.mean())

    result = data / factors[np.newaxis, :]
    return pd.DataFrame(result, index=expr.index, columns=expr.columns)


def normalize_hkg(expr: pd.DataFrame, housekeeping_genes: list[str]) -> pd.DataFrame:
    """Housekeeping gene normalization.

    Mirrors GCEN v1 norm_hkg():
      1. For each sample, compute the median expression of housekeeping genes.
      2. Scale factor = global_mean_of_medians / per_sample_median.
      3. Multiply each sample by its scale factor.
    """
    present = [g for g in housekeeping_genes if g in expr.index]
    if not present:
        raise GCENInputError(
            "None of the specified housekeeping genes are present in the expression matrix."
        )
    hkg_medians = expr.loc[present].median(axis=0)
    global_mean = hkg_medians.mean()
    factor = global_mean / hkg_medians.replace(0, np.nan).fillna(1.0)
    return expr.multiply(factor, axis=1)


def filter_genes(
    expr: Union[ExpressionMatrix, pd.DataFrame],
    cutoff_mean: float = 0.0,
    cutoff_sd: float = 0.0,
    percent_mean: float = 1.0,
    percent_sd: float = 1.0,
) -> Union[ExpressionMatrix, pd.DataFrame]:
    """Filter genes by expression mean and standard deviation.

    Mirrors GCEN v1 data_filter logic.

    Bug fix vs GCEN v1: standard deviation uses ddof=1 (sample std dev).

    Parameters
    ----------
    expr:
        Expression matrix (genes × samples).
    cutoff_mean:
        Minimum absolute mean expression. Genes below this are removed.
    cutoff_sd:
        Minimum absolute standard deviation. Genes below this are removed.
    percent_mean:
        Retain genes with mean expression in the top ``percent_mean`` fraction.
        E.g. 0.75 keeps the top 75% (genes above the 25th percentile of means).
    percent_sd:
        Retain genes with std dev in the top ``percent_sd`` fraction.
    """
    df, is_wrapped = _unwrap(expr)

    gene_mean = df.mean(axis=1)
    # GCEN2 fix: use sample std dev (ddof=1); v1 used population std dev (ddof=0)
    gene_sd = df.std(axis=1, ddof=1)

    mean_pct_thresh = gene_mean.quantile(1.0 - percent_mean) if percent_mean < 1.0 else 0.0
    sd_pct_thresh   = gene_sd.quantile(1.0 - percent_sd)     if percent_sd   < 1.0 else 0.0

    mean_thresh = max(cutoff_mean, mean_pct_thresh)
    sd_thresh   = max(cutoff_sd,   sd_pct_thresh)

    mask = (gene_mean >= mean_thresh) & (gene_sd >= sd_thresh)
    result_df = df.loc[mask]

    if is_wrapped:
        return ExpressionMatrix(
            data=result_df,
            gene_names=np.array(result_df.index, dtype=str),
            sample_names=np.array(result_df.columns, dtype=str),
            is_normalized=getattr(expr, "is_normalized", False),
            normalization_method=getattr(expr, "normalization_method", None),
        )
    return result_df


def log_transform(
    expr: Union[ExpressionMatrix, pd.DataFrame],
    base: str = "log2",
    pseudo: float = 1.0,
) -> Union[ExpressionMatrix, pd.DataFrame]:
    """Apply log transformation: log_base(x + pseudo).

    Mirrors GCEN v1 string_to_log / string_to_log2 / string_to_log10.

    Parameters
    ----------
    base:
        'log' (natural log), 'log2', or 'log10'.
    pseudo:
        Pseudo-count added before taking the log (default 1.0 = GCEN v1 default).
    """
    df, is_wrapped = _unwrap(expr)

    funcs = {"log": np.log, "log2": np.log2, "log10": np.log10}
    if base not in funcs:
        raise GCENInputError(
            f"Unknown log base '{base}'. Choose from: {', '.join(funcs.keys())}"
        )

    result_df = funcs[base](df + pseudo)

    if is_wrapped:
        return ExpressionMatrix(
            data=result_df,
            gene_names=np.array(result_df.index, dtype=str),
            sample_names=np.array(result_df.columns, dtype=str),
            is_normalized=getattr(expr, "is_normalized", False),
            normalization_method=getattr(expr, "normalization_method", None),
        )
    return result_df


# ---------------------------------------------------------------------------
# Private helpers
# ---------------------------------------------------------------------------

def _unwrap(expr: Union[ExpressionMatrix, pd.DataFrame]) -> tuple[pd.DataFrame, bool]:
    if isinstance(expr, ExpressionMatrix):
        return expr.data, True
    return expr, False


def _upper_quartile(arr: np.ndarray) -> float:
    """Upper quartile (75th percentile) — mirrors GCEN v1 upper_quartile()."""
    if len(arr) == 0:
        return 0.0
    return float(np.percentile(np.sort(arr), 75))


def _get_rank(arr: np.ndarray) -> np.ndarray:
    """Average-rank transform with tie handling — mirrors GCEN v1 get_rank()."""
    n = len(arr)
    index = np.argsort(arr)[::-1]   # descending (largest rank = 1)
    rank = np.zeros(n)
    i = 0
    while i < n:
        j = i
        # find extent of tied group
        while j + 1 < n and arr[index[j]] == arr[index[j + 1]]:
            j += 1
        avg = (sum(range(i, j + 1))) / (j - i + 1) + 1.0
        for k in range(i, j + 1):
            rank[index[k]] = avg
        i = j + 1
    return rank
