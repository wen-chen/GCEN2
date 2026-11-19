"""Tests for the dropout-aware correlation methods (pearson_nonzero,
pearson_pairwise_complete) and gcen2.single_cell.build_dropout_aware."""

import warnings

import numpy as np
import pandas as pd
import pytest
from scipy import stats

import gcen2
from gcen2.exceptions import GCENInputError
from gcen2.single_cell import build_dropout_aware

core = pytest.importorskip("gcen2.core._gcen_core")

SUBSETS = {
    "pearson_nonzero": lambda x, y: (x != 0) & (y != 0),
    "pearson_pairwise_complete": lambda x, y: (x != 0) | (y != 0),
}


def _sparse_matrix(G=40, N=300, seed=0):
    rng = np.random.default_rng(seed)
    X = rng.poisson(rng.uniform(0.05, 3, (G, 1)), (G, N)).astype(float)
    X = np.log1p(X * rng.uniform(0.5, 2, (1, N)))
    X[5] = 0                          # all zero
    X[6] = 0
    X[6, :2] = 1.0                    # two non-zeros
    X[7] = 0
    X[7, :20] = 2.0                   # constant where non-zero
    return np.ascontiguousarray(X)


@pytest.mark.parametrize("method", list(SUBSETS))
def test_matches_pearson_on_the_subset(method):
    X = _sparse_matrix()
    idx, r, p = core.correlation.pairwise_correlation(X, method, 4)
    assert len(r) == 40 * 39 // 2
    checked = 0
    for (i, j), rr, pp in zip(idx, r, p):
        m = SUBSETS[method](X[i], X[j])
        if m.sum() < 3 or np.ptp(X[i][m]) == 0 or np.ptp(X[j][m]) == 0:
            assert rr == 0 and pp == 1
            continue
        er, ep = stats.pearsonr(X[i][m], X[j][m])
        assert rr == pytest.approx(er, abs=1e-12)
        # v1's incomplete beta (relative accuracy ~1e-7), as for pearson
        assert pp == pytest.approx(ep, rel=1e-6, abs=1e-12)
        checked += 1
    assert checked > 600


@pytest.mark.parametrize("method", list(SUBSETS))
def test_no_zeros_is_pearson(method):
    Y = np.random.default_rng(1).normal(5, 1, (30, 50))
    a = core.correlation.pairwise_correlation(Y, "pearson", 2)
    b = core.correlation.pairwise_correlation(Y, method, 2)
    np.testing.assert_array_equal(a[0], b[0])
    np.testing.assert_allclose(b[1], a[1], rtol=0, atol=1e-13)
    np.testing.assert_allclose(b[2], a[2], rtol=1e-9)


@pytest.mark.parametrize("method", list(SUBSETS))
def test_filters_threads_and_blocks(method):
    X = _sparse_matrix(G=60, seed=3)
    full = core.correlation.pairwise_correlation(X, method, 1)
    for threads in (2, 7):
        got = core.correlation.pairwise_correlation(X, method, threads)
        for u, v in zip(full, got):
            np.testing.assert_array_equal(u, v)
    keep = (np.abs(full[1]) >= 0.2) & (full[2] <= 0.01) & (full[1] >= 0)
    got = core.correlation.pairwise_correlation(X, method, 3, min_abs_cor=0.2, max_p=0.01,
                                                positive_only=True)
    np.testing.assert_array_equal(got[0], full[0][keep])
    parts = [core.correlation.pairwise_correlation(X, method, 2, row_begin=b, row_end=e)
             for b, e in [(0, 10), (10, 45), (45, 60)]]
    np.testing.assert_array_equal(np.concatenate([q[1] for q in parts]), full[1])


def test_build_and_build_to_file(tmp_path):
    X = _sparse_matrix()
    genes = [f"g{i}" for i in range(len(X))]
    df = pd.DataFrame(X, index=genes, columns=[f"c{k}" for k in range(X.shape[1])])
    net = gcen2.network.build(df, method="pearson_nonzero", cor=0.3, pval=0.05)
    assert net.method == "pearson_nonzero" and net.edge_count > 0
    gcen2.network.build_to_file(df, tmp_path / "a.tsv", method="pearson_nonzero", cor=0.3,
                                pval=0.05, max_pairs_per_block=100)
    net.to_edgelist(tmp_path / "b.tsv")
    got = pd.read_csv(tmp_path / "a.tsv", sep="\t")
    assert list(zip(got.gene_a, got.gene_b)) == list(zip(net.edges.gene_a, net.edges.gene_b))
    # the TSV path of gcen2 build takes the numpy route for these methods
    df.to_csv(tmp_path / "x.tsv", sep="\t", lineterminator="\n")
    from gcen2._build_file import build_from_file
    info = build_from_file(tmp_path / "x.tsv", tmp_path / "c.tsv", method="pearson_nonzero",
                           cor=0.3, pval=0.05)
    assert info["edge_count"] == net.edge_count
    assert (tmp_path / "c.tsv").read_bytes() == (tmp_path / "a.tsv").read_bytes()


def test_unknown_method():
    with pytest.raises(GCENInputError, match="Unknown method"):
        gcen2.network.build(pd.DataFrame(np.ones((3, 4))), method="kendall")
    with pytest.raises(GCENInputError, match="method must be"):
        gcen2.trajectory.build_trajectory(pd.DataFrame(np.ones((3, 6))), np.arange(6.0),
                                          window_size=4, method="pearson_nonzero")


def _dropout_data(drop, lo, hi, seed=2, n_cells=1500, n_mod=3, per=15, extra=40):
    """Modules on latent factors, Poisson counts, then random technical dropout."""
    rng = np.random.default_rng(seed)
    G = n_mod * per + extra
    F = rng.normal(size=(n_cells, n_mod))
    lib = rng.lognormal(0, 0.3, n_cells)
    logmu = np.log(rng.uniform(lo, hi, G))[None, :].repeat(n_cells, 0)
    for m in range(n_mod):
        logmu[:, m * per:(m + 1) * per] += 0.6 * F[:, [m]]
    X = rng.poisson(np.exp(logmu) * lib[:, None]).astype(float)
    X[rng.random(X.shape) < drop] = 0
    genes = [f"g{i}" for i in range(G)]
    truth = {(genes[i], genes[j]) for m in range(n_mod)
             for i in range(m * per, (m + 1) * per) for j in range(i + 1, (m + 1) * per)}
    return pd.DataFrame(np.log1p(X.T / lib[None, :]), index=genes), truth


def _precision(net, truth):
    e = net.edges.assign(a=net.edges["correlation"].abs())
    top = e.sort_values("a", ascending=False).head(len(truth))
    return np.mean([(a, b) in truth for a, b in zip(top["gene_a"], top["gene_b"])])


def test_both_nonzero_helps_with_technical_dropout():
    df, truth = _dropout_data(drop=0.7, lo=1, hi=5)
    plain = _precision(gcen2.network.build(df, method="pearson", cor=0, pval=1), truth)
    nz = _precision(build_dropout_aware(df, cor=0, pval=1), truth)
    assert nz > 0.9 and nz > plain + 0.3


def test_pairwise_complete_bias_and_warning():
    df, truth = _dropout_data(drop=0.0, lo=0.2, hi=1.0)
    with pytest.warns(UserWarning, match="negative"):
        net = build_dropout_aware(df, zero_handling="pairwise_complete", cor=0, pval=1)
    t = np.array([(a, b) in truth for a, b in zip(net.edges.gene_a, net.edges.gene_b)])
    assert net.edges["correlation"][~t].mean() < -0.2   # the documented bias
    assert net.metadata["zero_handling"] == "pairwise_complete"


def test_build_dropout_aware_arguments():
    df, _ = _dropout_data(drop=0.5, lo=1, hi=5, n_cells=100)
    with pytest.raises(GCENInputError, match="pearson"):
        build_dropout_aware(df, method="spearman")
    with pytest.raises(GCENInputError, match="zero_handling"):
        build_dropout_aware(df, zero_handling="x")
    net = build_dropout_aware(df, method="pearson_nonzero", cor=0.1, pval=0.5)
    assert net.method == "pearson_nonzero" and net.metadata["zero_handling"] == "both_nonzero"


def test_build_dropout_aware_anndata():
    ad = pytest.importorskip("anndata")
    sp = pytest.importorskip("scipy.sparse")
    df, _ = _dropout_data(drop=0.5, lo=1, hi=5, n_cells=200)
    adata = ad.AnnData(X=sp.csr_matrix(df.to_numpy().T),
                       var=pd.DataFrame(index=df.index))
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        a = build_dropout_aware(adata, cor=0.1, pval=0.5)
    b = build_dropout_aware(df, cor=0.1, pval=0.5)
    pd.testing.assert_frame_equal(a.edges, b.edges)


def test_impute_warns_p_values_invalid(monkeypatch):
    """zero_handling='impute' warns that the p-values of imputed cells are not
    valid (a stand-in for MAGIC: magic-impute is an optional extra)."""
    import sys
    import types

    class FakeMAGIC:
        def fit_transform(self, x):
            return x

    monkeypatch.setitem(sys.modules, "magic", types.SimpleNamespace(MAGIC=FakeMAGIC))
    rng = np.random.default_rng(0)
    df = pd.DataFrame(rng.normal(size=(5, 40)), index=[f"g{i}" for i in range(5)])
    with pytest.warns(UserWarning, match="not independent samples"):
        net = build_dropout_aware(df, zero_handling="impute", cor=0, pval=1)
    assert net.metadata["zero_handling"] == "impute"

