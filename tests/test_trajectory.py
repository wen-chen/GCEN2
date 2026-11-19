"""Tests for gcen2.trajectory (sliding-window networks along pseudotime)."""

import numpy as np
import pandas as pd
import pytest

import gcen2
from gcen2.exceptions import GCENInputError
from gcen2.trajectory import PATTERNS, build_trajectory

try:
    import anndata as ad

    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False

requires_ad = pytest.mark.skipif(not HAS_ANNDATA, reason="anndata not installed")


def _cells(n=600, g=30, seed=0):
    """Cells x genes counts: g0-g1 co-expressed early, g2-g3 late, g4-g5 always."""
    rng = np.random.default_rng(seed)
    t = rng.uniform(0, 1, n)
    X = rng.poisson(3, (n, g)).astype(float)
    f = rng.normal(size=n) * 5
    h = rng.normal(size=n) * 5
    for a, b, on in ((0, 1, t < 0.4), (2, 3, t > 0.6), (4, 5, np.ones(n, bool))):
        sig = np.where(on, f if a != 4 else h, 0.0)
        X[:, a] += sig + 10
        X[:, b] += sig + 10
    return np.clip(X, 0, None), t


def _frame(X, t):
    cols = [f"c{i}" for i in range(X.shape[0])]
    expr = pd.DataFrame(X.T, index=[f"g{j}" for j in range(X.shape[1])], columns=cols)
    return expr, pd.Series(t, index=cols)


def _pattern(traj, a, b):
    d = traj.dynamics
    row = d[(d["gene_a"] == a) & (d["gene_b"] == b)]
    assert len(row) == 1
    return row.iloc[0]


class TestBuildTrajectory:
    def test_patterns_found(self):
        expr, pt = _frame(*_cells())
        traj = build_trajectory(expr, pt, window_size=0.2, cor=0.5)
        assert _pattern(traj, "g0", "g1")["pattern"] == "lost"
        assert _pattern(traj, "g2", "g3")["pattern"] == "gained"
        assert _pattern(traj, "g4", "g5")["pattern"] == "stable"
        assert set(traj.dynamics["pattern"]) <= set(PATTERNS)
        assert _pattern(traj, "g2", "g3")["trend"] > 0.5
        assert _pattern(traj, "g0", "g1")["cor_change"] < -0.5

    def test_windows_cover_all_cells(self):
        expr, pt = _frame(*_cells(n=103))
        traj = build_trajectory(expr, pt, window_size=20, step=15, cor=0.5)
        t = traj.window_table()
        assert (t["n_cells"] == 20).all()
        assert t["start"].iloc[0] == pt.min() and t["end"].iloc[-1] == pt.max()
        # starts 0, 15, ..., 75 and a last window at 83 so the tail is covered
        assert len(t) == 7
        assert t["center"].is_monotonic_increasing

    def test_window_networks_equal_network_build(self):
        X, t = _cells(n=200)
        expr, pt = _frame(X, t)
        traj = build_trajectory(expr, pt, window_size=50, step=50, method="pearson",
                                cor=0.3, pval=0.01)
        order = np.argsort(t, kind="stable")
        for w in traj.windows:
            cols = expr.columns[order[w.index * 50:(w.index + 1) * 50]]
            ref = gcen2.network.build(expr[cols], method="pearson", cor=0.3, pval=0.01)
            pd.testing.assert_frame_equal(w.network.edges, ref.edges, check_exact=False,
                                          rtol=1e-12, atol=1e-12)

    @pytest.mark.parametrize("method,log", [("spearman", None), ("pearson", "log2")])
    def test_dynamics_correlations_match_windows(self, method, log):
        expr, pt = _frame(*_cells(n=300))
        traj = build_trajectory(expr, pt, window_size=60, method=method, log=log, cor=0.2)
        assert len(traj.dynamics) > 0
        for w in traj.windows:
            m = traj.dynamics.merge(w.network.edges, on=["gene_a", "gene_b"])
            assert len(m) == w.network.edge_count
            np.testing.assert_allclose(m[f"cor_{w.index}"], m["correlation"], atol=1e-9)

    def test_n_windows_first_last(self):
        expr, pt = _frame(*_cells())
        traj = build_trajectory(expr, pt, window_size=0.2, cor=0.5)
        present = np.zeros((len(traj.dynamics), len(traj.windows)), bool)
        key = traj.dynamics["gene_a"] + "|" + traj.dynamics["gene_b"]
        for w in traj.windows:
            present[:, w.index] = key.isin(w.network.edges["gene_a"] + "|"
                                           + w.network.edges["gene_b"])
        assert (traj.dynamics["n_windows"] == present.sum(1)).all()
        assert (traj.dynamics["first"] == present.argmax(1)).all()

    def test_metacells(self):
        expr, pt = _frame(*_cells(n=400))
        traj = build_trajectory(expr, pt, window_size=100, metacell_size=8, cor=0.5)
        assert all(w.n_samples == 13 for w in traj.windows)  # ceil(100 / 8)

    def test_pool_averages_consecutive_cells(self):
        from gcen2.trajectory import _pool

        block = np.arange(20.0).reshape(10, 2)
        pooled = _pool(block, 4)
        np.testing.assert_array_equal(pooled, [block[0:4].mean(0), block[4:8].mean(0),
                                               block[8:10].mean(0)])

    def test_nan_pseudotime_dropped(self):
        expr, pt = _frame(*_cells(n=200))
        pt.iloc[:10] = np.nan
        with pytest.warns(UserWarning, match="10 cells"):
            traj = build_trajectory(expr, pt, window_size=50, cor=0.5)
        assert traj.params["n_cells"] == 190

    def test_degree_dynamics(self):
        expr, pt = _frame(*_cells())
        traj = build_trajectory(expr, pt, window_size=0.2, cor=0.5)
        deg = traj.degree_dynamics()
        assert list(deg.columns) == [f"degree_{w.index}" for w in traj.windows]
        for w in traj.windows:
            assert deg[f"degree_{w.index}"].sum() == 2 * w.network.edge_count

    def test_write(self, tmp_path):
        expr, pt = _frame(*_cells(n=200))
        traj = build_trajectory(expr, pt, window_size=0.5, cor=0.5)
        traj.write(tmp_path / "out")
        names = {p.name for p in (tmp_path / "out").iterdir()}
        assert {"windows.tsv", "edge_dynamics.tsv", "degree_dynamics.tsv"} <= names
        assert {f"window_{w.index}.network" for w in traj.windows} <= names
        back = gcen2.io.read_network(tmp_path / "out" / "window_0.network")
        assert len(back.edges) == traj.windows[0].network.edge_count

    @pytest.mark.parametrize("kwargs,match", [
        ({"window_size": 2.0}, "fraction"),
        ({"window_size": 0}, "positive"),
        ({"window_size": 1000}, "larger"),
        ({"window_size": 4, "metacell_size": 2}, "at least 3"),
        ({"metacell_size": 0}, "metacell_size"),
    ])
    def test_bad_arguments(self, kwargs, match):
        expr, pt = _frame(*_cells(n=100))
        with pytest.raises(GCENInputError, match=match):
            build_trajectory(expr, pt, **kwargs)

    def test_pseudotime_length_checked(self):
        expr, _ = _frame(*_cells(n=100))
        with pytest.raises(GCENInputError, match="pseudotime has"):
            build_trajectory(expr, np.arange(99.0))

    def test_series_missing_sample(self):
        expr, pt = _frame(*_cells(n=100))
        with pytest.raises(GCENInputError, match="No pseudotime"):
            build_trajectory(expr, pt.iloc[1:])


@requires_ad
class TestAnnData:
    def _adata(self):
        X, t = _cells()
        obs = pd.DataFrame({"pt": t}, index=[f"c{i}" for i in range(len(t))])
        var = pd.DataFrame(index=[f"g{j}" for j in range(X.shape[1])])
        import scipy.sparse as sp

        return ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var), X, t

    def test_same_as_dataframe(self):
        adata, X, t = self._adata()
        a = build_trajectory(adata, "pt", window_size=0.25, cor=0.5)
        expr, pt = _frame(X, t)
        b = build_trajectory(expr, pt, window_size=0.25, cor=0.5)
        pd.testing.assert_frame_equal(a.dynamics, b.dynamics)

    def test_missing_obs_column(self):
        adata, _, _ = self._adata()
        with pytest.raises(GCENInputError, match="not a column"):
            build_trajectory(adata, "nope")

    def test_highly_variable(self):
        adata, _, _ = self._adata()
        adata.var["highly_variable"] = [j < 6 for j in range(adata.n_vars)]
        traj = build_trajectory(adata, "pt", window_size=0.2, cor=0.5,
                                use_highly_variable=True)
        genes = set(traj.dynamics["gene_a"]) | set(traj.dynamics["gene_b"])
        assert genes <= {f"g{j}" for j in range(6)}


def test_plot(tmp_path):
    pytest.importorskip("matplotlib")
    import matplotlib

    matplotlib.use("Agg")
    expr, pt = _frame(*_cells())
    traj = build_trajectory(expr, pt, window_size=0.2, cor=0.5)
    fig = gcen2.plotting.trajectory_dynamics(traj, top=3, output=tmp_path / "t.png")
    assert (tmp_path / "t.png").stat().st_size > 0
    assert len(fig.axes[1].lines) >= 3


def test_cli(tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    expr, pt = _frame(*_cells(n=300))
    expr.to_csv(tmp_path / "expr.tsv", sep="\t", lineterminator="\n")
    pt.rename("pseudotime").to_csv(tmp_path / "pt.tsv", sep="\t", index_label="cell", lineterminator="\n")
    out = tmp_path / "traj"
    res = CliRunner().invoke(app, ["trajectory", "-i", str(tmp_path / "expr.tsv"),
                                   "--pseudotime", str(tmp_path / "pt.tsv"),
                                   "-o", str(out), "-w", "60", "-c", "0.5"])
    assert res.exit_code == 0, res.output
    assert "Edge patterns:" in res.output
    ref = build_trajectory(expr, pt, window_size=60, cor=0.5)
    dyn = pd.read_csv(out / "edge_dynamics.tsv", sep="\t")
    assert len(dyn) == len(ref.dynamics)
    assert list(dyn["pattern"]) == list(ref.dynamics["pattern"])
