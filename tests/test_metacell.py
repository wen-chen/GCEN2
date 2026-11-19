"""Tests for gcen2.single_cell.metacells and metacell cell-type networks."""

import warnings

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")
sp = pytest.importorskip("scipy.sparse")

from gcen2.exceptions import GCENInputError  # noqa: E402
from gcen2.single_cell import build_by_celltype, metacells  # noqa: E402


def _modules_adata(n_cells=900, n_mod=3, per=15, extra=30, seed=1, n_types=2):
    """Poisson counts whose gene modules follow latent factors."""
    rng = np.random.default_rng(seed)
    n_genes = n_mod * per + extra
    F = rng.normal(size=(n_cells, n_mod))
    logmu = np.log(rng.uniform(0.2, 1.0, n_genes))[None, :].repeat(n_cells, 0)
    for m in range(n_mod):
        logmu[:, m * per:(m + 1) * per] += 0.6 * F[:, [m]]
    lib = rng.lognormal(0, 0.3, n_cells)
    X = rng.poisson(np.exp(logmu) * lib[:, None]).astype(np.float32)
    genes = [f"g{i}" for i in range(n_genes)]
    truth = {(genes[i], genes[j]) for m in range(n_mod)
             for i in range(m * per, (m + 1) * per) for j in range(i + 1, (m + 1) * per)}
    obs = pd.DataFrame({"ct": [f"T{i % n_types}" for i in range(n_cells)],
                        "sample": [f"s{i % 3}" for i in range(n_cells)]},
                       index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame({"highly_variable": [True] * (n_mod * per) + [False] * extra},
                       index=genes)
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var), X, truth


def _overlaps(members):
    sets = [set(r) for r in members]
    return max((len(a & b) for i, a in enumerate(sets) for b in sets[i + 1:]), default=0)


class TestMetacells:
    def test_structure_and_constraints(self):
        adata, X, _ = _modules_adata()
        m = metacells(adata, groupby="ct", k=10, max_shared=3, max_metacells=40,
                      normalize=False)
        members = m.uns["metacell_members"]
        assert members.shape == (m.n_obs, 10)
        assert list(m.var_names) == list(adata.var_names)
        assert set(m.obs.columns) == {"ct", "n_cells", "center"}
        for ct in ("T0", "T1"):
            rows = members[(m.obs["ct"] == ct).to_numpy()]
            assert 0 < len(rows) <= 40
            # members from their own group only, distinct, centre first
            assert (adata.obs["ct"].to_numpy()[rows] == ct).all()
            assert all(len(set(r)) == 10 for r in rows)
            assert _overlaps(rows) <= 3
        np.testing.assert_array_equal(m.obs["center"].to_numpy(),
                                      adata.obs_names.to_numpy()[members[:, 0]])
        # mean of the member cells
        np.testing.assert_allclose(
            m.X, np.stack([X[r].astype(np.float64).mean(axis=0) for r in members]),
            rtol=1e-12)
        s = metacells(adata, groupby="ct", k=10, max_shared=3, max_metacells=40,
                      aggregate="sum", normalize=False)
        np.testing.assert_allclose(s.X, m.X * 10, rtol=1e-12)
        # counts: log-normalised by default (10,000 per metacell, log1p)
        n = metacells(adata, groupby="ct", k=10, max_shared=3, max_metacells=40)
        np.testing.assert_array_equal(n.uns["metacell_members"], members)
        ref = m.X / m.X.sum(axis=1, keepdims=True) * 1e4
        np.testing.assert_allclose(n.X, np.log1p(ref), rtol=1e-12)
        assert n.uns["metacell_params"]["normalize"] is True
        # values that are not counts are left as they are
        logged = adata.copy()
        logged.X = np.log1p(X)
        lm = metacells(logged, groupby="ct", k=10, max_shared=3, max_metacells=40)
        assert lm.uns["metacell_params"]["normalize"] is False

    def test_deterministic_and_dense_equals_sparse(self):
        adata, X, _ = _modules_adata(n_cells=300)
        a = metacells(adata, groupby="ct", k=8, max_metacells=20, seed=3)
        b = metacells(adata, groupby="ct", k=8, max_metacells=20, seed=3)
        np.testing.assert_array_equal(a.uns["metacell_members"], b.uns["metacell_members"])
        dense = adata.copy()
        dense.X = X
        c = metacells(dense, groupby="ct", k=8, max_metacells=20, seed=3)
        np.testing.assert_array_equal(a.uns["metacell_members"], c.uns["metacell_members"])
        d = metacells(adata, groupby="ct", k=8, max_metacells=20, seed=4)
        assert not np.array_equal(a.uns["metacell_members"], d.uns["metacell_members"])

    def test_uses_given_representation(self):
        adata, _, _ = _modules_adata(n_cells=200)
        # cells ordered on a line: neighbourhoods are runs of consecutive cells
        adata.obsm["X_line"] = np.arange(200, dtype=float)[:, None]
        m = metacells(adata, k=5, max_shared=0, use_rep="X_line", max_metacells=1000)
        for row in m.uns["metacell_members"]:
            assert np.ptp(row) == 4
        assert 20 <= m.n_obs <= 40  # disjoint runs of 5
        with pytest.raises(GCENInputError, match="obsm"):
            metacells(adata, use_rep="nope")

    def test_groupby_several_keys_and_small_groups(self):
        adata, _, _ = _modules_adata(n_cells=300)
        m = metacells(adata, groupby=["ct", "sample"], k=10, min_cells=10)
        pairs = set(zip(m.obs["ct"], m.obs["sample"]))
        assert len(pairs) == 6
        adata.obs["uneven"] = ["A"] * 250 + ["B"] * 50
        with pytest.warns(UserWarning, match="skipped"):
            m = metacells(adata, groupby="uneven", k=10, min_cells=100)
        assert set(m.obs["uneven"]) == {"A"}
        with pytest.raises(GCENInputError, match="No group"), warnings.catch_warnings():
            warnings.simplefilter("ignore")
            metacells(adata, groupby="ct", k=10, min_cells=1000)

    @pytest.mark.parametrize("kwargs, match", [
        ({"k": 1}, "k must"), ({"max_shared": 25}, "max_shared"),
        ({"aggregate": "median"}, "aggregate"), ({"groupby": "nope"}, "groupby"),
        ({"layer": "nope"}, "Layer"), ({"max_metacells": 0}, "max_metacells"),
    ])
    def test_bad_arguments(self, kwargs, match):
        with pytest.raises(GCENInputError, match=match):
            metacells(_modules_adata(n_cells=100)[0], **kwargs)


class TestCelltypeMetacells:
    def test_recovers_modules_better_than_pseudo_replicates(self):
        adata, _, truth = _modules_adata(n_cells=1200, n_types=1, extra=1000)

        def precision(strategy):
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                net = build_by_celltype(adata, "ct", strategy=strategy, method="pearson",
                                        cor=0.0, pval=1.0, threads=2)["T0"]
            e = net.edges.assign(a=net.edges["correlation"].abs())
            top = e.sort_values("a", ascending=False).head(len(truth))
            return np.mean([(a, b) in truth or (b, a) in truth
                            for a, b in zip(top["gene_a"], top["gene_b"])]), net

        p_meta, net = precision("metacell")
        p_rand, _ = precision("pseudo_replicates")
        assert p_meta > 0.9 and p_meta > p_rand + 0.5
        assert net.metadata["strategy"] == "metacell"
        assert net.metadata["n_cells"] == 1200 and net.metadata["n_samples"] > 100

    def test_pseudo_replicates_deprecated(self):
        adata, _, _ = _modules_adata(n_cells=200, n_types=1)
        with pytest.warns(FutureWarning, match="pseudo_replicates' is deprecated"):
            build_by_celltype(adata, "ct", strategy="pseudo_replicates", method="pearson",
                              cor=0.0, pval=1.0, threads=1)

    def test_hvg_layer_and_small_types(self):
        adata, X, _ = _modules_adata(n_cells=400)
        adata.layers["counts"] = sp.csr_matrix(X * 2)
        with pytest.warns(UserWarning, match="only"):
            nets = build_by_celltype(adata, "ct", use_highly_variable=True, layer="counts",
                                     k=20, max_metacells=8, cor=0.0, pval=1.0)
        assert set(nets) == {"T0", "T1"}
        genes = set(nets["T0"].edges["gene_a"]) | set(nets["T0"].edges["gene_b"])
        assert genes <= set(adata.var_names[adata.var["highly_variable"]])
        with pytest.warns(UserWarning, match="Skipping"):
            assert build_by_celltype(adata, "ct", min_cells=500, cor=0.0, pval=1.0) == {}

    def test_bad_strategy(self):
        with pytest.raises(GCENInputError, match="strategy"):
            build_by_celltype(_modules_adata(n_cells=100)[0], "ct", strategy="x")


def _two_type_h5ad(path, G=1500, seed=5):
    """T cells: modules g0-g19 and g20-g39; B cells: g0-g19 and g40-g59."""
    rng = np.random.default_rng(seed)

    def cells(n, modules):
        F = rng.normal(size=(n, len(modules)))
        lib = rng.lognormal(0, 0.3, n)
        logmu = np.log(rng.uniform(0.05, 1.5, G))[None, :].repeat(n, 0)
        for m, (lo, hi) in enumerate(modules):
            logmu[:, lo:hi] += 0.7 * F[:, [m]]
        return rng.poisson(np.exp(logmu) * lib[:, None])

    X = np.vstack([cells(700, [(0, 20), (20, 40)]),
                   cells(600, [(0, 20), (40, 60)])]).astype(np.float32)
    labels = ["T cell"] * 700 + ["B/cell"] * 600
    labels[:30] = ["rare"] * 30
    adata = ad.AnnData(X=sp.csr_matrix(X),
                       obs=pd.DataFrame({"cell_type": labels},
                                        index=[f"c{i}" for i in range(len(labels))]),
                       var=pd.DataFrame({"highly_variable": np.arange(G) < 100},
                                        index=[f"g{i}" for i in range(G)]))
    adata.write_h5ad(path)
    return adata


class TestCelltypeCommand:
    @pytest.fixture
    def runner(self):
        typer_testing = pytest.importorskip("typer.testing")
        return typer_testing.CliRunner()

    def test_networks_and_diff(self, tmp_path, runner):
        from gcen2.cli import app

        _two_type_h5ad(tmp_path / "x.h5ad")
        out = tmp_path / "out"
        r = runner.invoke(app, ["celltype", "-i", str(tmp_path / "x.h5ad"), "--groupby",
                                "cell_type", "-o", str(out), "-m", "pearson", "-c", "0.3",
                                "-p", "0.001", "--diff", "-t", "2", "-k", "20"])
        assert r.exit_code == 0, r.output
        assert "Cell type 'rare' has 30 cells" in r.output
        table = pd.read_csv(out / "celltypes.tsv", sep="\t")
        assert list(table["cell_type"]) == ["T cell", "B/cell"]
        assert list(table["file"]) == ["T_cell.network", "B_cell.network"]
        assert list(table["n_cells"]) == [670, 600]
        # the median gene's effective samples (per gene pair, from the overlap)
        assert all(4 <= n <= m for n, m in zip(table["effective_samples"], table["n_samples"]))
        nets = {ct: gcen2_read(out / f) for ct, f in zip(table["cell_type"], table["file"])}
        for ct, mods in {"T cell": [(0, 20), (20, 40)], "B/cell": [(0, 20), (40, 60)]}.items():
            e = nets[ct].edges
            a = e["gene_a"].str[1:].astype(int)
            b = e["gene_b"].str[1:].astype(int)
            true = np.zeros(len(e), dtype=bool)
            for lo, hi in mods:
                true |= (a >= lo) & (a < hi) & (b >= lo) & (b < hi)
            # p <= 0.001 on valid p-values: at most ~0.1% of the other pairs
            assert true.sum() >= 0.6 * 380 and (~true).sum() <= 0.001 * 1500 * 1499 / 2
        summary = pd.read_csv(out / "diff" / "summary.tsv", sep="\t")
        assert len(summary) == 1 and summary["file"][0] == "T_cell__vs__B_cell.tsv"
        d = pd.read_csv(out / "diff" / summary["file"][0], sep="\t")
        a = d["gene_a"].str[1:].astype(int)
        b = d["gene_b"].str[1:].astype(int)
        t_only = (a >= 20) & (a < 40) & (b >= 20) & (b < 40)
        b_only = (a >= 40) & (a < 60) & (b >= 40) & (b < 60)
        assert d.loc[t_only, "category"].isin(["A_only", "shared_changed"]).mean() > 0.9
        assert d.loc[b_only, "category"].isin(["B_only", "shared_changed"]).mean() > 0.9

    def test_groups_hvg_and_pseudo_replicates(self, tmp_path, runner):
        from gcen2.cli import app

        _two_type_h5ad(tmp_path / "x.h5ad")
        r = runner.invoke(app, ["celltype", "-i", str(tmp_path / "x.h5ad"), "--groupby",
                                "cell_type", "-o", str(tmp_path / "o"), "--group", "B/cell",
                                "--hvg", "--strategy", "pseudo_replicates", "--replicates", "8",
                                "-c", "0.5", "-p", "0.05"])
        assert r.exit_code == 0, r.output
        table = pd.read_csv(tmp_path / "o" / "celltypes.tsv", sep="\t")
        assert list(table["cell_type"]) == ["B/cell"]
        assert list(table["n_samples"]) == [8] and list(table["effective_samples"]) == [8]
        assert not (tmp_path / "o" / "diff").exists()

    @pytest.mark.parametrize("args, message", [
        (["--groupby", "nope"], "not a column"),
        (["--groupby", "cell_type", "--group", "NK"], "no cells"),
        (["--groupby", "cell_type", "--min-cells", "5000"], "no cell type"),
        (["--groupby", "cell_type", "--strategy", "x"], "strategy"),
    ])
    def test_errors(self, tmp_path, runner, args, message):
        from gcen2.cli import app

        _two_type_h5ad(tmp_path / "x.h5ad", G=80)
        r = runner.invoke(app, ["celltype", "-i", str(tmp_path / "x.h5ad"),
                                "-o", str(tmp_path / "o"), *args])
        assert r.exit_code == 1 and message in r.output


def gcen2_read(path):
    from gcen2.io import read_network

    return read_network(path)


def test_write_celltype_networks_names_and_small_diff(tmp_path):
    from gcen2.models import Network
    from gcen2.single_cell import effective_samples, write_celltype_networks

    edges = pd.DataFrame({"gene_a": ["A"], "gene_b": ["B"], "correlation": [0.9],
                          "p_value": [1e-5]})
    nets = {}
    for ct, n in [("a/b", 100), ("a b", 100), ("..", 30)]:
        nets[ct] = Network(edges=edges.copy(), metadata={
            "n_cells": n * 25, "n_samples": n, "strategy": "metacell", "k": 25})
    nets[".."].metadata["n_cells"] = 50  # 2 effective samples
    assert effective_samples(nets["a/b"]) == 100 and effective_samples(nets[".."]) == 2
    with pytest.warns(UserWarning, match="fewer than 4"):
        table = write_celltype_networks(nets, tmp_path, diff=True)
    assert list(table["file"]) == ["a_b.network", "a_b_2.network", "celltype.network"]
    summary = pd.read_csv(tmp_path / "diff" / "summary.tsv", sep="\t")
    assert len(summary) == 1 and summary["shared_stable"][0] == 1


@pytest.mark.parametrize("effective", ["overlap", "cells"])
def test_celltype_diff_uses_the_networks_sample_sizes(tmp_path, effective):
    """The diff between cell types tests every edge on the sample sizes of
    the networks' own p-values: per gene pair (smaller gene) for overlap,
    n_cells // k for cells."""
    from gcen2.differential import _fisher_z_arrays
    from gcen2.single_cell import effective_samples, write_celltype_networks

    adata = _two_type_h5ad(tmp_path / "x.h5ad", G=120)
    adata = adata[adata.obs["cell_type"] != "rare"].copy()
    nets = build_by_celltype(adata, "cell_type", k=20, method="pearson", cor=0.1, pval=0.5,
                             fdr=False, effective=effective, threads=2)
    write_celltype_networks(nets, tmp_path / "o", diff=True, pval=0.05)
    summary = pd.read_csv(tmp_path / "o" / "diff" / "summary.tsv", sep="\t")
    d = pd.read_csv(tmp_path / "o" / "diff" / summary["file"][0], sep="\t")
    d = d[d["category"].str.startswith("shared")]
    assert len(d) > 100
    n = {}
    for x, ct in zip("ab", ["T cell", "B/cell"]):
        md = nets[ct].metadata
        if effective == "overlap":
            g = md["gene_effective_samples"]
            n[x] = np.minimum(g[d.gene_a].to_numpy(), g[d.gene_b].to_numpy())
            assert summary[f"n_{x}"][0] == md["p_value_samples"]
            assert n[x].min() < n[x].max()
        else:
            n[x] = effective_samples(nets[ct])
            assert summary[f"n_{x}"][0] == n[x] == md["n_cells"] // 20
    _, p = _fisher_z_arrays(d.cor_a.to_numpy(), d.cor_b.to_numpy(), n["a"], n["b"])
    np.testing.assert_allclose(d.p_value, p, rtol=1e-6)


@pytest.mark.parametrize("fdr", [False, True])
def test_effective_pvalues_match_direct_computation(fdr):
    from scipy import stats

    import gcen2
    from gcen2.network import _benjamini_hochberg

    adata, _, _ = _modules_adata(n_cells=600, n_types=1, extra=60)
    kw = dict(method="pearson", cor=0.2, pval=0.01, fdr=fdr, threads=2)
    net = build_by_celltype(adata, "ct", k=20, effective="cells", **kw)["T0"]
    n_eff = net.metadata["p_value_samples"]
    assert n_eff == 600 // 20
    meta = metacells(adata, groupby="ct", k=20)
    full = gcen2.network.build(meta, method="pearson", all_edges=True).edges
    r = full["correlation"].to_numpy()
    p = 2 * stats.t.sf(np.abs(r) * np.sqrt((n_eff - 2) / (1 - r * r)), n_eff - 2)
    keep = (np.abs(r) >= 0.2) & (p <= 0.01)
    if fdr:
        q = _benjamini_hochberg(p, len(p))
        keep &= q <= 0.01
    ref = full[keep].reset_index(drop=True)
    pd.testing.assert_frame_equal(net.edges[["gene_a", "gene_b", "correlation"]],
                                  ref[["gene_a", "gene_b", "correlation"]])
    np.testing.assert_allclose(net.edges["p_value"], p[keep], rtol=1e-10)
    if fdr:
        np.testing.assert_allclose(net.edges["fdr"], q[keep], rtol=1e-10)


@pytest.mark.parametrize("fdr", [False, True])
def test_overlap_pvalues_per_pair(fdr):
    """effective="overlap": a pair is tested on the smaller of its genes'
    metacell_effective_samples."""
    from scipy import stats

    import gcen2
    from gcen2.network import _benjamini_hochberg
    from gcen2.single_cell import metacell_effective_samples

    adata, _, _ = _modules_adata(n_cells=600, n_types=1, extra=60)
    kw = dict(method="pearson", cor=0.2, pval=0.01, fdr=fdr, threads=2)
    net = build_by_celltype(adata, "ct", k=20, **kw)["T0"]
    meta = metacells(adata, groupby="ct", k=20)
    ne = metacell_effective_samples(pd.DataFrame(meta.X.T, index=meta.var_names),
                                    meta.uns["metacell_members"], 20)
    pd.testing.assert_series_equal(net.metadata["gene_effective_samples"], ne)
    assert net.metadata["p_value_samples"] == int(np.median(ne))
    assert (ne <= meta.n_obs).all() and (ne >= 3).all()
    full = gcen2.network.build(meta, method="pearson", all_edges=True).edges
    r = full["correlation"].to_numpy()
    n = np.minimum(ne[full["gene_a"]].to_numpy(), ne[full["gene_b"]].to_numpy())
    p = 2 * stats.t.sf(np.abs(r) * np.sqrt((n - 2) / (1 - r * r)), n - 2)
    keep = (np.abs(r) >= 0.2) & (p <= 0.01)
    if fdr:
        q = _benjamini_hochberg(p, len(p))
        keep &= q <= 0.01
    ref = full[keep].reset_index(drop=True)
    pd.testing.assert_frame_equal(net.edges[["gene_a", "gene_b", "correlation"]],
                                  ref[["gene_a", "gene_b", "correlation"]])
    np.testing.assert_allclose(net.edges["p_value"], p[keep], rtol=1e-10)


def test_overlap_effective_samples_calibrated():
    """Genes unrelated to everything: their correlations over overlapping
    metacells have the variance the per-gene effective samples predict, and
    the p-values are uniform (not too small)."""
    from scipy import stats

    from gcen2.single_cell import metacell_effective_samples

    rng = np.random.default_rng(0)
    n, g = 800, 120
    pcs = rng.normal(size=(n, 5))
    state = pcs[:, :1] @ rng.normal(size=(1, 20))  # 20 genes follow the state
    X = np.hstack([np.exp(0.5 * state), np.ones((n, g - 20))])
    X = rng.poisson(X).astype(np.float32)
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame({"ct": ["a"] * n},
                                                           index=[f"c{i}" for i in range(n)]),
                       obsm={"X_pca": pcs})
    adata.var_names = [f"g{i}" for i in range(g)]
    meta = metacells(adata, groupby="ct", k=10, use_rep="X_pca", n_pcs=5)
    x = np.asarray(meta.X, dtype=float).T
    ne = metacell_effective_samples(pd.DataFrame(x, index=adata.var_names),
                                    meta.uns["metacell_members"], 10)
    assert ne.iloc[:20].median() < ne.iloc[20:].median()  # state genes: fewer
    z = (x - x.mean(1, keepdims=True)) / x.std(1, keepdims=True)
    r = (z @ z.T) / x.shape[1]
    iu = np.triu_indices(g, 1)
    null = (iu[0] >= 20)  # pairs of two noise genes or a state and a noise gene
    rr = r[iu][null]
    m = np.minimum(ne.to_numpy()[iu[0]], ne.to_numpy()[iu[1]])[null]
    p = 2 * stats.t.sf(np.abs(rr) * np.sqrt((m - 2) / (1 - rr * rr)), m - 2)
    assert 0.02 < (p <= 0.05).mean() < 0.08
    assert meta.n_obs > n // 10 * 3  # heavy overlap: n // k would be far too strict


def test_too_few_effective_samples():
    adata, _, _ = _modules_adata(n_cells=200, n_types=1)
    with pytest.warns(UserWarning, match="effective samples"):
        assert build_by_celltype(adata, "ct", k=60, min_cells=10, cor=0.0, pval=1.0) == {}


@pytest.mark.parametrize("use_pca", [True, False])
def test_groups_same_as_subset(use_pca):
    """metacells(groups=...) equals metacells of the subset, bit for bit (the
    PCA, when computed, uses only those groups' cells)."""
    a, _, _ = _modules_adata(n_cells=900, n_types=3)
    if use_pca:
        a.obsm["X_pca"] = np.random.default_rng(0).normal(size=(a.n_obs, 10))
    labels = a.obs["ct"].astype(str)
    keep = sorted(labels.unique())[:2]
    sub = a[labels.isin(keep).to_numpy()].copy()
    want = metacells(sub, groupby="ct", k=10, min_cells=50)
    got = metacells(a, groupby="ct", k=10, min_cells=50, groups=keep)
    assert np.array_equal(got.X, want.X)
    pos = np.flatnonzero(labels.isin(keep).to_numpy())
    assert np.array_equal(got.uns["metacell_members"], pos[want.uns["metacell_members"]])
    assert list(got.obs["center"]) == list(want.obs["center"])
    with pytest.raises(GCENInputError, match="groupby column"):
        metacells(a, groupby=None, groups=keep)
