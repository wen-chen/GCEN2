"""Tests for gcen2.spatial.modules (hdWGCNA spatial workflow) and gcen2
spatial-modules."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")
sp = pytest.importorskip("scipy.sparse")

from gcen2.core import _HAVE_CORE  # noqa: E402
from gcen2.exceptions import GCENInputError  # noqa: E402

REF = Path(__file__).parent / "data" / "hdwgcna_spatial"


@pytest.fixture(scope="module")
def ref_adata():
    """The synthetic Visium sample of make_reference.R."""
    t = pd.read_csv(REF / "counts.tsv", sep="\t", index_col=0)
    genes = [c for c in t.columns if c[0] == "g" and c[1:].isdigit()]
    obs = t[["row", "col", "group"]].rename(columns={"row": "array_row", "col": "array_col"})
    a = ad.AnnData(X=sp.csr_matrix(t[genes].to_numpy(dtype=np.float32)), obs=obs,
                   var=pd.DataFrame(index=genes))
    a.obsm["spatial"] = obs[["array_col", "array_row"]].to_numpy(dtype=float) * 100
    return a


class TestMetaspots:
    def test_same_as_hdwgcna(self, ref_adata):
        from gcen2.spatial import metaspots

        with pytest.warns(UserWarning, match="Group 'D'"):
            ms = metaspots(ref_adata, "group")
        want = pd.read_csv(REF / "members.tsv", sep="\t")
        assert list(ms.obs_names) == list(want["metaspot"])
        names = np.asarray(ref_adata.obs_names)
        for members, (centre, nb) in zip(ms.uns["metaspot_members"],
                                          want.itertuples(index=False, name=None)):
            assert list(names[members]) == [centre] + nb.split(",")
        assert (ms.obs["n_spots"] == [len(m) for m in ms.uns["metaspot_members"]]).all()

    def test_normalized_as_hdwgcna(self, ref_adata):
        from gcen2.spatial import metaspots
        from gcen2.spatial.modules import log_normalize

        with pytest.warns(UserWarning):
            ms = metaspots(ref_adata, "group")
        want = pd.read_csv(REF / "metaspots.tsv", sep="\t", index_col=0)
        np.testing.assert_allclose(log_normalize(ms.X), want.loc[ms.obs_names].to_numpy(),
                                   rtol=1e-12, atol=1e-12)

    def test_sum_average_and_hexagon(self):
        from gcen2.spatial import metaspots

        rc = [(r, c) for r in range(9) for c in range(17) if (r + c) % 2 == 0]
        obs = pd.DataFrame(rc, columns=["array_row", "array_col"],
                           index=[f"s{i}" for i in range(len(rc))])
        rng = np.random.default_rng(0)
        X = rng.poisson(3, size=(len(rc), 4)).astype(float)
        a = ad.AnnData(X=X, obs=obs)
        ms = metaspots(a, min_spots=10)
        avg = metaspots(a, min_spots=10, mode="average")
        pos = {(r, c): i for i, (r, c) in enumerate(rc)}
        # centres: every 4th column from col 3, rows of the same parity
        centres = {rc[m[0]] for m in ms.uns["metaspot_members"]}
        assert centres == {(r, c) for (r, c) in rc if c in (3, 7, 11, 15) and r % 2 == 1}
        for m, s, v in zip(ms.uns["metaspot_members"], ms.X, avg.X):
            r, c = rc[m[0]]
            hexagon = {(r + dr, c + dc) for dr, dc in
                       ((-1, -1), (-1, 1), (1, -1), (1, 1), (0, -2), (0, 2))}
            assert {rc[i] for i in m[1:]} == {p for p in hexagon if p in pos}
            np.testing.assert_allclose(s, X[m].sum(axis=0))
            np.testing.assert_allclose(v, X[m].mean(axis=0))

    def test_errors(self, ref_adata):
        from gcen2.spatial import metaspots

        with pytest.raises(GCENInputError, match="array_row"):
            metaspots(ref_adata, "group", row_key="nope", method="grid")
        with pytest.raises(GCENInputError, match="method"):
            metaspots(ref_adata, "group", method="hex")
        with pytest.raises(GCENInputError, match="groupby"):
            metaspots(ref_adata, "nope")
        with pytest.raises(GCENInputError, match="mode"):
            metaspots(ref_adata, "group", mode="max")
        with pytest.raises(GCENInputError, match="no group"), pytest.warns(UserWarning):
            metaspots(ref_adata, "group", min_spots=10**6)


def test_reassign_negative_kme():
    from gcen2.spatial.modules import reassign_negative_kme

    labels = np.array([1, 1, 2, 0, 2])
    kme = np.array([[0.5, 0.1], [-0.2, 0.3], [0.4, -0.1], [-0.5, -0.5], [-0.3, -0.2]])
    new, moved = reassign_negative_kme(labels, kme, [1, 2])
    assert list(new) == [1, 2, 1, 0, 0] and moved == [1, 2, 4]


class TestSpotEigengenes:
    def test_same_as_hdwgcna(self, ref_adata):
        from gcen2.spatial import spot_eigengenes
        from gcen2.spatial.modules import log_normalize

        want = pd.read_csv(REF / "spot_me.tsv", sep="\t", index_col=0)
        labels = np.array([1] * 10 + [2] * 10 + [0] * 10)
        me = spot_eigengenes(log_normalize(ref_adata.X), labels, [1, 2])
        np.testing.assert_allclose(me, want.loc[ref_adata.obs_names].to_numpy(),
                                   rtol=1e-8, atol=1e-8)

    def test_wide_matrix(self):
        """More genes than spots: the Gram matrix of the spots."""
        from gcen2.spatial.modules import module_eigengene_scores

        x = np.random.default_rng(1).normal(size=(8, 30))
        z = (x - x.mean(axis=0)) / x.std(axis=0, ddof=1)  # ScaleData (nothing above 10)
        u, s, _ = np.linalg.svd(z, full_matrices=False)
        pc = module_eigengene_scores(x)
        np.testing.assert_allclose(np.abs(pc), np.abs(u[:, 0] * s[0]), rtol=1e-10)


@pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")
class TestSpatialModules:
    def test_modules_and_spot_eigengenes(self, ref_adata):
        from gcen2.spatial import group_trait_cor, spatial_modules
        from gcen2.spatial.modules import log_normalize
        from gcen2.tom_modules import HDWGCNA, tom_modules

        with pytest.warns(UserWarning):
            r = spatial_modules(ref_adata, "group", fraction=0.0, min_size=5, threads=2)
        md = r.metadata
        assert md["groups"] == ["A", "B", "C"] and md["n_spots"] == ref_adata.n_obs
        lab = pd.Series(r.labels, index=r.genes)
        for block in (range(1, 11), range(11, 21)):  # the two latent factors
            b = lab[[f"g{i}" for i in block]]
            assert (b > 0).all() and b.nunique() == 1
        # = tom_modules on the metaspots, then hdWGCNA's reassignment by spot kME
        again = tom_modules(md["metaspots"], threads=2, **{**HDWGCNA, "min_size": 5})
        assert again.power == r.power
        with pytest.warns(UserWarning):
            r0 = spatial_modules(ref_adata, "group", fraction=0.0, min_size=5, threads=2,
                                 reassign=False)
        assert np.array_equal(again.labels, r0.labels) and r0.metadata["reassigned"] == []
        kme = md["spot_kme"]
        np.testing.assert_allclose(kme.to_numpy(), r0.metadata["spot_kme"].to_numpy())
        levels = [int(c[3:]) for c in kme.columns]
        lab = pd.Series(r0.labels, index=r0.genes)
        moved = [g for g in r.genes if lab[g] > 0 and kme.loc[g, f"kME{lab[g]}"] < 0]
        assert md["reassigned"] == moved
        for g in r.genes:
            want = lab[g]
            if g in moved:
                best = kme.loc[g].to_numpy()
                want = levels[int(np.argmax(best))] if best.max() >= 0 else 0
            assert r.labels[r.genes.index(g)] == want
        me = md["spot_eigengenes"]
        assert list(me.index) == list(ref_adata.obs_names)
        assert list(me.columns) == [f"ME{k}" for k in r.me_labels if k > 0]
        tc = group_trait_cor(me, ref_adata.obs["group"])
        assert len(tc) == me.shape[1] * 4 and (tc["n"] == ref_adata.n_obs).all()
        # spot eigengene of a module: first PC of its scaled genes
        g = [i for i, k in enumerate(r.labels) if k == 1]
        x = log_normalize(ref_adata.X)[:, g]
        assert abs(np.corrcoef(me["ME1"], x.mean(axis=1))[0, 1]) > 0.5

    def test_gene_selection(self, ref_adata):
        from gcen2.single_cell import fraction_genes
        from gcen2.spatial import spatial_modules

        kw = dict(min_size=5, threads=2, spot_level=False)
        with pytest.warns(UserWarning):
            r = spatial_modules(ref_adata, "group", genes=["g1", "g2", "g3", "g12", "g13"],
                                **{**kw, "min_size": 2})
        assert r.genes == ["g1", "g2", "g3", "g12", "g13"]
        assert "spot_eigengenes" not in r.metadata
        with pytest.raises(GCENInputError, match="not in adata"):
            spatial_modules(ref_adata, "group", genes=["g1", "zz"], **kw)
        X = ref_adata.X
        keep = fraction_genes(X, 0.9)
        n = X.shape[0]
        pos = np.asarray((X > 0).sum(axis=0)).ravel()
        assert np.array_equal(keep, pos >= round(0.9 * n))
        dense = np.zeros((10, 2))
        dense[:3, 0] = 1  # 3 of 10: kept at fraction 0.25 (round(2.5) = 2)
        dense[:2, 1] = 1  # 2 of 10
        assert list(fraction_genes(dense, 0.25)) == [True, True]
        assert list(fraction_genes(dense, 0.3)) == [True, False]

    def test_needs_counts(self, ref_adata):
        from gcen2.spatial import spatial_modules

        b = ref_adata.copy()
        b.X = b.X * 0.5
        with pytest.raises(GCENInputError, match="raw counts"):
            spatial_modules(b, "group")


@pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")
def test_cli(tmp_path, ref_adata):
    from typer.testing import CliRunner

    from gcen2.cli import app

    ref_adata.write_h5ad(tmp_path / "x.h5ad")
    out = tmp_path / "out"
    r = CliRunner().invoke(app, [
        "spatial-modules", "-i", str(tmp_path / "x.h5ad"), "--groupby", "group", "-o",
        str(out), "--fraction", "0", "--min-size", "5", "--n-hubs", "3", "--plot",
        "--write-metaspots"])
    assert r.exit_code == 0, r.output
    for f in ("modules.tsv", "eigengenes.tsv", "kme.tsv", "hub_genes.tsv",
              "spot_eigengenes.tsv", "group_cor.tsv", "metaspots.tsv", "summary.tsv",
              "spatial_eigengenes.png"):
        assert (out / f).is_file(), f
    s = pd.read_csv(out / "summary.tsv", sep="\t")
    assert s.loc[0, "groups"] == "A,B,C" and s.loc[0, "modules"] >= 2
    me = pd.read_csv(out / "spot_eigengenes.tsv", sep="\t", index_col=0)
    assert list(me.columns) == [f"ME{i}" for i in range(s.loc[0, "modules"])]
    gc = pd.read_csv(out / "group_cor.tsv", sep="\t")
    assert set(gc["trait"]) == {"A", "B", "C", "D"}
    assert "WARNING" in r.output and "Group 'D'" in r.output
    assert s.loc[0, "metaspots"] == "grid"
    bad = CliRunner().invoke(app, ["spatial-modules", "-i", str(tmp_path / "x.h5ad"),
                                   "-o", str(out), "--row-key", "nope", "--metaspots", "grid"])
    assert bad.exit_code == 1 and "nope" in bad.output


# --------------------------------------------------------------------------
# knn metaspots (platforms without a Visium grid)
# --------------------------------------------------------------------------
def _scatter(n=600, seed=0):
    """Counts of spots at random positions (no grid): two halves with
    different gene programmes."""
    rng = np.random.default_rng(seed)
    xy = rng.uniform(0, 1000, size=(n, 2))
    left = xy[:, 0] < 500
    base = rng.gamma(2.0, 1.0, size=(n, 40))
    base[left, :15] *= 4
    base[~left, 15:30] *= 4
    X = rng.poisson(base).astype(np.float32)
    a = ad.AnnData(X=sp.csr_matrix(X), obs=pd.DataFrame(
        {"half": np.where(left, "L", "R")}, index=[f"s{i}" for i in range(n)]),
        var=pd.DataFrame(index=[f"g{j}" for j in range(40)]))
    a.obsm["spatial"] = xy
    return a


class TestKnnMetaspots:
    def test_members_and_overlap(self):
        from scipy.spatial import cKDTree

        from gcen2.spatial import metaspots

        a = _scatter()
        ms = metaspots(a, "half", k=7, max_shared=2, min_spots=10)
        assert ms.uns["metaspot_params"]["method"] == "knn"
        xy, half = a.obsm["spatial"], a.obs["half"].to_numpy()
        X = a.X.toarray()
        names = np.asarray(a.obs_names)
        members = ms.uns["metaspot_members"]
        assert len(members) > 50
        for m, s, centre in zip(members, ms.X.toarray(), ms.obs_names):
            assert len(m) == 7 and len(set(m)) == 7 and names[m[0]] == centre
            same = np.flatnonzero(half == half[m[0]])
            d, nn = cKDTree(xy[same]).query(xy[m[0]], k=7)
            assert set(same[nn]) == set(m)              # the 7 nearest of its group
            np.testing.assert_allclose(s, X[m].sum(axis=0))
        sets = [set(m.tolist()) for m in members]
        assert max(len(x & y) for i, x in enumerate(sets) for y in sets[:i]) <= 2
        again = metaspots(a, "half", min_spots=10)
        assert all(np.array_equal(x, y) for x, y in zip(members, again.uns["metaspot_members"]))
        other = metaspots(a, "half", min_spots=10, seed=1)
        assert list(other.obs_names) != list(ms.obs_names)
        avg = metaspots(a, "half", min_spots=10, mode="average")
        np.testing.assert_allclose(avg.X.toarray(), ms.X.toarray() / 7)

    def test_auto_and_errors(self, ref_adata):
        from gcen2.spatial import metaspots

        a = _scatter(n=200)
        with pytest.warns(UserWarning):
            assert metaspots(ref_adata, "group").uns["metaspot_params"]["method"] == "grid"
        # a square grid (Visium HD bins) is not a hexagonal one: knn
        sq = a.copy()
        sq.obs["array_row"] = np.repeat(np.arange(10), 20)
        sq.obs["array_col"] = np.tile(np.arange(20), 10)
        sq.obsm["spatial"] = sq.obs[["array_col", "array_row"]].to_numpy(dtype=float)
        assert metaspots(sq, min_spots=10).uns["metaspot_params"]["method"] == "knn"
        assert metaspots(a, max_shared=0, min_spots=10).n_obs <= 200 // 7
        no_xy = a.copy()
        del no_xy.obsm["spatial"]
        with pytest.raises(GCENInputError, match="spatial"):
            metaspots(no_xy)
        with pytest.raises(GCENInputError, match="max_shared"):
            metaspots(a, k=5, max_shared=5)
        with pytest.raises(GCENInputError, match="k must"):
            metaspots(a, k=1)
        bad = a.copy()
        bad.obsm["spatial"][0, 0] = np.nan
        with pytest.raises(GCENInputError, match="missing"):
            metaspots(bad)

    @pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")
    def test_spatial_modules_and_cli(self, tmp_path):
        from typer.testing import CliRunner

        from gcen2.cli import app
        from gcen2.spatial import spatial_modules

        a = _scatter()
        res = spatial_modules(a, "half", fraction=0, min_spots=10, min_size=5, threads=2)
        md = res.metadata
        assert md["metaspot_params"]["method"] == "knn" and md["n_metaspots"] > 50
        lab = dict(zip(res.genes, res.labels))
        prog_l = {lab[f"g{j}"] for j in range(15)} - {0}
        prog_r = {lab[f"g{j}"] for j in range(15, 30)} - {0}
        assert prog_l and prog_r and not prog_l & prog_r
        a.write_h5ad(tmp_path / "x.h5ad")
        out = tmp_path / "out"
        r = CliRunner().invoke(app, [
            "spatial-modules", "-i", str(tmp_path / "x.h5ad"), "--groupby", "half", "-o",
            str(out), "--fraction", "0", "--min-size", "5", "--min-spots", "10", "-k", "9",
            "--max-shared", "3", "--plot"])
        assert r.exit_code == 0, r.output
        assert "knn metaspots" in r.output and (out / "spatial_eigengenes.png").is_file()
        s = pd.read_csv(out / "summary.tsv", sep="\t")
        assert s.loc[0, "metaspots"] == "knn"
        bad = CliRunner().invoke(app, ["spatial-modules", "-i", str(tmp_path / "x.h5ad"),
                                       "-o", str(out), "--metaspots", "grid"])
        assert bad.exit_code == 1 and "array_row" in bad.output
