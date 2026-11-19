"""Tests for gcen2.single_cell.modules_by_celltype and gcen2 celltype-modules."""

import numpy as np
import pandas as pd
import pytest

ad = pytest.importorskip("anndata")
sp = pytest.importorskip("scipy.sparse")

from gcen2.core import _HAVE_CORE  # noqa: E402
from gcen2.exceptions import GCENInputError  # noqa: E402

pytestmark = pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")


def _adata(n_cells=1200, seed=3):
    """Two cell types; Poisson counts with three gene modules each (latent
    factors), plus noise genes and genes that are rarely expressed."""
    rng = np.random.default_rng(seed)
    n_mod, per, noise, rare = 3, 25, 40, 10
    n_genes = n_mod * per + noise + rare
    X = np.zeros((n_cells, n_genes), dtype=np.float32)
    types = np.array(["A", "B"])[np.arange(n_cells) % 2]
    for t in ("A", "B"):
        cells = np.flatnonzero(types == t)
        F = rng.normal(size=(len(cells), n_mod))
        logmu = np.log(rng.uniform(0.5, 2.0, n_genes))[None, :].repeat(len(cells), 0)
        for m in range(n_mod):
            logmu[:, m * per:(m + 1) * per] += 0.8 * F[:, [m]]
        logmu[:, -rare:] = np.log(0.01)
        X[cells] = rng.poisson(np.exp(logmu))
    genes = [f"g{i}" for i in range(n_genes)]
    obs = pd.DataFrame({"ct": types}, index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame({"highly_variable": [True] * (n_genes - rare) + [False] * rare},
                       index=genes)
    return ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var), per, n_mod, rare


@pytest.fixture(scope="module")
def data():
    return _adata()


KW = dict(k=10, max_shared=3, max_metacells=200, min_cells=50, min_size=10, threads=2)


class TestModulesByCelltype:
    def test_modules_found(self, data):
        from gcen2.single_cell import modules_by_celltype

        adata, per, n_mod, rare = data
        res = modules_by_celltype(adata, "ct", **KW)
        assert set(res) == {"A", "B"}
        for t, r in res.items():
            md = r.metadata
            assert md["cell_type"] == t and md["n_cells"] == 600
            assert md["n_metacells"] == r.eigengenes.shape[1] == len(r.samples)
            assert md["n_genes"] == len(r.genes) == adata.n_vars - rare  # rare genes dropped
            lab = pd.Series(r.labels, index=r.genes)
            for m in range(n_mod):  # each true module mostly in one found module
                block = lab[[f"g{i}" for i in range(m * per, (m + 1) * per)]]
                assert (block > 0).mean() > 0.8 and block[block > 0].nunique() == 1

    def test_same_as_module_on_metacells(self, data):
        from gcen2.single_cell import modules_by_celltype
        from gcen2.tom_modules import HDWGCNA, tom_modules

        adata, *_ = data
        r = modules_by_celltype(adata, "ct", groups=["B"], reassign=False, **KW)["B"]
        m = r.metadata["metacells"]
        again = tom_modules(m, threads=2, **{**HDWGCNA, "min_size": 10})
        assert np.array_equal(again.labels, r.labels) and again.power == r.power
        assert list(m.columns) == r.samples

    def test_fraction_and_hvg(self, data):
        from gcen2.single_cell import expressed_fraction, modules_by_celltype

        adata, _, _, rare = data
        f = expressed_fraction(adata, (adata.obs["ct"] == "A").to_numpy())
        dense = expressed_fraction(ad.AnnData(adata.X.toarray()),
                                   (adata.obs["ct"] == "A").to_numpy())
        np.testing.assert_allclose(f, dense)
        assert (f[-rare:] < 0.05).all() and (f[:-rare] >= 0.05).all()
        r = modules_by_celltype(adata, "ct", groups="A", fraction=0.0, **KW)["A"]
        assert len(r.genes) == adata.n_vars
        r = modules_by_celltype(adata, "ct", groups="A", fraction=0.0,
                                use_highly_variable=True, **KW)["A"]
        assert len(r.genes) == adata.n_vars - rare

    def test_hub_genes(self, data):
        from gcen2.single_cell import modules_by_celltype

        adata, *_ = data
        r = modules_by_celltype(adata, "ct", groups="A", **KW)["A"]
        hubs = r.hub_genes(5)
        kme, _ = r.kme()
        col = {k: i for i, k in enumerate(r.me_labels)}
        for mod, grp in hubs.groupby("module"):
            assert list(grp["rank"]) == list(range(1, len(grp) + 1)) and len(grp) == 5
            genes = np.flatnonzero(r.labels == mod + 1)
            best = np.sort(kme[genes, col[mod + 1]])[::-1][:5]
            np.testing.assert_allclose(grp["kME"], best)

    def test_small_types_and_errors(self, data):
        from gcen2.single_cell import modules_by_celltype

        adata, *_ = data
        with pytest.warns(UserWarning, match="Skipping"):
            assert modules_by_celltype(adata, "ct", **{**KW, "min_cells": 1000}) == {}
        with pytest.raises(GCENInputError):
            modules_by_celltype(adata, "ct", groups=["C"], **KW)
        with pytest.raises(GCENInputError):
            modules_by_celltype(adata, "nope", **KW)
        with pytest.raises(GCENInputError):
            modules_by_celltype(adata, "ct", preset="x", **KW)


class TestCli:
    def test_celltype_modules(self, tmp_path, data):
        from typer.testing import CliRunner

        from gcen2.cli import app
        from gcen2.single_cell import modules_by_celltype

        adata, *_ = data
        adata.write_h5ad(tmp_path / "x.h5ad")
        out = tmp_path / "out"
        r = CliRunner().invoke(app, [
            "celltype-modules", "-i", str(tmp_path / "x.h5ad"), "--groupby", "ct", "-o",
            str(out), "-k", "10", "--max-shared", "3", "--max-metacells", "200",
            "--min-cells", "50", "--min-size", "10", "--n-hubs", "3", "--write-metacells"])
        assert r.exit_code == 0, r.output
        table = pd.read_csv(out / "celltypes.tsv", sep="\t")
        assert list(table["cell_type"]) == ["A", "B"]
        want = modules_by_celltype(adata, "ct", **KW)
        for t in ("A", "B"):
            d = out / t
            for f in ("modules.tsv", "eigengenes.tsv", "kme.tsv", "hub_genes.tsv",
                      "metacells.tsv", "cell_eigengenes.tsv", "cell_kme.tsv"):
                assert (d / f).is_file()
            me = pd.read_csv(d / "cell_eigengenes.tsv", sep="\t", index_col=0)
            assert me.shape == (600, want[t].metadata["cell_eigengenes"].shape[1])
            np.testing.assert_allclose(
                me.to_numpy(), want[t].metadata["cell_eigengenes"].to_numpy(), rtol=1e-9)
            mods = pd.read_csv(d / "modules.tsv", sep="\t")
            assert [sorted(mods.loc[mods["module"] == i, "gene"])
                    for i in range(mods["module"].max() + 1)] == want[t].modules()
            hubs = pd.read_csv(d / "hub_genes.tsv", sep="\t")
            assert (hubs.groupby("module").size() == 3).all()
            mc = pd.read_csv(d / "metacells.tsv", sep="\t", index_col=0)
            assert mc.shape == (want[t].metadata["n_genes"], want[t].metadata["n_metacells"])
            mem = pd.read_csv(d / "metacell_members.tsv", sep="\t")
            assert list(mem.columns) == ["metacell", "cell", "center"]
            assert set(mem["metacell"]) == set(mc.columns)
            assert (mem.groupby("metacell").size() == 10).all()  # k cells each
            assert (mem.groupby("metacell")["center"].sum() == 1).all()
            assert set(adata.obs.loc[mem["cell"], "ct"]) == {t}  # cells of that type
        assert "power" in r.output


class TestCellLevel:
    def test_cell_eigengenes_as_spatial(self, data):
        """Eigengenes: the spatial (hdWGCNA ModuleEigengenes) code on the
        log-normalized cells, bit for bit; kME: dense Pearson to 1e-12."""
        from gcen2.preservation import _cor as cor_columns  # columns of x with those of y
        from gcen2.single_cell import cell_eigengenes
        from gcen2.spatial import spot_eigengenes
        from gcen2.spatial.modules import log_normalize

        adata, per, n_mod, rare = data
        cells = np.flatnonzero(adata.obs["ct"] == "A")
        idx = np.arange(n_mod * per + 20)
        labels = np.r_[np.repeat(np.arange(1, n_mod + 1), per), np.zeros(20, int)]
        me, kme = cell_eigengenes(adata.X[cells], idx, labels, [1, 2, 3])
        dense = log_normalize(adata.X[cells], idx)
        assert np.array_equal(me, spot_eigengenes(dense, labels, [1, 2, 3]))
        np.testing.assert_allclose(kme, cor_columns(dense, me), atol=1e-12)
        # a module's genes correlate with its eigengene
        assert (kme[np.arange(per), 0] > 0.3).mean() > 0.9
        raw, _ = cell_eigengenes(dense, idx, labels, [1, 2, 3], normalize=False)
        assert np.array_equal(raw, me)

    def test_metadata_and_reassign(self, data):
        from gcen2.single_cell import modules_by_celltype
        from gcen2.spatial.modules import reassign_negative_kme

        adata, *_ = data
        keep = modules_by_celltype(adata, "ct", groups="A", reassign=False, **KW)["A"]
        r = modules_by_celltype(adata, "ct", groups="A", **KW)["A"]
        md = r.metadata
        levels = [k for k in keep.me_labels if k > 0]
        me, kme = md["cell_eigengenes"], md["cell_kme"]
        assert list(me.index) == list(adata.obs_names[adata.obs["ct"] == "A"])
        assert list(me.columns) == [f"ME{k}" for k in levels]
        good = np.flatnonzero(keep._good)
        want, moved = reassign_negative_kme(keep.labels[good], kme.to_numpy(), levels)
        assert np.array_equal(r.labels[good], want)
        assert md["reassigned"] == [keep.genes[good[i]] for i in moved]
        off = modules_by_celltype(adata, "ct", groups="A", cell_level=False, **KW)["A"]
        assert "cell_eigengenes" not in off.metadata
        assert np.array_equal(off.labels, keep.labels)

    def test_cli_options(self, tmp_path, data):
        from typer.testing import CliRunner

        from gcen2.cli import app

        adata, *_ = data
        adata.write_h5ad(tmp_path / "x.h5ad")
        base = ["celltype-modules", "-i", str(tmp_path / "x.h5ad"), "--groupby", "ct",
                "-k", "10", "--max-shared", "3", "--max-metacells", "200", "--min-cells",
                "50", "--min-size", "10", "--group", "A"]
        r = CliRunner().invoke(app, [*base, "-o", str(tmp_path / "on")])
        assert r.exit_code == 0, r.output
        assert "reassigned" in r.output
        assert "reassigned" in pd.read_csv(tmp_path / "on/celltypes.tsv", sep="\t").columns
        kme = pd.read_csv(tmp_path / "on/A/cell_kme.tsv", sep="\t")
        assert list(kme.columns[:2]) == ["gene", "module"]
        r = CliRunner().invoke(app, [*base, "-o", str(tmp_path / "off"),
                                     "--no-cell-eigengenes"])
        assert r.exit_code == 0, r.output
        assert not (tmp_path / "off/A/cell_eigengenes.tsv").exists()


def test_backed_same_as_in_memory(tmp_path, data):
    """A backed AnnData (X on disk, one cell type's rows read at a time) gives
    the same results as one in memory."""
    from gcen2.single_cell import modules_by_celltype

    adata, *_ = data
    adata.write_h5ad(tmp_path / "x.h5ad")
    mem = modules_by_celltype(adata, "ct", **KW)
    backed = ad.read_h5ad(tmp_path / "x.h5ad", backed="r")
    disk = modules_by_celltype(backed, "ct", **KW)
    assert set(mem) == set(disk)
    for t in mem:
        assert np.array_equal(mem[t].labels, disk[t].labels)
        assert mem[t].metadata["metacells"].equals(disk[t].metadata["metacells"])
        assert mem[t].metadata["cell_eigengenes"].equals(disk[t].metadata["cell_eigengenes"])
    csc = adata.copy()
    csc.X = sp.csc_matrix(csc.X)
    csc.write_h5ad(tmp_path / "csc.h5ad")
    with pytest.raises(GCENInputError, match="CSR"):
        modules_by_celltype(ad.read_h5ad(tmp_path / "csc.h5ad", backed="r"), "ct", **KW)


def test_positive_counts_blocks():
    from gcen2.single_cell.modules import _positive_counts

    rng = np.random.default_rng(0)
    x = rng.poisson(0.3, size=(500, 40)).astype(np.float32)
    m = sp.csr_matrix(x)
    m.data[::7] = 0  # explicit zeros are not expressed
    want = (m.toarray() > 0).sum(axis=0)
    assert np.array_equal(_positive_counts(m), want)
    assert np.array_equal(_positive_counts(m.toarray()), want)
    assert np.array_equal(_positive_counts(sp.csc_matrix(m)), want)


class TestConsensusAndHarmony:
    def _with_samples(self, data):
        adata, *_ = data
        a = adata.copy()
        rng = np.random.default_rng(5)
        a.obs["sample"] = np.array(["s1", "s2", "s3"])[np.arange(a.n_obs) % 3]
        a.obs["donor"] = rng.choice(["d1", "d2", "d3", "d4"], a.n_obs)
        return a

    def test_consensus_by(self, data):
        from gcen2.single_cell import modules_by_celltype

        a = self._with_samples(data)
        _, per, n_mod, _ = data
        res = modules_by_celltype(a, "ct", consensus_by="sample", **KW)
        for r in res.values():
            md = r.metadata
            assert md["consensus"] and md["sets"] == ["s1", "s2", "s3"]
            assert sum(len(v) for v in md["set_samples"].values()) == md["n_metacells"]
            assert list(md["metacells"].columns) == r.samples
            lab = pd.Series(r.labels, index=r.genes)
            for m in range(n_mod):
                assert lab[[f"g{i}" for i in range(m * per, (m + 1) * per)]].value_counts() \
                    .iloc[0] >= per * 0.7

    def test_harmonize_by(self, data):
        from gcen2.single_cell import modules_by_celltype

        a = self._with_samples(data)
        res = modules_by_celltype(a, "ct", groups="A", harmonize_by="donor", **KW)
        md = res["A"].metadata
        h, me = md["cell_harmonized_eigengenes"], md["cell_eigengenes"]
        assert h.shape == me.shape and list(h.columns) == list(me.columns)
        for c in me.columns:  # same signal, little batch effect to remove
            assert np.corrcoef(h[c], me[c])[0, 1] > 0.9
        with pytest.raises(GCENInputError):
            modules_by_celltype(a, "ct", harmonize_by="donor", cell_level=False, **KW)
        with pytest.raises(GCENInputError):
            modules_by_celltype(a, "ct", consensus_by="nope", **KW)

    def test_cli(self, tmp_path, data):
        from typer.testing import CliRunner

        from gcen2.cli import app

        a = self._with_samples(data)
        path = tmp_path / "a.h5ad"
        a.write_h5ad(path)
        r = CliRunner().invoke(app, [
            "celltype-modules", "-i", str(path), "--groupby", "ct", "-o", str(tmp_path / "o"),
            "-k", "10", "--max-shared", "3", "--max-metacells", "200", "--min-cells", "50",
            "--min-size", "10", "--group", "A", "--consensus-by", "sample",
            "--harmonize-by", "donor"])
        assert r.exit_code == 0, r.output
        d = tmp_path / "o" / "A"
        sets = pd.read_csv(d / "sets.tsv", sep="\t")
        assert list(sets["set"]) == ["s1", "s2", "s3"]
        assert (d / "set_eigengenes.tsv").exists()
        assert (d / "cell_harmonized_eigengenes.tsv").exists()
        tab = pd.read_csv(tmp_path / "o" / "celltypes.tsv", sep="\t")
        assert tab.loc[0, "power"].count(",") == 2
        r = CliRunner().invoke(app, [
            "celltype-modules", "-i", str(path), "--groupby", "ct", "-o", str(tmp_path / "x"),
            "--consensus-quantile", "0.2"])
        assert r.exit_code == 1


def test_dmes_and_traits(tmp_path, data):
    from typer.testing import CliRunner

    from gcen2.cli import app
    from gcen2.single_cell import modules_by_celltype

    adata, *_ = data
    a = adata.copy()
    rng = np.random.default_rng(2)
    a.obs["cond"] = np.array(["ctl", "trt"])[np.arange(a.n_obs) % 4 // 2]
    a.obs["age"] = rng.normal(50, 10, a.n_obs)
    a.obs["sex"] = pd.Categorical(rng.choice(["f", "m"], a.n_obs))
    res = modules_by_celltype(a, "ct", groups="A", dme_by="cond", traits=["age", "sex"],
                              trait_groupby="cond", **KW)
    md = res["A"].metadata
    k = md["cell_eigengenes"].shape[1]
    assert set(md["dmes"]["group"]) == {"ctl", "trt"} and len(md["dmes"]) == 2 * k
    assert set(md["trait_cor"]["group"]) == {"all_cells", "ctl", "trt"}
    assert len(md["trait_cor"]) == 3 * 2 * k
    with pytest.raises(GCENInputError):
        modules_by_celltype(a, "ct", traits=["age"], cell_level=False, **KW)
    with pytest.raises(GCENInputError):
        modules_by_celltype(a, "ct", trait_groupby="cond", **KW)
    path = tmp_path / "a.h5ad"
    a.write_h5ad(path)
    r = CliRunner().invoke(app, [
        "celltype-modules", "-i", str(path), "--groupby", "ct", "-o", str(tmp_path / "o"),
        "-k", "10", "--max-shared", "3", "--max-metacells", "200", "--min-cells", "50",
        "--min-size", "10", "--group", "A", "--dme-by", "cond", "--trait", "age",
        "--trait", "sex", "--trait-groupby", "cond"])
    assert r.exit_code == 0, r.output
    d = pd.read_csv(tmp_path / "o" / "A" / "dmes.tsv", sep="\t")
    assert list(d.columns)[:3] == ["group", "module", "p_val"]
    assert d["module"].str.match(r"ME\d+$").all()
    t = pd.read_csv(tmp_path / "o" / "A" / "trait_cor.tsv", sep="\t")
    assert set(t["trait"]) == {"age", "sex"}
    r = CliRunner().invoke(app, [
        "celltype-modules", "-i", str(path), "--groupby", "ct", "-o", str(tmp_path / "x"),
        "--trait-groupby", "cond"])
    assert r.exit_code == 1


def test_consensus_modules_cli(tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    rng = np.random.default_rng(1)
    genes = [f"g{i}" for i in range(150)]
    paths = []
    for s in range(2):
        f = rng.normal(size=(3, 30))
        load = np.zeros((150, 3))
        for k in range(3):
            load[k * 40:(k + 1) * 40, k] = 1.0
        x = load @ f + 0.6 * rng.normal(size=(150, 30))
        p = tmp_path / f"set{s}.tsv"
        pd.DataFrame(x, index=genes, columns=[f"s{s}_{j}" for j in range(30)]).to_csv(p, sep="\t", lineterminator="\n")
        paths.append(p)
    r = CliRunner().invoke(app, ["consensus-modules", "-e", str(paths[0]), "-e", str(paths[1]),
                                 "--power", "6", "--min-size", "20", "--signed",
                                 "-o", str(tmp_path / "out")])
    assert r.exit_code == 0, r.output
    mods = pd.read_csv(tmp_path / "out" / "modules.tsv", sep="\t")
    assert mods["module"].nunique() == 3
    sets = pd.read_csv(tmp_path / "out" / "sets.tsv", sep="\t")
    assert list(sets["set"]) == ["set0", "set1"] and list(sets["samples"]) == [30, 30]
    me = pd.read_csv(tmp_path / "out" / "set_eigengenes.tsv", sep="\t")
    assert len(me) == 60
    r = CliRunner().invoke(app, ["consensus-modules", "-e", str(paths[0]), "-e", str(paths[1]),
                                 "--power", "6,7,8", "-o", str(tmp_path / "bad")])
    assert r.exit_code == 1
