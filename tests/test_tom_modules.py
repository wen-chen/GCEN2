"""WGCNA-style TOM modules (gcen2.tom_modules, src/hclust.cpp): each step
against WGCNA 1.74 / dynamicTreeCut / fastcluster output on
tests/data/wgcna_tom/modexpr.tsv (make_modules_reference.R)."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from gcen2.core import _HAVE_CORE
from gcen2.exceptions import GCENInputError

pytestmark = pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")

DATA = Path(__file__).parent / "data" / "wgcna_tom"
KINDS = ["signed", "unsigned"]


def _core():
    from gcen2.core import _gcen_core
    return _gcen_core


@pytest.fixture(scope="module")
def expr():
    return pd.read_csv(DATA / "modexpr.tsv", sep="\t", index_col=0)


@pytest.fixture(scope="module", params=KINDS)
def ref(request):
    kind = request.param
    tree = pd.read_csv(DATA / f"hclust_{kind}.tsv", sep="\t")
    return {
        "kind": kind,
        "diss": np.loadtxt(DATA / f"diss_{kind}.tsv.gz"),
        "merge": tree[["a", "b"]].to_numpy(),
        "height": tree["height"].to_numpy(),
        "cutree": pd.read_csv(DATA / f"cutree_{kind}.tsv", sep="\t")["label"].to_numpy(),
        "modules": pd.read_csv(DATA / f"modules_{kind}.tsv", sep="\t"),
    }


def _canonical(labels):
    """Partition as labels numbered by first appearance (0 stays 0)."""
    out, seen = np.zeros(len(labels), dtype=int), {}
    for i, v in enumerate(labels):
        if v != 0:
            out[i] = seen.setdefault(v, len(seen) + 1)
    return out


class TestSteps:
    def test_dissimilarity(self, expr, ref):
        diss, info = _core().module.tom_dissimilarity(
            np.ascontiguousarray(expr.to_numpy()), "pearson", ref["kind"] == "signed", 6.0)
        assert diss.dtype == np.float64 and info["power"] == 6.0
        np.testing.assert_allclose(diss, ref["diss"], atol=1e-12)

    def test_hclust_like_fastcluster(self, ref):
        merge, height = _core().module.hclust_average(ref["diss"])
        assert np.array_equal(merge, ref["merge"])
        np.testing.assert_allclose(height, ref["height"], rtol=1e-14)
        m32, _ = _core().module.hclust_average(ref["diss"].astype(np.float32))
        assert np.array_equal(m32, ref["merge"])

    def test_hclust_inplace_restores(self, ref):
        d = ref["diss"].copy()
        merge, height = _core().module.hclust_average_inplace(d)
        assert np.array_equal(merge, ref["merge"])
        assert np.array_equal(d, ref["diss"])

    def test_cutree_hybrid(self, ref):
        from gcen2.tom_modules import cutree_hybrid

        lab = cutree_hybrid(ref["merge"], ref["height"], ref["diss"], min_cluster_size=10,
                            deep_split=2, cut_height=0.995)
        assert np.array_equal(lab, ref["cutree"])

    def test_blockwise_modules(self, expr, ref):
        from gcen2.tom_modules import tom_modules

        res = tom_modules(expr, method="pearson", signed=ref["kind"] == "signed", power=6,
                          min_size=10, merge_cut=0.25)
        assert np.array_equal(res.labels, ref["modules"]["merged"].to_numpy())
        assert np.array_equal(_canonical(res.unmerged),
                              _canonical(ref["modules"]["unmerged"].to_numpy()))
        assert res.genes == list(expr.index)


class TestDownstream:
    """Eigengenes, kME and trait correlations against WGCNA
    (make_downstream_reference.R) for the merged modules."""

    @pytest.fixture(scope="class", params=KINDS)
    def res(self, request, expr):
        from gcen2.tom_modules import tom_modules

        kind = request.param
        return kind, tom_modules(expr, signed=kind == "signed", power=6, min_size=10,
                                 merge_cut=0.25)

    def test_eigengenes(self, res):
        kind, r = res
        want = pd.read_csv(DATA / f"me_{kind}.tsv", sep="\t", index_col=0)
        got = r.eigengene_frame()
        assert list(got.columns) == list(want.columns) and list(got.index) == list(want.index)
        np.testing.assert_allclose(got.to_numpy(), want.to_numpy(), atol=1e-12)
        var = pd.read_csv(DATA / f"varexp_{kind}.tsv", sep="\t")["var"].to_numpy()
        np.testing.assert_allclose(r.var_explained, var, atol=1e-12)

    def test_kme(self, res):
        kind, r = res
        want = pd.read_csv(DATA / f"kme_{kind}.tsv", sep="\t", index_col=0).to_numpy()
        kme, p = r.kme()
        k = kme.shape[1]
        np.testing.assert_allclose(kme, want[:, :k], atol=1e-12)
        np.testing.assert_allclose(p, want[:, k:], rtol=1e-9, atol=1e-300)

    def test_trait_cor(self, res):
        kind, r = res
        traits = pd.read_csv(DATA / "traits.tsv", sep="\t", index_col=0)
        want = pd.read_csv(DATA / f"trait_{kind}.tsv", sep="\t")
        got = r.trait_cor(traits)
        assert len(got) == len(want)
        m = got.merge(want, on=["module", "trait"], suffixes=("", "_r"))
        assert len(m) == len(want) and (m["n"] == m["n_r"]).all()
        np.testing.assert_allclose(m["cor"], m["cor_r"], atol=1e-12)
        np.testing.assert_allclose(m["p"], m["p_r"], rtol=1e-9)

    def test_trait_samples_matched_by_name(self, res):
        _, r = res
        traits = pd.read_csv(DATA / "traits.tsv", sep="\t", index_col=0)
        a = r.trait_cor(traits)
        b = r.trait_cor(traits.iloc[::-1].drop(index="s05"))
        assert (b["n"] <= a["n"]).all() and (b.loc[b.trait == "age", "n"] == 29).all()
        with pytest.raises(GCENInputError):
            r.trait_cor(traits.assign(group="x"))
        with pytest.raises(GCENInputError):
            r.trait_cor(traits.set_axis([f"z{i}" for i in range(len(traits))]))

    def test_left_out_genes_have_no_kme(self, expr):
        from gcen2.tom_modules import tom_modules

        x = expr.copy()
        x.iloc[3] = 1.0
        r = tom_modules(x, signed=True, power=6, min_size=10)
        kme, p = r.kme()
        assert np.isnan(kme[3]).all() and np.isfinite(np.delete(kme, 3, axis=0)).all()


class TestHdwgcna:
    """bicor power selection and the hdwgcna preset against WGCNA
    (make_bicor_reference.R, make_hdwgcna_reference.R)."""

    def test_bicor(self):
        x = pd.read_csv(DATA / "bicor_expr.tsv", sep="\t", index_col=0)
        b = pd.read_csv(DATA / "bicor.tsv", sep="\t", index_col=0).to_numpy()
        diss, _ = _core().module.tom_dissimilarity(np.ascontiguousarray(x.to_numpy()),
                                                   "bicor", True, 2.0)
        a = ((1 + b) / 2) ** 2
        np.fill_diagonal(a, 0)
        k = a.sum(axis=1)
        want = 1 - (a @ a + a) / (np.minimum.outer(k, k) + 1 - a)
        np.fill_diagonal(want, 0)
        np.testing.assert_allclose(diss, want, atol=1e-13)

    def test_bicor_soft_threshold(self):
        x = pd.read_csv(DATA / "bicor_expr.tsv", sep="\t", index_col=0)
        want = pd.read_csv(DATA / "sft_bicor_signed.tsv", sep="\t")
        fits = np.array(_core().correlation.soft_threshold(
            np.ascontiguousarray(x.to_numpy()), "bicor", True,
            [float(p) for p in want["Power"]], 2))
        np.testing.assert_allclose(fits[:, 1], want["SFT.R.sq"], atol=1e-6)
        np.testing.assert_allclose(fits[:, 2], want["slope"], rtol=1e-5)
        np.testing.assert_allclose(fits[:, 3], want["mean.k."], rtol=1e-6)

    def test_preset(self, expr):
        from gcen2.tom_modules import HDWGCNA, tom_modules

        want = pd.read_csv(DATA / "hdwgcna_modules.tsv", sep="\t")["module"].to_numpy()
        power = float((DATA / "hdwgcna_power.txt").read_text(encoding="utf-8"))
        res = tom_modules(expr, threads=3, **{**HDWGCNA, "min_size": 10})
        assert res.power == power and res.power_estimated
        assert "bicor" in res.power_rule
        assert np.array_equal(_canonical(res.labels), _canonical(want))

    def test_preset_cli(self, tmp_path, expr):
        from typer.testing import CliRunner

        from gcen2.cli import app

        src = tmp_path / "x.tsv"
        expr.to_csv(src, sep="\t", float_format="%.17g", lineterminator="\n")
        r = CliRunner().invoke(app, ["module", "-e", str(src), "-o", str(tmp_path / "m.tsv"),
                                     "--preset", "hdwgcna", "--min-size", "10"])
        assert r.exit_code == 0, r.output
        assert "TOM power 9 (bicor, lowest power above 3 with R^2 >= 0.8)" in r.output
        got = pd.read_csv(tmp_path / "m.tsv", sep="\t")
        want = pd.read_csv(DATA / "hdwgcna_modules.tsv", sep="\t")
        assert len(got) == int((want["module"] > 0).sum())


class TestBehaviour:
    def test_threads_do_not_matter(self, expr):
        from gcen2.tom_modules import tom_modules

        x = np.ascontiguousarray(expr.to_numpy())
        a, _ = _core().module.tom_dissimilarity(x, "spearman", False, 0.0, [], 1)
        b, _ = _core().module.tom_dissimilarity(x, "spearman", False, 0.0, [], 5)
        assert np.array_equal(a, b)
        r1 = tom_modules(expr, method="spearman", threads=1)
        r4 = tom_modules(expr, method="spearman", threads=4)
        assert np.array_equal(r1.labels, r4.labels) and r1.power == r4.power

    def test_dissimilarity_bands(self):
        """Many genes (several bands of the in-place 1 - TOM, partial last
        tiles) against a numpy TOM."""
        rng = np.random.default_rng(3)
        n, s = 4200, 12
        x = rng.normal(size=(n, s)) + np.repeat(rng.normal(size=(30, s)), n // 30, axis=0)
        diss, _ = _core().module.tom_dissimilarity(x, "pearson", True, 4.0, [], 3)
        a = ((1 + np.corrcoef(x)) / 2) ** 4
        np.fill_diagonal(a, 0)
        k = a.sum(axis=1)
        want = 1 - (a @ a + a) / (np.minimum.outer(k, k) + 1 - a)
        np.fill_diagonal(want, 0)
        np.testing.assert_allclose(diss, want, atol=1e-12)
        assert np.array_equal(diss, diss.T)

    def test_kernels_agree(self, tmp_path):
        """The AVX and baseline (SSE2 / plain C++) product kernels give the
        same bits (GCEN2_SIMD=baseline forces the latter), double and float."""
        import os
        import subprocess
        import sys

        code = (
            "import numpy as np, sys\n"
            "from gcen2.core import _gcen_core as c\n"
            "x = np.random.default_rng(1).normal(size=(700, 37))\n"
            "x[:350] += np.repeat(np.random.default_rng(2).normal(size=(7, 37)), 50, axis=0)\n"
            "d, _ = c.module.tom_dissimilarity(x, 'pearson', False, 5.0, [], 3)\n"
            "ij, *_, tom, info = c.correlation.network_edges(\n"
            "    x, 'pearson', 3, 0.0, 1.0, False, False, True, 5.0)\n"
            "np.save(sys.argv[1], np.concatenate([d.ravel(), np.asarray(tom, float)]))\n"
            "print(c.module.tom_kernel())\n")
        out = {}
        for mode in ("", "baseline"):
            env = dict(os.environ, GCEN2_SIMD=mode)
            f = tmp_path / f"t{mode}.npy"
            r = subprocess.run([sys.executable, "-c", code, str(f)], env=env,
                               capture_output=True, text=True, encoding="utf-8", check=True)
            out[r.stdout.strip()] = np.load(f)
        assert "baseline" in out
        if "avx" not in out:
            pytest.skip("no AVX kernel on this machine")
        assert np.array_equal(out["avx"], out["baseline"])

    def test_bad_genes_left_out(self, expr):
        from gcen2.tom_modules import tom_modules

        x = expr.copy()
        x.iloc[3] = 1.0
        x.iloc[8, 2] = np.nan
        res = tom_modules(x, signed=True, power=6, min_size=10, merge_cut=0.25)
        assert res.labels[3] == 0 and res.labels[8] == 0
        good = tom_modules(x.drop(index=x.index[[3, 8]]), signed=True, power=6, min_size=10,
                           merge_cut=0.25)
        assert np.array_equal(np.delete(res.labels, [3, 8]), good.labels)

    def test_modules_listing(self, expr):
        from gcen2.tom_modules import tom_modules

        res = tom_modules(expr, signed=True, power=6, min_size=10)
        mods = res.modules()
        assert [len(m) for m in mods] == sorted((len(m) for m in mods), reverse=True)
        assert sum(len(m) for m in mods) == int((res.labels > 0).sum())

    def test_small_trees(self):
        d = np.array([[0, 1, 4], [1, 0, 2], [4, 2, 0]], dtype=float)
        merge, height = _core().module.hclust_average(d)
        assert merge.tolist() == [[-1, -2], [-3, 1]]
        assert height.tolist() == [1.0, 3.0]
        merge, height = _core().module.hclust_average(np.array([[0, 0.5], [0.5, 0.0]]))
        assert merge.tolist() == [[-1, -2]] and height.tolist() == [0.5]

    @pytest.mark.parametrize("kw", [dict(method="pearson_nonzero"), dict(power=-1),
                                    dict(min_size=1), dict(merge_cut=2.0),
                                    dict(deep_split=7)])
    def test_rejected(self, expr, kw):
        from gcen2.tom_modules import tom_modules

        with pytest.raises(GCENInputError):
            tom_modules(expr, **kw)


class TestCli:
    def test_module_from_expression(self, tmp_path, expr):
        from typer.testing import CliRunner

        from gcen2.cli import app
        from gcen2.tom_modules import tom_modules

        src = tmp_path / "x.tsv"
        expr.to_csv(src, sep="\t", float_format="%.17g", lineterminator="\n")
        out = tmp_path / "m.tsv"
        r = CliRunner().invoke(app, ["module", "-e", str(src), "-o", str(out), "--signed",
                                     "--power", "6", "--min-size", "10", "--merge-cut", "0.25"])
        assert r.exit_code == 0, r.output
        assert "TOM power 6 (given)" in r.output
        got = pd.read_csv(out, sep="\t")
        want = tom_modules(expr, signed=True, power=6, min_size=10, merge_cut=0.25)
        assert [sorted(got.loc[got["module"] == i, "gene"]) for i in range(len(want.modules()))] \
            == want.modules()

    def test_downstream_files(self, tmp_path, expr):
        from typer.testing import CliRunner

        from gcen2.cli import app

        src = tmp_path / "x.tsv"
        expr.to_csv(src, sep="\t", float_format="%.17g", lineterminator="\n")
        o = {k: tmp_path / f"{k}.tsv" for k in ("m", "me", "kme", "tc")}
        r = CliRunner().invoke(app, [
            "module", "-e", str(src), "-o", str(o["m"]), "--signed", "--power", "6",
            "--min-size", "10", "--merge-cut", "0.25", "--me", str(o["me"]),
            "--kme", str(o["kme"]), "--traits", str(DATA / "traits.tsv"),
            "--trait-cor", str(o["tc"])])
        assert r.exit_code == 0, r.output
        mods = pd.read_csv(o["m"], sep="\t")
        k = mods["module"].max() + 1
        me = pd.read_csv(o["me"], sep="\t", index_col=0)
        assert list(me.columns) == ["ME_unassigned"] + [f"ME{i}" for i in range(k)]
        want = pd.read_csv(DATA / "me_signed.tsv", sep="\t", index_col=0)
        np.testing.assert_allclose(me.to_numpy(), want.to_numpy(), atol=1e-9)
        kme = pd.read_csv(o["kme"], sep="\t", index_col=0, keep_default_na=False)
        got = dict(zip(mods["gene"], mods["module"].astype(str)))
        assert all(kme.loc[g, "module"] == got.get(g, "NA") for g in kme.index)
        assert list(kme.columns[1:k + 2]) == ["kME_unassigned"] + [f"kME{i}" for i in range(k)]
        tc = pd.read_csv(o["tc"], sep="\t")
        assert list(tc.columns) == ["module", "trait", "cor", "p", "n"]
        assert set(tc["module"]) == set(me.columns)

    def test_h5ad_input(self, tmp_path, expr):
        ad = pytest.importorskip("anndata")
        from typer.testing import CliRunner

        from gcen2.cli import app

        a = ad.AnnData(expr.T.to_numpy())
        a.obs_names, a.var_names = list(expr.columns), list(expr.index)
        a.obs["g"] = ["a"] * 20 + ["b"] * 10
        a.write_h5ad(tmp_path / "x.h5ad")
        src = tmp_path / "x.tsv"
        expr.to_csv(src, sep="\t", float_format="%.17g", lineterminator="\n")
        args = ["--signed", "--power", "6", "--min-size", "10"]
        run = lambda *a: CliRunner().invoke(app, ["module", *a, *args])  # noqa: E731
        for inp, out in ((src, "t.tsv"), (tmp_path / "x.h5ad", "h.tsv")):
            assert run("-e", str(inp), "-o", str(tmp_path / out)).exit_code == 0
        assert (tmp_path / "t.tsv").read_text(encoding="utf-8") == (tmp_path / "h.tsv").read_text(encoding="utf-8")
        r = run("-e", str(tmp_path / "x.h5ad"), "-o", str(tmp_path / "g.tsv"),
                "--groupby", "g", "--group", "a", "--me", str(tmp_path / "me.tsv"))
        assert r.exit_code == 0, r.output
        me = pd.read_csv(tmp_path / "me.tsv", sep="\t", index_col=0)
        assert list(me.index) == list(expr.columns[:20])
        assert run("-e", str(src), "-o", str(tmp_path / "e.tsv"), "--hvg").exit_code == 1
        assert run("-i", str(src), "-o", str(tmp_path / "e.tsv"),
                   "--kme", str(tmp_path / "k.tsv")).exit_code == 1

    def test_needs_one_input(self, tmp_path):
        from typer.testing import CliRunner

        from gcen2.cli import app

        r = CliRunner().invoke(app, ["module", "-o", str(tmp_path / "m.tsv")])
        assert r.exit_code == 1


def test_blas_limit_covers_scipy():
    """scipy's own OpenBLAS is limited too, even when first imported inside
    the limit (it then ran on every CPU: eigengenes several times slower).
    Skipped where scipy's BLAS is not one threadpoolctl controls (Apple's
    Accelerate on macOS)."""
    import subprocess
    import sys

    probe = subprocess.run(
        [sys.executable, "-c", "import scipy.linalg\nfrom threadpoolctl import threadpool_info\n"
         "import re\nprint(any(i['user_api'] == 'blas' and "
         "re.search(r'[/\\\\]scipy([.]libs)?[/\\\\]', i['filepath']) is not None "
         "for i in threadpool_info()))"], capture_output=True, text=True, encoding="utf-8", check=True)
    if probe.stdout.strip() != "True":
        pytest.skip("scipy's BLAS is not controlled by threadpoolctl here")

    code = ("from threadpoolctl import threadpool_info\n"
            "from gcen2.tom_modules import _blas_threads\n"
            "with _blas_threads(2):\n"
            "    import scipy.linalg\n"
            "    print(sorted({i['num_threads'] for i in threadpool_info()"
            " if i['user_api'] == 'blas'}))\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8",
                         check=True).stdout
    assert out.strip() == "[2]"
