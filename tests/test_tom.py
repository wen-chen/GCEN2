"""TOM edge scores (src/tom.cpp): WGCNA 1.74 reference values
(tests/data/wgcna_tom, from make_reference.R), a dense numpy reference
across several tiles, and the same output from every build path."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.stats import rankdata

from gcen2.core import _HAVE_CORE
from gcen2.exceptions import GCENInputError

pytestmark = pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")

DATA = Path(__file__).parent / "data" / "wgcna_tom"


@pytest.fixture(scope="module")
def wgcna_expr():
    return pd.read_csv(DATA / "expr.tsv", sep="\t", index_col=0)


def _dense_tom(x, power, signed):
    """TOM of all gene pairs in float64 (genes with a constant row: adjacency 0)."""
    with np.errstate(invalid="ignore", divide="ignore"):
        r = np.corrcoef(x)
    r[~np.isfinite(r)] = np.nan
    base = (1 + r) / 2 if signed else np.abs(r)
    a = np.nan_to_num(base ** power, nan=0.0)
    np.fill_diagonal(a, 0.0)
    k = a.sum(1)
    return ((a @ a) + a) / (np.minimum.outer(k, k) + 1 - a)


class TestWgcnaReference:
    @pytest.mark.parametrize("kind", ["signed", "unsigned"])
    def test_pick_soft_threshold(self, wgcna_expr, kind):
        from gcen2.network import pick_soft_threshold

        got = pick_soft_threshold(wgcna_expr, method="pearson", signed=kind == "signed")
        ref = pd.read_csv(DATA / f"sft_{kind}.tsv", sep="\t")
        assert got["power"].tolist() == ref["power"].tolist()
        np.testing.assert_allclose(got["r2"], ref["r2"], atol=1e-6)
        np.testing.assert_allclose(got["slope"], ref["slope"], rtol=1e-5, atol=1e-6)
        np.testing.assert_allclose(got["mean_k"], ref["mean_k"], rtol=1e-6)
        assert got.loc[got["chosen"], "power"].tolist() == [ref["estimate"][0]]

    @pytest.mark.parametrize("kind", ["signed", "unsigned"])
    def test_tom_values(self, wgcna_expr, kind):
        from gcen2.network import build

        net = build(wgcna_expr, method="pearson", signed=kind == "signed", all_edges=True,
                    tom=True, power=6)
        got = net.edges.set_index(["gene_a", "gene_b"])["tom"]
        ref = pd.read_csv(DATA / f"tom_{kind}.tsv", sep="\t")
        g = got.loc[list(zip(ref["gene_a"], ref["gene_b"]))].to_numpy()
        np.testing.assert_allclose(g, ref["tom"], rtol=2e-6)
        assert net.metadata["tom"] == {"power": 6.0, "estimated": False, "fits": []}

    @pytest.mark.parametrize("kind", ["signed", "unsigned"])
    def test_auto_power(self, wgcna_expr, kind):
        from gcen2.network import build

        net = build(wgcna_expr, method="pearson", signed=kind == "signed", tom=True)
        ref = pd.read_csv(DATA / f"sft_{kind}.tsv", sep="\t")
        meta = net.metadata["tom"]
        assert meta["power"] == ref["estimate"][0] and meta["estimated"]
        assert [f["power"] for f in meta["fits"]] == ref["power"].tolist()


def _matrix(n=300, s=30, seed=3):
    """More genes than one tile (128 x 256): modules, a constant gene, a NaN."""
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(5, s))
    x = base[rng.integers(0, 5, n)] * rng.uniform(0.3, 1.5, (n, 1)) + rng.normal(size=(n, s))
    x[7] = 2.0
    x[11, 4] = np.nan
    return x


class TestDenseReference:
    @pytest.mark.parametrize("signed", [False, True])
    @pytest.mark.parametrize("method", ["pearson", "spearman"])
    def test_all_pairs(self, method, signed):
        from gcen2.core import _gcen_core

        x = _matrix()
        idx, r, p, q, t, info = _gcen_core.correlation.network_edges(
            x, method, 4, 0.0, 1.0, signed, False, True, 5.0)
        y = x.copy()
        if method == "spearman":  # ranks as the kernel computes them (NaN row: all NaN)
            y = rankdata(x, axis=1)
            y[np.isnan(x).any(axis=1)] = np.nan
        want = _dense_tom(y, 5.0, signed)[idx[:, 0], idx[:, 1]]
        np.testing.assert_allclose(t, want, rtol=1e-5, atol=1e-9)
        assert info["power"] == 5.0

    def test_invalid_genes_have_adjacency_zero(self):
        from gcen2.core import _gcen_core

        x = _matrix()
        idx, _, _, _, t, _ = _gcen_core.correlation.network_edges(
            x, "pearson", 2, 0.0, 1.0, False, False, True, 4.0)
        for g in (7, 11):  # TOM = l / (min k + 1) with a = 0
            sel = (idx[:, 0] == g) | (idx[:, 1] == g)
            assert np.all(t[sel] >= 0) and np.all(np.isfinite(t[sel]))


class TestPaths:
    OPTS = [
        dict(cor=0.3, pval=0.01, signed=True, tom_min=0.0),
        dict(cor=0.2, pval=0.05, signed=False, fdr=True, tom_min=0.02),
        dict(all_edges=True, tom_min=0.5),  # all_edges keeps every pair
    ]

    @pytest.mark.parametrize("opts", OPTS)
    def test_same_file_everywhere(self, tmp_path, opts):
        """build().to_edgelist, build_to_file (any block size and threads)
        and build_from_file (C++ reader) write the same bytes."""
        from gcen2._build_file import build_from_file
        from gcen2.network import build, build_to_file

        x = _matrix()
        genes = [f"g{i}" for i in range(len(x))]
        df = pd.DataFrame(x, index=genes, columns=[f"s{j}" for j in range(x.shape[1])])
        net = build(df, method="pearson", tom=True, power=6, threads=3, **opts)
        net.to_edgelist(tmp_path / "a.tsv")
        want = (tmp_path / "a.tsv").read_bytes()
        assert want.splitlines()[0].endswith(b"\ttom")
        for block, threads in [(1, 1), (997, 4), (10**7, 2)]:
            info = build_to_file(x, tmp_path / "b.tsv", method="pearson", gene_names=genes,
                                 tom=True, power=6, threads=threads,
                                 max_pairs_per_block=block, **opts)
            assert (tmp_path / "b.tsv").read_bytes() == want, (block, threads)
            assert info["edge_count"] == len(net.edges)
            assert info["tom"]["power"] == 6.0
        # the C++ file reader (no NaN: numpy path otherwise)
        df2 = df.fillna(0.0)
        df2.to_csv(tmp_path / "m.tsv", sep="\t", float_format="%.17g", lineterminator="\n")
        build(df2, method="pearson", tom=True, power=6, **opts).to_edgelist(tmp_path / "c.tsv")
        build_from_file(tmp_path / "m.tsv", tmp_path / "d.tsv", method="pearson", tom=True,
                        power=6, **opts)
        assert (tmp_path / "d.tsv").read_bytes() == (tmp_path / "c.tsv").read_bytes()

    def test_tom_min_filters(self):
        from gcen2.network import build

        x = _matrix()
        df = pd.DataFrame(x, index=[f"g{i}" for i in range(len(x))])
        full = build(df, method="pearson", cor=0.2, pval=0.05, tom=True, power=6)
        cut = build(df, method="pearson", cor=0.2, pval=0.05, tom=True, power=6, tom_min=0.05)
        keep = full.edges[full.edges["tom"] >= 0.05].reset_index(drop=True)
        pd.testing.assert_frame_equal(cut.edges, keep)
        assert 0 < len(cut.edges) < len(full.edges)

    def test_without_tom_unchanged(self):
        from gcen2.network import build

        df = pd.DataFrame(_matrix(), index=[f"g{i}" for i in range(300)])
        a = build(df, method="pearson", cor=0.3)
        b = build(df, method="pearson", cor=0.3, tom=True, power=6)
        pd.testing.assert_frame_equal(a.edges, b.edges.drop(columns="tom"))
        assert a.metadata == {}

    def test_read_back(self, tmp_path):
        from gcen2.io import read_network
        from gcen2.network import build

        df = pd.DataFrame(_matrix(), index=[f"g{i}" for i in range(300)])
        net = build(df, method="pearson", cor=0.3, tom=True, power=6)
        net.to_edgelist(tmp_path / "n.tsv")
        back = read_network(tmp_path / "n.tsv")
        assert list(back.edges.columns) == ["gene_a", "gene_b", "correlation", "p_value", "tom"]
        np.testing.assert_allclose(back.edges["tom"], net.edges["tom"], rtol=5e-7)


class TestFallbackAndErrors:
    def test_faq_power_when_no_fit(self):
        from gcen2.network import build

        rng = np.random.default_rng(0)
        df = pd.DataFrame(rng.normal(size=(60, 25)), index=[f"g{i}" for i in range(60)])
        with pytest.warns(UserWarning, match="FAQ power"):
            net = build(df, method="pearson", signed=True, tom=True)
        meta = net.metadata["tom"]
        assert not meta["estimated"] and meta["power"] == 16.0  # signed, 20-29 samples
        assert max(f["r2"] for f in meta["fits"]) <= 0.85

    @pytest.mark.parametrize("kw", [dict(method="pearson_nonzero"), dict(power=0),
                                    dict(power="x"), dict(power=-2)])
    def test_rejected(self, kw):
        from gcen2.network import build

        df = pd.DataFrame(_matrix(40), index=[f"g{i}" for i in range(40)])
        with pytest.raises(GCENInputError):
            build(df, **({"method": "pearson"} | kw), tom=True)
