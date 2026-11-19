"""Tests for gcen2.differential — Fisher's Z test differential co-expression."""

import numpy as np
import pandas as pd
import pytest

from gcen2.differential import _fisher_z_arrays as _fisher_z_test
from gcen2.differential import compare
from gcen2.exceptions import GCENInputError
from gcen2.models import Network


@pytest.fixture
def net_a():
    edges = pd.DataFrame({
        "gene_a": ["G1", "G1", "G2"],
        "gene_b": ["G2", "G3", "G3"],
        "correlation": [0.9, 0.7, 0.8],
    })
    return Network(edges=edges)


@pytest.fixture
def net_b():
    edges = pd.DataFrame({
        "gene_a": ["G1", "G2", "G3"],
        "gene_b": ["G2", "G4", "G4"],
        "correlation": [0.3, 0.85, 0.6],
    })
    return Network(edges=edges)


class TestFisherZTest:
    def test_identical_correlations(self):
        z, p = _fisher_z_test(0.8, 0.8, n_a=30, n_b=30)
        assert abs(z) < 1e-9
        assert abs(p - 1.0) < 0.01

    def test_very_different_correlations(self):
        z, p = _fisher_z_test(0.9, 0.1, n_a=100, n_b=100)
        assert abs(z) > 5
        assert p < 0.001

    def test_symmetric(self):
        z1, p1 = _fisher_z_test(0.8, 0.5, 50, 50)
        z2, p2 = _fisher_z_test(0.5, 0.8, 50, 50)
        assert abs(z1 + z2) < 1e-9
        assert abs(p1 - p2) < 1e-9


class TestCompare:
    def test_returns_diff_result(self, net_a, net_b):
        result = compare(net_a, net_b, n_a=30, n_b=30)
        assert not result.results.empty
        assert "category" in result.results.columns

    def test_categories_present(self, net_a, net_b):
        result = compare(net_a, net_b, n_a=30, n_b=30)
        categories = set(result.results["category"])
        valid = {"A_only", "B_only", "shared_changed", "shared_stable"}
        assert categories.issubset(valid)

    def test_shared_edge_detected(self, net_a, net_b):
        result = compare(net_a, net_b, n_a=30, n_b=30)
        shared = result.results[
            result.results["gene_a"].isin(["G1"]) & result.results["gene_b"].isin(["G2"])
        ]
        assert not shared.empty

    def test_invalid_method_raises(self, net_a, net_b):
        with pytest.raises(GCENInputError):
            compare(net_a, net_b, method="chi2", n_a=30, n_b=30)

    def test_too_few_samples_raises(self, net_a, net_b):
        with pytest.raises(GCENInputError):
            compare(net_a, net_b, n_a=3, n_b=3)

    def test_default_samples_from_metadata(self, net_a, net_b):
        """Without n_a / n_b: the per-gene effective samples of a metacell
        network, else p_value_samples, else n_samples."""
        genes = {"G1": 40.0, "G2": 12.5, "G3": 9.0, "G4": 50.0}
        a = Network(edges=net_a.edges, metadata={"n_samples": 1000, "p_value_samples": 20,
                                                 "gene_effective_samples": pd.Series(genes)})
        b = Network(edges=net_b.edges, metadata={"n_samples": 1000, "p_value_samples": 25})
        got = compare(a, b, pval=0.3).results
        want = compare(net_a, net_b, n_a=pd.Series(genes), n_b=25, pval=0.3).results
        pd.testing.assert_frame_equal(got, want)
        c = Network(edges=net_b.edges, metadata={"n_samples": 40})
        pd.testing.assert_frame_equal(compare(a, c).results,
                                      compare(net_a, net_b, n_a=pd.Series(genes), n_b=40).results)

    def test_per_gene_samples(self, net_a, net_b):
        """A mapping gene -> samples: every shared edge on the smaller of its
        genes' values; fewer than 4: no test; a missing gene is an error."""
        n_a = pd.Series({"G1": 40.0, "G2": 12.5, "G3": 3.0, "G4": 50.0})
        n_b = {"G1": 20, "G2": 100, "G3": 60, "G4": 8}
        res = compare(net_a, net_b, n_a=n_a, n_b=n_b, pval=0.3).results
        row = res[(res.gene_a == "G1") & (res.gene_b == "G2")].iloc[0]
        z, p = _fisher_z_test(0.9, 0.3, n_a=12.5, n_b=20)
        assert row.z_score == pytest.approx(z, rel=1e-12)
        assert row.p_value == pytest.approx(p, rel=1e-12)
        assert row.category == "shared_changed"
        scalar = compare(net_a, net_b, n_a=12.5, n_b=20, pval=0.3).results
        assert scalar.loc[row.name, "p_value"] == pytest.approx(row.p_value, rel=1e-12)
        res2 = compare(net_a, net_b, n_a=pd.Series({"G1": 40, "G2": 3, "G3": 9, "G4": 9}),
                       n_b=n_b).results
        row2 = res2[(res2.gene_a == "G1") & (res2.gene_b == "G2")].iloc[0]
        assert np.isnan(row2.p_value) and row2.category == "shared_stable"
        with pytest.raises(GCENInputError, match="G2"):
            compare(net_a, net_b, n_a={"G1": 40, "G3": 9}, n_b=n_b)

    def test_per_element_fisher_z(self):
        r_a, r_b = np.array([0.8, 0.2, 0.5]), np.array([0.1, 0.6, 0.5])
        n_a, n_b = np.array([30.0, 3.0, 10.0]), np.array([20.0, 40.0, 4.0])
        z, p = _fisher_z_test(r_a, r_b, n_a, n_b)
        for i in (0, 2):
            zi, pi = _fisher_z_test(r_a[i], r_b[i], n_a[i], n_b[i])
            assert z[i] == pytest.approx(zi, rel=1e-12) and p[i] == pytest.approx(pi, rel=1e-12)
        assert np.isnan(z[1]) and np.isnan(p[1])


def _reference(net_a, net_b, pval, n_a, n_b):
    """The original row-by-row implementation, rows sorted by gene pair."""
    import numpy as np

    def to_dict(net):
        d = {}
        for _, row in net.edges.iterrows():
            d[tuple(sorted([row["gene_a"], row["gene_b"]]))] = float(row["correlation"])
        return d

    ea, eb = to_dict(net_a), to_dict(net_b)
    rows = []
    for pair in sorted(set(ea) | set(eb)):
        ca, cb = ea.get(pair), eb.get(pair)
        if ca is None:
            z, p, cat = np.nan, np.nan, "B_only"
        elif cb is None:
            z, p, cat = np.nan, np.nan, "A_only"
        else:
            z, p = _fisher_z_test(ca, cb, n_a, n_b)
            cat = "shared_changed" if p <= pval else "shared_stable"
        rows.append((pair[0], pair[1], ca, cb, z, p, cat))
    df = pd.DataFrame(rows, columns=["gene_a", "gene_b", "cor_a", "cor_b", "z_score",
                                     "p_value", "category"])
    df.insert(6, "fdr", _bh_reference(df["p_value"].to_numpy()))
    return df


def _bh_reference(p):
    """Benjamini-Hochberg from its definition: q_i = min over p_j >= p_i of
    p_j m / rank_j (finite p only), capped at 1."""
    import numpy as np

    q = np.full(len(p), np.nan)
    ok = [i for i in range(len(p)) if np.isfinite(p[i])]
    m = len(ok)
    ranked = sorted(ok, key=lambda i: p[i])
    for i in ok:
        q[i] = min(1.0, min(p[j] * m / (r + 1) for r, j in enumerate(ranked) if p[j] >= p[i]))
    return q


class TestVectorisedCompare:
    def _net(self, rng, n_genes, n_edges, nan=False):
        import numpy as np

        g = [f"G{i:02d}" for i in range(n_genes)]
        pairs = rng.integers(0, n_genes, size=(n_edges, 2))
        cor = rng.uniform(-1, 1, n_edges).round(4)
        cor[:3] = [1.0, -1.0, 0.99995]  # clipped values
        if nan:
            cor[5] = np.nan
        return Network(edges=pd.DataFrame({
            "gene_a": [g[a] for a, _ in pairs], "gene_b": [g[b] for _, b in pairs],
            "correlation": cor}))

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_same_as_row_by_row(self, seed):
        import numpy as np

        rng = np.random.default_rng(seed)
        a = self._net(rng, 30, 200, nan=seed == 1)  # repeats/reversed pairs occur
        b = self._net(rng, 30, 150)
        got = compare(a, b, pval=0.05, n_a=12, n_b=20).results
        ref = _reference(a, b, 0.05, 12, 20)
        pd.testing.assert_frame_equal(got, ref, check_exact=True)
        assert set(got["category"]) == {"A_only", "B_only", "shared_changed", "shared_stable"}

    def test_bh_as_r_p_adjust(self):
        import numpy as np

        from gcen2.differential import _bh

        # p.adjust(c(0.01, 0.04, 0.03, 0.2, NA, 0.03), "BH") in R
        got = _bh(np.array([0.01, 0.04, 0.03, 0.2, np.nan, 0.03]))
        np.testing.assert_allclose(got, [0.05, 0.05, 0.05, 0.2, np.nan, 0.05])
        np.testing.assert_allclose(_bh(np.array([0.01, 0.04, 0.03, 0.2])),
                                   [0.04, 0.16 / 3, 0.16 / 3, 0.2])

    def test_fisher_z_threads_same_values(self):
        import numpy as np

        rng = np.random.default_rng(6)
        a, b = rng.uniform(-1, 1, (2, (1 << 20) * 2 + 123))
        one = _fisher_z_test(a, b, 30, 40, "spearman")
        four = _fisher_z_test(a, b, 30, 40, "spearman", threads=4)
        for x, y in zip(one, four):
            assert np.array_equal(x, y)

    def test_bh_extension_same_as_numpy(self, monkeypatch):
        import numpy as np

        import gcen2.core
        from gcen2.differential import _bh

        if gcen2.core.extension() is None:
            pytest.skip("compiled extension not available")
        rng = np.random.default_rng(5)
        p = rng.uniform(0, 1, 300_000) ** 3
        p[::7] = p[::7].round(3)  # ties
        p[::1001] = np.nan
        fast = [_bh(p, threads) for threads in (1, 4)]
        monkeypatch.setattr(gcen2.core, "_gcen_core", None)
        ref = _bh(p)
        for q in fast:
            assert np.array_equal(q, ref, equal_nan=True)
        np.testing.assert_array_equal(_bh(np.array([])), np.array([]))

    def test_fdr_category_and_spearman_variance(self):
        import numpy as np

        rng = np.random.default_rng(3)
        a, b = self._net(rng, 30, 200), self._net(rng, 30, 150)
        p_res = compare(a, b, pval=0.05, n_a=12, n_b=20).results
        q_res = compare(a, b, pval=0.05, n_a=12, n_b=20, fdr=True).results
        shared = p_res["z_score"].notna()
        assert (q_res.loc[shared, "category"] == np.where(
            q_res.loc[shared, "fdr"] <= 0.05, "shared_changed", "shared_stable")).all()
        assert (q_res["fdr"][shared] >= p_res["p_value"][shared]).all()
        s_res = compare(a, b, pval=0.05, n_a=12, n_b=20, correlation="spearman").results
        np.testing.assert_allclose(s_res["z_score"][shared],
                                   p_res["z_score"][shared] / np.sqrt(1.06))
        with pytest.raises(GCENInputError, match="correlation"):
            compare(a, b, n_a=12, n_b=20, correlation="kendall")

    def test_empty_network(self, net_a):
        empty = Network(edges=pd.DataFrame({"gene_a": [], "gene_b": [], "correlation": []}))
        got = compare(net_a, empty, n_a=10, n_b=10).results
        assert list(got["category"]) == ["A_only"] * 3
        assert len(compare(empty, empty, n_a=10, n_b=10).results) == 0


class TestCompareFiles:
    """compare_files (edges joined and written by the C++ extension) writes
    the same file as compare(read_network(a), read_network(b)).results.to_csv."""

    _HEADERS = {
        "v1": "#Node1\tNode2\tCorrelation\tPvalue\tFDR\n",
        "gcen2": "gene_a\tgene_b\tcorrelation\tp_value\n",
        "none": "",
    }

    def _write(self, path, seed, header, n_rows=300):
        rng = np.random.default_rng(seed)
        genes = [f"g{i:02d}" for i in range(25)] + ["G10", "a b", "101", 'q"t']
        lines = [self._HEADERS[header]]
        for k in range(n_rows):
            a, b = rng.choice(len(genes), 2)
            ga, gb = genes[a], genes[b]
            ga = f'"{ga}"' if ga == 'q"t' else ga  # quoted in the file
            gb = '"q""t"' if gb == 'q"t' else gb
            if k % 7 == 0:
                ga = f'"{ga}"' if not ga.startswith('"') else ga
            if k % 31 == 0:
                gb = ga  # self-loop
            cor = rng.uniform(-1, 1) * 10.0 ** -int(rng.integers(0, 3))
            c = "NA" if k % 13 == 0 else ("1.0" if k % 17 == 0 else repr(cor))
            extra = {"v1": "" if k % 5 == 0 else "\t1e-3", "gcen2": "", "none": "\t1e-3"}[header]
            lines.append(f"{ga}\t{gb}\t{c}\t2.5E-07{extra}\n")
            if k % 9 == 0:  # repeated / reversed edge: the last row wins
                lines.append(f"{gb}\t{ga}\t{repr(cor / 3)}\t1e-2{extra}\n")
        lines.append("\n")
        path.write_text("".join(lines), encoding="utf-8")
        return path

    def _both(self, tmp_path, monkeypatch, a, b, **kw):
        import gcen2.core
        from gcen2.core import _gcen_core
        from gcen2.differential import compare_files

        if _gcen_core is None:
            pytest.skip("compiled extension not available")
        got = compare_files(a, b, tmp_path / "fast.tsv", **kw)
        monkeypatch.setattr(gcen2.core, "_gcen_core", None)
        ref = compare_files(a, b, tmp_path / "ref.tsv", **kw)
        monkeypatch.undo()
        from tests.golden.golden import assert_same_file

        assert_same_file(tmp_path / "fast.tsv", tmp_path / "ref.tsv")
        assert got == ref
        return got

    @pytest.mark.parametrize("ha,hb", [("v1", "gcen2"), ("none", "v1"), ("gcen2", "none")])
    @pytest.mark.parametrize("pval,n_a,n_b", [(0.05, 30, 30), (0.5, 8, 100)])
    def test_same_as_pandas(self, tmp_path, monkeypatch, ha, hb, pval, n_a, n_b):
        a = self._write(tmp_path / "a.tsv", 1, ha)
        b = self._write(tmp_path / "b.tsv", 2, hb)
        got = self._both(tmp_path, monkeypatch, a, b, pval=pval, n_a=n_a, n_b=n_b)
        assert set(got) == {"A_only", "B_only", "shared_changed", "shared_stable"}

    @pytest.mark.parametrize("fdr,correlation", [(True, "pearson"), (True, "spearman"),
                                                 (False, "spearman")])
    def test_fdr_and_spearman_same_as_pandas(self, tmp_path, monkeypatch, fdr, correlation):
        a = self._write(tmp_path / "a.tsv", 4, "gcen2")
        b = self._write(tmp_path / "b.tsv", 5, "v1")
        self._both(tmp_path, monkeypatch, a, b, pval=0.3, n_a=10, n_b=40, fdr=fdr,
                   correlation=correlation)

    @pytest.mark.parametrize("correlation", ["pearson", "spearman"])
    def test_per_gene_samples_same_as_pandas(self, tmp_path, monkeypatch, correlation):
        a = self._write(tmp_path / "a.tsv", 6, "gcen2")
        b = self._write(tmp_path / "b.tsv", 7, "v1")
        from gcen2.io import read_network

        genes = sorted({str(g) for f in (a, b) for c in ("gene_a", "gene_b")
                        for g in read_network(f).edges[c]})  # names as parsed
        rng = np.random.default_rng(8)
        n_a = pd.Series(rng.uniform(2, 60, len(genes)), index=genes)
        n_b = dict(zip(genes, rng.integers(3, 80, len(genes)).tolist()))
        got = self._both(tmp_path, monkeypatch, a, b, pval=0.2, n_a=n_a, n_b=n_b,
                         correlation=correlation)
        assert got.get("shared_changed", 0) > 0
        res = pd.read_csv(tmp_path / "fast.tsv", sep="\t", keep_default_na=False,
                          na_values=[""])
        res = res[res["category"].str.startswith("shared")]
        res = res[np.isfinite(res["cor_a"]) & np.isfinite(res["cor_b"])]
        na = np.minimum(n_a[res.gene_a.astype(str)].to_numpy(),
                        n_a[res.gene_b.astype(str)].to_numpy())
        nb = np.minimum([n_b[g] for g in res.gene_a.astype(str)],
                        [n_b[g] for g in res.gene_b.astype(str)])
        _, p = _fisher_z_test(res.cor_a.to_numpy(), res.cor_b.to_numpy(), na, nb, correlation)
        np.testing.assert_allclose(res.p_value, p, rtol=1e-6)  # file: 6-7 digits

    def test_same_file_twice_and_empty(self, tmp_path, monkeypatch):
        a = self._write(tmp_path / "a.tsv", 3, "gcen2")
        empty = tmp_path / "empty.tsv"
        empty.write_text(self._HEADERS["gcen2"], encoding="utf-8")
        got = self._both(tmp_path, monkeypatch, a, a)
        assert set(got) <= {"shared_stable"}  # NaN p (NA correlation) is stable too
        assert self._both(tmp_path, monkeypatch, empty, empty) == {}
        assert set(self._both(tmp_path, monkeypatch, a, empty)) == {"A_only"}

    @staticmethod
    def _write_large(path, seed, n_rows=60_000, crlf=False):
        rng = np.random.default_rng(seed)
        genes = np.array([f"g{i:03d}" for i in range(400)] + ['"q""t"', "NA_x"])
        a = genes[rng.integers(0, len(genes), n_rows)]
        b = genes[rng.integers(0, len(genes), n_rows)]
        cor = rng.uniform(-1, 1, n_rows)
        nl = "\r\n" if crlf else "\n"
        lines = [f"{x}\t{y}\t{'NA' if k % 101 == 0 else repr(float(c))}\t1e-3"
                 for k, (x, y, c) in enumerate(zip(a, b, cor))]
        path.write_text("gene_a\tgene_b\tcorrelation\tp_value" + nl + nl.join(lines), encoding="utf-8")  # no final newline
        return path

    def test_threads_and_blocks_give_the_same_result(self, tmp_path, monkeypatch):
        from gcen2._matrix_io import correlation_specs
        from gcen2.core import _gcen_core
        from gcen2.differential import compare_files
        from tests.golden.golden import assert_same_file

        if _gcen_core is None:
            pytest.skip("compiled extension not available")
        a = self._write_large(tmp_path / "a.tsv", 1)
        b = self._write_large(tmp_path / "b.tsv", 2, crlf=True)
        specs = correlation_specs([a, b])
        ref = _gcen_core.network.diff_edge_files(specs[0], specs[1], 1)
        assert ref.size > 1 << 15  # the writer formats several blocks
        for threads, block in ((4, 1 << 26), (4, 4096), (3, 1000)):
            got = _gcen_core.network.diff_edge_files(specs[0], specs[1], threads, block)
            for col in ("cor_a", "cor_b", "flags"):
                np.testing.assert_array_equal(getattr(got, col)(), getattr(ref, col)())
        compare_files(a, b, tmp_path / "t1.tsv", threads=1, fdr=True)
        compare_files(a, b, tmp_path / "t4.tsv", threads=4, fdr=True)
        assert_same_file(tmp_path / "t4.tsv", tmp_path / "t1.tsv")
        assert (tmp_path / "t4.tsv").read_bytes() == (tmp_path / "t1.tsv").read_bytes()
        self._both(tmp_path, monkeypatch, a, b, threads=4)  # same as pandas

    def test_error_line_number_with_threads(self, tmp_path):
        from gcen2._matrix_io import correlation_specs
        from gcen2.core import _gcen_core

        if _gcen_core is None:
            pytest.skip("compiled extension not available")
        a = self._write_large(tmp_path / "a.tsv", 3, n_rows=5000)
        lines = a.read_text(encoding="utf-8").split("\n")
        lines[4321] = "g001\tNA\t0.5\t1e-3"  # line 4322 of the file
        a.write_text("\n".join(lines), encoding="utf-8")
        specs = correlation_specs([a, a])
        for threads, block in ((1, 1 << 26), (4, 2048)):
            with pytest.raises(RuntimeError, match="line 4322: missing gene name"):
                _gcen_core.network.diff_edge_files(specs[0], specs[1], threads, block)

    def test_missing_gene_name_raises(self, tmp_path):
        from gcen2.core import _gcen_core
        from gcen2.differential import compare_files
        from gcen2.exceptions import GCENIOError

        if _gcen_core is None:
            pytest.skip("compiled extension not available")
        a = tmp_path / "a.tsv"
        a.write_text("gene_a\tgene_b\tcorrelation\nA\tNA\t0.5\n", encoding="utf-8")
        with pytest.raises(GCENIOError, match="missing gene name"):
            compare_files(a, a, tmp_path / "out.tsv")

    def test_too_few_samples_raises(self, tmp_path):
        from gcen2.differential import compare_files

        with pytest.raises(GCENInputError):
            compare_files(tmp_path / "x", tmp_path / "y", tmp_path / "o", n_a=3)
