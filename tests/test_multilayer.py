"""Tests for gcen2.multilayer (weighted fusion and multiplex RWR)."""

import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import gcen2
from gcen2.exceptions import GCENInputError
from gcen2.models import Network
from gcen2.multilayer import as_layer, fuse, read_layer, rwr_multiplex

DATA = Path(__file__).parent / "data"


def _random_net(m=150, n_genes=60, seed=1):
    rng = np.random.default_rng(seed)
    genes = np.array([f"G{i}" for i in range(n_genes)])
    a, b = rng.choice(genes, m), rng.choice(genes, m)
    k = a != b
    return Network(edges=pd.DataFrame({
        "gene_a": a[k], "gene_b": b[k],
        "correlation": rng.uniform(-1, 1, k.sum()), "p_value": 0.0}))


def _edges(rows, cols=("gene_a", "gene_b", "weight")):
    return pd.DataFrame(rows, columns=list(cols))


class TestLayers:
    def test_network_weight_is_abs_correlation(self):
        net = Network(edges=_edges([("A", "B", -0.8, 0.0)],
                                   ("gene_a", "gene_b", "correlation", "p_value")))
        assert as_layer(net)["weight"].tolist() == [0.8]

    def test_unweighted_table(self):
        assert as_layer(_edges([("A", "B")], ("gene_a", "gene_b")))["weight"].tolist() == [1.0]

    def test_read_gcen_network(self, tmp_path):
        net = _random_net()
        net.to_edgelist(tmp_path / "n.network")
        layer = read_layer(tmp_path / "n.network")
        np.testing.assert_allclose(layer["weight"], net.edges["correlation"].abs(), atol=1e-6)
        assert layer["gene_a"].tolist() == net.edges["gene_a"].tolist()

    def test_read_string_style_file(self, tmp_path):
        p = tmp_path / "ppi.txt"
        p.write_text("protein1 protein2 combined_score\nA B 900\nB C 400\n", encoding="utf-8")
        layer = read_layer(p, sep=" ")
        assert layer.values.tolist() == [["A", "B", 900.0], ["B", "C", 400.0]]

    def test_read_two_column_file(self, tmp_path):
        p = tmp_path / "reg.tsv"
        p.write_text("TF1\tG1\nTF1\tG2\n", encoding="utf-8")
        assert read_layer(p)["weight"].tolist() == [1.0, 1.0]

    def test_read_weight_column_by_name_and_position(self, tmp_path):
        p = tmp_path / "x.tsv"
        p.write_text("a\tb\tscore\tother\nA\tB\t0.5\t7\n", encoding="utf-8")
        assert read_layer(p, weight_column="other")["weight"].tolist() == [7.0]
        assert read_layer(p, weight_column=2)["weight"].tolist() == [0.5]
        with pytest.raises(GCENInputError, match="no column"):
            read_layer(p, weight_column="nope")

    def test_numeric_gene_ids_kept_as_text(self, tmp_path):
        p = tmp_path / "e.tsv"
        p.write_text("0007\t12\t1\n", encoding="utf-8")
        assert read_layer(p)[["gene_a", "gene_b"]].values.tolist() == [["0007", "12"]]


class TestFuse:
    def test_weighted_mean(self):
        a = _edges([("A", "B", 1.0), ("B", "C", 0.5)])
        b = _edges([("B", "A", 2.0), ("C", "D", 1.0)])
        e = fuse({"coexp": a, "ppi": b}, weights={"coexp": 3, "ppi": 1},
                 normalize="none").edges
        row = e.set_index(["gene_a", "gene_b"])
        assert row.loc[("A", "B"), "score"] == pytest.approx((3 * 1.0 + 1 * 2.0) / 4)
        assert row.loc[("A", "B"), "n_layers"] == 2
        assert row.loc[("A", "B"), "layers"] == "coexp,ppi"
        assert row.loc[("B", "C"), "score"] == pytest.approx(3 * 0.5 / 4)
        assert row.loc[("C", "D"), "w_coexp"] == 0.0
        assert list(e.columns) == ["gene_a", "gene_b", "score", "n_layers", "layers",
                                   "w_coexp", "w_ppi"]
        assert e["score"].is_monotonic_decreasing

    def test_normalize_max_and_rank(self):
        a = _edges([("A", "B", 900.0), ("B", "C", 300.0)])
        e = fuse([a], normalize="max").edges.set_index(["gene_a", "gene_b"])
        assert e.loc[("B", "C"), "score"] == pytest.approx(1 / 3)
        e = fuse([a], normalize="rank").edges.set_index(["gene_a", "gene_b"])
        assert e.loc[("B", "C"), "score"] == pytest.approx(0.5)
        assert e.loc[("A", "B"), "score"] == pytest.approx(1.0)

    def test_duplicates_keep_max_and_self_loops_dropped(self):
        a = _edges([("A", "B", 0.2), ("B", "A", 0.7), ("A", "A", 1.0)])
        e = fuse([a], normalize="none").edges
        assert e[["gene_a", "gene_b", "score"]].values.tolist() == [["A", "B", 0.7]]

    def test_min_layers_and_min_score(self):
        a = _edges([("A", "B", 1.0), ("B", "C", 1.0)])
        b = _edges([("A", "B", 1.0), ("C", "D", 0.2)])
        e = fuse({"a": a, "b": b}, min_layers=2).edges
        assert e[["gene_a", "gene_b"]].values.tolist() == [["A", "B"]]
        e = fuse({"a": a, "b": b}, normalize="none", min_score=0.5).edges
        assert set(map(tuple, e[["gene_a", "gene_b"]].values)) == {("A", "B"), ("B", "C")}

    def test_result_is_a_network(self):
        net = fuse({"a": _random_net(seed=1), "b": _random_net(seed=2)})
        assert net.method == "fused"
        assert net.metadata["layers"] == ["a", "b"]
        assert net.stats().node_count == net.node_count

    @pytest.mark.parametrize("kwargs,match", [
        ({"normalize": "zscore"}, "normalize"),
        ({"weights": {"x": 1}}, "unknown"),
        ({"weights": {"a": -1}}, ">= 0"),
        ({"min_layers": 3}, "min_layers"),
    ])
    def test_bad_arguments(self, kwargs, match):
        with pytest.raises(GCENInputError, match=match):
            fuse({"a": _edges([("A", "B", 1.0)]), "b": _edges([("A", "B", 1.0)])}, **kwargs)

    def test_no_layers(self):
        with pytest.raises(GCENInputError, match="At least one"):
            fuse({})


class TestRwrMultiplex:
    SEEDS = ["G1", "G2", "G3"]

    def test_one_layer_equals_rwr_predict(self):
        net = _random_net()
        ref = gcen2.rwr.predict(net, self.SEEDS, prob_cutoff=0).ranked_genes
        res = rwr_multiplex({"a": net}, self.SEEDS, aggregate="sum").ranked_genes
        m = ref.merge(res, on="gene")
        assert len(m) == len(ref) == len(res)
        np.testing.assert_allclose(m["probability_x"], m["probability_y"], atol=1e-14)

    def test_one_layer_weighted_equals_rwr_predict(self):
        net = _random_net()
        ref = gcen2.rwr.predict(net, self.SEEDS, prob_cutoff=0,
                                weighted_network=True).ranked_genes
        res = rwr_multiplex([net], self.SEEDS, weighted=True, aggregate="sum").ranked_genes
        m = ref.merge(res, on="gene")
        np.testing.assert_allclose(m["probability_x"], m["probability_y"], atol=1e-14)

    @pytest.mark.parametrize("gamma,delta", [(0.5, 0.5), (0.3, 0.8), (0.7, 0.1)])
    def test_identical_layers_equal_lazy_single_layer(self, gamma, delta):
        # transition="layer": a jump to the same gene in the other (identical)
        # layer is a lazy step: single-layer RWR with restart gamma / (1 - (1 - gamma) delta)
        net = _random_net()
        res = rwr_multiplex({"a": net, "b": net}, self.SEEDS, gamma=gamma, delta=delta,
                            aggregate="sum", transition="layer", tol=1e-14).ranked_genes
        np.testing.assert_allclose(res["p_a"], res["p_b"], atol=1e-15)
        ref = gcen2.rwr.predict(net, self.SEEDS, gamma=gamma / (1 - (1 - gamma) * delta),
                                prob_cutoff=0, tol=1e-14).ranked_genes
        m = ref.merge(res, on="gene")
        np.testing.assert_allclose(m["probability_x"], m["probability_y"], atol=1e-12)

    @pytest.mark.parametrize("weighted", [False, True])
    @pytest.mark.parametrize("delta", [0.5, 0.2])
    def test_matches_randomwalkrestartmh(self, weighted, delta):
        # tests/data/rwrmh/make_reference.R: RandomWalkRestartMH, r = 0.7
        d = DATA / "rwrmh"
        layers = {k: read_layer(d / f"{k}.tsv") for k in "abc"}
        res = rwr_multiplex(layers, (d / "seeds.txt").read_text(encoding="utf-8").split(), gamma=0.7,
                            delta=delta, weighted=weighted, tol=1e-14).ranked_genes
        ref = pd.read_csv(d / f"{'weighted' if weighted else 'unweighted'}_delta{delta}.tsv",
                          sep="\t")
        m = ref.merge(res, left_on="NodeNames", right_on="gene")
        assert len(m) == len(res) == len(ref) - 3  # seeds excluded
        # RandomWalkRestartMH stops at an L2 change below 1e-10
        np.testing.assert_allclose(m["probability"], m["Score"], rtol=0, atol=1e-11)

    def test_transitions_differ_and_delta_zero(self):
        a, b = _random_net(seed=1), _random_net(seed=2)
        sup = rwr_multiplex({"a": a, "b": b}, self.SEEDS).ranked_genes.set_index("gene")
        lay = rwr_multiplex({"a": a, "b": b}, self.SEEDS,
                            transition="layer").ranked_genes.set_index("gene")
        assert not np.allclose(sup["probability"], lay.loc[sup.index, "probability"])
        for t in ("supra", "layer"):  # genes without edges in a layer still jump
            res = rwr_multiplex({"a": a, "b": b}, self.SEEDS, delta=0.0,
                                transition=t, aggregate="sum").ranked_genes
            assert np.isfinite(res["probability"]).all()

    def test_layer_only_gene_reached_through_jump(self):
        # D is only in layer b; it is reached from the seed A in layer a by
        # jumping to b, so it gets probability in both layers
        a = _edges([("A", "B", 1.0)])
        b = _edges([("A", "D", 1.0)])
        res = rwr_multiplex({"a": a, "b": b}, ["A"]).ranked_genes.set_index("gene")
        assert res.loc["D", "p_a"] > 0 and res.loc["D", "p_b"] > res.loc["D", "p_a"]
        assert res.loc["D", "probability"] == pytest.approx(
            np.sqrt(res.loc["D", "p_a"] * res.loc["D", "p_b"]))

    def test_aggregates(self):
        a, b = _random_net(seed=1), _random_net(seed=2)
        g = rwr_multiplex({"a": a, "b": b}, self.SEEDS).ranked_genes
        np.testing.assert_allclose(g["probability"], np.sqrt(g["p_a"] * g["p_b"]))
        m = rwr_multiplex({"a": a, "b": b}, self.SEEDS, aggregate="arithmetic").ranked_genes
        np.testing.assert_allclose(m["probability"], (m["p_a"] + m["p_b"]) / 2)
        assert m["probability"].is_monotonic_decreasing

    def test_seeds_excluded_and_missing_seeds_ignored(self):
        res = rwr_multiplex([_random_net()], ["G1", "NOT_A_GENE"]).ranked_genes
        assert "G1" not in set(res["gene"])
        with pytest.raises(GCENInputError, match="None of the seed"):
            rwr_multiplex([_random_net()], ["NOT_A_GENE"])

    @pytest.mark.parametrize("kwargs,match", [
        ({"gamma": 1.0}, "gamma"),
        ({"delta": 1.5}, "delta"),
        ({"aggregate": "median"}, "aggregate"),
        ({"transition": "row"}, "transition"),
    ])
    def test_bad_arguments(self, kwargs, match):
        with pytest.raises(GCENInputError, match=match):
            rwr_multiplex([_random_net()], self.SEEDS, **kwargs)

    def test_python_fallback_matches(self, monkeypatch):
        a, b = _random_net(seed=1), _random_net(seed=2)
        ref = rwr_multiplex({"a": a, "b": b}, self.SEEDS).ranked_genes
        monkeypatch.setattr("gcen2.core._gcen_core", None)
        res = rwr_multiplex({"a": a, "b": b}, self.SEEDS).ranked_genes
        pd.testing.assert_frame_equal(ref, res, check_exact=False, atol=1e-12)


class TestCli:
    def _files(self, tmp_path):
        _random_net(seed=1).to_edgelist(tmp_path / "coexp.network")
        (tmp_path / "ppi.tsv").write_text(
            "p1\tp2\tcombined_score\nG1\tG5\t900\nG5\tG7\t500\nG2\tG9\t300\n", encoding="utf-8")
        (tmp_path / "seeds.txt").write_text("G1\nG2\n", encoding="utf-8")
        return [tmp_path / "coexp.network", tmp_path / "ppi.tsv"]

    def test_fuse(self, tmp_path):
        from typer.testing import CliRunner

        from gcen2.cli import app

        a, b = self._files(tmp_path)
        out = tmp_path / "fused.tsv"
        res = CliRunner().invoke(app, ["fuse", "-i", str(a), "-i", str(b), "-w", "2", "-w",
                                       "1", "-o", str(out)])
        assert res.exit_code == 0, res.output
        got = pd.read_csv(out, sep="\t", dtype={"gene_a": str, "gene_b": str})
        ref = fuse({"coexp": read_layer(a), "ppi": read_layer(b)},
                   weights={"coexp": 2, "ppi": 1}).edges
        assert list(got.columns) == list(ref.columns)
        np.testing.assert_allclose(got["score"], ref["score"])

    def test_fuse_name_count_mismatch(self, tmp_path):
        from typer.testing import CliRunner

        from gcen2.cli import app

        a, b = self._files(tmp_path)
        res = CliRunner().invoke(app, ["fuse", "-i", str(a), "-i", str(b), "--name", "x",
                                       "-o", str(tmp_path / "f.tsv")])
        assert res.exit_code == 1

    def test_rwr_multiplex(self, tmp_path):
        from typer.testing import CliRunner

        from gcen2.cli import app

        a, b = self._files(tmp_path)
        out = tmp_path / "rwr.tsv"
        res = CliRunner().invoke(app, ["rwr-multiplex", "-n", str(a), "-n", str(b),
                                       "--name", "coexp", "--name", "ppi",
                                       "-g", str(tmp_path / "seeds.txt"), "-o", str(out)])
        assert res.exit_code == 0, res.output
        got = pd.read_csv(out, sep="\t")
        assert list(got.columns) == ["gene", "probability", "p_coexp", "p_ppi"]
        ref = rwr_multiplex({"coexp": read_layer(a), "ppi": read_layer(b)},
                            ["G1", "G2"]).ranked_genes
        np.testing.assert_allclose(got["probability"], ref["probability"])


class TestFuseToFile:
    def _layers(self):
        a = _edges([("A", "B", 1.0), ("B", "C", 1e-7), ('q"t', "tab\there", 0.3),
                    ("nl\nx", "é", 2.5e20), ("C", "D", 123456789.123)])
        b = _edges([("B", "A", 2.0), ("C", "D", 1.0), ("é", "Z", 0.1), ("1", "2", 0.0)])
        return {"co,exp": a, 'p\tpi"': b}

    @pytest.mark.parametrize("normalize", ["max", "rank", "none"])
    def test_same_bytes_as_pandas(self, tmp_path, normalize):
        from gcen2.multilayer import fuse_to_file

        layers = self._layers()
        fuse(layers, normalize=normalize).to_edgelist(tmp_path / "ref.tsv")
        r = fuse_to_file(layers, tmp_path / "out.tsv", normalize=normalize, threads=3)
        assert (tmp_path / "out.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()
        ref = fuse(layers, normalize=normalize)
        assert r["edge_count"] == ref.edge_count and r["node_count"] == ref.node_count

    def test_random_networks_many_blocks(self, tmp_path):
        from gcen2.multilayer import fuse_to_file

        layers = {"a": _random_net(m=60000, n_genes=3000, seed=1),
                  "b": _random_net(m=40000, n_genes=3000, seed=2),
                  "c": _random_net(m=500, n_genes=3000, seed=3)}
        fuse(layers, weights={"a": 3}, min_score=0.05).to_edgelist(tmp_path / "ref.tsv")
        r = fuse_to_file(layers, tmp_path / "out.tsv", weights={"a": 3}, min_score=0.05,
                         threads=4)
        assert (tmp_path / "out.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()
        counts = fuse(layers, weights={"a": 3}, min_score=0.05).edges["n_layers"]
        assert r["layer_counts"] == {int(k): int(v) for k, v in
                                     counts.value_counts().sort_index().items()}

    def test_empty_result(self, tmp_path):
        from gcen2.multilayer import fuse_to_file

        layers = self._layers()
        fuse(layers, min_score=10).to_edgelist(tmp_path / "ref.tsv")
        r = fuse_to_file(layers, tmp_path / "out.tsv", min_score=10)
        assert r["edge_count"] == 0 and r["layer_counts"] == {}
        assert (tmp_path / "out.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()

    def test_without_extension(self, tmp_path, monkeypatch):
        from gcen2.multilayer import fuse_to_file

        layers = self._layers()
        fuse(layers).to_edgelist(tmp_path / "ref.tsv")
        monkeypatch.setattr("gcen2.core._gcen_core", None)
        r = fuse_to_file(layers, tmp_path / "out.tsv")
        assert (tmp_path / "out.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()
        assert sum(r["layer_counts"].values()) == r["edge_count"]

    def test_bad_code_raises(self, tmp_path):
        from gcen2.core import _gcen_core

        with pytest.raises(IndexError):
            _gcen_core.network.write_coded_table(
                str(tmp_path / "x.tsv"), "a\tb\n", ["A"], np.array([0]), np.array([1]), [])


class TestLayerReaderCpp:
    """read_layer (C++ when certain) against the pandas reader."""

    CASES = {
        "gcen": "gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t-0.5\t1e-3\nB\tC\t0.25\t0.1\n",
        "no_header": "A\tB\t0.5\nB\tC\t-1.5e-7\n",
        "two_cols": "TF1\tG1\nTF1\tG2\n\nTF2\tG1\n",
        "numbers": "a\tb\tw\nA\tB\t1.\nB\tC\t.5\nC\tD\t-0\nD\tE\t1E5\nE\tF\t00.10\n"
                   "F\tG\t3778353374068124117\nG\tH\t0.324682034176848787\n",
        "crlf": "a\tb\tw\r\nA\tB\t1\r\nB\tC\t2\r\n",
        "short_row": "a\tb\tw\textra\nA\tB\t1\nB\tC\t2\tx\n",
        "unicode": "a\tb\tw\né\t基因\t0.5\n",
        "numeric_ids": "0007\t12\t1\n12\t0007\t2\n",
        "blank_lead": "\n\nA\tB\t1\n",
        "header_only": "a\tb\tweight\n",
        "weight_col": "a\tb\tscore\tother\nA\tB\t0.5\t7\nB\tC\t0.25\t-8\n",
        "nan_header": "A\tB\tNA\nB\tC\t1\n",
        # declined by C++, read by pandas
        "plus": "a\tb\tw\nA\tB\t+0.5\n",
        "quoted": 'a\tb\tw\n"A x"\tB\t1\n',
        "bom": "﻿a\tb\tw\nA\tB\t1\n",
        "inf": "a\tb\tw\nA\tB\tinf\n",
        "padded": "a\tb\tw\nA\tB\t 0.5\n",
    }
    DECLINED = {"plus", "quoted", "bom", "inf", "padded"}

    def _compare(self, path, **kw):
        from gcen2.multilayer import _read_layer_pandas

        ref = _read_layer_pandas(path, kw.get("weight_column"), kw.get("sep", "\t"))
        got = read_layer(path, **kw)
        # pandas < 3 parses some numbers 1 ulp off (the C++ reader rounds
        # correctly, as pandas 3 does)
        exact = int(pd.__version__.split(".")[0]) >= 3
        pd.testing.assert_frame_equal(got, ref, check_exact=exact, rtol=1e-15, atol=0)

    @pytest.mark.parametrize("name", sorted(CASES))
    def test_same_table_as_pandas(self, tmp_path, name):
        p = tmp_path / f"{name}.tsv"
        p.write_bytes(self.CASES[name].encode("utf-8"))
        self._compare(p)
        if name == "weight_col":
            self._compare(p, weight_column="other")
            self._compare(p, weight_column=3)

    @pytest.mark.parametrize("name", sorted(CASES))
    def test_cpp_handles_or_declines(self, tmp_path, name):
        from gcen2.core import _gcen_core
        from gcen2.multilayer import _layout

        p = tmp_path / f"{name}.tsv"
        p.write_bytes(self.CASES[name].encode("utf-8"))
        first = pd.read_csv(p, sep="\t", header=None, dtype=str, nrows=1)
        has_header, header, col = _layout(first, None, p)
        r = _gcen_core.network.read_layer_file(str(p), "\t", has_header, first.shape[1],
                                               -1 if col is None else header.index(col))
        assert r["ok"] == (name not in self.DECLINED), r["reason"]

    def test_space_separated(self, tmp_path):
        p = tmp_path / "ppi.txt"
        p.write_text("protein1 protein2 combined_score\nA B 900\nB C 400\n", encoding="utf-8")
        self._compare(p, sep=" ")

    def test_errors_match_pandas(self, tmp_path):
        p = tmp_path / "bad.tsv"
        p.write_text("a\tb\tw\nA\tB\tx\n", encoding="utf-8")
        with pytest.raises(GCENInputError, match="non-numeric"):
            read_layer(p)
        p.write_text("A\tB\t1\nB\tC\t1\t2\n", encoding="utf-8")  # longer row: pandas tokenizer error
        with pytest.raises(gcen2.GCENIOError):
            read_layer(p)
        p.write_text("A\n", encoding="utf-8")
        with pytest.raises(GCENInputError, match="two columns"):
            read_layer(p)

    def test_random_network_file(self, tmp_path):
        net = _random_net(m=20000, n_genes=500, seed=4)
        net.to_edgelist(tmp_path / "n.network")
        self._compare(tmp_path / "n.network")

    def test_codes_same_as_files_as_tables(self, tmp_path):
        a, b = _random_net(seed=1), _random_net(seed=2, n_genes=90)
        a.to_edgelist(tmp_path / "a.network")
        b.to_edgelist(tmp_path / "b.network")
        from_files = fuse({"a": tmp_path / "a.network", "b": str(tmp_path / "b.network")})
        from_tables = fuse({"a": read_layer(tmp_path / "a.network"),
                            "b": read_layer(tmp_path / "b.network")})
        pd.testing.assert_frame_equal(from_files.edges, from_tables.edges)
        r1 = rwr_multiplex({"a": tmp_path / "a.network", "b": tmp_path / "b.network"},
                           ["G1", "G2"]).ranked_genes
        r2 = rwr_multiplex({"a": read_layer(tmp_path / "a.network"),
                            "b": read_layer(tmp_path / "b.network")}, ["G1", "G2"]).ranked_genes
        pd.testing.assert_frame_equal(r1, r2)

    def test_unique_ints(self):
        from gcen2.multilayer import _unique_ints

        x = np.random.default_rng(0).integers(0, 50, 1000)
        u, inv = _unique_ints(x, return_inverse=True)
        ru, rinv = np.unique(x, return_inverse=True)
        np.testing.assert_array_equal(u, ru)
        np.testing.assert_array_equal(inv, rinv)
        np.testing.assert_array_equal(_unique_ints(x), ru)
        assert len(_unique_ints(np.empty(0, np.int64))) == 0


class TestFuseCpp:
    """C++ fuse_layer_codes against the numpy reference _fuse_arrays_numpy."""

    @staticmethod
    def _coded(n_layers, n_genes, m, seed, levels=None):
        rng = np.random.default_rng(seed)
        out = []
        for _ in range(n_layers):
            mk = int(m * rng.uniform(0.2, 1.0))
            a = rng.integers(0, n_genes, mk)
            b = rng.integers(0, n_genes, mk)  # self-loops and repeats included
            w = (rng.integers(0, levels, mk) / 7.0 if levels else rng.uniform(0, 3, mk))
            out.append((a.astype(np.int64), b.astype(np.int64), w))
        return out

    @staticmethod
    def _check(coded, n_genes, normalize, alpha, min_layers=1, min_score=0.0, threads=4):
        from gcen2.core import _gcen_core
        from gcen2.multilayer import _fuse_arrays_numpy

        names = [f"l{k}" for k in range(len(coded))]
        r = _gcen_core.network.fuse_layer_codes(coded, n_genes, normalize, alpha, min_layers,
                                                min_score, threads)
        assert r["ok"], r["reason"]
        f = _fuse_arrays_numpy(names, n_genes, coded, normalize, dict(zip(names, alpha)),
                               min_layers, min_score)
        for key in ("a", "b", "score", "n_layers", "label_idx"):
            np.testing.assert_array_equal(r[key], f[key], err_msg=key)
        assert len(r["wide"]) == len(f["wide"])
        for x, y in zip(r["wide"], f["wide"]):
            np.testing.assert_array_equal(x, y)
        labels = [",".join(n for k, n in enumerate(names) if v >> k & 1) for v in r["patterns"]]
        assert labels == f["labels"] and r["node_count"] == f["node_count"]

    @pytest.mark.parametrize("normalize", ["max", "rank", "none"])
    @pytest.mark.parametrize("levels", [None, 5])  # 5: many tied weights
    def test_small(self, normalize, levels):
        coded = self._coded(3, 40, 400, seed=7, levels=levels)
        self._check(coded, 40, normalize, [1.0, 0.3, 2.5])
        self._check(coded, 40, normalize, [1.0, 0.0, 1.0], min_layers=2, min_score=0.2)

    @pytest.mark.parametrize("normalize", ["max", "rank"])
    @pytest.mark.parametrize("threads", [1, 3, 8])
    def test_large_parallel(self, normalize, threads):
        # enough rows for the parallel sorts and merges to split
        coded = self._coded(2, 3000, 400_000, seed=11, levels=50)
        self._check(coded, 3000, normalize, [1.0, 2.0], threads=threads)

    def test_many_layers(self):
        # more than 16 layers: patterns indexed through a hash map
        coded = self._coded(18, 30, 60, seed=3)
        self._check(coded, 30, "max", [1.0] * 18)

    def test_empty_and_all_zero(self):
        empty = (np.empty(0, np.int64), np.empty(0, np.int64), np.empty(0))
        zero = (np.array([0, 1]), np.array([1, 2]), np.zeros(2))
        self._check([empty, zero], 3, "max", [1.0, 1.0])
        self._check([empty, empty], 3, "rank", [1.0, 1.0])
        self._check([zero], 0 + 3, "none", [1.0], min_score=1.0)

    def test_nan_declines_and_fuse_falls_back(self, monkeypatch):
        from gcen2.core import _gcen_core

        lay = (np.array([0, 1]), np.array([1, 2]), np.array([1.0, np.nan]))
        r = _gcen_core.network.fuse_layer_codes([lay], 3, "max", [1.0], 1, 0.0)
        assert not r["ok"] and "NaN" in r["reason"]
        layers = {"a": _edges([("A", "B", 1.0), ("B", "C", np.nan)]),
                  "b": _edges([("A", "B", 0.5)])}
        # numpy's NaN in np.maximum.at warns on x86 but not on ARM: not checked
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", RuntimeWarning)
            got = fuse(layers).edges
            monkeypatch.setattr("gcen2.core._gcen_core", None)
            ref = fuse(layers).edges
        pd.testing.assert_frame_equal(got, ref)

    def test_bad_code_raises(self):
        from gcen2.core import _gcen_core

        lay = (np.array([0, 5]), np.array([1, 2]), np.array([1.0, 1.0]))
        with pytest.raises(ValueError, match="outside"):
            _gcen_core.network.fuse_layer_codes([lay], 3, "max", [1.0], 1, 0.0)
        with pytest.raises(ValueError, match="normalize"):
            _gcen_core.network.fuse_layer_codes([lay], 6, "x", [1.0], 1, 0.0)
