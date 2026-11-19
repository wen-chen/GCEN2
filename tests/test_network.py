"""Tests for gcen2.network — utility functions (C++-independent)."""

import numpy as np
import pandas as pd
import pytest

from gcen2 import network as net_module
from gcen2.exceptions import GCENInputError
from gcen2.models import Network


class TestNetworkMerge:
    def _make_net(self, edges: list[tuple]) -> Network:
        df = pd.DataFrame(edges, columns=["gene_a", "gene_b", "correlation"])
        return Network(edges=df)

    def test_merge_sums_correlations(self):
        a = self._make_net([("G1", "G2", 0.8), ("G1", "G3", 0.7)])
        b = self._make_net([("G1", "G2", 0.6), ("G2", "G3", 0.5)])
        merged = net_module.merge([a, b], cor_cutoff=0.0)
        row = merged.edges[
            (merged.edges["gene_a"] == "G1") & (merged.edges["gene_b"] == "G2")
        ]
        assert len(row) == 1
        assert abs(row.iloc[0]["correlation"] - 1.4) < 1e-9

    def test_merge_direction_is_separate_column(self):
        a = self._make_net([("G1", "G2", 0.8)])
        b = self._make_net([("G1", "G2", -0.6)])
        merged = net_module.merge([a, b], cor_cutoff=0.0)
        assert "direction" in merged.edges.columns
        # correlation is not a string — GCEN v1 bug fixed
        assert merged.edges["correlation"].dtype != object

    def test_merge_cutoff_filters(self):
        a = self._make_net([("G1", "G2", 0.3)])
        merged = net_module.merge([a], cor_cutoff=0.5)
        assert merged.edges.empty

    def test_merge_empty_raises(self):
        from gcen2.exceptions import GCENInputError
        with pytest.raises(GCENInputError):
            net_module.merge([])


class TestNetworkExtract:
    def test_extract_subnetwork(self, small_network):
        sub = net_module.extract(small_network, ["GeneA", "GeneB"])
        genes_in_sub = set(sub.edges["gene_a"]) | set(sub.edges["gene_b"])
        assert "GeneA" in genes_in_sub or "GeneB" in genes_in_sub


class TestNetworkShuffle:
    def test_shuffle_preserves_degree(self, small_network):
        pytest.importorskip("networkx")

        original_G = small_network.to_networkx()
        shuffled = net_module.shuffle(small_network, seed=42)
        shuffled_G = shuffled.to_networkx()

        orig_degrees = sorted(d for _, d in original_G.degree())
        shuf_degrees = sorted(d for _, d in shuffled_G.degree())
        assert orig_degrees == shuf_degrees

    def test_shuffle_is_reproducible(self, small_network):
        s1 = net_module.shuffle(small_network, seed=123)
        s2 = net_module.shuffle(small_network, seed=123)
        pd.testing.assert_frame_equal(
            s1.edges.reset_index(drop=True),
            s2.edges.reset_index(drop=True),
        )


class TestNetworkStats:
    def test_stats_fields(self, small_network):
        s = net_module.stats(small_network)
        assert s.node_count == 5
        assert s.edge_count == 5
        assert 0.0 <= s.density <= 1.0
        assert s.average_degree > 0
        assert s.connected_components >= 1


class TestBenjaminiHochberg:
    def test_fdr_monotonic(self):
        p = np.array([0.01, 0.04, 0.1, 0.002, 0.3])
        fdr = net_module._benjamini_hochberg(p, n_total=100)
        # FDR values for sorted p-values should be non-decreasing
        sorted_p_idx = np.argsort(p)
        fdr_sorted = fdr[sorted_p_idx]
        assert np.all(np.diff(fdr_sorted) >= -1e-12)

    def test_fdr_leq_1(self):
        p = np.array([0.001, 0.01, 0.05])
        fdr = net_module._benjamini_hochberg(p, n_total=1000)
        assert np.all(fdr <= 1.0)


# ---------------------------------------------------------------------------
# Regression: in-kernel edge filtering must not change results
# ---------------------------------------------------------------------------

def _reference_bh(p_values, n_total):
    """The original loop-based BH implementation (full set of tests)."""
    m = n_total
    order = np.argsort(p_values)[::-1]
    fdr = np.empty_like(p_values)
    prev = 1.0
    for rank, idx in enumerate(order):
        prev = min(p_values[idx] * m / (m - rank), prev)
        fdr[idx] = prev
    return fdr


def _reference_build(expr, method, cor, pval, signed, fdr):
    """Old behaviour: compute every pair, then filter in Python."""
    from gcen2.network import build
    full = build(expr, method=method, all_edges=True, threads=1).edges
    if fdr:
        n = len(expr)
        full["fdr"] = _reference_bh(full["p_value"].to_numpy(), n * (n - 1) // 2)
    mask = (full["correlation"].abs() >= cor) & (full["p_value"] <= pval)
    if fdr:
        mask &= full["fdr"] <= pval
    if signed:
        mask &= full["correlation"] >= 0
    return full[mask]


def _sorted(edges):
    return edges.sort_values(["gene_a", "gene_b"]).reset_index(drop=True)


class TestKernelFilterEquivalence:
    @pytest.fixture(scope="class")
    def expr(self):
        rng = np.random.default_rng(123)
        n_genes, n_samples = 80, 12
        # shared latent factors -> a realistic mix of strong +/- correlations
        latent = rng.standard_normal((3, n_samples))
        load = rng.standard_normal((n_genes, 3))
        data = load @ latent + 0.8 * rng.standard_normal((n_genes, n_samples))
        return pd.DataFrame(np.exp(data),
                            index=[f"g{i:03d}" for i in range(n_genes)],
                            columns=[f"s{j}" for j in range(n_samples)])

    @pytest.mark.parametrize("method", ["pearson", "spearman"])
    @pytest.mark.parametrize("signed", [False, True])
    @pytest.mark.parametrize("fdr", [False, True])
    def test_matches_unfiltered_reference(self, expr, method, signed, fdr):
        from gcen2.network import build
        cor, pval = 0.5, 0.01
        got = _sorted(build(expr, method=method, cor=cor, pval=pval,
                            signed=signed, fdr=fdr, threads=3).edges)
        ref = _sorted(_reference_build(expr, method, cor, pval, signed, fdr))
        assert len(ref) > 0, "fixture should produce some edges"
        pd.testing.assert_frame_equal(got[ref.columns], ref, check_exact=False,
                                      rtol=0, atol=1e-12)

    def test_signed_semantics_match_v1(self, expr):
        """v1 ``-s``: positive edges only; default: both signs by |r|."""
        from gcen2.network import build
        unsigned = build(expr, cor=0.3, pval=0.05, threads=2)
        signed = build(expr, cor=0.3, pval=0.05, signed=True, threads=2)
        assert (unsigned.edges["correlation"] < 0).any()
        assert (signed.edges["correlation"] >= 0).all()
        pos = unsigned.edges[unsigned.edges["correlation"] >= 0]
        assert len(pos) == len(signed.edges)
        assert signed.is_signed and not unsigned.is_signed

    def test_kernel_defaults_return_all_pairs(self, expr):
        from gcen2.core import _gcen_core
        m = np.ascontiguousarray(expr.to_numpy(dtype=np.float64))
        idx, _, _ = _gcen_core.correlation.pairwise_correlation(m, "pearson", 2)
        n = len(expr)
        assert len(idx) == n * (n - 1) // 2


class TestBenjaminiHochbergSubset:
    def test_full_set_matches_reference(self):
        from gcen2.network import _benjamini_hochberg
        p = np.random.default_rng(0).uniform(size=500) ** 3
        p[::50] = p[1]  # ties
        np.testing.assert_allclose(_benjamini_hochberg(p, p.size),
                                   _reference_bh(p, p.size), rtol=0, atol=1e-15)

    def test_subset_matches_full_below_cutoff(self):
        from gcen2.network import _benjamini_hochberg
        p = np.random.default_rng(1).uniform(size=2000) ** 4
        cutoff = 0.01
        full = _reference_bh(p, p.size)
        keep = p <= cutoff
        sub = _benjamini_hochberg(p[keep], p.size)
        passed = full[keep] <= cutoff
        np.testing.assert_array_equal(sub <= cutoff, passed)
        np.testing.assert_allclose(sub[passed], full[keep][passed], rtol=0, atol=1e-15)

    def test_empty(self):
        from gcen2.network import _benjamini_hochberg
        assert _benjamini_hochberg(np.array([]), 10).size == 0


class TestStatsMatchNetworkx:
    """Regression: sparse NetworkStats must equal the networkx definitions."""

    def _nx_stats(self, net):
        nx = pytest.importorskip("networkx")
        G = nx.from_pandas_edgelist(net.edges, "gene_a", "gene_b")
        deg = np.array([d for _, d in G.degree()])
        comps = list(nx.connected_components(G))
        n, e = G.number_of_nodes(), G.number_of_edges()
        return dict(node_count=n, edge_count=e, density=2 * e / (n * (n - 1)),
                    average_degree=deg.mean(), max_degree=deg.max(),
                    connected_components=len(comps),
                    giant_component_size=max(len(c) for c in comps),
                    heterogeneity=deg.std() / deg.mean(),
                    clustering_coefficient=nx.average_clustering(G),
                    transitivity=nx.transitivity(G))

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_random_graph(self, seed):
        rng = np.random.default_rng(seed)
        genes = [f"g{i}" for i in range(70)]
        rows = [(genes[a], genes[b], 0.5) for a, b in rng.integers(0, 70, size=(400, 2))]
        rows += [(b, a, c) for a, b, c in rows[:30]]  # reverse duplicates
        net = Network(edges=pd.DataFrame(rows, columns=["gene_a", "gene_b", "correlation"]))
        got = net_module.stats(net)
        for field, expected in self._nx_stats(net).items():
            assert getattr(got, field) == pytest.approx(expected, rel=1e-12, abs=1e-12), field

    def test_small_chunks_same_result(self, small_network):
        from gcen2.models import NetworkStats
        assert NetworkStats.compute(small_network, chunk_size=2) == \
            NetworkStats.compute(small_network)


class TestStatsKernel:
    """C++ graph_stats (NetworkStats.compute, stats_from_file) against the
    scipy implementation, on graphs with self-loops, duplicates and hubs."""

    def _net(self, seed, n=80, m=600, loops=10, hub=True):
        rng = np.random.default_rng(seed)
        genes = [f"g{i:03d}" for i in range(n)]
        pairs = [tuple(p) for p in rng.integers(0, n, size=(m, 2))]
        pairs += [(i, i) for i in rng.integers(0, n, loops)]
        if hub:  # a hub linked to most nodes, plus a clique
            pairs += [(0, j) for j in range(1, n, 2)] + [
                (a, b) for a in range(5, 15) for b in range(a + 1, 15)]
        pairs += [(b, a) for a, b in pairs[:40]]
        rows = [(genes[a], genes[b], 0.5, 0.01) for a, b in pairs]
        return Network(edges=pd.DataFrame(
            rows, columns=["gene_a", "gene_b", "correlation", "p_value"]))

    @pytest.mark.parametrize("seed", [0, 1, 2])
    def test_matches_sparse(self, seed):
        from gcen2.models import NetworkStats
        net = self._net(seed, hub=seed != 2)
        got, ref = NetworkStats.compute(net), NetworkStats._compute_sparse(net)
        for f in ("node_count", "edge_count", "max_degree", "connected_components",
                  "giant_component_size"):
            assert getattr(got, f) == getattr(ref, f), f
        for f in ("density", "average_degree", "heterogeneity", "clustering_coefficient",
                  "transitivity"):
            assert getattr(got, f) == pytest.approx(getattr(ref, f), rel=1e-12), f

    def test_threads_and_file_path_identical(self, tmp_path):
        from gcen2.io import read_network
        from gcen2.models import NetworkStats
        path = tmp_path / "n.tsv"
        self._net(5).to_edgelist(path)
        net = read_network(path)
        ref = NetworkStats.compute(net, threads=1)
        assert NetworkStats.compute(net, threads=7) == ref
        assert net_module.stats_from_file(path, threads=3) == ref

    def test_disconnected(self):
        edges = pd.DataFrame({"gene_a": ["a", "c", "e", "e"], "gene_b": ["b", "d", "f", "g"],
                              "correlation": 0.9, "p_value": 0.0})
        s = net_module.stats(Network(edges=edges))
        assert (s.node_count, s.edge_count, s.connected_components,
                s.giant_component_size) == (7, 4, 3, 3)
        assert s.clustering_coefficient == 0.0 and s.transitivity == 0.0


class TestShuffleNoDuplicates:
    def test_no_reverse_duplicates_or_self_loops(self):
        rng = np.random.default_rng(4)
        genes = [f"g{i}" for i in range(25)]
        pairs = {tuple(sorted(rng.choice(25, 2, replace=False))) for _ in range(150)}
        net = Network(edges=pd.DataFrame(
            [(genes[a], genes[b], 0.9) for a, b in pairs],
            columns=["gene_a", "gene_b", "correlation"]))
        out = net_module.shuffle(net, swap_factor=50, seed=1).edges
        undirected = {frozenset(e) for e in zip(out["gene_a"], out["gene_b"])}
        assert len(undirected) == len(out)                      # no (a,b)/(b,a) dupes
        assert (out["gene_a"] != out["gene_b"]).all()           # no self-loops
        changed = set(zip(out["gene_a"], out["gene_b"])) != set(zip(net.edges["gene_a"], net.edges["gene_b"]))
        assert changed


class TestShuffleKernel:
    """The C++ EdgeSwapper applies the same numpy random stream as the
    pure-Python loop, so a seed gives the same network either way."""

    @pytest.mark.parametrize("seed,dup", [(0, False), (3, True), (8, True)])
    def test_same_as_python_loop(self, monkeypatch, seed, dup):
        import gcen2.core
        rng = np.random.default_rng(seed)
        genes = [f"g{i}" for i in range(40)]
        rows = [(genes[a], genes[b], 0.5 + i * 1e-3)
                for i, (a, b) in enumerate(rng.integers(0, 40, size=(160, 2))) if a != b]
        if dup:  # repeated and reversed edges in the input
            rows += rows[:10] + [(b, a, c) for a, b, c in rows[10:20]]
        net = Network(edges=pd.DataFrame(rows, columns=["gene_a", "gene_b", "correlation"]))
        fast = net_module.shuffle(net, swap_factor=30, seed=seed).edges
        monkeypatch.setattr(gcen2.core, "_gcen_core", None)
        slow = net_module.shuffle(net, swap_factor=30, seed=seed).edges
        pd.testing.assert_frame_equal(fast, slow)
        assert not fast["gene_b"].equals(net.edges["gene_b"])

    def test_key_set_survives_many_swaps(self):
        from gcen2.core import _gcen_core
        if _gcen_core is None:
            pytest.skip("compiled extension not available")
        a = np.arange(0, 2000, 2)
        b = a + 1
        s = _gcen_core.network.EdgeSwapper(a, b, 2000)
        rng = np.random.default_rng(1)
        applied = s.run(rng.integers(0, 1000, 200_000), rng.integers(0, 1000, 200_000))
        assert applied > 1000
        nb = s.b()
        assert sorted(nb.tolist()) == sorted(b.tolist())          # degrees kept
        assert len({(min(x, y), max(x, y)) for x, y in zip(a, nb)}) == 1000


class TestFileMergeExtract:
    """merge_files / extract_file (streamed by the C++ extension) give the
    same file as the pandas pipeline: read_network -> merge/extract ->
    to_edgelist."""

    _HEADERS = {
        "v1": "#Node1\tNode2\tCorrelation\tPvalue\tFDR\n",
        "gcen2": "gene_a\tgene_b\tcorrelation\tp_value\n",
        "none": "",
    }

    def _write(self, path, seed, header, n_rows=400, quoted=False, nan=False):
        rng = np.random.default_rng(seed)
        genes = [f"g{i:02d}" for i in range(30)] + ["G10", "a b"]
        lines = [self._HEADERS[header]]
        for k in range(n_rows):
            a, b = rng.choice(len(genes), 2)
            cor = rng.uniform(-1, 1) * 10.0 ** -int(rng.integers(0, 4))
            ga, gb = genes[a], genes[b]
            if quoted and k % 7 == 0:
                ga = f'"{ga}"'
            c = "NA" if nan and k % 11 == 0 else repr(cor)
            # v1 header: FDR column sometimes missing; header-less: always 5 columns
            extra = {"v1": "" if k % 5 == 0 else "\t1e-3", "gcen2": "", "none": "\t1e-3"}[header]
            lines.append(f"{ga}\t{gb}\t{c}\t2.5E-07{extra}\n")
            if k % 9 == 0:  # repeated and reversed pairs
                lines.append(f"{gb}\t{ga}\t{repr(cor / 3)}\t1e-2{extra}\n")
        lines.append("\n")
        path.write_text("".join(lines), encoding="utf-8")
        return path

    def _reference(self, monkeypatch, fn, *args):
        import gcen2.core
        monkeypatch.setattr(gcen2.core, "_gcen_core", None)
        fn(*args)
        monkeypatch.undo()

    @pytest.mark.parametrize("cutoff", [0.0, 0.5])
    def test_merge_same_as_pandas(self, tmp_path, monkeypatch, cutoff):
        from gcen2.core import _gcen_core
        if _gcen_core is None:
            pytest.skip("compiled extension not available")
        files = [self._write(tmp_path / "a.tsv", 0, "v1", nan=True),
                 self._write(tmp_path / "b.tsv", 1, "gcen2", quoted=True),
                 self._write(tmp_path / "c.tsv", 2, "none"),
                 self._write(tmp_path / "d.tsv", 3, "v1")]
        got = net_module.merge_files(files, tmp_path / "fast.tsv", cutoff)
        self._reference(monkeypatch, net_module.merge_files, files, tmp_path / "ref.tsv", cutoff)
        from tests.golden.golden import assert_same_file

        assert_same_file(tmp_path / "fast.tsv", tmp_path / "ref.tsv")
        ref = pd.read_csv(tmp_path / "ref.tsv", sep="\t", dtype={"gene_a": str, "gene_b": str})
        assert got == {"node_count": len(set(ref["gene_a"]) | set(ref["gene_b"])),
                       "edge_count": len(ref)}

    def test_merge_quote_in_gene_name(self, tmp_path, monkeypatch):
        # a name holding '"' is written back quoted with the quote doubled,
        # as pandas' to_csv does
        from gcen2.core import _gcen_core
        if _gcen_core is None:
            pytest.skip("compiled extension not available")
        f = tmp_path / "q.tsv"
        # "q"t" is read as qt" (text after the closing quote taken literally)
        f.write_text('gene_a\tgene_b\tcorrelation\n"C""x"\tA\t0.7\nA\tB"\t-0.6\n'
                     '"q"t"\tA\t0.9\n', encoding="utf-8")
        net_module.merge_files([f], tmp_path / "fast.tsv", 0.5)
        self._reference(monkeypatch, net_module.merge_files, [f], tmp_path / "ref.tsv", 0.5)
        assert (tmp_path / "fast.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()
        assert b'"C""x"' in (tmp_path / "fast.tsv").read_bytes()

    def test_merge_kahan_sum(self, tmp_path, monkeypatch):
        # many values for one pair: the compensated sum differs from a plain
        # one. Rows not in pair order, as in real merges (pandas >= 3.1 sums
        # pre-sorted rows with np.add.reduceat instead of its Kahan kernel).
        vals = [1.0] + [1e-16] * 10  # plain sum 1.0, compensated 1.000000000000001
        files = []
        for i, v in enumerate(vals):
            files.append(tmp_path / f"n{i}.tsv")
            files[-1].write_text(f"x\ty\t{v!r}\t0.01\na\tb\t{v!r}\t0.01\n", encoding="utf-8")
        net_module.merge_files(files, tmp_path / "fast.tsv", 0.0)
        self._reference(monkeypatch, net_module.merge_files, files, tmp_path / "ref.tsv", 0.0)
        assert (tmp_path / "fast.tsv").read_text(encoding="utf-8") == (tmp_path / "ref.tsv").read_text(encoding="utf-8")
        assert "\t1.000000000000001\t" in (tmp_path / "fast.tsv").read_text(encoding="utf-8")

    @pytest.mark.parametrize("header", ["v1", "gcen2", "none"])
    def test_extract_same_as_pandas(self, tmp_path, monkeypatch, header):
        path = self._write(tmp_path / "n.tsv", 4, header, quoted=True)
        genes = ["g03", "g17", "a b", "G10", "missing"]
        got = net_module.extract_file(path, genes, tmp_path / "fast.tsv")
        self._reference(monkeypatch, net_module.extract_file, path, genes, tmp_path / "ref.tsv")
        assert (tmp_path / "fast.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()
        assert got["edge_count"] == len((tmp_path / "ref.tsv").read_text(encoding="utf-8").splitlines()) - 1 > 0

    def test_extract_no_match_and_numeric_names(self, tmp_path):
        path = tmp_path / "n.tsv"
        path.write_text("gene_a\tgene_b\tcorrelation\tp_value\n"
                        "1017\t672\t0.9\t0.001\n672\t5925\t-0.8\t0.002\n", encoding="utf-8")
        assert net_module.extract_file(path, ["999"], tmp_path / "e.tsv")["edge_count"] == 0
        assert (tmp_path / "e.tsv").read_text(encoding="utf-8") == "gene_a\tgene_b\tcorrelation\tp_value\n"
        # numeric gene IDs are names: they match a gene list (read as text)
        r = net_module.extract_file(path, ["1017"], tmp_path / "e.tsv")
        assert r == {"node_count": 2, "edge_count": 1}
        assert (tmp_path / "e.tsv").read_text(encoding="utf-8").splitlines()[1] == "1017\t672\t0.900000\t1.000000e-03"


class TestStreamingBuild:
    """build_to_file must write exactly what build().to_edgelist() writes."""

    @pytest.fixture
    def expr(self):
        rng = np.random.default_rng(21)
        n_genes, n_samples = 90, 12
        latent = rng.standard_normal((3, n_samples))
        load = rng.standard_normal((n_genes, 3))
        data = load @ latent + 0.8 * rng.standard_normal((n_genes, n_samples))
        data[5] = data[6]  # ties in p-values
        return pd.DataFrame(np.exp(data), index=[f"g{i:03d}" for i in range(n_genes)],
                            columns=[f"s{j}" for j in range(n_samples)])

    @pytest.mark.parametrize("kw", [
        dict(method="spearman", cor=0.5, pval=0.01),
        dict(method="pearson", cor=0.5, pval=0.01, signed=True),
        dict(method="spearman", cor=0.3, pval=0.05, fdr=True),
        dict(method="pearson", cor=0.3, pval=0.05, fdr=True, signed=True),
        dict(method="pearson", all_edges=True),
        dict(method="spearman", cor=0.5, pval=0.01, log="log2"),
    ])
    @pytest.mark.parametrize("max_pairs", [1, 37, 500, 10**9])
    def test_byte_identical(self, expr, tmp_path, kw, max_pairs):
        from gcen2.network import build, build_to_file
        ref = tmp_path / "ref.tsv"
        got = tmp_path / "got.tsv"
        net = build(expr, threads=3, **kw)
        net.to_edgelist(ref)
        info = build_to_file(expr, got, threads=3, max_pairs_per_block=max_pairs, **kw)
        assert net.edge_count > 0
        assert got.read_bytes() == ref.read_bytes()
        assert info == {"node_count": net.node_count, "edge_count": net.edge_count}

    @pytest.mark.parametrize("kw", [
        dict(method="spearman", cor=0.5, pval=0.01),
        dict(method="pearson", cor=0.3, pval=0.05, fdr=True),
        dict(method="spearman", cor=0.5, pval=0.01, log="log10"),
    ])
    def test_ndarray_input(self, expr, tmp_path, kw):
        """A plain array plus gene names (the pandas-free CLI path) gives the same file."""
        from gcen2.network import build_to_file
        ref, got = tmp_path / "ref.tsv", tmp_path / "got.tsv"
        a = build_to_file(expr, ref, threads=2, **kw)
        b = build_to_file(expr.to_numpy(), got, threads=2, gene_names=list(expr.index), **kw)
        assert got.read_bytes() == ref.read_bytes()
        assert a == b

    def test_ndarray_input_needs_gene_names(self, expr, tmp_path):
        from gcen2.network import build_to_file
        with pytest.raises(GCENInputError):
            build_to_file(expr.to_numpy(), tmp_path / "x.tsv")
        with pytest.raises(GCENInputError):
            build_to_file(expr.to_numpy(), tmp_path / "x.tsv", gene_names=["a"])

    @pytest.mark.parametrize("fdr", [False, True])
    def test_v1_number_format(self, expr, tmp_path, fdr):
        """Numbers are written like GCEN v1: r as %f, p / FDR as %e."""
        from gcen2.network import build, build_to_file
        out = tmp_path / "net.tsv"
        build_to_file(expr, out, cor=0.3, pval=0.05, fdr=fdr, threads=2)
        net = build(expr, cor=0.3, pval=0.05, fdr=fdr, threads=2)
        lines = out.read_text(encoding="utf-8").splitlines()
        assert len(lines) == net.edge_count + 1
        for line, row in zip(lines[1:], net.edges.itertuples(index=False)):
            want = [row.gene_a, row.gene_b, f"{row.correlation:f}", f"{row.p_value:e}"]
            if fdr:
                want.append(f"{row.fdr:e}")
            assert line.split("\t") == want

    def test_python_fallback_matches_compiled(self, expr, tmp_path, monkeypatch):
        import sys

        import gcen2.core
        from gcen2.network import build
        net = build(expr, cor=0.3, pval=0.05, fdr=True, threads=2)
        net.to_edgelist(tmp_path / "c.tsv")
        monkeypatch.delattr(gcen2.core, "_gcen_core", raising=False)
        monkeypatch.setitem(sys.modules, "gcen2.core._gcen_core", None)
        with pytest.raises(ImportError):
            from gcen2.core import _gcen_core  # noqa: F401
        net.to_edgelist(tmp_path / "py.tsv")
        assert (tmp_path / "py.tsv").read_bytes() == (tmp_path / "c.tsv").read_bytes()

    def test_roundtrip_read(self, expr, tmp_path):
        from gcen2.io import read_network
        from gcen2.network import build
        net = build(expr, cor=0.3, pval=0.05, threads=2)
        net.to_edgelist(tmp_path / "n.tsv")
        back = read_network(tmp_path / "n.tsv")
        assert list(back.edges["gene_a"]) == list(net.edges["gene_a"])
        np.testing.assert_allclose(back.edges["correlation"], net.edges["correlation"], atol=5e-7)
        np.testing.assert_allclose(back.edges["p_value"], net.edges["p_value"], rtol=5e-7)

    def test_no_edges_writes_header_only(self, expr, tmp_path):
        from gcen2.network import build_to_file
        out = tmp_path / "empty.tsv"
        info = build_to_file(expr, out, cor=1.01, pval=1e-30)
        assert info == {"node_count": 0, "edge_count": 0}
        assert out.read_text(encoding="utf-8") == "gene_a\tgene_b\tcorrelation\tp_value\n"



def _numpy_edges(matrix, method, cor, pval, signed, fdr, all_edges):
    """The edges of build(): kernel, then BH and filters in numpy (the
    Python implementation the C++ stream_network replaced)."""
    from gcen2.core import _gcen_core
    from gcen2.network import _benjamini_hochberg, _final_mask

    kw = ({} if all_edges else {"max_p": pval} if fdr else
          {"min_abs_cor": cor, "max_p": pval, "positive_only": signed})
    idx, r, p = _gcen_core.correlation.pairwise_correlation(matrix, method, 2, **kw)
    n = len(matrix)
    q = _benjamini_hochberg(p, n * (n - 1) // 2) if fdr else None
    if not all_edges:
        keep = _final_mask(r, p, q, cor, pval, signed)
        idx, r, p = idx[keep], r[keep], p[keep]
        q = q[keep] if fdr else None
    return idx, r, p, q


class TestCppStreaming:
    """network_edges / build_network_array (src/network_stream.cpp) give the
    edges and BH values of the numpy implementation, bit for bit, whatever
    the block size, the threads and the number of edges kept between the
    FDR passes."""

    @staticmethod
    def _matrix():
        rng = np.random.default_rng(1)
        x = rng.normal(size=(90, 25))
        x[10:20] += x[0]                       # a correlated block
        x[30:40] = np.round(rng.poisson(0.7, size=(10, 25)).astype(float))  # ties, zeros
        x[41] = x[42]                          # duplicate genes: tied p-values
        x[50, 3] = np.nan                      # NaN r, infinite p
        x[51] = 1.0                            # constant
        return x

    @pytest.mark.parametrize("method", ["pearson", "spearman", "pearson_nonzero",
                                        "pearson_pairwise_complete"])
    @pytest.mark.parametrize("opts", [
        dict(cor=0.3, pval=0.05, signed=False, fdr=True, all_edges=False),
        dict(cor=0.0, pval=1.0, signed=True, fdr=True, all_edges=False),
        dict(cor=0.3, pval=0.05, signed=False, fdr=True, all_edges=True),
        dict(cor=0.3, pval=0.05, signed=True, fdr=False, all_edges=False),
        dict(cor=0.3, pval=0.05, signed=False, fdr=False, all_edges=True),
    ])
    def test_same_as_numpy(self, tmp_path, method, opts):
        from gcen2.core import _gcen_core
        from gcen2.models import Network

        m = self._matrix()
        genes = [f"g{i}" for i in range(len(m))]
        idx, r, p, q = _numpy_edges(m, method, **opts)
        got = _gcen_core.correlation.network_edges(
            m, method, 3, opts["cor"], opts["pval"], opts["signed"], opts["fdr"],
            opts["all_edges"])
        assert np.array_equal(got[0], idx)
        for a, b in zip(got[1:], (r, p, q)):
            assert (a is None) == (b is None)
            if a is not None:
                assert np.array_equal(a, b, equal_nan=True)
        g = np.array(genes)
        want_net = Network(edges=pd.DataFrame(
            {"gene_a": g[idx[:, 0]], "gene_b": g[idx[:, 1]], "correlation": r,
             "p_value": p} | ({"fdr": q} if opts["fdr"] else {})))
        want = tmp_path / "want.tsv"
        want_net.to_edgelist(want)
        for block, limit, threads in [(1, 0, 1), (97, 3, 4), (500_000, 1 << 24, 2)]:
            out = tmp_path / "got.tsv"
            _gcen_core.correlation.build_network_array(
                m, genes, str(out), method, opts["cor"], opts["pval"], opts["signed"],
                opts["fdr"], opts["all_edges"], threads, block, limit)
            assert out.read_bytes() == want.read_bytes(), (block, limit, threads)

    @pytest.mark.parametrize("method", ["pearson", "spearman", "pearson_nonzero",
                                        "pearson_pairwise_complete"])
    @pytest.mark.parametrize("n_samples", [3, 5, 12, 60])
    def test_p_filter_skips_no_pair(self, method, n_samples):
        """The kernels skip the p-value of pairs below the |r| bound implied
        by max_p; the pairs kept must be exactly those with p <= max_p."""
        from gcen2.core import _gcen_core
        rng = np.random.default_rng(n_samples)
        x = rng.normal(size=(120, n_samples))
        x[60:] = np.maximum(x[60:], 0.0)   # zeros for the dropout methods
        x[1:40] += 2.0 * x[0]
        idx, r, p = _gcen_core.correlation.pairwise_correlation(x, method, 2)
        for max_p in [1e-12, 1e-6, 1e-3, 0.0123, 0.5, 0.999999, 1.0, 0.0]:
            want = (p <= max_p)
            for cut in (0.0, 0.05):
                got = _gcen_core.correlation.pairwise_correlation(
                    x, method, 2, min_abs_cor=cut, max_p=max_p)
                w = want & (np.abs(r) >= cut)
                assert np.array_equal(got[0], idx[w]), (max_p, cut)
                assert np.array_equal(got[2], p[w])
            # pairs just at the cutoff: max_p equal to an observed p-value
            for q in p[(p > 0) & (p < 1)][:20]:
                got = _gcen_core.correlation.pairwise_correlation(x, method, 2, max_p=q)
                assert np.array_equal(got[0], idx[p <= q])

    def test_unwritable_path_raises_oserror(self, tmp_path):
        from gcen2.network import build_to_file
        with pytest.raises(OSError):
            build_to_file(self._matrix(), tmp_path / "no" / "such.tsv",
                          gene_names=[f"g{i}" for i in range(90)])


def _core():
    from gcen2.core import _gcen_core
    if _gcen_core is None:
        pytest.skip("compiled extension not available")
    return _gcen_core


class TestPcg64:
    """_gcen_core.network.Pcg64 draws Generator.integers(0, high, n) bit for
    bit, from numpy's bit_generator.state."""

    @pytest.mark.parametrize("seed", [0, 42, 2**63 + 5])
    def test_integers_match_numpy(self, seed):
        core = _core()
        for high in [1, 2, 3, 10, 1000, 2**16 + 1, 3_609_722, 2**31, 2**32 - 1, 2**32,
                     2**32 + 1, 2**40 + 3, 2**63 - 1]:
            g = np.random.default_rng(seed)
            g.integers(0, 7, 1)  # leaves a buffered uint32 (has_uint32 = 1)
            c = core.network.Pcg64(g.bit_generator.state)
            for size in (1, 5, 333):
                np.testing.assert_array_equal(c.integers(high, size), g.integers(0, high, size))
            assert c.state == g.bit_generator.state


class TestShuffleStream:
    """EdgeSwapper.shuffle (numpy's stream reproduced in C++) applies the
    same swaps as run() fed by numpy, per batch of 2**20, for 32-bit and
    64-bit edge keys."""

    @pytest.mark.parametrize("n_nodes", [3000, 70_000])  # keys: 32-bit, 64-bit
    def test_same_as_numpy_batches(self, n_nodes):
        core = _core()
        rng = np.random.default_rng(5)
        n = 3000
        a = rng.integers(0, n_nodes, n)
        b = (a + 1 + rng.integers(0, n_nodes - 1, n)) % n_nodes  # no self-loops
        n_swaps = 2 * (1 << 20) + 12345                     # three batches
        fast = core.network.EdgeSwapper(a, b, n_nodes)
        applied, state = fast.shuffle(np.random.default_rng(9).bit_generator.state, n_swaps)
        ref = core.network.EdgeSwapper(a, b, n_nodes)
        g = np.random.default_rng(9)
        done = ref_applied = 0
        while done < n_swaps:
            m = min(1 << 20, n_swaps - done)
            ref_applied += ref.run(g.integers(0, n, m), g.integers(0, n, m))
            done += m
        np.testing.assert_array_equal(fast.b(), ref.b())
        assert applied == ref_applied > 0
        assert state == g.bit_generator.state


class TestShuffleFile:
    """shuffle_file / extract_file (read and written by the extension) write
    the same file as shuffle(read_network(path)).to_edgelist(output) /
    extract(read_network(path), genes).to_edgelist(output)."""

    _GENES = ["g01", "g05", "g17", "101", "a b", "nope"]

    def _ref(self, path, out, sf, seed):
        from gcen2.io import read_network
        net_module.shuffle(read_network(path), swap_factor=sf, seed=seed).to_edgelist(out)

    def _rows(self, seed, n_rows, fmt):
        rng = np.random.default_rng(seed)
        genes = [f"g{i:02d}" for i in range(40)] + ["101", "7", "a b"]
        out = []
        for k in range(n_rows):
            a, b = rng.choice(len(genes), 2, replace=False)
            out.append(fmt(k, genes[a], genes[b], rng))
        return out

    def _check(self, tmp_path, text, sf=20, seed=3, expect_cpp=True, genes=None):
        core = _core()
        from gcen2._matrix_io import network_file_columns
        from gcen2.io import read_network

        p = tmp_path / "in.network"
        p.write_bytes(text.encode())
        genes = self._GENES if genes is None else genes
        if expect_cpp is not None:  # which path the file functions will take
            table = net_module._edge_table_args(*network_file_columns(p))
            r = core.network.shuffle_edge_file(
                str(p), *table, str(tmp_path / "probe"), 1,
                np.random.default_rng(0).bit_generator.state)
            assert r["ok"] is expect_cpp, r["reason"]
            r = core.network.extract_edge_file(str(p), *table, str(tmp_path / "probe"), genes)
            assert r["ok"] is expect_cpp, r["reason"]
        got = net_module.shuffle_file(p, tmp_path / "fast.tsv", swap_factor=sf, seed=seed)
        self._ref(p, tmp_path / "ref.tsv", sf, seed)
        assert (tmp_path / "fast.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()
        ref = read_network(tmp_path / "ref.tsv")
        assert got == {"node_count": ref.node_count, "edge_count": ref.edge_count}
        # extract: same file as the pandas pipeline on the whole file
        got = net_module.extract_file(p, genes, tmp_path / "fast_e.tsv")
        sub = net_module.extract(read_network(p), genes)
        sub.to_edgelist(tmp_path / "ref_e.tsv")
        assert (tmp_path / "fast_e.tsv").read_bytes() == (tmp_path / "ref_e.tsv").read_bytes()
        assert got == {"node_count": sub.node_count, "edge_count": sub.edge_count}

    def test_extract_dtypes_from_every_row(self, tmp_path):
        # "count" has a missing value only in a row that is not extracted:
        # read_network reads the whole file, so the column is float (1.0)
        rows = "".join(f"g{i:02d}\tg{i + 1:02d}\t{i}\tx\n" for i in range(30))
        text = "gene_a\tgene_b\tcount\tnote\n" + rows + "g90\tg91\t\tx\n"
        self._check(tmp_path, text, genes=["g03"])
        assert (tmp_path / "fast_e.tsv").read_text(encoding="utf-8").splitlines()[1:] == [
            "g02\tg03\t2.0\tx", "g03\tg04\t3.0\tx"]

    def test_extract_missing_gene_in_other_rows(self, tmp_path):
        # a missing gene name outside the extracted rows does not matter
        core = _core()
        p = tmp_path / "x.tsv"
        p.write_text("gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t0.5\t0.1\n"
                     "C\tNA\t0.4\t0.1\nB\tC\t0.3\t0.1\n", encoding="utf-8")
        r = core.network.extract_edge_file(str(p), True, 4, True,
                                           "gene_a\tgene_b\tcorrelation\tp_value\n",
                                           str(tmp_path / "o.tsv"), ["A"])
        assert r["ok"] and r["edge_count"] == 1
        from gcen2.io import read_network
        net_module.extract(read_network(p), ["A"]).to_edgelist(tmp_path / "ref.tsv")
        assert (tmp_path / "o.tsv").read_bytes() == (tmp_path / "ref.tsv").read_bytes()

    def test_v1_header_missing_fdr_na_quotes_crlf(self, tmp_path):
        def fmt(k, a, b, rng):
            a = f'"{a}"' if k % 7 == 0 else a
            c = "NA" if k % 13 == 0 else f"{rng.uniform(-1, 1):.6f}"
            fdr = "" if k % 5 == 0 else f"\t{rng.uniform():.4E}"
            return f"{a}\t{b}\t{c}\t{rng.uniform():.5e}{fdr}\r\n"
        text = "#Node1\tNode2\tCorrelation\tPvalue\tFDR\r\n" + "".join(self._rows(1, 300, fmt))
        self._check(tmp_path, text + "\r\n")

    def test_gcen2_header_four_columns(self, tmp_path):
        def fmt(k, a, b, rng):
            return f"{a}\t{b}\t{rng.uniform(-1, 1)!r}\t{rng.uniform()!r}\n"
        self._check(tmp_path, "gene_a\tgene_b\tcorrelation\tp_value\n" + "".join(self._rows(2, 300, fmt)))

    def test_headerless_six_columns(self, tmp_path):
        # GCEN v1 network_build output: not the standard layout, written by to_csv
        def fmt(k, a, b, rng):
            r = rng.uniform(-1, 1)
            return f"{a}\t{b}\t{r:.6f}\t{rng.uniform():.5E}\t{rng.uniform():.5E}\t{1 if r > 0 else -1}\n"
        self._check(tmp_path, "".join(self._rows(3, 300, fmt)), sf=50, seed=11)

    def test_custom_columns_dtypes(self, tmp_path):
        # int column, int column with NA (-> float), float, all-NA, strings
        # (with a quote: written quoted), numbers inside a string column
        def fmt(k, a, b, rng):
            n_na = "" if k % 9 == 0 else str(int(rng.integers(-5, 50)))
            label = ['x', 'q"t', '0.5', 'NA', '007'][k % 5]
            return f"{a}\t{b}\t{int(rng.integers(0, 9))}\t{n_na}\t{rng.uniform():.3f}\tNA\t{label}\n"
        text = "gene_a\tgene_b\tcount\tmaybe\tweight\tempty\tlabel\n" + "".join(self._rows(4, 300, fmt))
        self._check(tmp_path, text)

    def test_short_rows_empty_and_single_edge(self, tmp_path):
        rows = "gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t0.5\n" + "".join(
            f"g{i}\tg{i + 1}\t0.{i + 1}\t1e-3\n" for i in range(30))
        self._check(tmp_path, rows)
        self._check(tmp_path, "gene_a\tgene_b\tcorrelation\tp_value\n")
        self._check(tmp_path, "gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t0.5\t0.1\n")

    @pytest.mark.parametrize("bad", ["True", " 0.5", "-inf", "1e", "0.12345678901234567891"])
    def test_uncertain_values_fall_back(self, tmp_path, bad):
        body = "".join(f"g{i}\tg{(i * 7 + 3) % 40}\t0.{i + 1}\tx{i}\n" for i in range(40))
        text = "gene_a\tgene_b\tcorrelation\tnote\n" + body + f"g1\tg9\t{bad}\tx\n"
        self._check(tmp_path, text, expect_cpp=False)

    def test_extra_fields_fall_back(self, tmp_path):
        core = _core()
        p = tmp_path / "x.tsv"
        p.write_text("gene_a\tgene_b\tcorrelation\nA\tB\t0.5\nB\tC\t0.4\textra\n", encoding="utf-8")
        r = core.network.shuffle_edge_file(str(p), True, 3, False, "h\n", str(tmp_path / "o"),
                                           1, np.random.default_rng(0).bit_generator.state)
        assert not r["ok"] and "more fields" in r["reason"]
        assert not (tmp_path / "o").exists()  # nothing written when not handled
