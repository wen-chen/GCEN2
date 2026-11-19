"""Tests for gcen2.rwr — Random Walk with Restart (Python fallback)."""

import numpy as np
import pandas as pd
import pytest

from gcen2.exceptions import GCENInputError
from gcen2.models import Network
from gcen2.rwr import _rwr_python


class TestRWRPython:
    def _column_normalize(self, W):
        col_sums = W.sum(axis=0)
        col_sums[col_sums == 0] = 1.0
        return W / col_sums

    def test_converges_simple(self):
        W = self._column_normalize(np.array([
            [0, 1, 0],
            [1, 0, 1],
            [0, 1, 0],
        ], dtype=float))
        P0 = np.array([1.0, 0.0, 0.0])
        PT, steps = _rwr_python(W, P0, gamma=0.5, max_iter=10000, tol=1e-10)
        assert steps < 10000
        assert abs(PT.sum() - 1.0) < 1e-9

    def test_probability_sums_to_one(self):
        W = self._column_normalize(np.eye(4, dtype=float))
        P0 = np.array([0.25, 0.25, 0.25, 0.25])
        PT, _ = _rwr_python(W, P0, gamma=0.5, max_iter=1000, tol=1e-10)
        assert abs(PT.sum() - 1.0) < 1e-9

    def test_seed_genes_get_higher_probability(self):
        W = self._column_normalize(np.ones((5, 5), dtype=float))
        P0 = np.zeros(5)
        P0[0] = 1.0
        PT, _ = _rwr_python(W, P0, gamma=0.9, max_iter=1000, tol=1e-10)
        assert PT[0] > PT[1]


class TestRWRPredict:
    def test_missing_seed_raises(self, small_network):
        from gcen2 import rwr
        with pytest.raises(GCENInputError):
            rwr.predict(small_network, seed_genes=["NONEXISTENT"])

    def test_empty_seed_raises(self, small_network):
        from gcen2 import rwr
        with pytest.raises(GCENInputError):
            rwr.predict(small_network, seed_genes=[])

    def test_invalid_gamma_raises(self, small_network):
        from gcen2 import rwr
        with pytest.raises(GCENInputError):
            rwr.predict(small_network, seed_genes=["GeneA"], gamma=1.5)


class TestRWRCoreMatchesPython:
    """Regression: the sparse C++ RWR must agree with the Python reference.

    An earlier dense binding silently iterated with W.T instead of W; that only
    shows up when node degrees differ (W.T is then not column-stochastic).
    """

    def test_irregular_degree_graph(self):
        core = pytest.importorskip("gcen2.core._gcen_core")
        A = np.array([
            [0, 1, 1, 1, 0],
            [1, 0, 0, 0, 0],
            [1, 0, 0, 1, 0],
            [1, 0, 1, 0, 1],
            [0, 0, 0, 1, 0],
        ], dtype=float)
        W = A / A.sum(axis=0)
        P0 = np.array([0.0, 1.0, 0.0, 0.0, 0.0])
        PT_py, _ = _rwr_python(W, P0, gamma=0.3, max_iter=10000, tol=1e-12)
        from scipy.sparse import csr_matrix
        Ws = csr_matrix(W)
        PT_cpp, steps_cpp = core.rwr.random_walk_restart(
            Ws.indptr.astype(np.int64), Ws.indices.astype(np.int32), Ws.data, P0,
            0.3, 10000, 1e-12)
        _, steps_py = _rwr_python(W, P0, gamma=0.3, max_iter=10000, tol=1e-12)
        np.testing.assert_allclose(PT_cpp, PT_py, atol=1e-12)
        assert steps_cpp == steps_py
        assert np.isclose(PT_cpp.sum(), 1.0)

    def test_python_fallback_matches(self, monkeypatch):
        """predict() without the extension (scipy path) gives the same ranking."""
        import gcen2.core
        from gcen2.rwr import predict
        rng = np.random.default_rng(9)
        genes = [f"g{i}" for i in range(40)]
        rows = [(genes[a], genes[b], 0.8) for a, b in rng.integers(0, 40, size=(120, 2)) if a != b]
        net = Network(edges=pd.DataFrame(rows, columns=["gene_a", "gene_b", "correlation"]))
        with_cpp = predict(net, ["g0", "g1"], prob_cutoff=0.0, tol=1e-13).ranked_genes
        monkeypatch.setattr(gcen2.core, "_gcen_core", None)
        without = predict(net, ["g0", "g1"], prob_cutoff=0.0, tol=1e-13).ranked_genes
        pd.testing.assert_frame_equal(with_cpp, without, check_exact=False, atol=1e-12)


def _dense_reference_predict(network, seeds, gamma, directed, weighted_network):
    """The previous dense, iterrows-based implementation (kept as a reference)."""
    all_genes = sorted(set(network.edges["gene_a"]) | set(network.edges["gene_b"]))
    idx = {g: i for i, g in enumerate(all_genes)}
    n = len(all_genes)
    W = np.zeros((n, n))
    for _, row in network.edges.iterrows():
        i, j = idx[row["gene_a"]], idx[row["gene_b"]]
        w = abs(row["correlation"]) if weighted_network else 1.0
        W[i, j] += w
        if not directed:
            W[j, i] += w
    cs = W.sum(axis=0)
    cs[cs == 0] = 1.0
    W = W / cs
    P0 = np.zeros(n)
    for g in seeds:
        P0[idx[g]] = 1.0
    P0 /= P0.sum()
    PT, _ = _rwr_python(W, P0, gamma, 1_000_000, 1e-12)
    return dict(zip(all_genes, PT))


class TestSparseRWRMatchesDense:
    @pytest.fixture
    def random_network(self):
        rng = np.random.default_rng(5)
        genes = [f"g{i}" for i in range(60)]
        pairs = {tuple(sorted(rng.choice(60, 2, replace=False))) for _ in range(200)}
        rows = [(genes[a], genes[b], rng.uniform(-1, 1)) for a, b in pairs]
        rows.append((genes[0], genes[1], 0.5))  # duplicate edge -> weights add up
        return Network(edges=pd.DataFrame(rows, columns=["gene_a", "gene_b", "correlation"]))

    @pytest.mark.parametrize("directed", [False, True])
    @pytest.mark.parametrize("weighted", [False, True])
    def test_probabilities_match(self, random_network, directed, weighted):
        from gcen2.rwr import predict
        seeds = ["g0", "g7"]
        ref = _dense_reference_predict(random_network, seeds, 0.4, directed, weighted)
        res = predict(random_network, seeds, gamma=0.4, prob_cutoff=0.0, directed=directed,
                      weighted_network=weighted, tol=1e-12)
        got = dict(zip(res.ranked_genes["gene"], res.ranked_genes["probability"]))
        assert set(got) == set(ref) - set(seeds)
        for g, p in got.items():
            assert p == pytest.approx(ref[g], abs=1e-12)


class TestNeighborSets:
    def test_matches_iterrows(self):
        rng = np.random.default_rng(3)
        genes = [f"g{i}" for i in range(30)]
        rows = [(genes[a], genes[b], 0.5) for a, b in rng.integers(0, 30, size=(120, 2))]
        net = Network(edges=pd.DataFrame(rows, columns=["gene_a", "gene_b", "correlation"]))
        ref = {}
        for _, r in net.edges.iterrows():
            ref.setdefault(r["gene_a"], set()).add(r["gene_b"])
            ref.setdefault(r["gene_b"], set()).add(r["gene_a"])
        assert net.neighbor_sets() == ref


class TestPredictFromFile:
    """predict_from_file (C++ reader, used by ``gcen2 rwr``) must give
    bit-identical probabilities, in the same order, as predict()."""

    @pytest.fixture
    def net_file(self, tmp_path):
        rng = np.random.default_rng(5)
        genes = [f"G{i:02d}" for i in range(40)]
        pairs = rng.integers(0, 40, size=(160, 2))
        rows = [(genes[i], genes[j]) for i, j in pairs]
        rows += rows[:15] + [(b, a) for a, b in rows[15:30]] + [("G03", "G03")]
        # two symmetric leaves -> tied probabilities
        rows += [("G00", "L1"), ("G00", "L2")]
        p = tmp_path / "net.tsv"
        p.write_text("gene_a\tgene_b\tcorrelation\tp_value\n"
                     + "".join(f"{a}\t{b}\t0.9\t1e-5\n" for a, b in rows), encoding="utf-8")
        return p

    def test_counts_kernel_matches_data_kernel(self, net_file):
        from gcen2.core import _HAVE_CORE
        if not _HAVE_CORE:
            pytest.skip("compiled extension not available")
        import gc

        from gcen2.core import _gcen_core as core
        from gcen2.module import read_file_adjacency
        indptr, indices, names, counts = read_file_adjacency(net_file, with_counts=True,
                                                             descending=True)
        gc.collect()  # arrays own their (zero-copy) buffers
        assert counts.max() > 1 and indices.flags.owndata is False
        n = len(names)
        col = np.bincount(indices, weights=counts, minlength=n)
        col[col == 0] = 1.0
        data = (1.0 / col)[indices] * counts
        P0 = np.zeros(n)
        P0[[0, 5]] = 0.5
        a, sa = core.rwr.random_walk_restart_counts(indptr, indices, counts, P0, 0.3)
        b, sb = core.rwr.random_walk_restart(indptr, indices, data, P0, 0.3)
        assert sa == sb
        np.testing.assert_array_equal(a, b)  # bit-identical
        with pytest.raises(ValueError):
            core.rwr.random_walk_restart_counts(indptr, indices, counts[:-1], P0, 0.3)
        bad = counts.copy()
        bad[0] = -1
        with pytest.raises(ValueError):
            core.rwr.random_walk_restart_counts(indptr, indices, bad, P0, 0.3)

    @pytest.mark.parametrize("cutoff", [0.0, 0.01])
    def test_same_as_predict(self, net_file, cutoff):
        from gcen2.core import _HAVE_CORE
        from gcen2.io import read_network
        from gcen2.rwr import predict, predict_from_file
        if not _HAVE_CORE:
            pytest.skip("compiled extension not available")
        seeds = ["G00", "G07", "missing"]
        ref = predict(read_network(net_file), seeds, gamma=0.3, prob_cutoff=cutoff).ranked_genes
        genes, prob = predict_from_file(net_file, seeds, gamma=0.3, prob_cutoff=cutoff)
        assert len(genes) > 3
        assert genes == ref["gene"].tolist()
        np.testing.assert_array_equal(prob, ref["probability"].to_numpy())

    def test_rank_desc_ties_by_gene_name(self):
        """Ties ordered by gene name (not by how a sort algorithm happens to
        leave them): the same ranking with every numpy version."""
        from gcen2.rwr import _rank_desc
        rng = np.random.default_rng(0)
        v = rng.integers(0, 5, 200).astype(float) / 7
        genes = [f"g{i:03d}" for i in rng.permutation(200)]
        want = (pd.DataFrame({"p": v, "g": genes})
                .sort_values(["p", "g"], ascending=[False, True]).index.to_numpy())
        np.testing.assert_array_equal(_rank_desc(v, genes), want)

    def test_errors(self, net_file):
        from gcen2.rwr import predict_from_file
        with pytest.raises(GCENInputError):
            predict_from_file(net_file, ["nope"])
        with pytest.raises(GCENInputError):
            predict_from_file(net_file, [])
        with pytest.raises(GCENInputError):
            predict_from_file(net_file, ["G00"], gamma=1.5)
