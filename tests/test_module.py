"""Tests for gene module identification (Jaccard + DFS).

The pure-Python integration path and the C++ extension path are both tested.
C++ tests are skipped when the extension is not compiled.
"""

import pandas as pd
import pytest

from gcen2 import module as mod_module
from gcen2.exceptions import GCENInputError
from gcen2.models import ModuleSet, Network

try:
    from gcen2.core import _gcen_core  # type: ignore[import]
    HAS_CORE = True
except ImportError:
    HAS_CORE = False

skip_no_core = pytest.mark.skipif(not HAS_CORE, reason="_gcen_core not compiled")


def _make_clique_network(genes: list[str]) -> Network:
    """Create a fully connected clique network among ``genes``."""
    edges = []
    for i, a in enumerate(genes):
        for b in genes[i + 1:]:
            edges.append({"gene_a": a, "gene_b": b, "correlation": 0.9, "p_value": 0.001})
    return Network(edges=pd.DataFrame(edges))


@skip_no_core
class TestIdentifyModules:
    def test_single_clique_forms_one_module(self):
        genes = ["G1", "G2", "G3", "G4"]
        net = _make_clique_network(genes)
        result = mod_module.identify(net, similarity=0.3, threads=1)
        assert isinstance(result, ModuleSet)
        # All 4 genes should form a single module
        assert result.module_count == 1
        assert set(result.modules[0]) == set(genes)

    def test_two_cliques_form_two_modules(self):
        clique_a = ["A1", "A2", "A3"]
        clique_b = ["B1", "B2", "B3"]
        net_a = _make_clique_network(clique_a)
        net_b = _make_clique_network(clique_b)
        combined = pd.concat([net_a.edges, net_b.edges], ignore_index=True)
        net = Network(edges=combined)
        result = mod_module.identify(net, similarity=0.3, threads=1)
        assert result.module_count == 2

    def test_singletons_excluded(self):
        # Genes not connected via Jaccard > cutoff should not appear in any module.
        # G1 is isolated; G2-G3 share neighbor G4 giving Jaccard=1 via the formula.
        edges = pd.DataFrame({
            "gene_a": ["G2", "G3", "G2"],
            "gene_b": ["G3", "G4", "G4"],
            "correlation": [0.9, 0.9, 0.9],
            "p_value": [0.001, 0.001, 0.001],
        })
        net = Network(edges=edges)
        result = mod_module.identify(net, similarity=0.5, threads=1)
        # G2, G3, G4 form a module; isolated G1 should not appear
        all_module_genes = {g for m in result.modules for g in m}
        assert "G1" not in all_module_genes

    def test_multithreaded_matches_single_threaded(self):
        genes = [f"G{i}" for i in range(10)]
        net = _make_clique_network(genes)
        r1 = mod_module.identify(net, similarity=0.3, threads=1)
        r2 = mod_module.identify(net, similarity=0.3, threads=2)
        # Same total gene coverage
        genes_1 = {g for m in r1.modules for g in m}
        genes_2 = {g for m in r2.modules for g in m}
        assert genes_1 == genes_2

    def test_empty_network_returns_empty_moduleset(self):
        net = Network(edges=pd.DataFrame(
            columns=["gene_a", "gene_b", "correlation", "p_value"]
        ))
        result = mod_module.identify(net, similarity=0.5)
        assert result.module_count == 0

    def test_invalid_similarity_raises(self):
        net = _make_clique_network(["A", "B", "C"])
        with pytest.raises(GCENInputError):
            mod_module.identify(net, similarity=1.5)


def _csr(adj_pairs, n):
    import numpy as np
    from scipy.sparse import coo_matrix
    pairs = np.array(sorted({(a, b) for a, b in adj_pairs} | {(b, a) for a, b in adj_pairs}))
    A = coo_matrix((np.ones(len(pairs)), (pairs[:, 0], pairs[:, 1])), shape=(n, n)).tocsr()
    A.sort_indices()
    return A.indptr.astype(np.int64), A.indices.astype(np.int32)


def _bruteforce_modules(net, cutoff):
    """All-pairs GCEN v1 algorithm (module_identify.cpp), pure Python."""
    adj = net.neighbor_sets()
    genes = sorted(adj)
    links = {g: set() for g in genes}
    for i, a in enumerate(genes):
        for b in genes[i + 1:]:
            shared = len(adj[a] & adj[b])
            ab = 1 if b in adj[a] else 0
            denom = len(adj[a]) + len(adj[b]) - shared - ab
            if denom > 0 and (shared + ab) / denom > cutoff:
                links[a].add(b)
                links[b].add(a)
    seen, modules = set(), []
    for g in genes:
        if g in seen or not links[g]:
            continue
        stack, comp = [g], set()
        while stack:
            x = stack.pop()
            if x in comp:
                continue
            comp.add(x)
            stack.extend(links[x] - comp)
        seen |= comp
        modules.append(frozenset(comp))
    return set(modules)


def _random_network(seed):
    import numpy as np
    import pandas as pd

    from gcen2.models import Network
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(60)]
    # dense communities + sparse noise + a few self-loops and duplicate edges
    rows = []
    for c in range(4):
        members = genes[c * 12:(c + 1) * 12]
        for a in range(12):
            for b in range(a + 1, 12):
                if rng.uniform() < 0.7:
                    rows.append((members[a], members[b]))
    rows += [(genes[a], genes[b]) for a, b in rng.integers(0, 60, size=(40, 2)) if a != b]
    rows += [(genes[i], genes[i]) for i in (3, 20)]
    rows += [(b, a) for a, b in rows[:5]]
    return Network(edges=pd.DataFrame(
        [(a, b, 0.9) for a, b in rows], columns=["gene_a", "gene_b", "correlation"]))


@skip_no_core
class TestCppModuleKernel:
    """Direct tests of the sparse _gcen_core.module kernel."""

    def test_triangle_is_one_module(self):
        indptr, indices = _csr([(0, 1), (1, 2), (0, 2)], 3)
        labels = _gcen_core.module.jaccard_modules(indptr, indices, 0.3, 1)
        assert list(labels) == [0, 0, 0]

    def test_unlinked_nodes_get_minus_one(self):
        # star 0-1, 0-2: J(1,2) = 1/(1+1-1) = 1 ; J(0,1) = 1/(2+1-0-1) = 0.5
        indptr, indices = _csr([(0, 1), (0, 2)], 4)  # node 3 isolated
        labels = _gcen_core.module.jaccard_modules(indptr, indices, 0.6, 1)
        assert list(labels) == [-1, 0, 0, -1]

    def test_rejects_bad_csr(self):
        import numpy as np
        with pytest.raises(ValueError):
            _gcen_core.module.jaccard_modules(np.array([0, 1], dtype=np.int64),
                                              np.array([5], dtype=np.int32), 0.5, 1)


class TestModulesMatchV1Algorithm:
    """The C++ sparse kernel and the scipy fallback must both reproduce the
    all-pairs GCEN v1 algorithm exactly (also checked against a compiled
    v1 0.6.3 on the 3.8M-edge sample network: 1708/1708 modules identical)."""

    @pytest.mark.parametrize("seed", [0, 1, 2])
    @pytest.mark.parametrize("cutoff", [0.2, 0.5, 0.8])
    def test_same_modules(self, seed, cutoff):
        from gcen2.module import _jaccard_modules, identify
        net = _random_network(seed)
        expected = _bruteforce_modules(net, cutoff)
        assert len(expected) > 0
        scipy_mods = {frozenset(m) for m in _jaccard_modules(net, cutoff)}
        assert scipy_mods == expected
        if _gcen_core is not None:
            for threads in (1, 4):
                got = identify(net, similarity=cutoff, threads=threads)
                assert {frozenset(m) for m in got.modules} == expected
                assert len(got.modules) == len(expected)

    def test_cpp_and_scipy_same_order(self, small_network):
        if _gcen_core is None:
            pytest.skip("_gcen_core not compiled")
        from gcen2.module import _jaccard_modules, _jaccard_modules_cpp
        assert _jaccard_modules_cpp(small_network, 0.3, 2, _gcen_core) == \
            _jaccard_modules(small_network, 0.3)

    def test_chunking_does_not_change_result(self, small_network):
        from gcen2.module import _jaccard_modules
        assert _jaccard_modules(small_network, 0.3, chunk_size=1) == \
            _jaccard_modules(small_network, 0.3)


class TestAdjacencyInputs:
    def test_categorical_and_duplicate_edges(self):
        import numpy as np

        from gcen2.models import Network
        from gcen2.module import identify
        rng = np.random.default_rng(7)
        genes = [f"g{i}" for i in range(60)]
        pairs = rng.integers(0, 60, size=(400, 2))
        df = pd.DataFrame({"gene_a": [genes[i] for i in pairs[:, 0]],
                           "gene_b": [genes[j] for j in pairs[:, 1]]})
        # duplicates, reversed duplicates and a self-loop must not change anything
        noisy = pd.concat([df, df.iloc[:50], df.iloc[:50][["gene_b", "gene_a"]].set_axis(
            ["gene_a", "gene_b"], axis=1), pd.DataFrame({"gene_a": ["g1"], "gene_b": ["g1"]})],
            ignore_index=True)
        clean = pd.concat([df, pd.DataFrame({"gene_a": ["g1"], "gene_b": ["g1"]})],
                          ignore_index=True)
        ref = identify(Network(edges=clean), similarity=0.3)
        for edges in (noisy, noisy.astype("category")):
            got = identify(Network(edges=edges), similarity=0.3)
            assert got.modules == ref.modules
        assert ref.module_count > 0


@skip_no_core
class TestModulesFromFile:
    """modules_from_file (C++ reader, used by ``gcen2 module``) must match
    identify(read_network(path))."""

    @pytest.fixture
    def edges(self):
        import numpy as np
        rng = np.random.default_rng(11)
        genes = [f"G{i:03d}" for i in range(80)] + ["Ä-gene", "b", "B", "a1"]
        pairs = rng.integers(0, len(genes), size=(600, 2))
        a = [genes[i] for i in pairs[:, 0]]
        b = [genes[j] for j in pairs[:, 1]]
        # reversed duplicates and a self-loop
        a += b[:40] + ["G001"]
        b += a[:40] + ["G001"]
        return a, b

    @pytest.mark.parametrize("header", [
        "gene_a\tgene_b\tcorrelation\tp_value\n",
        "#node1\tnode2\tcorrelation\tp-value\n",
        "",  # header-less v1
    ])
    @pytest.mark.parametrize("eol", ["\n", "\r\n"])
    def test_same_as_identify(self, tmp_path, edges, header, eol):
        from gcen2.io import read_network
        from gcen2.module import identify, modules_from_file
        a, b = edges
        text = header.replace("\n", eol) + "".join(
            f"{x}\t{y}\t0.95\t1e-5{eol}" for x, y in zip(a, b))
        p = tmp_path / "net.tsv"
        p.write_bytes(text.encode())
        ref = identify(read_network(p, genes_only=True), similarity=0.3, threads=2).modules
        got = modules_from_file(p, similarity=0.3, threads=3)
        assert len(ref) > 0
        assert got == ref

    def test_two_columns_blank_lines_and_quotes(self, tmp_path):
        from gcen2.io import read_network
        from gcen2.module import identify, modules_from_file
        p = tmp_path / "net.tsv"
        p.write_text('gene_a\tgene_b\n"A"\tB\n\nA\tC\nB\tC\nC\t"D""x"\nB\tD"x\n', encoding="utf-8")
        mods = modules_from_file(p, similarity=0.3)
        assert mods == identify(read_network(p, genes_only=True), similarity=0.3).modules
        assert {"A", "B", "C", 'D"x'} <= set().union(*mods)

    def test_empty_and_header_only(self, tmp_path):
        from gcen2.module import modules_from_file
        p = tmp_path / "net.tsv"
        p.write_text("gene_a\tgene_b\tcorrelation\tp_value\n", encoding="utf-8")
        assert modules_from_file(p) == []
        p.write_text("", encoding="utf-8")
        assert modules_from_file(p) == []

    def test_errors(self, tmp_path):
        from gcen2.exceptions import GCENIOError
        from gcen2.module import modules_from_file
        p = tmp_path / "net.tsv"
        p.write_text("gene_a\tgene_b\nA\tB\nlonely\n", encoding="utf-8")
        with pytest.raises(GCENIOError, match="line 3"):
            modules_from_file(p)
        with pytest.raises(GCENIOError):
            modules_from_file(tmp_path / "missing.tsv")
        with pytest.raises(GCENInputError):
            modules_from_file(p, similarity=1.5)

    def test_no_pandas_import(self):
        import subprocess
        import sys
        code = ("import sys; from gcen2.module import modules_from_file; import gcen2.cli; "
                "print(sorted(m for m in ('pandas', 'scipy') if m in sys.modules))")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8",
                             check=True)
        assert out.stdout.strip() == "[]"


class TestReaderAndKernelMemoryPaths:
    """The reader's uint16 -> int32 node-id switch and the per-thread
    union-find of jaccard_modules give the same results as before."""

    def test_reader_more_than_65536_nodes(self, tmp_path):
        import numpy as np

        from gcen2.core import _HAVE_CORE
        from gcen2.io import read_network
        from gcen2.module import _unique_adjacency, read_file_adjacency
        if not _HAVE_CORE:
            pytest.skip("compiled extension not available")
        rng = np.random.default_rng(2)
        n = 70_000  # node ids pass 65535 midway through the file
        a = np.arange(0, n, 2)
        b = a + 1
        # links back to early nodes after the switch, plus duplicates
        extra = rng.integers(0, n, size=(3000, 2))
        a = np.r_[a, extra[:, 0], b[:50]]
        b = np.r_[b, extra[:, 1], a[:50]]
        p = tmp_path / "big.tsv"
        p.write_text("".join(f"n{x}\tn{y}\t0.9\t1e-5\n" for x, y in zip(a, b)), encoding="utf-8")
        indptr, indices, names, counts = read_file_adjacency(p, with_counts=True)
        A, genes = _unique_adjacency(read_network(p, genes_only=True))
        assert len(names) == n and names == list(genes)
        np.testing.assert_array_equal(indptr, A.indptr)
        np.testing.assert_array_equal(indices, A.indices)
        assert counts.sum() == 2 * len(a)

    def test_jaccard_threads_agree(self):
        import numpy as np

        from gcen2.core import _HAVE_CORE
        if not _HAVE_CORE:
            pytest.skip("compiled extension not available")
        from gcen2.core import _gcen_core as core
        # 12 dense groups of 30 nodes with sparse noise between them
        rng = np.random.default_rng(4)
        grp = np.repeat(np.arange(12), 30)
        P = np.where(grp[:, None] == grp[None, :], 0.5, 0.01)
        M = np.triu(rng.random((360, 360)) < P, 1)
        from scipy.sparse import csr_matrix
        A = csr_matrix(M | M.T)
        A.sort_indices()
        ip, ix = A.indptr.astype(np.int64), A.indices.astype(np.int32)
        for cutoff in (0.2, 0.3):
            ref = np.asarray(core.module.jaccard_modules(ip, ix, cutoff, 1))
            assert ref.max() >= 5  # several modules
            for threads in (2, 4, 7):
                np.testing.assert_array_equal(
                    np.asarray(core.module.jaccard_modules(ip, ix, cutoff, threads)), ref)
