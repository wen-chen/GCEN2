"""Tests for gcen2.annotation — GO parsing, enrichment analysis."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

import gcen2.core
from gcen2.annotation import (
    enrich,
    load_association,
    load_go,
)
from gcen2.exceptions import GCENIOError
from tests.reference import hypergeometric_pvalue

SAMPLE_OBO = """\
format-version: 1.2
[Term]
id: GO:0001
name: test_process
namespace: biological_process
[Term]
id: GO:0002
name: child_process
namespace: biological_process
is_a: GO:0001 ! test_process
[Term]
id: GO:0003
name: obsolete_term
namespace: biological_process
is_obsolete: true
"""

SAMPLE_ASSOC = """\
GeneA\tGO:0002
GeneB\tGO:0002
GeneC\tGO:0001
GeneD\tGO:0001
"""


@pytest.fixture
def obo_file(tmp_path):
    p = tmp_path / "test.obo"
    p.write_text(SAMPLE_OBO, encoding="utf-8")
    return p


@pytest.fixture
def assoc_file(tmp_path):
    p = tmp_path / "test.assoc"
    p.write_text(SAMPLE_ASSOC, encoding="utf-8")
    return p


class TestLoadGO:
    def test_parses_terms(self, obo_file):
        onto = load_go(obo_file)
        assert "GO:0001" in onto.terms
        assert "GO:0002" in onto.terms

    def test_parent_relationship(self, obo_file):
        onto = load_go(obo_file)
        assert "GO:0001" in onto.terms["GO:0002"].parents

    def test_obsolete_flag(self, obo_file):
        onto = load_go(obo_file)
        assert onto.terms["GO:0003"].is_obsolete

    def test_missing_file_raises(self):
        with pytest.raises(GCENIOError):
            load_go("/nonexistent/path/test.obo")


class TestGOPropagation:
    def test_propagation_adds_ancestors(self, obo_file):
        onto = load_go(obo_file)
        gene_terms = {"GeneX": {"GO:0002"}}
        propagated = onto.propagate(gene_terms)
        assert "GO:0001" in propagated["GeneX"]
        assert "GO:0002" in propagated["GeneX"]


class TestLoadAssociation:
    def test_parses_associations(self, assoc_file):
        assoc = load_association(assoc_file)
        assert "GeneA" in assoc
        assert "GO:0002" in assoc["GeneA"]

    def test_propagation_applied_when_ontology_provided(self, obo_file, assoc_file):
        onto = load_go(obo_file)
        assoc = load_association(assoc_file, ontology=onto)
        # GeneA annotated with GO:0002 should also have GO:0001 after propagation
        assert "GO:0001" in assoc.get("GeneA", set())


class TestHypergeometric:
    def test_pvalue_range(self):
        p = hypergeometric_pvalue(100, 20, 10, 5)
        assert 0.0 <= p <= 1.0

    def test_perfect_enrichment(self):
        # All study genes in the term → low p-value
        p = hypergeometric_pvalue(1000, 20, 20, 20)
        assert p < 0.05

    def test_no_enrichment(self):
        # study count matches background proportion → high p-value
        p = hypergeometric_pvalue(100, 10, 10, 1)
        assert p > 0.05


class TestEnrich:
    def test_enrichment_returns_result(self, obo_file, assoc_file):
        onto = load_go(obo_file)
        assoc = load_association(assoc_file, ontology=onto)
        result = enrich(
            query_genes=["GeneA", "GeneB"],
            background_genes=["GeneA", "GeneB", "GeneC", "GeneD"],
            ontology=onto,
            assoc=assoc,
            pval=1.0,
        )
        assert result.ontology_type == "GO"


def _reference_enrichment(query, background, gene_terms, pval):
    """The original per-query / per-term loop implementation."""
    term_pop = {}
    for g in background:
        for t in gene_terms.get(g, set()):
            term_pop[t] = term_pop.get(t, 0) + 1
    pop_n, study_n = len(background), len(query)
    rows = []
    for term_id, pop_count in term_pop.items():
        sc = sum(1 for g in query if term_id in gene_terms.get(g, set()))
        if sc == 0:
            continue
        rows.append({"term_id": term_id, "study_count": sc, "pop_count": pop_count,
                     "study_n": study_n, "pop_n": pop_n,
                     "enrichment": (sc / study_n) / (pop_count / pop_n),
                     "p_value": hypergeometric_pvalue(pop_n, pop_count, study_n, sc)})
    df = pd.DataFrame(rows)
    if not df.empty:
        df = df[df["p_value"] <= pval]
    return df


class TestBatchEnrichmentMatchesReference:
    @pytest.fixture
    def data(self):
        rng = np.random.default_rng(11)
        genes = [f"g{i}" for i in range(300)]
        terms = [f"T{i}" for i in range(40)]
        # skewed term sizes so that some terms are genuinely enriched
        weights = np.linspace(3, 0.2, len(terms))
        weights /= weights.sum()
        gene_terms = {g: set(rng.choice(terms, size=rng.integers(1, 6), p=weights))
                      for g in genes[:250]}            # g250..g299 unannotated
        queries = {f"q{k}": set(rng.choice(genes, size=rng.integers(1, 60), replace=False))
                   for k in range(25)}
        queries["empty_hits"] = {"g299", "g298"}       # only unannotated genes
        return gene_terms, queries

    @pytest.mark.parametrize("bg_kind", ["annotated", "subset"])
    def test_matches(self, data, bg_kind):
        from gcen2.annotation import _batch_enrichment
        gene_terms, queries = data
        background = set(gene_terms) if bg_kind == "annotated" else set(list(gene_terms)[:150])
        pval = 0.3
        got = _batch_enrichment(queries, background, gene_terms, "GO", pval, chunk_size=7)
        cols = ["term_id", "study_count", "pop_count", "study_n", "pop_n", "enrichment", "p_value"]
        n_nonempty = 0
        for k, q in queries.items():
            ref = _reference_enrichment(q, background, gene_terms, pval)
            res = got[k].results
            assert res.empty == ref.empty, k
            if ref.empty:
                continue
            n_nonempty += 1
            ref = ref.sort_values("term_id").reset_index(drop=True)[cols]
            res = res.sort_values("term_id").reset_index(drop=True)[cols]
            pd.testing.assert_frame_equal(res, ref, check_dtype=False, rtol=1e-12)
            assert got[k].results["p_value"].is_monotonic_increasing
        assert n_nonempty > 5


class TestVectorizedHypergeomSF:
    def test_matches_scipy_random_grid(self):
        from scipy.stats import hypergeom

        from gcen2.annotation import _hypergeom_sf_vec
        rng = np.random.default_rng(0)
        size = 20000
        M = rng.integers(2, 30000, size)
        n = (rng.uniform(size=size) ** 3 * M).astype(int)        # successes, skewed small
        N = (rng.uniform(size=size) ** 2 * M).astype(int)        # draws
        hi = np.minimum(n, N)
        k = (rng.uniform(size=size) * (hi + 2)).astype(int)      # includes 0 and > max
        got = _hypergeom_sf_vec(k, M, n, N)
        ref = hypergeom.sf(k - 1, M, n, N)
        np.testing.assert_allclose(got, ref, rtol=1e-9, atol=1e-300)

    def test_extreme_left_of_mode_is_one(self):
        from scipy.stats import hypergeom

        from gcen2.annotation import _hypergeom_sf_vec
        # mean ~475: P(X >= 1) == 1 even though pmf(1) underflows
        got = _hypergeom_sf_vec(np.array([1, 400, 600]), 21056, 5000, 2000)
        ref = hypergeom.sf(np.array([0, 399, 599]), 21056, 5000, 2000)
        np.testing.assert_allclose(got, ref, rtol=1e-9, atol=1e-300)

    def test_scalar_helper_consistent(self):
        from gcen2.annotation import _hypergeom_sf_vec
        for args in [(100, 20, 10, 5), (1000, 20, 20, 20), (100, 10, 10, 1), (100, 10, 10, 0)]:
            pop_n, pop_k, study_n, study_k = args
            assert _hypergeom_sf_vec([study_k], pop_n, [pop_k], [study_n])[0] == \
                pytest.approx(hypergeometric_pvalue(*args), rel=1e-9)


# ---------------------------------------------------------------------------
# GCEN v1 annotation semantics
# ---------------------------------------------------------------------------

FIXTURE = Path(__file__).parent / "data" / "v1_annotate"


class TestLoadAssociationFormats:
    def test_tab_and_comma_separated_terms(self, tmp_path):
        p = tmp_path / "a.assoc"
        p.write_text("G1\tGO:1\tGO:2\tGO:3\nG2\tGO:4,GO:5\nG1\tGO:6\nG3\n# comment\n", encoding="utf-8")
        assoc = load_association(p)
        assert assoc["G1"] == {"GO:1", "GO:2", "GO:3", "GO:6"}
        assert assoc["G2"] == {"GO:4", "GO:5"}
        assert assoc["G3"] == set()  # listed without terms: still "annotated" (v1)


def _v1_network_queries(edges, annotated):
    """Literal port of GCEN v1 enrich.hpp ``load_network``."""
    network, background = {}, set()
    for a, b in edges:
        if a not in annotated and b in annotated:
            network.setdefault(a, set()).add(b)
            background.add(b)
        elif a in annotated and b not in annotated:
            network.setdefault(b, set()).add(a)
            background.add(a)
    return network, background


class TestAnnotationQueries:
    def test_network_queries_match_v1(self):
        from gcen2.annotation import build_annotation_queries
        from gcen2.models import Network
        rng = np.random.default_rng(4)
        genes = [f"g{i}" for i in range(80)]
        annotated = {g: {"T"} for g in genes[:45]}
        rows = [(genes[a], genes[b]) for a, b in rng.integers(0, 80, size=(400, 2))]
        rows += [("g1", "g1"), ("g70", "g70")]  # self-loops are ignored
        net = Network(edges=pd.DataFrame([(a, b, 0.9) for a, b in rows],
                                         columns=["gene_a", "gene_b", "correlation"]))
        q, bg = build_annotation_queries(net, annotated)
        q_ref, bg_ref = _v1_network_queries(rows, annotated)
        assert q == q_ref and bg == bg_ref
        assert set(q).isdisjoint(annotated)

    def test_module_queries_use_all_module_genes(self):
        from gcen2.annotation import build_annotation_queries
        from gcen2.models import ModuleSet
        ms = ModuleSet(modules=[{"a", "b", "x"}, {"c", "y"}])
        q, bg = build_annotation_queries(ms, {"a": {"T"}, "b": {"T"}, "c": {"U"}})
        assert q == {"module_0": {"a", "b", "x"}, "module_1": {"c", "y"}}
        assert bg == {"a", "b", "c", "x", "y"}


class TestNamesAndOrder:
    def test_name_namespace_and_v1_order(self, obo_file, assoc_file):
        onto = load_go(obo_file)
        assoc = load_association(assoc_file, ontology=onto)
        res = enrich(["GeneA", "GeneB"], ["GeneA", "GeneB", "GeneC", "GeneD"],
                     onto, assoc, pval=1.0).results
        assert dict(zip(res.term_id, res.name)) == {"GO:0001": "test_process",
                                                    "GO:0002": "child_process"}
        assert set(res.namespace) == {"biological_process"}

    def test_enrich_counts_only_annotated_genes(self, obo_file, assoc_file):
        onto = load_go(obo_file)
        assoc = load_association(assoc_file, ontology=onto)
        res = enrich(["GeneA", "NotAnnotated"], ["GeneA", "GeneB", "GeneC", "GeneD", "Other"],
                     onto, assoc, pval=1.0).results
        assert set(res.study_n) == {1} and set(res.pop_n) == {4}


class TestKEGG:
    def test_k2ko_mapping_and_names(self, tmp_path):
        from gcen2.annotation import annotate_kegg, load_kegg
        from gcen2.models import ModuleSet
        k2ko = tmp_path / "K2ko.tsv"
        k2ko.write_text("K1\tko001\tGlycolysis\nK1\tko002\tTCA cycle\nK2\tko002\tTCA cycle\n", encoding="utf-8")
        assoc = tmp_path / "kegg.assoc"
        assoc.write_text("A\tK1\nB\tK2\nC\tK9\n", encoding="utf-8")
        kegg = load_kegg(k2ko)
        assert kegg.k_to_pathways == {"K1": {"ko001", "ko002"}, "K2": {"ko002"}}
        res = annotate_kegg(ModuleSet(modules=[{"A", "B", "C", "D"}]), k2ko, assoc, pval=1.0)
        r = res["module_0"].results.set_index("term_id")
        assert r.loc["ko002", "name"] == "TCA cycle"
        assert r.loc["ko002", "study_count"] == 2 and r.loc["ko001", "study_count"] == 1
        assert set(r.pop_n) == {4} and set(r.study_n) == {4}


def _two_sided(df):
    from scipy.stats import fisher_exact
    a = df.study_count.to_numpy()
    b = df.study_n.to_numpy() - a
    c = df.pop_count.to_numpy() - a
    d = df.pop_n.to_numpy() - df.pop_count.to_numpy() - b
    return np.array([fisher_exact([[w, x], [y, z]])[1] for w, x, y, z in zip(a, b, c, d)])


class TestV1AnnotateFixture:
    """GCEN2 vs output of a compiled GCEN v1 0.6.3 ``annotate`` (see FIXTURE/README.md).

    Queries, background, study/pop counts, names and namespaces must be
    identical; only the test differs (GCEN2 one-sided, v1 two-sided), so the
    two-sided p computed from GCEN2's counts must reproduce v1's p-values.
    """

    @pytest.mark.parametrize("mode", ["network", "module"])
    def test_matches_v1(self, mode):
        from gcen2.annotation import annotate_go
        from gcen2.io import read_modules, read_network
        target = (read_network(FIXTURE / "net.network") if mode == "network"
                  else read_modules(FIXTURE / "modules.txt"))
        res = annotate_go(target, FIXTURE / "go.obo", FIXTURE / "gene_go.assoc", pval=1.0)
        v1_files = sorted((FIXTURE / f"v1_{mode}").glob("*.go"))
        assert len(v1_files) >= 4
        nonempty = {k for k, r in res.items() if not r.results.empty}
        assert nonempty == {f.stem for f in v1_files}
        rank = {"biological_process": 0, "molecular_function": 1, "cellular_component": 2}
        for f in v1_files:
            v1 = pd.read_csv(f, sep="\t", keep_default_na=False)
            v1 = v1.rename(columns={v1.columns[0]: "term_id"}).set_index("term_id")
            got = res[f.stem].results
            got = got.assign(p_two=_two_sided(got)).set_index("term_id")
            # v1 keeps p < 1: every v1 row is in GCEN2, the rest have p == 1
            assert set(v1.index) <= set(got.index), f.stem
            extra = got.loc[got.index.difference(v1.index)]
            np.testing.assert_allclose(extra.p_two, 1.0, rtol=1e-9)
            got = got.loc[v1.index]
            for col in ["study_count", "study_n", "pop_count", "pop_n"]:
                assert list(got[col]) == list(v1[col]), (f.stem, col)
            assert list(got.name) == list(v1.name)
            assert list(got.namespace) == list(v1.name_space)
            np.testing.assert_allclose(got.p_two, v1.p_val, rtol=1e-5)
            # GCEN2 order: namespace (BP, MF, CC), then p-value
            r = res[f.stem].results
            key = list(zip(r.namespace.map(rank).fillna(3), r.p_value))
            assert key == sorted(key)


# ---------------------------------------------------------------------------
# Compiled kernels and the streaming writer (annotate_to_dir)
# ---------------------------------------------------------------------------
class TestCompiledHypergeom:
    def _grid(self, size=20000, seed=3):
        rng = np.random.default_rng(seed)
        M = rng.integers(2, 30000, size)
        n = (rng.uniform(size=size) ** 3 * M).astype(int)
        N = (rng.uniform(size=size) ** 2 * M).astype(int)
        k = (rng.uniform(size=size) * (np.minimum(n, N) + 2)).astype(int)
        return k, M, n, N

    def test_numpy_fallback_matches_scipy(self):
        from scipy.stats import hypergeom

        from gcen2.annotation import _hypergeom_sf_numpy
        k, M, n, N = self._grid()
        np.testing.assert_allclose(_hypergeom_sf_numpy(k, M, n, N),
                                   hypergeom.sf(k - 1, M, n, N), rtol=1e-9, atol=1e-300)

    def test_compiled_matches_numpy(self):
        from gcen2.annotation import _hypergeom_sf_numpy, _hypergeom_sf_vec
        from gcen2.core import extension
        if extension() is None:
            pytest.skip("compiled extension not available")
        k, M, n, N = self._grid(seed=4)
        np.testing.assert_allclose(_hypergeom_sf_vec(k, M, n, N, threads=3),
                                   _hypergeom_sf_numpy(k, M, n, N), rtol=1e-9, atol=1e-300)

    def test_rejects_invalid(self):
        from gcen2.core import extension
        core = extension()
        if core is None:
            pytest.skip("compiled extension not available")
        with pytest.raises(ValueError):
            core.enrichment.hypergeom_sf(np.array([1.0]), np.array([10.0]),
                                         np.array([20.0]), np.array([5.0]))

    def test_py_float_repr(self):
        from gcen2.core import extension
        core = extension()
        if core is None:
            pytest.skip("compiled extension not available")
        rng = np.random.default_rng(5)
        with np.errstate(all="ignore"):
            vals = list(rng.standard_normal(5000) * 10.0 ** rng.integers(-320, 300, 5000))
        vals += [0.0, -0.0, 1.0, 0.1, 1e-4, 1e-5, 1e15, 1e16, 2.5, 1 / 3, 5e-324,
                 1.7976931348623157e308, float("inf")]
        assert [core.enrichment.format_py_float(v) for v in vals] == [repr(float(v)) for v in vals]


class TestAnnotateToDir:
    """annotate_to_dir writes exactly what to_csv of the in-memory results writes."""

    def _expected(self, res, out, suffix, modules=None):
        for k, r in res.items():
            if r.results.empty:
                continue
            r.results.to_csv(out / f"{k}.{suffix}", sep="\t", index=False, lineterminator="\n")
            if modules is not None:
                idx = int(k.split("_", 1)[1])
                with open(out / f"{k}.txt", "w", encoding="utf-8", newline="\n") as fh:
                    fh.write("\n".join(sorted(modules[idx])) + "\n")

    def _assert_same_dir(self, a, b):
        fa = sorted(p.name for p in a.iterdir())
        assert fa == sorted(p.name for p in b.iterdir()) and fa
        for name in fa:
            assert (a / name).read_bytes() == (b / name).read_bytes(), name

    @pytest.mark.parametrize("mode", ["network", "module"])
    @pytest.mark.parametrize("compiled", [True, False])
    @pytest.mark.parametrize("chunk_size", [1, 3, 500])
    def test_matches_in_memory(self, tmp_path, monkeypatch, mode, compiled, chunk_size):
        import gcen2.annotation as A
        import gcen2.core
        from gcen2.io import read_modules, read_network
        if not compiled:
            monkeypatch.setattr(gcen2.core, "_gcen_core", None)
        elif gcen2.core.extension() is None:
            pytest.skip("compiled extension not available")
        target = (read_network(FIXTURE / "net.network") if mode == "network"
                  else read_modules(FIXTURE / "modules.txt"))
        res = A.annotate_go(target, FIXTURE / "go.obo", FIXTURE / "gene_go.assoc", pval=0.5)
        ref, got = tmp_path / "ref", tmp_path / "got"
        ref.mkdir()
        self._expected(res, ref, "go", None if mode == "network" else target.modules)
        n = A.annotate_to_dir(target, got, FIXTURE / "gene_go.assoc", go=FIXTURE / "go.obo",
                              pval=0.5, threads=2, chunk_size=chunk_size)
        assert n == sum(not r.results.empty for r in res.values())
        self._assert_same_dir(ref, got)

    def test_genes_only_network_same_output(self, tmp_path):
        from gcen2.annotation import annotate_to_dir
        from gcen2.io import read_network
        a, b = tmp_path / "a", tmp_path / "b"
        for out, genes_only in [(a, False), (b, True)]:
            annotate_to_dir(read_network(FIXTURE / "net.network", genes_only=genes_only), out,
                            FIXTURE / "gene_go.assoc", go=FIXTURE / "go.obo")
        self._assert_same_dir(a, b)

    def test_quoted_names_like_pandas(self, tmp_path, monkeypatch):
        import gcen2.annotation as A
        if gcen2.core.extension() is None:
            pytest.skip("compiled extension not available")
        ontology = A.GOOntology()
        for i, name in enumerate(['plain', 'has "quote"', 'has\ttab']):
            ontology.terms[f"GO:{i}"] = A.GOTerm(id=f"GO:{i}", name=name,
                                                 namespace="biological_process")
        gene_terms = {f"a{j}": {f"GO:{j % 3}"} for j in range(30)}
        from gcen2.models import ModuleSet
        mods = ModuleSet(modules=[{f"a{j}" for j in range(0, 30, 3)},
                                  {f"a{j}" for j in range(1, 30, 3)}, {"a2", "a5", "a8"}])
        res = A.annotate_go(mods, ontology, gene_terms, pval=1.0)
        ref, got = tmp_path / "ref", tmp_path / "got"
        ref.mkdir()
        self._expected(res, ref, "go", mods.modules)
        A.annotate_to_dir(mods, got, gene_terms, go=ontology, pval=1.0)
        self._assert_same_dir(ref, got)


class TestPandasFreeAnnotate:
    """NetworkAdjacency input, the C++ study_counts kernel and the
    pandas-free loaders give the same results as before."""

    def test_adjacency_queries_equal_network_queries(self, tmp_path):
        import gcen2.annotation as A
        from gcen2.io import read_network
        if gcen2.core.extension() is None:
            pytest.skip("compiled extension not available")
        rng = np.random.default_rng(3)
        genes = [f"g{i}" for i in range(50)]
        pairs = rng.integers(0, 50, size=(300, 2))
        lines = [f"{genes[i]}\t{genes[j]}\t0.9\t1e-5" for i, j in pairs]
        lines += [lines[0], "\t".join(lines[1].split("\t")[1::-1] + ["0.9", "1e-5"]),
                  "g1\tg1\t1\t0"]
        p = tmp_path / "net.tsv"
        p.write_text("gene_a\tgene_b\tcorrelation\tp_value\n" + "\n".join(lines) + "\n", encoding="utf-8")
        gene_terms = {g: {"GO:1"} for g in genes[::3]}
        ref = A.build_annotation_queries(read_network(p), gene_terms)
        got = A.build_annotation_queries(A.NetworkAdjacency.from_file(p), gene_terms)
        assert ref[0] and got == ref

    def test_adjacency_annotate_to_dir_same_files(self, tmp_path):
        import gcen2.annotation as A
        from gcen2.io import read_network
        if gcen2.core.extension() is None:
            pytest.skip("compiled extension not available")
        a, b = tmp_path / "a", tmp_path / "b"
        A.annotate_to_dir(read_network(FIXTURE / "net.network"), a, FIXTURE / "gene_go.assoc",
                          go=FIXTURE / "go.obo", pval=0.5)
        A.annotate_to_dir(A.NetworkAdjacency.from_file(FIXTURE / "net.network"), b,
                          FIXTURE / "gene_go.assoc", go=FIXTURE / "go.obo", pval=0.5)
        TestAnnotateToDir()._assert_same_dir(a, b)

    def test_study_counts_matches_sparse_product(self):
        from scipy.sparse import random as sprandom

        core = gcen2.core.extension()
        if core is None:
            pytest.skip("compiled extension not available")
        G = (sprandom(200, 40, density=0.1, random_state=1, format="csr") > 0).astype(np.int32)
        Q = (sprandom(30, 200, density=0.05, random_state=2, format="csr") > 0).astype(np.int32)
        for threads in (1, 4):
            r, t, c = core.enrichment.study_counts(
                Q.indptr.astype(np.int64), Q.indices.astype(np.int64),
                G.indptr.astype(np.int64), G.indices.astype(np.int32), 40, threads)
            assert np.all(np.diff(r) >= 0)
            S = (Q @ G).tocoo()
            assert sorted(zip(r, t, c)) == sorted(zip(S.row, S.col, S.data))
        with pytest.raises(ValueError):
            core.enrichment.study_counts(Q.indptr.astype(np.int64), Q.indices.astype(np.int64),
                                         G.indptr.astype(np.int64),
                                         G.indices.astype(np.int32), 5, 1)

    def test_duplicate_terms_in_dict_assoc_counted_once(self, tmp_path):
        import gcen2.annotation as A
        from gcen2.models import ModuleSet
        mods = ModuleSet(modules=[{"a", "b"}, {"c", "d", "e"}])
        as_sets = {"a": {"T1", "T2"}, "b": {"T1"}, "c": {"T2"}, "d": {"T2"}, "e": {"T1"}}
        as_lists = {g: sorted(ts) * 2 for g, ts in as_sets.items()}
        r1 = A.annotate_go(mods, A.GOOntology(), as_sets, pval=1.0)
        r2 = A.annotate_go(mods, A.GOOntology(), as_lists, pval=1.0)
        for k in r1:
            pd.testing.assert_frame_equal(r1[k].results, r2[k].results)

    def test_propagate_matches_per_gene_walk(self, obo_file, assoc_file):
        ontology = load_go(obo_file)
        gene_terms = load_association(assoc_file)
        gene_terms["x"] = {"GO:unknown"}  # term not in the ontology is kept
        got = ontology.propagate(gene_terms)
        for gene, terms in gene_terms.items():
            want, queue = set(), list(terms)
            while queue:
                t = queue.pop()
                if t not in want:
                    want.add(t)
                    if t in ontology.terms:
                        queue.extend(ontology.terms[t].parents)
            assert got[gene] == want

    def test_load_go_crlf(self, obo_file, tmp_path):
        crlf = tmp_path / "crlf.obo"
        crlf.write_bytes(Path(obo_file).read_bytes().replace(b"\n", b"\r\n"))
        def key(o):
            return {k: (v.name, v.namespace, v.parents, v.alt_ids, v.is_obsolete)
                    for k, v in o.terms.items()}
        assert key(load_go(crlf)) == key(load_go(obo_file))

    def test_cli_paths_do_not_import_pandas(self):
        import subprocess
        import sys
        code = ("import sys; import gcen2.cli; from gcen2 import annotation, rwr; "
                "from gcen2.io import read_gene_list, read_modules; "
                "print(sorted(m for m in ('pandas', 'scipy') if m in sys.modules))")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8",
                             check=True)
        assert out.stdout.strip() == "[]"


class TestEncodedAnnotate:
    """Integer CSR gene terms / queries and the fused significant_counts
    kernel used by annotate_to_dir."""

    @staticmethod
    def _as_dict(gt):
        return {g: set(gt.term_ids[gt.g_idx[gt.row_ptr[i]:gt.row_ptr[i + 1]]].tolist())
                for i, g in enumerate(gt.genes)}

    def test_gene_terms_with_closure_equal_propagate(self, obo_file, assoc_file):
        import gcen2.annotation as A
        ontology = load_go(obo_file)
        listed = load_association(assoc_file)
        listed["x"] = {"GO:unknown"}
        listed["empty"] = set()
        gt = A._GeneTerms._from_rows(listed.keys(), listed.values(), ontology._closure_fn())
        assert self._as_dict(gt) == ontology.propagate(listed)
        assert list(gt.term_ids) == sorted(gt.term_ids)
        for i in range(len(gt.genes)):  # terms unique per gene
            row = gt.g_idx[gt.row_ptr[i]:gt.row_ptr[i + 1]]
            assert len(set(row.tolist())) == len(row)

    def test_read_association_matches_load_association(self, assoc_file):
        import gcen2.annotation as A
        assert A._read_association(Path(assoc_file)) == load_association(assoc_file)

    def test_queries_from_adjacency_equal_dict(self, tmp_path):
        import gcen2.annotation as A
        if gcen2.core.extension() is None:
            pytest.skip("compiled extension not available")
        rng = np.random.default_rng(5)
        genes = [f"g{i}" for i in range(60)]
        pairs = rng.integers(0, 60, size=(400, 2))
        p = tmp_path / "net.tsv"
        p.write_text("".join(f"{genes[i]}\t{genes[j]}\t0.9\t1e-5\n" for i, j in pairs), encoding="utf-8")
        gene_terms = {g: {"T"} for g in genes[::4]}
        gene_terms["absent"] = {"T"}  # annotated but not in the network
        gt = A._GeneTerms.from_dict(gene_terms)
        adj = A.NetworkAdjacency.from_file(p)
        got = A._Queries.from_adjacency(adj, gt)
        want = A._Queries.from_dict(*A.build_annotation_queries(adj, gene_terms), gt)
        assert got.keys and got.keys == want.keys and got.pop_n == want.pop_n
        for a in ("ptr", "study_n", "in_bg"):
            np.testing.assert_array_equal(getattr(got, a), getattr(want, a))
        for b, e in zip(got.ptr[:-1], got.ptr[1:]):  # same genes (order is free)
            assert sorted(got.idx[b:e]) == sorted(want.idx[b:e])
        # no annotated / no unannotated genes: empty queries
        for gtx in (A._GeneTerms.from_dict({}), A._GeneTerms.from_dict(
                {g: {"T"} for g in genes})):
            q = A._Queries.from_adjacency(adj, gtx)
            assert q.keys == [] and list(q.ptr) == [0] and q.pop_n == 0

    def test_significant_counts_matches_unfused(self):
        from scipy.sparse import random as sprandom

        import gcen2.annotation as A
        core = gcen2.core.extension()
        if core is None:
            pytest.skip("compiled extension not available")
        G = (sprandom(200, 40, density=0.15, random_state=1, format="csr") > 0).astype(np.int32)
        Q = (sprandom(30, 200, density=0.08, random_state=2, format="csr") > 0).astype(np.int32)
        args = (Q.indptr.astype(np.int64), Q.indices.astype(np.int64),
                G.indptr.astype(np.int64), G.indices.astype(np.int32), 40)
        in_bg = np.arange(200) % 3 != 0
        pop_counts = np.asarray(G[in_bg].sum(axis=0)).ravel().astype(np.int64)
        pop_counts[:3] = 0  # terms absent from the background
        pop_n = 200
        study_n = np.diff(Q.indptr).astype(np.int64) + 2
        r, t, c = core.enrichment.study_counts(*args, 1)
        ok = pop_counts[t] > 0
        hits_ref = np.bincount(r[ok], minlength=30) > 0
        r, t, c = r[ok], t[ok], c[ok]
        p = A._hypergeom_sf_vec(c, pop_n, pop_counts[t], study_n[r])
        fdr = A._bh_by_row(r, p, int(np.count_nonzero(pop_counts)))
        for pval in (1.0, 0.05):
            keep = p <= pval
            ref = sorted(zip(r[keep], t[keep], c[keep], p[keep], fdr[keep]))
            for threads in (1, 4):
                rr, tt, cc, pp, ff, hh = core.enrichment.significant_counts(
                    *args, study_n, pop_counts, pop_n, pval, threads)
                assert np.all(np.diff(rr) >= 0)
                assert sorted(zip(rr, tt, cc, pp, ff)) == ref  # p and fdr bit-identical
                np.testing.assert_array_equal(hh, hits_ref)
        with pytest.raises(ValueError):
            core.enrichment.significant_counts(*args, study_n, pop_counts, 5, 0.05, 1)
        with pytest.raises(ValueError):
            core.enrichment.significant_counts(*args, study_n[:-1], pop_counts, pop_n, 0.05, 1)

    def test_enrich_query_larger_than_background(self):
        # study_n > pop_n: annotate never produces this, enrich() may. The
        # fused kernel is skipped and hypergeom_sf rejects it, as before.
        import gcen2.annotation as A
        if gcen2.core.extension() is None:
            pytest.skip("compiled extension not available")
        assoc = {"a": {"T"}, "b": {"T"}, "c": {"T", "U"}, "d": {"U"}}
        with pytest.raises(ValueError, match="hypergeom_sf"):
            A.enrich(["a", "b", "c", "d"], ["a", "b"], None, assoc, pval=1.0)
