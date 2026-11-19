"""GO/KEGG functional annotation and enrichment analysis."""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, NamedTuple, Union

import numpy as np

from gcen2._matrix_io import require_file, tsv_field
from gcen2.core import extension
from gcen2.models import EnrichmentResult, ModuleSet, Network

# pandas and scipy are imported only where needed: ``gcen2 annotate``
# (NetworkAdjacency or ModuleSet input, compiled extension) needs neither.
if TYPE_CHECKING:
    import pandas as pd


class NetworkAdjacency(NamedTuple):
    """A network as symmetric CSR adjacency over sorted gene names, as
    returned by :func:`gcen2.module.read_file_adjacency`. Accepted wherever
    a :class:`Network` is (annotation only needs who is adjacent to whom)."""

    indptr: np.ndarray
    indices: np.ndarray
    names: list

    @classmethod
    def from_file(cls, path) -> NetworkAdjacency:
        from gcen2.module import read_file_adjacency

        indptr, indices, names = read_file_adjacency(path)
        return cls(indptr, indices, names)


@dataclass
class GOTerm:
    id: str
    name: str = ""
    namespace: str = ""
    parents: set[str] = field(default_factory=set)
    alt_ids: set[str] = field(default_factory=set)
    is_obsolete: bool = False


@dataclass
class GOOntology:
    terms: dict[str, GOTerm] = field(default_factory=dict)

    def propagate(self, gene_terms: dict[str, set[str]]) -> dict[str, set[str]]:
        """Add all ancestor terms to each gene's annotation set.

        The ancestor closure of each distinct term is computed once and
        shared, instead of walking the DAG again for every gene.
        """
        ancestors = self._closure_fn()
        return {gene: set().union(*(ancestors(t) for t in term_set))
                for gene, term_set in gene_terms.items()}

    def _closure_fn(self):
        """A memoised ``term -> frozenset(term and its is_a ancestors)``."""
        closure: dict[str, frozenset] = {}
        terms = self.terms

        def ancestors(tid: str) -> frozenset:
            got = closure.get(tid)
            if got is None:
                expanded: set[str] = set()
                queue = [tid]
                while queue:
                    t = queue.pop()
                    if t in expanded:
                        continue
                    expanded.add(t)
                    if t in terms:
                        queue.extend(terms[t].parents)
                got = closure[tid] = frozenset(expanded)
            return got

        return ancestors


def load_go(obo_file: Union[str, Path]) -> GOOntology:
    """Parse a GO OBO file and return a GOOntology object."""
    obo_file = require_file(obo_file)

    ontology = GOOntology()
    terms = ontology.terms
    current: GOTerm | None = None

    from gcen2._matrix_io import iter_lines

    for line in iter_lines(obo_file):
        # most lines (def, synonym, xref, ...) are skipped by first character
        c = line[:1]
        if c not in _OBO_FIRST_CHARS:
            c = line.lstrip()[:1]
            if c not in _OBO_FIRST_CHARS:
                continue
        line = line.strip()
        if line == "[Term]":
            if current is not None:
                terms[current.id] = current
            current = GOTerm(id="")
        elif line == "[Typedef]":
            if current is not None:
                terms[current.id] = current
            current = None
        elif current is None:
            continue
        elif line.startswith("id:"):
            current.id = line[3:].strip()
        elif line.startswith("name:"):
            current.name = line[5:].strip()
        elif line.startswith("namespace:"):
            current.namespace = line[10:].strip()
        elif line.startswith("alt_id:"):
            current.alt_ids.add(line[7:].strip())
        elif line.startswith("is_a:"):
            parent_id = line[5:].split("!")[0].strip()
            current.parents.add(parent_id)
        elif line.startswith("is_obsolete: true"):
            current.is_obsolete = True

    if current is not None:
        terms[current.id] = current

    return ontology


# first characters of the OBO lines load_go reads ([Term], id, name,
# namespace, alt_id, is_a, is_obsolete)
_OBO_FIRST_CHARS = frozenset("[ina")


def load_association(
    assoc_file: Union[str, Path],
    ontology: GOOntology | None = None,
) -> dict[str, set[str]]:
    """Load a gene-term association file (GCEN v1 format).

    Each line is ``gene <tab> term [<tab> term ...]``; a field may also hold
    several comma-separated terms. Repeated genes are merged. Genes listed
    without any term are kept (with an empty set) because GCEN v1 counts them
    as annotated.

    If a GOOntology is provided, terms are propagated to all ``is_a``
    ancestors (as GCEN v1 ``propagate``).
    """
    assoc_file = require_file(assoc_file)

    gene_terms = _read_association(assoc_file)
    if ontology is not None:
        gene_terms = ontology.propagate(gene_terms)

    return gene_terms


def _read_association(assoc_file: Path) -> dict[str, set[str]]:
    """Gene -> terms as listed in an association file (no propagation)."""
    from gcen2._matrix_io import iter_lines

    gene_terms: dict[str, set[str]] = {}
    # one string object per distinct term (the same few thousand ids repeat
    # across genes; separate copies took most of this dict's memory)
    shared: dict[str, str] = {}
    for line in iter_lines(assoc_file):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        terms = gene_terms.setdefault(parts[0].strip(), set())
        for fld in parts[1:]:
            for t in fld.split(","):
                t = t.strip()
                if t:
                    terms.add(shared.setdefault(t, t))
    return gene_terms


@dataclass
class KEGGMapping:
    """KO -> pathway mapping from a GCEN v1 ``K2ko`` file."""

    k_to_pathways: dict[str, set[str]] = field(default_factory=dict)
    pathway_names: dict[str, str] = field(default_factory=dict)


def load_kegg(k2ko_file: Union[str, Path]) -> KEGGMapping:
    """Parse a ``K2ko`` file: ``K number <tab> pathway id <tab> description``."""
    k2ko_file = require_file(k2ko_file)
    mapping = KEGGMapping()
    from gcen2._matrix_io import iter_lines

    for line in iter_lines(k2ko_file):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        k, ko = parts[0].strip(), parts[1].strip()
        mapping.k_to_pathways.setdefault(k, set()).add(ko)
        if ko not in mapping.pathway_names:
            mapping.pathway_names[ko] = parts[2].strip() if len(parts) > 2 else ""
    return mapping


def load_kegg_association(
    assoc_file: Union[str, Path], kegg: KEGGMapping
) -> dict[str, set[str]]:
    """Gene -> KEGG pathways, from a gene-to-K-number association file."""
    gene_k = load_association(assoc_file)
    return {
        gene: set().union(*(kegg.k_to_pathways.get(k, set()) for k in ks)) if ks else set()
        for gene, ks in gene_k.items()
    }


def annotate_go(
    network_or_modules: Union[Network, ModuleSet],
    obo: Union[str, Path, GOOntology],
    assoc: Union[str, Path, dict],
    pval: float = 0.05,
    threads: int = 2,
) -> dict[str, EnrichmentResult]:
    """GO annotation of a network or module set (GCEN v1 ``annotate``).

    Network: every gene *without* GO terms is annotated from its annotated
    neighbours (guilt by association); the background is the set of
    annotated genes adjacent to an unannotated gene. Modules: each module is
    tested against all module genes. See :func:`build_annotation_queries`.
    The test is one-sided (over-representation), keeping ``p <= pval``.

    Returns a dict mapping gene / ``module_<i>`` -> EnrichmentResult; queries
    without any hit map to an empty result.
    """
    queries, background, gene_terms, term_info = _go_inputs(network_or_modules, obo, assoc)
    return _batch_enrichment(queries, background, gene_terms, ontology_type="GO",
                             pval=pval, term_info=term_info, threads=threads)


def _go_inputs(network_or_modules, obo, assoc):
    ontology = load_go(obo) if isinstance(obo, (str, Path)) else obo
    gene_terms = load_association(assoc, ontology) if isinstance(assoc, (str, Path)) else assoc
    queries, background = build_annotation_queries(network_or_modules, gene_terms)
    term_info = {t.id: (t.name, t.namespace) for t in ontology.terms.values()}
    return queries, background, gene_terms, term_info


def _kegg_inputs(network_or_modules, k2ko, assoc):
    kegg = load_kegg(k2ko) if isinstance(k2ko, (str, Path)) else k2ko
    gene_terms = (load_kegg_association(assoc, kegg) if isinstance(assoc, (str, Path))
                  else assoc)
    queries, background = build_annotation_queries(network_or_modules, gene_terms)
    term_info = {ko: (name, "") for ko, name in kegg.pathway_names.items()}
    return queries, background, gene_terms, term_info


def annotate_kegg(
    network_or_modules: Union[Network, ModuleSet],
    k2ko: Union[str, Path, KEGGMapping],
    assoc: Union[str, Path, dict],
    pval: float = 0.05,
    threads: int = 2,
) -> dict[str, EnrichmentResult]:
    """KEGG pathway annotation; same query/background rules as :func:`annotate_go`.

    ``assoc`` maps genes to KO (K) numbers; ``k2ko`` maps K numbers to
    pathways. A dict ``assoc`` is taken as gene -> pathways directly.
    """
    queries, background, gene_terms, term_info = _kegg_inputs(network_or_modules, k2ko, assoc)
    return _batch_enrichment(queries, background, gene_terms, ontology_type="KEGG",
                             pval=pval, term_info=term_info, threads=threads)


def annotate_to_dir(
    network_or_modules: Union[Network, ModuleSet],
    out_dir: Union[str, Path],
    assoc: Union[str, Path, dict],
    go: Union[str, Path, GOOntology, None] = None,
    kegg: Union[str, Path, KEGGMapping, None] = None,
    pval: float = 0.05,
    threads: int = 2,
    chunk_size: int = 500,
) -> int:
    """Run :func:`annotate_go` (``go``) or :func:`annotate_kegg` (``kegg``) and
    write the results like GCEN v1 ``annotate``, chunk by chunk.

    Writes ``<gene>.go`` / ``module_<i>.go`` (``.kegg`` for KEGG) for every
    query with significant terms, plus ``module_<i>.txt`` with the module's
    genes in module mode. Each file has the same content as
    ``result.results.to_csv(path, sep="\t", index=False)`` of the
    corresponding in-memory result, but results are written as each chunk
    of ``chunk_size`` queries is finished, so memory does not grow with the
    number of results. Returns the number of annotation files written.
    """
    # Gene terms and queries are kept as integer CSR arrays (see _GeneTerms):
    # the propagated GO sets as Python string sets took most of the memory.
    if go is not None:
        ontology = load_go(go) if isinstance(go, (str, Path)) else go
        if isinstance(assoc, (str, Path)):
            listed = _read_association(require_file(assoc))
            gt = _GeneTerms._from_rows(listed.keys(), listed.values(), ontology._closure_fn())
            del listed
        else:
            gt = _GeneTerms.from_dict(assoc)
        term_info = {t.id: (t.name, t.namespace) for t in ontology.terms.values()}
        del ontology
        suffix = "go"
    elif kegg is not None:
        kegg = load_kegg(kegg) if isinstance(kegg, (str, Path)) else kegg
        gt = _GeneTerms.from_dict(load_kegg_association(assoc, kegg)
                                  if isinstance(assoc, (str, Path)) else assoc)
        term_info = {ko: (name, "") for ko, name in kegg.pathway_names.items()}
        suffix = "kegg"
    else:
        raise ValueError("one of go or kegg is required")
    if isinstance(network_or_modules, NetworkAdjacency):
        q = _Queries.from_adjacency(network_or_modules, gt)
    else:
        gene_sets = dict.fromkeys(gt.genes)  # membership is all queries need
        q = _Queries.from_dict(*build_annotation_queries(network_or_modules, gene_sets), gt)
    modules = (None if isinstance(network_or_modules, (Network, NetworkAdjacency))
               else network_or_modules.modules)
    out_dir = Path(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    header = ("\t".join(_RESULT_COLUMNS) + "\n").encode()
    written = 0
    core = extension()
    prefix = None
    for c in _iter_encoded(q, gt, pval, chunk_size, term_info, threads):
        if len(c.terms) == 0:
            continue
        if core is not None:
            if prefix is None:
                prefix = ["\t".join(tsv_field(v) for v in row)
                          for row in zip(*c.terms_meta)]
            blobs = core.enrichment.format_annotation_rows(
                prefix, c.terms, c.enrichment, c.study_count, c.study_n, c.pop_count,
                c.pop_n, c.p_value, c.fdr, c.bounds, max(1, int(threads)))
        else:
            # one to_csv call per chunk; rows of one query are consecutive lines
            lines = c.table().to_csv(sep="\t", index=False, header=False,
                                     lineterminator="\n").split("\n")
            if len(lines) != len(c.terms) + 1:  # a field with a line break
                raise ValueError("annotation fields must not contain line breaks")
            blobs = [("\n".join(lines[b:e]) + "\n").encode() if e > b else b""
                     for b, e in zip(c.bounds[:-1], c.bounds[1:])]
        for key, blob in zip(c.keys, blobs):
            if not blob:
                continue
            with open(out_dir / f"{key}.{suffix}", "wb") as fh:
                fh.write(header)
                fh.write(blob)
            if modules is not None:
                idx = int(key.split("_", 1)[1])
                with open(out_dir / f"{key}.txt", "w", newline="\n", encoding="utf-8") as fh:
                    fh.write("\n".join(sorted(modules[idx])) + "\n")
            written += 1
    return written


def enrich(
    query_genes: list[str],
    background_genes: list[str],
    ontology: Union[GOOntology, None],
    assoc: dict[str, set[str]],
    pval: float = 0.05,
) -> EnrichmentResult:
    """Enrichment of a gene list against a background list (GCEN v1 ``enrich``).

    As in v1, only genes present in the association are counted, both in the
    query (``study_n``) and in the background (``pop_n``).
    """
    query = {g for g in query_genes if g in assoc}
    background = {g for g in background_genes if g in assoc}
    term_info = ({t.id: (t.name, t.namespace) for t in ontology.terms.values()}
                 if ontology is not None else None)
    return _batch_enrichment({"query": query}, background, assoc, ontology_type="GO",
                             pval=pval, term_info=term_info)["query"]


def build_annotation_queries(
    network_or_modules: Union[Network, ModuleSet],
    gene_terms: dict[str, set[str]],
) -> tuple[dict[str, set[str]], set[str]]:
    """Query gene sets and background, following GCEN v1 ``annotate``.

    Network (``load_network`` in v1 ``enrich.hpp``): for every edge between an
    unannotated gene ``u`` and an annotated gene ``a`` (annotated = listed in
    the association), ``a`` joins the query of ``u`` and the background.
    Edges between two annotated or two unannotated genes are ignored, so only
    unannotated genes are annotated.

    Modules (``load_module``): each module (all its genes) is a query and the
    background is the union of all module genes.

    A :class:`NetworkAdjacency` gives the same queries and background as the
    :class:`Network` it was read from.
    """
    if isinstance(network_or_modules, NetworkAdjacency):
        indptr, indices, names = network_or_modules
        n = len(names)
        if n == 0:
            return {}, set()
        names = np.array(names, dtype=object)
        is_ann = np.fromiter((g in gene_terms for g in names), dtype=bool, count=n)
        rows = np.repeat(np.arange(n), np.diff(indptr))
        # unannotated gene -> its annotated neighbours (each edge is stored
        # in both rows, so every (u, a) pair appears once with u as the row)
        m = ~is_ann[rows] & is_ann[indices]
        rows, nbr = rows[m], indices[m]
        if rows.size == 0:
            return {}, set()
        bounds = np.flatnonzero(np.diff(rows)) + 1
        starts = np.r_[0, bounds]
        queries = {names[rows[b]]: set(names[g].tolist())
                   for b, g in zip(starts, np.split(nbr, bounds))}
        return queries, set(names[np.unique(nbr)].tolist())

    if isinstance(network_or_modules, Network):
        import pandas as pd

        edges = network_or_modules.edges
        if edges.empty:
            return {}, set()
        annotated = set(gene_terms)
        # integer codes over a shared vocabulary of gene names (cheap for
        # categorical columns, e.g. read_network(..., genes_only=True))
        cat = pd.Categorical(edges["gene_a"])
        vocab = cat.categories.union(pd.Index(pd.Categorical(edges["gene_b"]).categories))
        ca = pd.Categorical(edges["gene_a"], categories=vocab).codes
        cb = pd.Categorical(edges["gene_b"], categories=vocab).codes
        vocab = vocab.to_numpy(dtype=object)
        is_ann = np.fromiter((g in annotated for g in vocab), dtype=bool, count=len(vocab))
        ina, inb = is_ann[ca], is_ann[cb]
        m1, m2 = ~ina & inb, ina & ~inb
        target = vocab[np.concatenate([ca[m1], cb[m2]])]
        neighbour = vocab[np.concatenate([cb[m1], ca[m2]])]
        if target.size == 0:
            return {}, set()
        codes, uniques = pd.factorize(target)
        order = np.argsort(codes, kind="stable")
        bounds = np.flatnonzero(np.diff(codes[order])) + 1
        groups = np.split(neighbour[order], bounds)
        keys = uniques[codes[order][np.r_[0, bounds]]]
        queries = {k: set(g.tolist()) for k, g in zip(keys.tolist(), groups)}
        return queries, set(neighbour.tolist())

    queries = {f"module_{i}": set(m) for i, m in enumerate(network_or_modules.modules)}
    background = set().union(*queries.values()) if queries else set()
    return queries, background


_NAMESPACE_RANK = {"biological_process": 0, "molecular_function": 1, "cellular_component": 2}

_RESULT_COLUMNS = [
    "term_id", "name", "namespace", "enrichment", "study_count",
    "study_n", "pop_count", "pop_n", "p_value", "fdr",
]


def _batch_enrichment(
    queries: dict[str, set[str]],
    background: set[str],
    gene_terms: dict[str, set[str]],
    ontology_type: str,
    pval: float,
    chunk_size: int = 2000,
    term_info: dict[str, tuple[str, str]] | None = None,
    threads: int = 1,
) -> dict[str, EnrichmentResult]:
    """Hypergeometric enrichment for many query gene sets at once.

    For every query and every term annotated in ``background``:

    * ``study_n``     = number of query genes
    * ``study_count`` = query genes annotated with the term
    * ``pop_n``       = number of background genes
    * ``pop_count``   = background genes annotated with the term
    * ``p_value``     = P(X >= study_count), X ~ Hypergeom(pop_n, pop_count, study_n)
    * ``fdr``         = Benjamini-Hochberg adjusted p over every term annotated in
      the background (terms without a query gene count with p = 1), as R's
      ``p.adjust(p, "BH")`` on all of them

    Terms with ``study_count == 0`` are skipped and only ``p_value <= pval`` is
    reported. ``term_info`` maps term id -> (name, namespace); results are
    ordered as in GCEN v1: biological_process, molecular_function,
    cellular_component, other, then by p-value (then term id).

    Study counts for all queries come from one sparse product
    (queries x genes) @ (genes x terms), processed in chunks of
    ``chunk_size`` queries, and p-values are evaluated vectorised (see
    :func:`_iter_enrichment`, which yields the chunks without keeping them).
    """
    import pandas as pd

    results: dict[str, EnrichmentResult] = {}
    for c in _iter_enrichment(queries, background, gene_terms, pval, chunk_size,
                              term_info, threads):
        table, bounds = c.table(), c.bounds
        for r, k in enumerate(c.keys):
            if not c.has_hits[r]:
                part = pd.DataFrame()
            else:
                part = table.iloc[bounds[r]:bounds[r + 1]].reset_index(drop=True)
            results[k] = EnrichmentResult(
                results=part, ontology_type=ontology_type, gene_or_module=k, p_cutoff=pval
            )
    return results


@dataclass
class _Chunk:
    """Significant results of one chunk of queries, sorted by query then in
    v1 order; rows of ``keys[r]`` are ``bounds[r]:bounds[r + 1]``.
    ``has_hits[r]`` is False when the query shares no term with the
    background at all."""

    keys: list
    bounds: np.ndarray
    has_hits: np.ndarray
    terms: np.ndarray        # term codes into terms_meta
    study_count: np.ndarray
    study_n: np.ndarray
    pop_count: np.ndarray
    pop_n: int
    p_value: np.ndarray
    fdr: np.ndarray
    terms_meta: tuple        # (term_ids, names, namespaces) by term code

    @property
    def enrichment(self) -> np.ndarray:
        return (self.study_count / self.study_n) / (self.pop_count / self.pop_n)

    def table(self) -> pd.DataFrame:
        import pandas as pd

        term_ids, term_name, term_ns = self.terms_meta
        t = self.terms
        return pd.DataFrame({
            "term_id": term_ids[t],
            "name": term_name[t],
            "namespace": term_ns[t],
            "enrichment": self.enrichment,
            "study_count": self.study_count,
            "study_n": self.study_n,
            "pop_count": self.pop_count,
            "pop_n": np.full(len(t), self.pop_n, dtype=np.int64),
            "p_value": self.p_value,
            "fdr": self.fdr,
        }, columns=_RESULT_COLUMNS)


class _GeneTerms(NamedTuple):
    """Gene -> term incidence as CSR: the terms of ``genes[i]`` are
    ``term_ids[g_idx[row_ptr[i]:row_ptr[i + 1]]]``, each listed once.
    ``term_ids`` is sorted (any coding works: output order is set later).
    Integer arrays need a fraction of the memory of a dict of string sets."""

    genes: list
    row_ptr: np.ndarray   # int64
    g_idx: np.ndarray     # int32 term codes
    term_ids: np.ndarray  # object

    @classmethod
    def _from_rows(cls, genes, rows, closure=None) -> _GeneTerms:
        """``rows``: an iterable of term sets, one per gene. With
        ``closure`` (term -> set of terms), each term is replaced by its
        closure (GO propagation) without building the propagated sets."""
        from array import array

        rows = list(rows)
        if closure is not None:
            # code the closures of the distinct listed terms once
            listed = set().union(*rows) if rows else set()
            closed = {t: closure(t) for t in listed}
            vocab = set().union(*closed.values()) if closed else set()
        else:
            closed = None
            vocab = set().union(*rows) if rows else set()
        term_list = sorted(vocab, key=str)
        code = {t: i for i, t in enumerate(term_list)}
        if closed is not None:
            closed = {t: frozenset(map(code.__getitem__, c)) for t, c in closed.items()}
            sets = (set().union(*map(closed.__getitem__, ts)) for ts in rows)
        else:
            # each gene's terms must be unique (they are sets when loaded from files)
            sets = ({code[t] for t in ts} for ts in rows)
        g_idx = array("i")
        lens = np.empty(len(rows), dtype=np.int64)
        for i, codes in enumerate(sets):
            lens[i] = len(codes)
            g_idx.extend(codes)
        row_ptr = np.zeros(len(rows) + 1, dtype=np.int64)
        np.cumsum(lens, out=row_ptr[1:])
        term_ids = np.empty(len(term_list), dtype=object)
        term_ids[:] = term_list
        return cls(list(genes), row_ptr, np.frombuffer(g_idx, dtype=np.int32).copy()
                   if len(g_idx) else np.zeros(0, dtype=np.int32), term_ids)

    @classmethod
    def from_dict(cls, gene_terms: dict) -> _GeneTerms:
        return cls._from_rows(gene_terms.keys(), gene_terms.values())


class _Queries(NamedTuple):
    """Query gene sets as CSR over ``_GeneTerms.genes`` (annotated genes
    only); ``study_n`` counts every query gene. ``in_bg`` marks annotated
    background genes, ``pop_n`` is the background size."""

    keys: list
    ptr: np.ndarray       # int64
    idx: np.ndarray       # int64 rows of _GeneTerms
    study_n: np.ndarray   # int64
    in_bg: np.ndarray     # bool, per _GeneTerms gene
    pop_n: int

    @classmethod
    def from_dict(cls, queries: dict, background, gt: _GeneTerms) -> _Queries:
        gene_idx = {g: i for i, g in enumerate(gt.genes)}
        keys = list(queries)
        study_n = np.fromiter((len(queries[k]) for k in keys), dtype=np.int64, count=len(keys))
        cols = [[gene_idx[g] for g in queries[k] if g in gene_idx] for k in keys]
        ptr = np.zeros(len(keys) + 1, dtype=np.int64)
        ptr[1:] = np.cumsum([len(c) for c in cols])
        idx = np.fromiter((c for cs in cols for c in cs), dtype=np.int64, count=ptr[-1])
        in_bg = np.fromiter((g in background for g in gt.genes), dtype=bool,
                            count=len(gt.genes))
        return cls(keys, ptr, idx, study_n, in_bg, len(background))

    @classmethod
    def from_adjacency(cls, adj: NetworkAdjacency, gt: _GeneTerms) -> _Queries:
        """Same queries and background as :func:`build_annotation_queries`,
        without building Python sets."""
        indptr, indices, names = adj
        gene_idx = {g: i for i, g in enumerate(gt.genes)}
        row_of = np.fromiter((gene_idx.get(g, -1) for g in names), dtype=np.int64,
                             count=len(names))
        is_ann = row_of >= 0
        # edges from an unannotated row to an annotated neighbour; row numbers
        # are found only for those (a per-edge row array is the largest temporary)
        m = is_ann[indices]
        m &= np.repeat(~is_ann, np.diff(indptr))
        pos = np.flatnonzero(m)
        del m
        rows = np.searchsorted(indptr, pos, side="right") - 1
        nbr = indices[pos]
        del pos
        starts = np.flatnonzero(np.r_[True, np.diff(rows) != 0]) if rows.size else rows
        ptr = np.r_[starts, rows.size].astype(np.int64)
        idx = row_of[nbr]  # adjacency rows hold distinct neighbours
        in_bg = np.zeros(len(gt.genes), dtype=bool)
        in_bg[idx] = True
        keys = [names[r] for r in rows[starts].tolist()]
        return cls(keys, ptr, idx, np.diff(ptr), in_bg, int(in_bg.sum()))


def _iter_enrichment(queries, background, gene_terms, pval, chunk_size=2000,
                     term_info=None, threads=1):
    """:func:`_iter_encoded` for dict queries and gene terms."""
    gt = _GeneTerms.from_dict(gene_terms)
    yield from _iter_encoded(_Queries.from_dict(queries, background, gt), gt, pval,
                             chunk_size, term_info, threads)


def _iter_encoded(q: _Queries, gt: _GeneTerms, pval, chunk_size=2000,
                  term_info=None, threads=1):
    """Yield one :class:`_Chunk` per ``chunk_size`` queries.

    Study counts come from the C++ ``study_counts`` kernel (scipy sparse
    product without the extension); nothing here needs pandas.
    """
    core = extension()
    row_ptr, g_idx, term_ids = gt.row_ptr, gt.g_idx, gt.term_ids
    n_terms = len(term_ids)
    pop_counts = np.bincount(g_idx[np.repeat(q.in_bg, np.diff(row_ptr))], minlength=n_terms
                             ).astype(np.int64)
    pop_n = q.pop_n
    n_tests = int(np.count_nonzero(pop_counts))  # BH over the terms of the background

    # per-term name, namespace and v1 namespace rank, looked up once
    info = term_info or {}
    term_name = np.array([info.get(t, ("", ""))[0] for t in term_ids], dtype=object)
    term_ns = np.array([info.get(t, ("", ""))[1] for t in term_ids], dtype=object)
    term_rank = np.array([_NAMESPACE_RANK.get(ns, 3) for ns in term_ns], dtype=np.int64)
    term_order = np.empty(n_terms, dtype=np.int64)  # rank of term id (tie-break)
    term_order[np.argsort(term_ids.astype(str), kind="stable")] = np.arange(n_terms)
    terms_meta = (term_ids, term_name, term_ns)

    # the fused kernel needs study_n <= pop_n (true for annotate queries; a
    # gene list passed to enrich() may be larger than its background)
    fused = (core is not None and hasattr(core.enrichment, "significant_counts")
             and int(q.study_n.max(initial=0)) <= pop_n and pop_n < 2_000_000_000)
    G = None
    if core is None:
        from scipy.sparse import csr_matrix

        G = csr_matrix((np.ones(len(g_idx), dtype=np.int32), g_idx, row_ptr),
                       shape=(len(gt.genes), n_terms))

    for start in range(0, len(q.keys), chunk_size):
        chunk = q.keys[start:start + chunk_size]
        stop = start + len(chunk)
        study_n = q.study_n[start:stop]
        q_ptr = q.ptr[start:stop + 1] - q.ptr[start]
        q_idx = q.idx[q.ptr[start]:q.ptr[stop]]
        if fused:
            # counts, p-values and the cutoff in one pass: only significant
            # entries are materialised (a few % of all query-term pairs)
            rows, terms, counts, p, fdr, has_hits = core.enrichment.significant_counts(
                q_ptr, q_idx, row_ptr, g_idx, n_terms, study_n, pop_counts, pop_n,
                float(pval), max(1, int(threads)))
            sn, pc = study_n[rows], pop_counts[terms]
        else:
            if core is not None:
                rows, terms, counts = core.enrichment.study_counts(
                    q_ptr, q_idx, row_ptr, g_idx, n_terms, max(1, int(threads)))
            else:
                from scipy.sparse import csr_matrix

                Q = csr_matrix((np.ones(q_ptr[-1], dtype=np.int32), q_idx, q_ptr),
                               shape=(len(chunk), len(gt.genes)))
                S = (Q @ G).tocoo()
                rows, terms, counts = S.row, S.col, S.data.astype(np.int64)
            hit = (counts > 0) & (pop_counts[terms] > 0)
            has_hits = np.bincount(rows[hit], minlength=len(chunk)) > 0
            rows, terms, counts = rows[hit], terms[hit], counts[hit]

            sn = study_n[rows]
            pc = pop_counts[terms]
            p = _hypergeom_sf_vec(counts, pop_n, pc, sn, threads)
            fdr = _bh_by_row(rows, p, n_tests)
            keep = p <= pval
            rows, terms, counts, sn, pc, p, fdr = (a[keep] for a in
                                                   (rows, terms, counts, sn, pc, p, fdr))
        # sort by query, namespace rank, p-value, term id (v1 order)
        order = np.lexsort((term_order[terms], p, term_rank[terms], rows))
        rows, terms, counts, sn, pc, p, fdr = (a[order] for a in
                                               (rows, terms, counts, sn, pc, p, fdr))
        bounds = np.searchsorted(rows, np.arange(len(chunk) + 1))
        yield _Chunk(chunk, bounds, has_hits, terms, counts, sn, pc, pop_n, p, fdr,
                     terms_meta)


def _bh_by_row(rows: np.ndarray, p: np.ndarray, m: int) -> np.ndarray:
    """Benjamini-Hochberg adjusted p of each entry within its query (``rows``),
    over ``m`` tests: the entries hold every p < 1 of a query, the other
    terms of the m have p = 1 (as the C++ ``significant_counts``)."""
    fdr = np.ones(len(p))
    if not len(p):
        return fdr
    order = np.lexsort((np.arange(len(p)), p, rows))
    r, ps = rows[order], p[order]
    starts = np.flatnonzero(np.r_[True, r[1:] != r[:-1]])
    ends = np.r_[starts[1:], len(r)]
    for a, b in zip(starts, ends):
        q = ps[a:b] * m / np.arange(1, b - a + 1)
        fdr[order[a:b]] = np.minimum(np.minimum.accumulate(q[::-1])[::-1], 1.0)
    return fdr


def _log_binom(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    from scipy.special import gammaln

    return gammaln(a + 1) - gammaln(b + 1) - gammaln(a - b + 1)


def _hypergeom_sf_vec(k, M, n, N, threads: int = 1) -> np.ndarray:
    """Vectorised P(X >= k) for X ~ Hypergeom(M total, n successes, N draws).

    Uses the compiled kernel ``_gcen_core.enrichment.hypergeom_sf`` (same
    algorithm, multithreaded) when available, else :func:`_hypergeom_sf_numpy`.
    """
    k, M, n, N = np.broadcast_arrays(*(np.asarray(v, dtype=np.float64) for v in (k, M, n, N)))
    core = extension()
    if core is None or k.ndim != 1:
        return _hypergeom_sf_numpy(k, M, n, N)
    return core.enrichment.hypergeom_sf(
        np.ascontiguousarray(k), np.ascontiguousarray(M), np.ascontiguousarray(n),
        np.ascontiguousarray(N), max(1, int(threads)))


def _hypergeom_sf_numpy(k, M, n, N) -> np.ndarray:
    """Vectorised P(X >= k) for X ~ Hypergeom(M total, n successes, N draws).

    Same quantity as ``scipy.stats.hypergeom.sf(k - 1, M, n, N)`` (which sums
    the pmf in a Python loop per element and dominates runtime for millions
    of tests). Each tail is summed away from the mode with the pmf ratio
    recurrence, so the starting term only underflows when the whole tail is
    negligible:

    * k > mode : p = sum_{x=k}^{xmax} pmf(x)
    * k <= mode: p = 1 - sum_{x=xmin}^{k-1} pmf(x)
    """
    k, n, N = (np.asarray(v, dtype=np.float64) for v in (k, n, N))
    M = np.broadcast_to(np.asarray(M, dtype=np.float64), k.shape)
    k, n, N = np.broadcast_arrays(k, n, N)
    xmin = np.maximum(0.0, N - (M - n))
    xmax = np.minimum(n, N)
    mode = np.floor((N + 1) * (n + 1) / (M + 2))

    out = np.zeros(k.shape, dtype=np.float64)
    out[k <= xmin] = 1.0
    work = (k > xmin) & (k <= xmax)

    def _tail(idx, start, upward):
        x = start[idx].copy()
        nn, NN, MM = n[idx], N[idx], M[idx]
        term = np.exp(_log_binom(nn, x) + _log_binom(MM - nn, NN - x) - _log_binom(MM, NN))
        total = term.copy()
        stop = xmax[idx] if upward else xmin[idx]
        live = (x < stop) if upward else (x > stop)
        while live.any():
            s = np.flatnonzero(live)
            xs = x[s]
            if upward:
                ratio = (nn[s] - xs) * (NN[s] - xs) / ((xs + 1) * (MM[s] - nn[s] - NN[s] + xs + 1))
                x[s] = xs + 1
            else:
                ratio = xs * (MM[s] - nn[s] - NN[s] + xs) / ((nn[s] - xs + 1) * (NN[s] - xs + 1))
                x[s] = xs - 1
            term[s] *= ratio
            total[s] += term[s]
            more = (x[s] < stop[s]) if upward else (x[s] > stop[s])
            live[s] = more & (term[s] > total[s] * 1e-17)
        return total

    upper = np.flatnonzero(work & (k > mode))
    if upper.size:
        out[upper] = np.minimum(_tail(upper, k, upward=True), 1.0)
    lower = np.flatnonzero(work & (k <= mode))
    if lower.size:
        out[lower] = np.clip(1.0 - _tail(lower, k - 1, upward=False), 0.0, 1.0)
    return out
