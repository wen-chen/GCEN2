"""Core domain data models for GCEN2."""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING, Optional

import numpy as np

# pandas is imported inside the methods that need it, so that code paths
# working on numpy arrays (gcen2 build/module/rwr/annotate) do not load it.
if TYPE_CHECKING:
    import networkx as nx
    import pandas as pd
    import scipy.sparse


@dataclass
class ExpressionMatrix:
    """Container for gene expression data (bulk or single-cell)."""

    data: pd.DataFrame | any  # pd.DataFrame for bulk, AnnData for single-cell
    gene_names: np.ndarray
    sample_names: np.ndarray
    is_single_cell: bool = False
    is_normalized: bool = False
    normalization_method: Optional[str] = None


@dataclass
class Network:
    """Gene co-expression network represented as an edge list."""

    edges: pd.DataFrame  # columns: gene_a, gene_b, correlation, p_value[, fdr][, tom]
    node_count: int = 0
    edge_count: int = 0
    is_signed: bool = False  # True: negative correlations excluded (GCEN v1 -s)
    method: str = "spearman"  # pearson | spearman
    metadata: dict = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.edge_count == 0:
            self.edge_count = len(self.edges)
        if self.node_count == 0 and not self.edges.empty:
            genes = set(self.edges["gene_a"]) | set(self.edges["gene_b"])
            self.node_count = len(genes)

    def to_networkx(self) -> nx.Graph:
        import networkx as nx

        G = nx.from_pandas_edgelist(
            self.edges,
            source="gene_a",
            target="gene_b",
            edge_attr=True,
        )
        return G

    def to_adjacency(self, weight: Optional[str] = "correlation", absolute: bool = False,
                     genes: Optional[list[str]] = None) -> scipy.sparse.csr_array:
        """Weighted adjacency matrix (symmetric, scipy ``csr_array``).

        weight:
            Edge column used as the weight: ``"correlation"`` (default; the
            sign is kept), another numeric column such as ``"tom"``, or None
            for an unweighted (0/1) matrix.
        absolute:
            Use the absolute values (``|correlation|``, the weights of
            ``gcen2.rwr.predict(weighted_network=True)``).
        genes:
            Order of the rows and columns (genes not in the network give empty
            rows, edges of other genes are left out); default: the genes in
            order of first appearance in the edge table (the node order of
            :meth:`to_networkx`).

        Repeated edges between two genes are summed.
        """
        import pandas as pd
        from scipy.sparse import coo_array

        from gcen2.exceptions import GCENInputError

        a = self.edges["gene_a"].to_numpy()
        b = self.edges["gene_b"].to_numpy()
        if genes is None:
            genes = pd.unique(np.column_stack([a, b]).ravel()) if len(a) else []
        if weight is None:
            w = np.ones(len(a), dtype=np.float64)
        else:
            if weight not in self.edges.columns:
                raise GCENInputError(f"no edge column '{weight}' for the weights "
                                     f"(columns: {', '.join(map(str, self.edges.columns))}).")
            w = pd.to_numeric(self.edges[weight], errors="coerce").to_numpy(dtype=np.float64)
            if absolute:
                w = np.abs(w)
        cat = pd.CategoricalDtype(categories=list(genes))
        i = pd.Series(a).astype(cat).cat.codes.to_numpy()
        j = pd.Series(b).astype(cat).cat.codes.to_numpy()
        keep = (i >= 0) & (j >= 0)
        i, j, w = i[keep], j[keep], w[keep]
        off = i != j  # a self-loop is one entry, not two
        n = len(genes)
        return coo_array((np.concatenate([w, w[off]]), (np.concatenate([i, j[off]]),
                                                         np.concatenate([j, i[off]]))),
                         shape=(n, n)).tocsr()

    def to_edgelist(self, path: str) -> None:
        """Write the edges as a tab-separated file with a header line.

        Standard networks (``gene_a gene_b correlation p_value [fdr] [tom]``)
        are written like GCEN v1 (correlation ``%f``, p-value/FDR/TOM ``%e``, i.e. 6
        decimals) by the compiled formatter (pure-Python formatting if the
        extension is missing); other column layouts are written by pandas.
        """
        cols = list(self.edges.columns)
        std = ["gene_a", "gene_b", "correlation", "p_value"]
        if cols in (std, std + ["fdr"], std + ["tom"], std + ["fdr", "tom"]):
            import pandas as pd

            from gcen2.core import extension

            core = extension()
            if core is not None:
                e = self.edges
                codes, names = pd.factorize(
                    pd.concat([e["gene_a"], e["gene_b"]], ignore_index=True).astype(str))
                m = len(e)
                idx = np.stack([codes[:m], codes[m:]], axis=1).astype(np.int32)
                fdr = e["fdr"].to_numpy(dtype=float) if "fdr" in cols else None
                tom = e["tom"].to_numpy(dtype=float) if "tom" in cols else None
                writer = core.correlation.EdgeFormatter(list(names))
                with open(path, "wb") as fh:
                    fh.write(("\t".join(cols) + "\n").encode())
                    if m:
                        fh.write(writer.format(idx, e["correlation"].to_numpy(dtype=float),
                                               e["p_value"].to_numpy(dtype=float), fdr, 4, tom))
                return
            fmt = {"correlation": "{:f}".format, "p_value": "{:e}".format, "fdr": "{:e}".format,
                   "tom": "{:e}".format}
            self.edges.assign(**{c: self.edges[c].map(f) for c, f in fmt.items() if c in cols}
                              ).to_csv(path, sep="\t", index=False, lineterminator="\n")
            return
        self.edges.to_csv(path, sep="\t", index=False, lineterminator="\n")

    def neighbor_sets(self) -> dict[str, set[str]]:
        """Undirected adjacency as ``{gene: set(neighbors)}``.

        Vectorised replacement for row-wise ``iterrows`` loops; a 3.6M-edge
        network takes seconds instead of minutes.
        """
        if self.edges.empty:
            return {}
        a = self.edges["gene_a"].to_numpy()
        b = self.edges["gene_b"].to_numpy()
        src = np.concatenate([a, b])
        dst = np.concatenate([b, a])
        import pandas as pd

        codes, uniques = pd.factorize(src)
        order = np.argsort(codes, kind="stable")
        bounds = np.flatnonzero(np.diff(codes[order])) + 1
        groups = np.split(dst[order], bounds)
        keys = uniques[codes[order][np.r_[0, bounds]]]
        return {k: set(g.tolist()) for k, g in zip(keys.tolist(), groups)}

    def sparse_adjacency(
        self,
        genes: Optional[list[str]] = None,
        weighted: bool = False,
        directed: bool = False,
    ) -> scipy.sparse.csr_matrix:
        """Sparse adjacency matrix ``W`` with ``W[i, j]`` for edge gene_a(i) -> gene_b(j).

        Rows/columns follow ``genes`` (default: sorted node names). Weights are
        ``|correlation|`` when ``weighted`` else 1; parallel edges are summed.
        Undirected networks get both ``W[i, j]`` and ``W[j, i]``.
        """
        from scipy.sparse import coo_matrix

        if genes is None:
            genes = sorted(set(self.edges["gene_a"]) | set(self.edges["gene_b"]))
        import pandas as pd

        cat = pd.CategoricalDtype(categories=genes)
        i = pd.Series(self.edges["gene_a"]).astype(cat).cat.codes.to_numpy()
        j = pd.Series(self.edges["gene_b"]).astype(cat).cat.codes.to_numpy()
        w = (
            self.edges["correlation"].abs().to_numpy(dtype=np.float64)
            if weighted
            else np.ones(len(i), dtype=np.float64)
        )
        keep = (i >= 0) & (j >= 0)
        i, j, w = i[keep], j[keep], w[keep]
        if not directed:
            i, j, w = np.concatenate([i, j]), np.concatenate([j, i]), np.concatenate([w, w])
        n = len(genes)
        return coo_matrix((w, (i, j)), shape=(n, n)).tocsr()

    def subnetwork(self, genes: list[str]) -> Network:
        gene_set = set(genes)
        mask = self.edges["gene_a"].isin(gene_set) | self.edges["gene_b"].isin(gene_set)
        sub_edges = self.edges[mask].reset_index(drop=True)
        return Network(edges=sub_edges, is_signed=self.is_signed, method=self.method)

    def stats(self) -> NetworkStats:
        return NetworkStats.compute(self)


@dataclass
class NetworkStats:
    """Computed statistics for a Network."""

    node_count: int
    edge_count: int
    density: float
    average_degree: float
    max_degree: int
    connected_components: int
    giant_component_size: int
    heterogeneity: float
    clustering_coefficient: float
    transitivity: float

    @classmethod
    def compute(cls, network: Network, chunk_size: int = 2000,
                threads: int = 2) -> NetworkStats:
        """Compute statistics with sparse-matrix operations.

        Numerically identical to the networkx definitions (undirected simple
        graph; parallel/reverse duplicate edges collapse; a self-loop adds 2 to
        the degree but is ignored for clustering/transitivity), but scales to
        millions of edges. With the compiled extension, degrees, triangles and
        components come from ``_gcen_core.network.graph_stats`` on the
        adjacency of ``module._unique_adjacency`` (nodes in sorted order, as
        :func:`gcen2.network.stats_from_file`, so both give the same values);
        otherwise triangles are counted in row chunks of the sparse adjacency.
        """
        edges = network.edges
        if edges.empty:
            return cls(0, 0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0, 0.0)
        from gcen2.core import extension

        if extension() is not None:
            from gcen2.module import _unique_adjacency

            A, _ = _unique_adjacency(network)
            return cls.from_adjacency(A.indptr.astype(np.int64), A.indices.astype(np.int32),
                                      threads)
        return cls._compute_sparse(network, chunk_size)

    @classmethod
    def from_adjacency(cls, indptr, indices, threads: int = 2) -> NetworkStats:
        """Statistics of the symmetric CSR adjacency ``(indptr, indices)``
        (rows sorted and unique, a self-loop is a diagonal entry), as given by
        ``module._unique_adjacency`` or ``module.read_file_adjacency``."""
        return cls.with_degrees(indptr, indices, threads)[0]

    @classmethod
    def with_degrees(cls, indptr, indices, threads: int = 2):
        """:meth:`from_adjacency` and the node degrees (networkx ``G.degree``
        values: a self-loop adds 2), in node order, as an int64 array."""
        from gcen2.core import require_extension

        n = len(indptr) - 1
        if n <= 0:
            return cls(0, 0, 0.0, 0.0, 0, 0, 0, 0.0, 0.0, 0.0), np.zeros(0, dtype=np.int64)
        g = require_extension().network.graph_stats(indptr, indices, max(1, int(threads)))
        k = np.asarray(g["degree_no_loop"])
        degrees = k + 2 * np.asarray(g["self_loop"])
        tri = np.asarray(g["triangles"]).astype(np.float64)
        stats = cls._finish(n, int(g["edge_count"]), k, degrees, tri, int(g["components"]),
                            int(g["giant_component"]))
        return stats, degrees

    @classmethod
    def _finish(cls, n, e, k, degrees, tri, n_comp, giant) -> NetworkStats:
        pairs_at_node = k * (k - 1) / 2.0
        with np.errstate(divide="ignore", invalid="ignore"):
            local_cc = np.where(k >= 2, tri / pairs_at_node, 0.0)
        triads = pairs_at_node.sum()
        transitivity = float(tri.sum() / triads) if tri.sum() > 0 and triads > 0 else 0.0
        avg_deg = float(degrees.mean())
        std_deg = float(degrees.std())
        return cls(
            node_count=n,
            edge_count=e,
            density=(2 * e) / (n * (n - 1)) if n > 1 else 0.0,
            average_degree=avg_deg,
            max_degree=int(degrees.max()),
            connected_components=n_comp,
            giant_component_size=giant,
            heterogeneity=std_deg / avg_deg if avg_deg > 0 else 0.0,
            clustering_coefficient=float(local_cc.mean()),
            transitivity=transitivity,
        )

    @classmethod
    def _compute_sparse(cls, network: Network, chunk_size: int = 2000) -> NetworkStats:
        """scipy implementation (no compiled extension)."""
        from scipy.sparse import coo_matrix
        from scipy.sparse.csgraph import connected_components

        edges = network.edges

        import pandas as pd

        codes, uniques = pd.factorize(
            np.concatenate([edges["gene_a"].to_numpy(), edges["gene_b"].to_numpy()])
        )
        m = len(edges)
        a, b = codes[:m], codes[m:]
        n = len(uniques)
        pairs = np.unique(np.stack([np.minimum(a, b), np.maximum(a, b)], axis=1), axis=0)
        u, v = pairs[:, 0], pairs[:, 1]
        e = len(pairs)
        loop = u == v
        u_off, v_off = u[~loop], v[~loop]

        A = coo_matrix(
            (np.ones(2 * len(u_off), dtype=np.float64),
             (np.concatenate([u_off, v_off]), np.concatenate([v_off, u_off]))),
            shape=(n, n),
        ).tocsr()

        k = np.asarray(A.sum(axis=1)).ravel().astype(np.int64)  # degree w/o self-loops
        degrees = k + 2 * np.bincount(u[loop], minlength=n)

        # triangles through each node: T_i = (A^3)_ii / 2, computed row-chunked
        tri = np.zeros(n, dtype=np.float64)
        for start in range(0, n, chunk_size):
            rows = A[start:start + chunk_size]
            tri[start:start + chunk_size] = np.asarray(
                (rows @ A).multiply(rows).sum(axis=1)
            ).ravel() / 2.0

        n_comp, labels = connected_components(A, directed=False)
        return cls._finish(n, e, k, degrees, tri, int(n_comp), int(np.bincount(labels).max()))

@dataclass
class ModuleSet:
    """Collection of gene modules identified from a network."""

    modules: list[set[str]]
    module_count: int = 0
    source_network: Optional[Network] = None
    similarity_cutoff: float = 0.5

    def __post_init__(self) -> None:
        if self.module_count == 0:
            self.module_count = len(self.modules)


@dataclass
class EnrichmentResult:
    """Results of GO or KEGG enrichment analysis."""

    # columns: term_id, name, namespace, enrichment, study_count,
    #          study_n, pop_count, pop_n, p_value
    results: pd.DataFrame
    ontology_type: str  # "GO" or "KEGG"
    gene_or_module: str
    p_cutoff: float = 0.05


@dataclass
class RWRResult:
    """Results of Random Walk with Restart gene function prediction."""

    ranked_genes: pd.DataFrame  # columns: gene, probability
    seed_genes: list[str]
    gamma: float = 0.5
    steps: int = 0


@dataclass
class DiffCoExResult:
    """Results of differential co-expression analysis between two conditions."""

    # columns: gene_a, gene_b, cor_a, cor_b, z_score, p_value, fdr, category
    results: pd.DataFrame
    network_a_name: str = "A"
    network_b_name: str = "B"
    method: str = "fisher_z"
    p_cutoff: float = 0.05
