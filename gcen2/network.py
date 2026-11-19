"""Network construction and utility operations."""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Union

import numpy as np

from gcen2._build_file import (  # noqa: F401  (METHODS, KERNELS, build_from_file re-exported)
    _FIT_COLUMNS,
    KERNELS,
    METHODS,
    _build_result,
    _check_kernel,
    _check_method,
    _tom_metadata,
    _tom_power,
    build_from_file,
)
from gcen2._matrix_io import require_file, tsv_header
from gcen2.core import extension, require_extension
from gcen2.exceptions import GCENInputError, GCENIOError

# pandas and gcen2.models (which imports pandas) are imported inside the
# functions that need them: ``gcen2 build`` streams a numpy matrix to disk and
# would otherwise pay ~50 MB and seconds of start-up for pandas alone.
if TYPE_CHECKING:
    import pandas as pd

    from gcen2.models import ExpressionMatrix, Network, NetworkStats


#: Samples from which ``kernel="auto"`` uses the matrix-product kernel.
GEMM_MIN_SAMPLES = 64


def build(
    expr: Union[ExpressionMatrix, pd.DataFrame],
    method: str = "spearman",
    cor: float = 0.1,
    pval: float = 0.001,
    signed: bool = False,
    fdr: bool = False,
    log: str | None = None,
    all_edges: bool = False,
    threads: int = 2,
    layer: str | None = None,
    use_highly_variable: bool = False,
    tom: bool = False,
    power: float | str = "auto",
    tom_min: float = 0.0,
    kernel: str = "auto",
) -> Network:
    """Build a gene co-expression network.

    Parameters
    ----------
    expr:
        Expression matrix (genes × samples), or an AnnData (cells × genes,
        dense or sparse; see :func:`gcen2.single_cell.adata_matrix`).
    method:
        Correlation method: 'pearson', 'spearman', or a dropout-aware Pearson
        (see :data:`METHODS`): 'pearson_nonzero' (per gene pair, only the
        samples where both genes are non-zero) or
        'pearson_pairwise_complete' (where at least one is non-zero).
    cor:
        Minimum absolute correlation coefficient for an edge to be retained.
    pval:
        Maximum p-value cutoff.
    signed:
        If True, build a signed network: only positive correlations are kept
        (same as GCEN v1 ``network_build -s`` and WGCNA's "signed" network).
        If False (default), edges are selected by ``|r|`` and both signs are
        kept.
    fdr:
        If True, apply Benjamini-Hochberg FDR correction.
    log:
        Optional log transformation before correlation: 'log', 'log2', 'log10'.
    all_edges:
        If True, return all edges (ignoring cor and pval cutoffs).
    threads:
        Number of threads for the C++ correlation kernel.
    layer, use_highly_variable:
        AnnData input only: the layer instead of ``X``, and keep the highly
        variable genes.
    tom, power, tom_min:
        Add WGCNA's topological overlap as a ``tom`` column (pearson or
        spearman only; see :data:`TOM_DOC`): ``power`` is the soft
        threshold, ``"auto"`` picks it as WGCNA's pickSoftThreshold. Edges
        with TOM below ``tom_min`` are dropped (unless ``all_edges``). The
        power and the scale-free fit table go to ``network.metadata["tom"]``.
    kernel:
        Pearson / Spearman kernel: ``"pairwise"`` (each pair summed in turn,
        exactly as GCEN v1), ``"gemm"`` (rows standardised once, r from a
        tiled matrix product with a fixed summation order; about 4x faster at
        128 samples and 10x at 512 and more; r agrees with ``"pairwise"`` to
        rounding, about 1e-15 to 1e-14, so an edge exactly on a cutoff can
        differ) or
        ``"auto"`` (``"gemm"`` from :data:`GEMM_MIN_SAMPLES` samples). Both
        give the same results for any number of threads. Ignored by the
        dropout-aware methods.

    Returns
    -------
    Network
    """
    import pandas as pd

    from gcen2.models import Network

    matrix, gene_names = _prepare_matrix(expr, method, log, layer=layer,
                                         use_highly_variable=use_highly_variable)
    gene_names = np.asarray(gene_names, dtype=str)
    tom_power = _tom_power(tom, power, method)
    _check_kernel(kernel)
    indices, correlations, p_values, q, t, info = require_extension().correlation.network_edges(
        matrix, method, threads, cor, pval, signed, fdr, all_edges, tom_power, tom_min,
        kernel=kernel)
    edges = pd.DataFrame(
        {
            "gene_a": gene_names[indices[:, 0]],
            "gene_b": gene_names[indices[:, 1]],
            "correlation": correlations,
            "p_value": p_values,
        }
    )
    if fdr:
        edges["fdr"] = q
    metadata = {}
    if tom:
        edges["tom"] = t
        metadata["tom"] = _tom_metadata(info)

    return Network(edges=edges, is_signed=signed, method=method, metadata=metadata)


def build_to_file(
    expr: Union[ExpressionMatrix, pd.DataFrame, np.ndarray],
    path,
    method: str = "spearman",
    cor: float = 0.1,
    pval: float = 0.001,
    signed: bool = False,
    fdr: bool = False,
    log: str | None = None,
    all_edges: bool = False,
    threads: int = 2,
    max_pairs_per_block: int = 500_000,
    gene_names=None,
    layer: str | None = None,
    use_highly_variable: bool = False,
    tom: bool = False,
    power: float | str = "auto",
    tom_min: float = 0.0,
    kernel: str = "auto",
) -> dict:
    """Build a network and stream it to ``path`` block by block.

    Same edges, columns, order and number formatting as
    ``build(...).to_edgelist(path)`` (the files are byte-identical; numbers
    are written like GCEN v1: correlation ``%f``, p-value/FDR ``%e``), but the
    gene pairs are computed in blocks of rows of at most
    ``max_pairs_per_block`` pairs and written immediately, so peak memory is
    bounded by one block instead of the whole network (GCEN v1 also streams
    its output). This matters for large or ``all_edges`` networks.

    With ``fdr=True`` the Benjamini-Hochberg ranks need every p-value first:
    one pass keeps the p-values of all tests with p <= ``pval`` (8 bytes
    each) and the edges that pass the other filters (24 bytes each, up to
    2^24 of them; beyond that the blocks are correlated a second time to
    write them). All of it runs in the compiled extension
    (``_gcen_core.correlation.build_network_array``, src/network_stream.cpp).

    ``expr`` may also be a plain genes x samples array, with the gene names
    in ``gene_names`` (this path does not need pandas), or an AnnData (with
    ``layer`` / ``use_highly_variable`` as in :func:`build`).

    ``tom`` / ``power`` / ``tom_min`` and ``kernel`` as in :func:`build`
    (TOM needs the whole adjacency matrix in memory: 4 bytes x genes^2).

    Returns ``{"node_count": ..., "edge_count": ...}``, plus ``"tom"`` (power
    and fit table, see :func:`build`) with ``tom=True``.
    """
    matrix, gene_names = _prepare_matrix(expr, method, log, gene_names, layer=layer,
                                         use_highly_variable=use_highly_variable)
    tom_power = _tom_power(tom, power, method)
    _check_kernel(kernel)
    nodes, edges, info = require_extension().correlation.build_network_array(
        matrix, gene_names, str(path), method, cor, pval, signed, fdr, all_edges, threads,
        max_pairs_per_block, tom_power=tom_power, tom_min=tom_min, kernel=kernel)
    return _build_result(nodes, edges, info)


#: How the ``tom`` edge score is computed (src/tom.hpp).
TOM_DOC = """\
Adjacency a_ij = ((1 + r_ij) / 2)^power for a signed network, |r_ij|^power
otherwise (r: Pearson, or Pearson of the ranks for spearman), connectivity
k_i = sum_u a_iu, and TOM_ij = (sum_u a_iu a_uj + a_ij) / (min(k_i, k_j) + 1 - a_ij),
as WGCNA's TOMsimilarityFromExpr (TOMType = networkType, TOMDenom = "min").
power="auto": the first of 1..10, 12..20 whose scale-free fit R^2 exceeds
0.85 (WGCNA 1.74 pickSoftThreshold), else the WGCNA FAQ power for the number
of samples (signed 18/16/14/12, unsigned 9/8/7/6 below 20/30/40 samples and
above). Single precision for the adjacency: relative error ~1e-6."""


def pick_soft_threshold(
    expr: Union[ExpressionMatrix, pd.DataFrame, np.ndarray],
    method: str = "pearson",
    signed: bool = False,
    powers=None,
    log: str | None = None,
    threads: int = 2,
    gene_names=None,
) -> pd.DataFrame:
    """WGCNA's pickSoftThreshold table: for each power the scale-free fit
    R^2 and slope and the mean / median / max connectivity (see
    :data:`TOM_DOC`), with the power ``power="auto"`` would pick marked in
    the ``chosen`` column."""
    import pandas as pd

    if method not in ("pearson", "spearman"):
        raise GCENInputError("pick_soft_threshold needs method 'pearson' or 'spearman'.")
    if gene_names is None and isinstance(expr, np.ndarray):
        gene_names = [str(i) for i in range(len(expr))]
    matrix, _ = _prepare_matrix(expr, method, log, gene_names)
    fits = require_extension().correlation.soft_threshold(
        matrix, method, signed, [float(p) for p in (powers or [])], threads)
    df = pd.DataFrame(fits, columns=list(_FIT_COLUMNS))
    above = np.flatnonzero(df["r2"].to_numpy() > 0.85)
    df["chosen"] = False
    if len(above):
        df.loc[above[0], "chosen"] = True
    return df


_LOG_FUNCS = {"log": np.log, "log2": np.log2, "log10": np.log10}


def _prepare_matrix(expr, method: str, log: str | None, gene_names=None, layer=None,
                    use_highly_variable=False):
    """Validate input; return a C-contiguous float64 matrix and gene names."""
    _check_method(method)

    from gcen2._anndata import adata_matrix, is_anndata

    inner = getattr(expr, "data", None) if type(expr).__name__ == "ExpressionMatrix" else expr
    if is_anndata(inner):
        expr, gene_names = adata_matrix(inner, layer=layer,
                                        use_highly_variable=use_highly_variable)
    elif layer is not None or use_highly_variable:
        raise GCENInputError("layer and use_highly_variable apply to AnnData input only.")

    if isinstance(expr, np.ndarray):
        if expr.ndim != 2:
            raise GCENInputError("expr must be a 2-D genes x samples array.")
        if gene_names is None or len(gene_names) != expr.shape[0]:
            raise GCENInputError("gene_names must give one name per row of expr.")
        if expr.shape[1] < 3:
            raise GCENInputError("At least 3 samples are required for correlation computation.")
        matrix = np.asarray(expr, dtype=np.float64)
        if log is not None:
            if log not in _LOG_FUNCS:
                raise GCENInputError(
                    f"Unknown log base '{log}'. Choose from: {', '.join(_LOG_FUNCS)}")
            matrix = _LOG_FUNCS[log](matrix + 1.0)  # as preprocessing.log_transform
        return np.ascontiguousarray(matrix), list(map(str, gene_names))

    import pandas as pd

    from gcen2.models import ExpressionMatrix
    from gcen2.preprocessing import log_transform

    df = expr.data if isinstance(expr, ExpressionMatrix) else expr

    if not isinstance(df, pd.DataFrame):
        raise GCENInputError(
            "expr must be a pd.DataFrame or ExpressionMatrix wrapping a DataFrame.")
    if df.shape[1] < 3:
        raise GCENInputError("At least 3 samples are required for correlation computation.")

    if log is not None:
        df = log_transform(df, base=log)

    # The C++ kernel reads rows (genes) as contiguous blocks. DataFrames built
    # column-wise (e.g. pd.concat of Series) yield Fortran-ordered arrays, so
    # force C order explicitly.
    matrix = np.ascontiguousarray(df.to_numpy(dtype=np.float64))
    return matrix, [str(g) for g in df.index]


def _final_mask(r, p, q, cor, pval, signed) -> np.ndarray:
    mask = (np.abs(r) >= cor) & (p <= pval)
    if q is not None:
        mask &= q <= pval
    if signed:
        mask &= r >= 0
    return mask


def merge(networks: list[Network], cor_cutoff: float = 0.5) -> Network:
    """Merge multiple networks.

    Edges present in multiple networks are aggregated by summing their correlations.
    The direction (sign) is output as a separate column, fixing the GCEN v1 bug.

    Parameters
    ----------
    networks:
        List of Network objects to merge.
    cor_cutoff:
        Minimum absolute aggregated correlation for an edge to appear in the output.
    """
    import pandas as pd

    from gcen2.models import Network

    if not networks:
        raise GCENInputError("At least one network is required for merging.")

    all_edges = pd.concat([n.edges[["gene_a", "gene_b", "correlation"]] for n in networks])
    grouped = all_edges.groupby(["gene_a", "gene_b"])["correlation"].sum().reset_index()
    grouped["direction"] = np.sign(grouped["correlation"]).astype(int)
    mask = grouped["correlation"].abs() >= cor_cutoff
    result = grouped[mask].reset_index(drop=True)
    return Network(edges=result, is_signed=all(n.is_signed for n in networks), method="merged")


def extract(network: Network, genes: list[str]) -> Network:
    """Extract a subnetwork containing edges where at least one node is in ``genes``."""
    return network.subnetwork(genes)


def extract_file(path, genes: list[str], output) -> dict:
    """:func:`extract` from an edge-list file to an edge-list file.

    Writes the same file as ``extract(read_network(path),
    genes).to_edgelist(output)``. The compiled extension streams the file,
    keeps the rows with a listed gene and writes them
    (``_gcen_core.network.extract_edge_file``), so pandas is not imported.
    It falls back to that pipeline without the extension, for gene columns
    that are not the first two, and when it cannot be sure how pandas would
    read the file (it then writes nothing). Returns ``{"node_count",
    "edge_count"}``.
    """
    path = require_file(path)
    core = extension()
    if core is not None:
        from gcen2._matrix_io import network_file_columns

        table = _edge_table_args(*network_file_columns(path))
        if table is not None:
            try:
                r = core.network.extract_edge_file(str(path), *table, str(output),
                                                   [str(g) for g in genes])
            except (RuntimeError, UnicodeDecodeError) as exc:
                raise GCENIOError(f"Failed to read network file {path}: {exc}") from exc
            if r["ok"]:
                return {"node_count": r["node_count"], "edge_count": r["edge_count"]}
    from gcen2.io import read_network

    sub = extract(read_network(path), genes)
    sub.to_edgelist(output)
    return {"node_count": sub.node_count, "edge_count": sub.edge_count}


def _edge_table_args(names: list[str], has_header: bool):
    """``(skip_header, n_columns, v1_layout, header)`` for the extension's
    edge-table writers (``shuffle_edge_file`` / ``extract_edge_file``), or
    None when the table is not theirs to write: gene columns not first, or
    repeated column names (pandas refuses those)."""
    if names[:2] != ["gene_a", "gene_b"] or len(set(names)) != len(names):
        return None
    std = ["gene_a", "gene_b", "correlation", "p_value"]
    v1 = names in (std, std + ["fdr"])
    header = "\t".join(names) + "\n" if v1 else tsv_header(names)
    return has_header, len(names), v1, header


def merge_files(paths: list, output, cor_cutoff: float = 0.5) -> dict:
    """:func:`merge` from edge-list files to an edge-list file.

    Same output as ``merge([read_network(p) for p in paths], cor_cutoff)
    .to_edgelist(output)``; the compiled extension reads, sums, sorts and
    writes the edges (``_gcen_core.network.merge_edge_files``) without
    pandas. Falls back to that call without the extension or for files
    whose gene columns are not the first two. Returns
    ``{"node_count", "edge_count"}``.
    """
    from gcen2._matrix_io import correlation_specs

    if not paths:
        raise GCENInputError("At least one network is required for merging.")
    paths = [require_file(p) for p in paths]
    core = extension()
    specs = correlation_specs(paths) if core is not None else None
    if specs is None:
        from gcen2.io import read_network

        merged = merge([read_network(p) for p in paths], cor_cutoff=cor_cutoff)
        merged.to_edgelist(output)
        return {"node_count": merged.node_count, "edge_count": merged.edge_count}
    try:
        body, node_count, edge_count = core.network.merge_edge_files(specs, float(cor_cutoff))
    except (RuntimeError, UnicodeDecodeError) as exc:
        raise GCENIOError(f"Failed to merge network files: {exc}") from exc
    with open(output, "wb") as fh:
        fh.write(b"gene_a\tgene_b\tcorrelation\tdirection\n")
        fh.write(body)
    return {"node_count": node_count, "edge_count": edge_count}


def shuffle_file(path, output, swap_factor: int = 100, seed: int = 42) -> dict:
    """:func:`shuffle` from an edge-list file to an edge-list file.

    Writes the same file as ``shuffle(read_network(path), swap_factor,
    seed).to_edgelist(output)``: the compiled extension reads the file,
    applies the swaps with numpy's random stream reproduced in C++ and
    writes the rows (``_gcen_core.network.shuffle_edge_file``), so pandas
    is not imported. Falls back to that pipeline without the extension, for
    gene columns that are not the first two, and for files whose pandas
    reading the extension cannot be sure of (it then writes nothing).
    Returns ``{"node_count", "edge_count"}``.
    """
    path = require_file(path)
    core = extension()
    if core is not None:
        from gcen2._matrix_io import network_file_columns

        table = _edge_table_args(*network_file_columns(path))
        if table is not None:
            state = np.random.default_rng(seed).bit_generator.state
            try:
                r = core.network.shuffle_edge_file(
                    str(path), *table, str(output), int(swap_factor), state)
            except (RuntimeError, UnicodeDecodeError) as exc:
                raise GCENIOError(f"Failed to shuffle network file {path}: {exc}") from exc
            if r["ok"]:
                return {"node_count": r["node_count"], "edge_count": r["edge_count"]}
    from gcen2.io import read_network

    shuffled = shuffle(read_network(path), swap_factor=swap_factor, seed=seed)
    shuffled.to_edgelist(output)
    return {"node_count": shuffled.node_count, "edge_count": shuffled.edge_count}


def shuffle(network: Network, swap_factor: int = 100, seed: int = 42) -> Network:
    """Degree-preserving edge swap randomization.

    Repeatedly picks two edges (a1, b1), (a2, b2) and rewires them to
    (a1, b2), (a2, b1) unless that would create a self-loop or an edge that
    already exists in either orientation.

    Parameters
    ----------
    swap_factor:
        Number of swaps = swap_factor × number_of_edges.
    seed:
        Random seed for reproducibility.
    """
    import pandas as pd

    from gcen2.models import Network

    rng = np.random.default_rng(seed)
    edges = network.edges.copy()
    n = len(edges)
    n_swaps = swap_factor * n
    if n < 2:
        return Network(edges=edges, is_signed=network.is_signed, method=network.method)

    # Integer gene codes and an undirected int64 key per edge: much cheaper
    # than (str, str) tuples, and (a, b) / (b, a) count as the same edge so a
    # swap can no longer create a reversed duplicate.
    codes, genes = pd.factorize(
        np.concatenate([edges["gene_a"].to_numpy(), edges["gene_b"].to_numpy()])
    )
    A = codes[:n].tolist()
    B = codes[n:].tolist()
    base = len(genes)

    def key(u: int, v: int) -> int:
        return u * base + v if u < v else v * base + u

    core = extension()
    if core is not None:
        # the swaps and numpy's random stream (PCG64, Generator.integers)
        # both run in C++ (EdgeSwapper.shuffle): same draws as the loop below
        swapper = core.network.EdgeSwapper(codes[:n], codes[n:], base)
        swapper.shuffle(rng.bit_generator.state, n_swaps)
        edges["gene_a"] = genes[codes[:n]]
        edges["gene_b"] = genes[swapper.b()]
        return Network(edges=edges, is_signed=network.is_signed, method=network.method)

    edge_set = {key(u, v) for u, v in zip(A, B)}

    batch = 1 << 20
    done = 0
    while done < n_swaps:
        m = min(batch, n_swaps - done)
        for i, j in zip(rng.integers(0, n, m).tolist(), rng.integers(0, n, m).tolist()):
            if i == j:
                continue
            a1, b1, a2, b2 = A[i], B[i], A[j], B[j]
            if a1 == b2 or a2 == b1:
                continue
            k1, k2 = key(a1, b2), key(a2, b1)
            if k1 == k2 or k1 in edge_set or k2 in edge_set:
                continue
            edge_set.discard(key(a1, b1))
            edge_set.discard(key(a2, b2))
            B[i], B[j] = b2, b1
            edge_set.add(k1)
            edge_set.add(k2)
        done += m

    edges["gene_a"] = genes[np.asarray(A)]
    edges["gene_b"] = genes[np.asarray(B)]
    return Network(edges=edges, is_signed=network.is_signed, method=network.method)


def stats(network: Network) -> NetworkStats:
    """Compute network statistics."""
    return network.stats()


def stats_from_file(path, threads: int = 2) -> NetworkStats:
    """Statistics of the network in an edge-list file, read without pandas.

    Same values as ``stats(read_network(path))``: the compiled extension
    reads the two gene columns into the symmetric adjacency
    (``module.read_file_adjacency``) and counts degrees, triangles and
    components on it. Without the extension it falls back to that call.
    """
    from gcen2.models import NetworkStats

    if extension() is None:
        from gcen2.io import read_network

        return read_network(path).stats()
    from gcen2.module import read_file_adjacency

    indptr, indices, _ = read_file_adjacency(path)
    return NetworkStats.from_adjacency(indptr, indices, threads)


@dataclass
class NetworkSummary:
    """What a report needs of a network, without its edge table: the
    statistics, the node degrees (any order) and the build settings
    (``Network`` defaults for a network read from a file)."""

    stats: NetworkStats
    degrees: np.ndarray
    method: str = "spearman"
    is_signed: bool = False


def summary_from_file(path, threads: int = 2) -> NetworkSummary | Network:
    """:class:`NetworkSummary` of an edge-list file, read without pandas.

    Statistics as :func:`stats_from_file` and degrees as networkx
    ``G.degree`` gives them, from one ``module.read_file_adjacency`` /
    ``graph_stats`` pass. Without the compiled extension it returns
    ``read_network(path)`` instead, which :mod:`gcen2.report` accepts too.
    """
    from gcen2.models import Network, NetworkStats

    if extension() is None:
        from gcen2.io import read_network

        return read_network(path)
    from gcen2.module import read_file_adjacency

    indptr, indices, _ = read_file_adjacency(path)
    st, degrees = NetworkStats.with_degrees(indptr, indices, threads)
    return NetworkSummary(stats=st, degrees=degrees, method=Network.method,
                          is_signed=Network.is_signed)


def _benjamini_hochberg(p_values: np.ndarray, n_total: int) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values with ``n_total`` tests in total.

    ``p_values`` may be a subset of all tests, as long as it contains *every*
    test whose p-value is <= max(p_values) (e.g. all tests with p <= cutoff).
    The rank of each p-value among all ``n_total`` tests is then its rank
    within the subset. For any adjusted value <= max(p_values) the result is
    identical to running BH on the full set; the omitted terms are all larger.
    """
    p = np.asarray(p_values, dtype=np.float64)
    k = p.size
    if k == 0:
        return p.copy()
    order = np.argsort(p, kind="stable")
    ranks = np.arange(1, k + 1, dtype=np.float64)
    adjusted = p[order] * n_total / ranks
    # step-up: q_(i) = min_{j >= i} p_(j) * m / j, capped at 1
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.minimum(adjusted, 1.0)
    fdr = np.empty(k, dtype=np.float64)
    fdr[order] = adjusted
    return fdr
