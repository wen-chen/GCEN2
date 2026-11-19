"""Differential co-expression analysis between two conditions."""

from __future__ import annotations

import numpy as np

from gcen2.exceptions import GCENInputError, GCENIOError
from gcen2.models import DiffCoExResult, Network

# pandas is imported inside compare(): compare_files writes the table with the
# compiled extension and never needs it.

CATEGORIES = ("A_only", "B_only", "shared_changed", "shared_stable")
CORRELATIONS = ("pearson", "spearman")


def compare(
    network_a: Network,
    network_b: Network,
    method: str = "fisher_z",
    pval: float = 0.05,
    n_a=None,
    n_b=None,
    fdr: bool = False,
    correlation: str = "pearson",
) -> DiffCoExResult:
    """Compare two co-expression networks to identify rewired edges.

    Uses Fisher's Z-transformation test to determine whether correlation
    coefficients are significantly different between conditions: z is the
    difference of the two ``atanh(r)`` divided by its standard error,
    ``sqrt(1/(n_a - 3) + 1/(n_b - 3))`` for Pearson and
    ``sqrt(1.06/(n_a - 3) + 1.06/(n_b - 3))`` for Spearman correlations
    (Fieller et al., Biometrika 1957; as DGCA). The p-values of the edges
    shared by both networks (the tested pairs) are adjusted by
    Benjamini-Hochberg into the ``fdr`` column.

    Parameters
    ----------
    network_a:
        Network from condition A.
    network_b:
        Network from condition B.
    method:
        Statistical test method. Currently only 'fisher_z' is supported.
    pval:
        P-value cutoff for significance.
    n_a:
        Number of samples used to build network_a. By default from its
        metadata: the per-gene effective samples of a metacell network
        (``gene_effective_samples``, as ``gcen2 celltype --diff`` uses
        them), else ``p_value_samples``, else ``n_samples``, else 30. Or a mapping (e.g. a
        ``pandas.Series``) gene → effective number of samples, such as
        ``metadata["gene_effective_samples"]`` of a metacell network
        (:func:`gcen2.single_cell.metacell_effective_samples`): an edge is
        then tested on the smaller of its two genes' values, and edges with
        fewer than 4 get no z or p-value (``shared_stable``).
    n_b:
        Number of samples used to build network_b (or a mapping, as ``n_a``).
    fdr:
        Call a shared edge ``shared_changed`` when its FDR (instead of its
        p-value) is at most ``pval``.
    correlation:
        How the networks' correlations were computed: ``pearson`` or
        ``spearman`` (sets the standard error of z).

    Returns
    -------
    DiffCoExResult
        ``results`` has one row per undirected edge in either network
        (gene names ordered within the pair), sorted by ``gene_a, gene_b``:
        ``gene_a, gene_b, cor_a, cor_b, z_score, p_value, fdr, category``
        (z, p and FDR only for edges in both networks).
    """
    import pandas as pd

    if method != "fisher_z":
        raise GCENInputError(f"Unknown method '{method}'. Only 'fisher_z' is supported.")
    _check_correlation(correlation)

    if n_a is None:
        n_a = _default_samples(network_a)
    if n_b is None:
        n_b = _default_samples(network_b)
    _check_samples(n_a, n_b)

    # Genes as integer codes in sorted-name order, so that comparing codes
    # compares names; each undirected edge is the key lo * n + hi.
    cols = [net.edges[c].to_numpy(dtype=object)
            for net in (network_a, network_b) for c in ("gene_a", "gene_b")]
    names = np.sort(pd.unique(np.concatenate(cols)))
    codes = [pd.Categorical(c, categories=names).codes.astype(np.int64) for c in cols]
    n = len(names)
    key_a, cor_a_all = _edge_keys(codes[0], codes[1], n, network_a)
    key_b, cor_b_all = _edge_keys(codes[2], codes[3], n, network_b)
    keys = np.union1d(key_a, key_b)
    pos_a = np.searchsorted(key_a, keys)
    pos_b = np.searchsorted(key_b, keys)
    in_a = (pos_a < len(key_a)) & (key_a[np.minimum(pos_a, len(key_a) - 1)] == keys) \
        if len(key_a) else np.zeros(len(keys), bool)
    in_b = (pos_b < len(key_b)) & (key_b[np.minimum(pos_b, len(key_b) - 1)] == keys) \
        if len(key_b) else np.zeros(len(keys), bool)
    cor_a = np.full(len(keys), np.nan)
    cor_b = np.full(len(keys), np.nan)
    cor_a[in_a] = cor_a_all[pos_a[in_a]]
    cor_b[in_b] = cor_b_all[pos_b[in_b]]
    shared = in_a & in_b  # (a NaN correlation still counts as an edge)
    result = pd.DataFrame({"gene_a": names[keys // max(n, 1)], "gene_b": names[keys % max(n, 1)],
                           "cor_a": cor_a, "cor_b": cor_b})

    z_score = np.full(len(result), np.nan)
    p = np.full(len(result), np.nan)
    q = np.full(len(result), np.nan)
    lo, hi = keys[shared] // max(n, 1), keys[shared] % max(n, 1)
    z_score[shared], p[shared] = _fisher_z_arrays(
        cor_a[shared], cor_b[shared], _pair_samples(n_a, names, lo, hi),
        _pair_samples(n_b, names, lo, hi), correlation)
    q[shared] = _bh(p[shared])
    category = np.where(in_a, "A_only", "B_only").astype(object)
    category[shared] = np.where((q if fdr else p)[shared] <= pval, "shared_changed",
                                "shared_stable")
    result["z_score"] = z_score
    result["p_value"] = p
    result["fdr"] = q
    result["category"] = category
    return DiffCoExResult(
        results=result,
        network_a_name=network_a.metadata.get("name", "A"),
        network_b_name=network_b.metadata.get("name", "B"),
        method=method,
        p_cutoff=pval,
    )


def compare_files(path_a, path_b, output, pval: float = 0.05,
                  n_a=None, n_b=None, fdr: bool = False,
                  correlation: str = "pearson", threads: int = 2) -> dict[str, int]:
    """:func:`compare` from two edge-list files to a result TSV.

    Writes the same file as ``compare(read_network(path_a),
    read_network(path_b), ...).results.to_csv(output, sep="\\t",
    index=False)``. The compiled extension reads and joins the edges
    (``_gcen_core.network.diff_edge_files``) and writes the table; only the
    Fisher z test and the FDR of the shared edges run in numpy/scipy
    (as in :func:`compare`), so pandas is never
    imported. Falls back to that pipeline without the extension or for
    files whose gene columns are not the first two. The extension reads,
    sorts and writes with ``threads`` threads (the file does not depend on
    the number). ``n_a`` / ``n_b`` may be mappings gene → effective
    samples, as in :func:`compare`. Returns the row count
    of each category present (``A_only``, ``B_only``, ``shared_changed``,
    ``shared_stable``).
    """
    from gcen2._matrix_io import correlation_specs, require_file
    from gcen2.core import extension

    n_a = 30 if n_a is None else n_a  # read_network gives no n_samples metadata
    n_b = 30 if n_b is None else n_b
    _check_samples(n_a, n_b)
    _check_correlation(correlation)
    paths = [require_file(path_a), require_file(path_b)]
    core = extension()
    specs = correlation_specs(paths) if core is not None else None
    if specs is None:
        from gcen2.io import read_network

        res = compare(read_network(paths[0]), read_network(paths[1]), pval=pval,
                      n_a=n_a, n_b=n_b, fdr=fdr, correlation=correlation).results
        res.to_csv(output, sep="\t", index=False, lineterminator="\n")
        counts = res["category"].value_counts()
        return {c: int(counts[c]) for c in CATEGORIES if c in counts}
    try:
        table = core.network.diff_edge_files(specs[0], specs[1], int(threads))
    except (RuntimeError, UnicodeDecodeError) as exc:
        raise GCENIOError(f"Failed to read network files: {exc}") from exc
    shared = table.flags() == 3
    if _is_scalar(n_a) and _is_scalar(n_b):
        pa, pb = n_a, n_b
    else:
        names = np.asarray(table.names(), dtype=object)
        lo, hi = table.lo()[shared], table.hi()[shared]
        pa, pb = _pair_samples(n_a, names, lo, hi), _pair_samples(n_b, names, lo, hi)
    z, p = _fisher_z_arrays(table.cor_a()[shared], table.cor_b()[shared], pa, pb,
                            correlation, threads)
    del shared
    q = _bh(p, threads)
    try:
        counts = table.write(str(output), z, p, q, float(pval), bool(fdr), int(threads))
    except RuntimeError as exc:
        raise GCENIOError(f"Failed to write {output}: {exc}") from exc
    return {c: int(n) for c, n in zip(CATEGORIES, counts) if n}


def _default_samples(network: Network):
    """The sample size of a network's own p-values (see :func:`compare`)."""
    md = network.metadata
    for key in ("gene_effective_samples", "p_value_samples", "n_samples"):
        if md.get(key) is not None:
            return md[key]
    return 30


def _is_scalar(n) -> bool:
    return np.ndim(n) == 0 and not hasattr(n, "keys")


def _check_samples(n_a, n_b) -> None:
    """Scalars: at least 4 each. Mappings gene -> samples are checked per
    edge (:func:`_pair_samples`; fewer than 4: no test)."""
    for n in (n_a, n_b):
        if _is_scalar(n):
            if n < 4:
                raise GCENInputError(
                    "At least 4 samples are required per condition for Fisher's Z test.")
        elif not hasattr(n, "keys"):
            raise GCENInputError("n_a / n_b must be a number or a mapping gene -> samples.")


def _pair_samples(n, names: np.ndarray, lo: np.ndarray, hi: np.ndarray):
    """The number of samples of each edge (genes ``names[lo]``, ``names[hi]``):
    ``n`` itself for a number, else the smaller of the two genes' values of
    the mapping ``n`` (gene -> effective samples)."""
    if _is_scalar(n):
        return n
    import pandas as pd

    ser = pd.Series(n, dtype=np.float64) if not isinstance(n, pd.Series) \
        else n.astype(np.float64)
    pos = pd.Index(ser.index.astype(str)).get_indexer(np.asarray(names, dtype=str))
    missing = np.unique(np.concatenate([lo, hi]))
    missing = missing[pos[missing] < 0]
    if len(missing):
        raise GCENInputError(
            f"No effective sample size for {len(missing)} gene(s) of the edges, e.g. "
            f"'{names[missing[0]]}'.")
    vals = ser.to_numpy()[pos] if len(pos) else np.empty(0)
    return np.minimum(vals[lo], vals[hi])


def _edge_keys(ca: np.ndarray, cb: np.ndarray, n: int, network: Network):
    """Sorted unique undirected edge keys ``min * n + max`` of the code pairs
    and the correlation of each (the last row wins for a repeated edge)."""
    key = np.minimum(ca, cb) * n + np.maximum(ca, cb)
    cor = network.edges["correlation"].to_numpy(dtype=np.float64)
    # last occurrence of each key: unique on the reversed array
    rev_key, rev_idx = np.unique(key[::-1], return_index=True)
    return rev_key, cor[::-1][rev_idx]


def _check_correlation(correlation: str) -> None:
    if correlation not in CORRELATIONS:
        raise GCENInputError(
            f"correlation must be one of {CORRELATIONS}, got '{correlation}'.")


def _bh(p: np.ndarray, threads: int = 1) -> np.ndarray:
    """Benjamini-Hochberg adjusted p-values over the finite entries of ``p``
    (NaN stays NaN), as R's ``p.adjust(p, "BH")``. The compiled extension
    (``_gcen_core.network.bh_adjust``, sorting with ``threads`` threads)
    gives the same values bit for bit."""
    from gcen2.core import extension

    core = extension()
    if core is not None:
        return np.asarray(core.network.bh_adjust(p, int(threads)))
    q = np.full(len(p), np.nan)
    ok = np.flatnonzero(np.isfinite(p))
    m = len(ok)
    if m:
        order = ok[np.argsort(p[ok], kind="stable")]
        adj = p[order] * m / np.arange(1, m + 1)
        q[order] = np.minimum(np.minimum.accumulate(adj[::-1])[::-1], 1.0)
    return q


def _fisher_z_arrays(r_a: np.ndarray, r_b: np.ndarray, n_a, n_b,
                     correlation: str = "pearson", threads: int = 1):
    """Fisher's z test of two correlations, element-wise: z and the two-sided p.
    ``n_a`` / ``n_b``: numbers of samples, or arrays of one per element
    (elements with fewer than 4 get NaN).
    Spearman correlations have the variance 1.06 / (n - 3) (Fieller et al.
    1957), Pearson 1 / (n - 3). Large inputs are split into chunks computed
    by ``threads`` threads (numpy releases the GIL; the values do not change)."""
    from scipy.stats import norm

    v = 1.06 if correlation == "spearman" else 1.0
    per_edge = np.ndim(n_a) > 0 or np.ndim(n_b) > 0
    if per_edge:  # per-edge sample sizes: fewer than 4 -> no test (NaN)
        na = np.broadcast_to(np.asarray(n_a, dtype=np.float64), np.shape(r_a))
        nb = np.broadcast_to(np.asarray(n_b, dtype=np.float64), np.shape(r_a))
        ok = (na >= 4) & (nb >= 4)
        with np.errstate(divide="ignore", invalid="ignore"):
            se_all = np.where(ok, np.sqrt(v / (na - 3) + v / (nb - 3)), np.nan)
    else:
        se_all = np.sqrt(v / (n_a - 3) + v / (n_b - 3))

    def test(a, b, se):
        a = np.clip(a, -0.9999, 0.9999)
        b = np.clip(b, -0.9999, 0.9999)
        z = (0.5 * np.log((1 + a) / (1 - a)) - 0.5 * np.log((1 + b) / (1 - b))) / se
        return z, 2.0 * norm.sf(np.abs(z))

    def se_of(lo, hi):
        return se_all[lo:hi] if per_edge else se_all

    n = np.size(r_a)
    step = 1 << 20
    if threads <= 1 or n <= step:
        return test(r_a, r_b, se_all)
    from concurrent.futures import ThreadPoolExecutor

    z_stat, p = np.empty(n), np.empty(n)

    def chunk(lo: int) -> None:
        hi = min(n, lo + step)
        z_stat[lo:hi], p[lo:hi] = test(r_a[lo:hi], r_b[lo:hi], se_of(lo, hi))

    with ThreadPoolExecutor(max_workers=int(threads)) as pool:
        list(pool.map(chunk, range(0, n, step)))
    return z_stat, p
