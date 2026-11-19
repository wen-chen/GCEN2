"""``gcen2 build`` from an expression file without importing numpy.

Importing numpy costs about 17 MB and up to a couple of seconds of start-up,
more than the rest of a typical build needs. For plain matrix files the
compiled extension reads the matrix, computes the correlations (and the FDR)
and writes the edge list in one call
(``_gcen_core.correlation.build_network_file``); other inputs and log
transforms go through :func:`gcen2._matrix_io.read_matrix`
and :func:`gcen2.network.build_to_file`. Both paths write the same file.
"""

from __future__ import annotations

from pathlib import Path

from gcen2.core import extension
from gcen2.exceptions import GCENInputError

#: Correlation methods of the C++ kernel. Dropout-aware ones, per gene pair
#: on a subset of the samples, p-value on its size (src/correlation.hpp):
#: ``pearson_nonzero``, the samples where both genes are non-zero (zeros as
#: missing values); ``pearson_pairwise_complete``, the samples where at least
#: one is non-zero (biased towards negative r on sparse data).
METHODS = ("pearson", "spearman", "pearson_nonzero", "pearson_pairwise_complete")

#: Default gene pairs per block. The C++ kernels reuse their threads across
#: blocks and rank the matrix once, so small blocks cost little time and cap
#: the memory of ``all_edges`` builds (5k genes: about 120 MB instead of
#: 180 MB at 2M).
MAX_PAIRS = 500_000

#: Pearson / Spearman kernels (:func:`gcen2.network.build`).
KERNELS = ("auto", "pairwise", "gemm")


def _check_method(method: str) -> None:
    if method not in METHODS:
        raise GCENInputError(f"Unknown method '{method}'. Choose from: {', '.join(METHODS)}.")


def _check_kernel(kernel: str) -> None:
    if kernel not in KERNELS:
        raise GCENInputError(f"kernel must be one of {', '.join(KERNELS)}, not '{kernel}'.")


def build_from_file(
    input_path,
    path,
    method: str = "spearman",
    cor: float = 0.1,
    pval: float = 0.001,
    signed: bool = False,
    fdr: bool = False,
    log: str | None = None,
    all_edges: bool = False,
    threads: int = 2,
    max_pairs_per_block: int | None = None,
    header: bool | None = None,
    tom: bool = False,
    power: float | str = "auto",
    tom_min: float = 0.0,
    kernel: str = "auto",
) -> dict:
    """Build a network from the tab-separated matrix ``input_path`` and
    write it to ``path``.

    Same file and return value as ``read_matrix(input_path, header=header)``
    followed by ``build_to_file(...)`` with these options. ``header`` works
    as in :func:`gcen2.io.read_tsv` (``None`` auto-detects). ``kernel`` as in
    :func:`gcen2.network.build`.
    ``max_pairs_per_block=None`` uses :data:`MAX_PAIRS` (the output does not
    depend on it).
    """
    _check_method(method)
    tom_power = _tom_power(tom, power, method)
    _check_kernel(kernel)
    if max_pairs_per_block is None:
        max_pairs_per_block = MAX_PAIRS
    input_path = Path(input_path)
    has_header = header
    core = extension()
    if (core is not None and log is None and method in ("pearson", "spearman")
            and input_path.is_file()):
        from gcen2._matrix_io import _detect_header

        if has_header is None:
            has_header = _detect_header(input_path, "\t")
        try:
            res = core.correlation.build_network_file(
                str(input_path), str(path), bool(has_header), method, float(cor), float(pval),
                bool(signed), bool(all_edges), int(threads), int(max_pairs_per_block),
                bool(fdr), tom_power, float(tom_min), kernel=kernel)
        except RuntimeError as exc:  # writing failed
            raise OSError(str(exc)) from exc
        if res is not None:
            return _build_result(*res)

    # anything outside the plain subset: the numpy path handles or reports it
    from gcen2._matrix_io import read_matrix
    from gcen2.network import build_to_file

    matrix, genes, _ = read_matrix(input_path, header=has_header)
    return build_to_file(
        matrix, path, method=method, cor=cor, pval=pval, fdr=fdr, signed=signed, log=log,
        all_edges=all_edges, threads=threads, max_pairs_per_block=max_pairs_per_block,
        gene_names=genes, tom=tom, power=power, tom_min=tom_min, kernel=kernel,
    )


def _tom_power(tom: bool, power, method: str):
    """The kernels' tom_power argument: None (no TOM), 0 (auto) or the power."""
    if not tom:
        return None
    if method not in ("pearson", "spearman"):
        raise GCENInputError("TOM needs method 'pearson' or 'spearman'.")
    if power == "auto" or power is None:
        return 0.0
    try:
        p = float(power)
    except (TypeError, ValueError):
        raise GCENInputError(f"power must be a positive number or 'auto', not {power!r}.") from None
    if not p > 0:
        raise GCENInputError(f"power must be a positive number or 'auto', not {power!r}.")
    return p


_FIT_COLUMNS = ("power", "r2", "slope", "mean_k", "median_k", "max_k")


def _tom_metadata(info) -> dict:
    """Power chosen for the TOM and pickSoftThreshold's table (list of dicts,
    empty when the power was given); warns when no power reached R^2 0.85."""
    out = {"power": info["power"], "estimated": info["estimated"],
           "fits": [dict(zip(_FIT_COLUMNS, f)) for f in info["fits"]]}
    if out["fits"] and not out["estimated"]:
        import warnings

        warnings.warn(f"no soft-threshold power reached scale-free fit R^2 > 0.85; "
                      f"using the WGCNA FAQ power {out['power']:g}", stacklevel=3)
    return out


def _build_result(nodes, edges, info) -> dict:
    res = {"node_count": int(nodes), "edge_count": int(edges)}
    if info is not None:
        res["tom"] = _tom_metadata(info)
    return res
