"""Trajectory co-expression networks: sliding windows along pseudotime.

Cells (or samples of a time course) are ordered by pseudotime and cut into
overlapping windows of equal size. One co-expression network is built per
window with :func:`gcen2.network.build`, and every edge seen in any window is
followed across all windows (its correlation in each window, the windows
where it passes the cutoffs, and a change pattern).

Requirement §3.6.1 of docs/design/GCEN2_Requirements_Specification.md.
"""

from __future__ import annotations

import warnings
from collections.abc import Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import TYPE_CHECKING, Any, Union

import numpy as np

from gcen2.exceptions import GCENInputError

if TYPE_CHECKING:
    import pandas as pd

    from gcen2.models import Network

#: edge patterns of :attr:`TrajectoryNetworks.dynamics`
PATTERNS = ("stable", "gained", "lost", "transient", "intermittent")


@dataclass
class TrajectoryWindow:
    """One pseudotime window and its network."""

    index: int
    start: float  # smallest pseudotime in the window
    end: float  # largest pseudotime in the window
    center: float  # median pseudotime
    n_cells: int
    n_samples: int  # correlation samples (metacells when metacell_size > 1)
    network: Network


@dataclass
class TrajectoryNetworks:
    """Result of :func:`build_trajectory`.

    ``dynamics`` has one row per edge present in at least one window:
    ``gene_a, gene_b``, the correlation in every window (``cor_0`` ...),
    ``n_windows`` (windows where the edge passes the cutoffs), ``first`` and
    ``last`` (first / last such window), ``pattern`` (see below),
    ``cor_change`` (correlation in the last window minus the first) and
    ``trend`` (Pearson correlation of the edge correlation with the window
    centre; NaN with fewer than 3 windows or a constant series).

    Patterns: ``stable`` (all windows), ``gained`` (absent from the first
    window, present in the last), ``lost`` (present first, absent last),
    ``transient`` (only in inner windows), ``intermittent`` (first and last
    but not all).
    """

    windows: list[TrajectoryWindow]
    dynamics: pd.DataFrame
    params: dict = field(default_factory=dict)

    def window_table(self) -> pd.DataFrame:
        """One row per window: range, sizes and network size."""
        import pandas as pd

        return pd.DataFrame(
            {
                "window": [w.index for w in self.windows],
                "start": [w.start for w in self.windows],
                "end": [w.end for w in self.windows],
                "center": [w.center for w in self.windows],
                "n_cells": [w.n_cells for w in self.windows],
                "n_samples": [w.n_samples for w in self.windows],
                "node_count": [w.network.node_count for w in self.windows],
                "edge_count": [w.network.edge_count for w in self.windows],
            }
        )

    def degree_dynamics(self) -> pd.DataFrame:
        """Gene degree in every window (genes x windows, 0 when absent).

        Rows are sorted by the largest degree reached in any window, so the
        hubs of any stage come first.
        """
        import pandas as pd

        cols = {}
        for w in self.windows:
            e = w.network.edges
            cols[w.index] = pd.concat([e["gene_a"], e["gene_b"]]).value_counts()
        deg = pd.DataFrame(cols).fillna(0).astype(np.int64)
        deg.columns = [f"degree_{c}" for c in deg.columns]
        deg.index.name = "gene"
        order = np.lexsort((deg.index.to_numpy(dtype=str), -deg.max(axis=1).to_numpy()))
        return deg.iloc[order]

    def write(self, directory: Union[str, Path]) -> None:
        """Write ``window_<i>.network``, ``windows.tsv``, ``edge_dynamics.tsv``
        and ``degree_dynamics.tsv`` into ``directory`` (created if needed)."""
        out = Path(directory)
        out.mkdir(parents=True, exist_ok=True)
        for w in self.windows:
            w.network.to_edgelist(out / f"window_{w.index}.network")
        self.window_table().to_csv(out / "windows.tsv", sep="\t", index=False, lineterminator="\n")
        self.dynamics.to_csv(out / "edge_dynamics.tsv", sep="\t", index=False, lineterminator="\n")
        self.degree_dynamics().to_csv(out / "degree_dynamics.tsv", sep="\t", lineterminator="\n")


def build_trajectory(
    data: Any,
    pseudotime: Union[str, Sequence[float], np.ndarray, pd.Series],
    window_size: Union[int, float] = 0.2,
    step: Union[int, float, None] = None,
    metacell_size: int = 1,
    layer: str | None = None,
    use_highly_variable: bool = False,
    method: str = "spearman",
    cor: float = 0.5,
    pval: float = 0.001,
    signed: bool = False,
    fdr: bool = False,
    log: str | None = None,
    threads: int = 2,
) -> TrajectoryNetworks:
    """Build co-expression networks in sliding windows along pseudotime.

    Parameters
    ----------
    data:
        AnnData (cells x genes), or a genes x samples ``pd.DataFrame``
        (e.g. a bulk time course; ``ExpressionMatrix`` is accepted too).
    pseudotime:
        For AnnData, the name of an ``adata.obs`` column, or one value per
        cell. For a DataFrame, one value per column, or a ``pd.Series``
        indexed by column name. Cells with a missing or infinite pseudotime
        are left out (with a warning).
    window_size:
        Cells per window: an int (>= 3), or a float in (0, 1] for a fraction
        of the ordered cells.
    step:
        Shift between window starts, int or fraction as ``window_size``.
        Default: half the window size. A last window is added at the end when
        the steps do not reach the last cell, so every cell is in a window.
    metacell_size:
        If > 1, consecutive cells (in pseudotime order) of a window are
        averaged into metacells of this many cells before correlation,
        which reduces single-cell dropout noise.
    layer:
        AnnData layer to use instead of ``X``.
    use_highly_variable:
        Subset to ``adata.var['highly_variable']`` genes (AnnData only).
    method, cor, pval, signed, fdr, log, threads:
        Passed to :func:`gcen2.network.build` for every window.

    Returns
    -------
    TrajectoryNetworks
    """
    import pandas as pd

    if method not in ("pearson", "spearman"):
        raise GCENInputError(f"method must be 'pearson' or 'spearman', got '{method}'.")

    from gcen2 import network as net_module

    matrix, genes, times = _expression_and_time(data, pseudotime, layer, use_highly_variable)
    finite = np.isfinite(times)
    if not finite.all():
        warnings.warn(f"{int((~finite).sum())} cells without a finite pseudotime are left out.",
                      stacklevel=2)
        matrix, times = matrix[finite], times[finite]
    order = np.argsort(times, kind="stable")
    matrix, times = matrix[order], times[order]
    n = len(times)

    if not isinstance(metacell_size, (int, np.integer)) or metacell_size < 1:
        raise GCENInputError(f"metacell_size must be an integer >= 1, got {metacell_size!r}.")
    w = _count(window_size, n, "window_size")
    s = max(1, w // 2) if step is None else _count(step, n, "step")
    if w > n:
        raise GCENInputError(f"window_size ({w}) is larger than the number of cells ({n}).")
    starts = list(range(0, n - w + 1, s))
    if starts[-1] + w < n:
        starts.append(n - w)
    n_samples = -(-w // metacell_size)
    if n_samples < 3:
        raise GCENInputError(
            f"A window of {w} cells with metacell_size {metacell_size} gives {n_samples} "
            "samples; at least 3 are needed for correlation.")

    windows: list[TrajectoryWindow] = []
    window_mats: list[np.ndarray] = []
    for k, a in enumerate(starts):
        block = _pool(matrix[a:a + w], metacell_size).T  # genes x samples
        expr = pd.DataFrame(block, index=genes)
        net = net_module.build(expr, method=method, cor=cor, pval=pval, signed=signed,
                               fdr=fdr, log=log, threads=threads)
        tw = times[a:a + w]
        net.metadata.update({"window": k, "pseudotime_start": float(tw[0]),
                             "pseudotime_end": float(tw[-1])})
        windows.append(TrajectoryWindow(index=k, start=float(tw[0]), end=float(tw[-1]),
                                        center=float(np.median(tw)), n_cells=w,
                                        n_samples=block.shape[1], network=net))
        window_mats.append(block)

    dynamics = _edge_dynamics(windows, window_mats, genes, method, log)
    params = {"window_size": w, "step": s, "metacell_size": metacell_size, "n_cells": n,
              "method": method, "cor": cor, "pval": pval, "signed": signed, "fdr": fdr,
              "log": log}
    return TrajectoryNetworks(windows=windows, dynamics=dynamics, params=params)


def _count(value, n: int, name: str) -> int:
    """An int count, or a fraction of n (float in (0, 1])."""
    if isinstance(value, bool):
        raise GCENInputError(f"{name} must be an int or a float, got {value!r}.")
    if isinstance(value, (float, np.floating)):
        if not 0.0 < value <= 1.0:
            raise GCENInputError(f"{name} as a fraction must be in (0, 1], got {value}.")
        return max(1, int(round(value * n)))
    if isinstance(value, (int, np.integer)) and value >= 1:
        return int(value)
    raise GCENInputError(f"{name} must be a positive int or a fraction, got {value!r}.")


def _pool(block: np.ndarray, size: int) -> np.ndarray:
    """Average consecutive rows in groups of ``size`` (the last may be smaller)."""
    if size == 1:
        return block
    groups = np.arange(block.shape[0]) // size
    sums = np.zeros((groups[-1] + 1, block.shape[1]))
    np.add.at(sums, groups, block)
    return sums / np.bincount(groups)[:, None]


def _expression_and_time(data, pseudotime, layer, use_highly_variable):
    """Cells x genes float64 matrix, gene names and pseudotime per cell."""
    import pandas as pd

    from gcen2.models import ExpressionMatrix

    if isinstance(data, ExpressionMatrix):
        data = data.data
    if isinstance(data, pd.DataFrame):
        if layer is not None or use_highly_variable:
            raise GCENInputError("layer and use_highly_variable apply to AnnData input only.")
        if isinstance(pseudotime, str):
            raise GCENInputError(
                "For a DataFrame, pseudotime must be values (one per column), not a name.")
        if isinstance(pseudotime, pd.Series):
            missing = set(map(str, data.columns)) - set(map(str, pseudotime.index))
            if missing:
                raise GCENInputError(
                    f"No pseudotime for {len(missing)} samples, e.g. {sorted(missing)[:3]}.")
            pt = pseudotime.copy()
            pt.index = pt.index.map(str)
            times = pt.loc[[str(c) for c in data.columns]].to_numpy(dtype=np.float64)
        else:
            times = np.asarray(pseudotime, dtype=np.float64)
        matrix = data.to_numpy(dtype=np.float64).T
        genes = [str(g) for g in data.index]
    elif hasattr(data, "obs") and hasattr(data, "var_names"):
        if use_highly_variable:
            if "highly_variable" not in data.var.columns:
                raise GCENInputError(
                    "use_highly_variable=True but 'highly_variable' not in adata.var.")
            data = data[:, data.var["highly_variable"].to_numpy(dtype=bool)]
        if isinstance(pseudotime, str):
            if pseudotime not in data.obs.columns:
                raise GCENInputError(f"'{pseudotime}' is not a column of adata.obs.")
            times = pd.to_numeric(data.obs[pseudotime], errors="coerce").to_numpy(np.float64)
        else:
            times = np.asarray(pseudotime, dtype=np.float64)
        X = data.layers[layer] if layer is not None else data.X
        if hasattr(X, "toarray"):
            X = X.toarray()
        matrix = np.asarray(X, dtype=np.float64)
        genes = [str(g) for g in data.var_names]
    else:
        raise GCENInputError("data must be an AnnData object or a genes x samples DataFrame.")

    if times.ndim != 1 or len(times) != matrix.shape[0]:
        raise GCENInputError(
            f"pseudotime has {times.size} values for {matrix.shape[0]} cells/samples.")
    return matrix, genes, times


def _standardize(block: np.ndarray, method: str, log: str | None) -> np.ndarray:
    """Rows scaled so that a dot product of two rows is their correlation
    (Spearman: of the ranks, ties averaged); constant rows become 0."""
    from scipy.stats import rankdata

    from gcen2.network import _LOG_FUNCS

    x = block
    if log is not None:
        x = _LOG_FUNCS[log](x + 1.0)
    if method == "spearman":
        x = rankdata(x, axis=1)
    x = x - x.mean(axis=1, keepdims=True)
    norm = np.sqrt((x * x).sum(axis=1, keepdims=True))
    with np.errstate(invalid="ignore", divide="ignore"):
        z = np.where(norm > 0, x / np.where(norm > 0, norm, 1.0), 0.0)
    return z


def _edge_dynamics(windows, window_mats, genes, method, log) -> pd.DataFrame:
    import pandas as pd

    n_win = len(windows)
    idx = {g: i for i, g in enumerate(genes)}
    n_genes = len(genes)
    codes_per_window = []
    for w in windows:
        e = w.network.edges
        a = e["gene_a"].map(idx).to_numpy(np.int64)
        b = e["gene_b"].map(idx).to_numpy(np.int64)
        lo, hi = np.minimum(a, b), np.maximum(a, b)
        codes_per_window.append(lo * n_genes + hi)
    union = (np.unique(np.concatenate(codes_per_window)) if codes_per_window
             else np.empty(0, np.int64))
    ga, gb = union // n_genes, union % n_genes

    cors = np.empty((len(union), n_win))
    present = np.zeros((len(union), n_win), dtype=bool)
    for k, (block, codes) in enumerate(zip(window_mats, codes_per_window)):
        z = _standardize(block, method, log)
        chunk = max(1, 4_000_000 // max(1, z.shape[1]))
        for s in range(0, len(union), chunk):
            r = np.einsum("ij,ij->i", z[ga[s:s + chunk]], z[gb[s:s + chunk]])
            cors[s:s + chunk, k] = np.clip(r, -1.0, 1.0)
        present[:, k] = np.isin(union, codes, assume_unique=False)

    n_present = present.sum(axis=1)
    first = np.argmax(present, axis=1)
    last = n_win - 1 - np.argmax(present[:, ::-1], axis=1)
    pattern = np.full(len(union), "intermittent", dtype=object)
    p0, pl = present[:, 0], present[:, -1]
    pattern[~p0 & ~pl] = "transient"
    pattern[p0 & ~pl] = "lost"
    pattern[~p0 & pl] = "gained"
    pattern[n_present == n_win] = "stable"

    gene_arr = np.asarray(genes, dtype=object)
    df = pd.DataFrame({"gene_a": gene_arr[ga], "gene_b": gene_arr[gb]})
    for k in range(n_win):
        df[f"cor_{k}"] = cors[:, k]
    df["n_windows"] = n_present
    df["first"] = first
    df["last"] = last
    df["pattern"] = pattern
    df["cor_change"] = cors[:, -1] - cors[:, 0]
    df["trend"] = _trend(cors, np.array([w.center for w in windows]))
    return df.sort_values(["pattern", "cor_change"], key=_sort_key,
                          ascending=[True, False]).reset_index(drop=True)


def _sort_key(col):
    if col.name == "pattern":
        return col.map({p: i for i, p in enumerate(PATTERNS)})
    return col.abs()


def _trend(cors: np.ndarray, centers: np.ndarray) -> np.ndarray:
    """Pearson correlation of each row of ``cors`` with ``centers``."""
    out = np.full(cors.shape[0], np.nan)
    if cors.shape[1] < 3 or np.ptp(centers) == 0:
        return out
    t = centers - centers.mean()
    x = cors - cors.mean(axis=1, keepdims=True)
    den = np.sqrt((x * x).sum(axis=1)) * np.sqrt((t * t).sum())
    ok = den > 0
    out[ok] = (x[ok] @ t) / den[ok]
    return out
