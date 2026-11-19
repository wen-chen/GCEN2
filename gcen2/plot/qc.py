"""Expression quality-control figures: sample distributions before and
after normalization, the mean-variance relation with the genes kept by
filtering, and sample-sample correlation."""

from __future__ import annotations

import numpy as np
import pandas as pd

from gcen2.plot._common import canvas, figure, finish
from gcen2.plot._style import DIVERGING, OKABE_ITO, figsize, style


def _expr(x) -> pd.DataFrame:
    """genes x samples DataFrame from a DataFrame, ExpressionMatrix or TSV."""
    if isinstance(x, pd.DataFrame):
        return x
    if type(x).__name__ == "ExpressionMatrix":
        return x.data if isinstance(x.data, pd.DataFrame) else pd.DataFrame(x.data)
    from gcen2._matrix_io import read_matrix

    m, g, s = read_matrix(x)
    return pd.DataFrame(m, index=g, columns=s)


def expression_distribution(expr, normalized=None, log: bool = True, max_samples: int = 60,
                            size="double", height=None, output=None, formats=None):
    """Per-sample expression distributions as boxes (log10(x + 1) with
    ``log``), and, with ``normalized``, the same after normalization next to
    them (``gcen2 normalize`` input and output). At most ``max_samples``
    samples are shown (evenly spaced)."""
    panels = [("Raw" if normalized is not None else "Expression", _expr(expr))]
    if normalized is not None:
        panels.append(("Normalized", _expr(normalized)))
    with style():
        fig = figure(size, height or 55 * len(panels))
        axes = np.atleast_1d(fig.subplots(len(panels), 1))
        for ax, (name, x) in zip(axes, panels):
            cols = x.columns
            if len(cols) > max_samples:
                cols = cols[np.linspace(0, len(cols) - 1, max_samples).astype(int)]
            v = x[cols].to_numpy(float)
            if log:
                v = np.log10(np.clip(v, 0, None) + 1)
            bp = ax.boxplot([c[np.isfinite(c)] for c in v.T], widths=0.6, patch_artist=True,
                            showfliers=False)
            for b in bp["boxes"]:
                b.set_facecolor(OKABE_ITO[1])
                b.set_linewidth(0.3)
            for k in ("medians", "whiskers", "caps"):
                for line in bp[k]:
                    line.set_linewidth(0.4)
                    line.set_color("black")
            ax.set_xticks(range(1, len(cols) + 1),
                          [str(c) for c in cols] if len(cols) <= 40 else [""] * len(cols),
                          rotation=90, fontsize=4)
            ax.set_ylabel("log10(x + 1)" if log else "Expression")
            ax.set_title(f"{name}: {x.shape[0]:,} genes x {x.shape[1]} samples")
    return finish(fig, output, formats)


def mean_variance(expr, kept=None, log: bool = True, size="single", height=None, output=None,
                  formats=None):
    """Mean against coefficient of variation of every gene (density of
    points), the genes kept by ``gcen2 filter`` (``kept``: its output or a
    gene list) in colour."""
    x = _expr(expr)
    v = x.to_numpy(float)
    if log:
        v = np.log10(np.clip(v, 0, None) + 1)
    m = np.nanmean(v, axis=1)
    sd = np.nanstd(v, axis=1, ddof=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        cv = sd / m
    ok = np.isfinite(m) & np.isfinite(cv) & (m > 0)
    keep = None
    if kept is not None:
        if isinstance(kept, (list, set, tuple, pd.Index)):
            ks = set(map(str, kept))
        else:
            ks = set(map(str, _expr(kept).index))
        keep = np.array([str(g) in ks for g in x.index])
    with canvas(None, size, height or 60, 0.8) as (fig, ax):
        if keep is None:
            hb = ax.hexbin(m[ok], cv[ok], gridsize=60, bins="log", cmap="viridis", mincnt=1,
                           linewidths=0)
            fig.colorbar(hb, ax=ax, shrink=0.7).set_label("Genes (log10)")
        else:
            out = ok & ~keep
            sc1 = ax.scatter(m[out], cv[out], s=1, color="#BBBBBB", edgecolors="none",
                             label=f"filtered out ({out.sum():,})")
            sc2 = ax.scatter(m[ok & keep], cv[ok & keep], s=1, color=OKABE_ITO[5],
                             edgecolors="none", label=f"kept ({(ok & keep).sum():,})")
            for s in (sc1, sc2):
                s.set_rasterized(True)
            ax.legend(loc="upper right", markerscale=4)
        ax.set_xlabel("Mean" + (" log10(x + 1)" if log else ""))
        ax.set_ylabel("Coefficient of variation")
    return finish(fig, output, formats)


def sample_correlation(expr, method: str = "pearson", log: bool = True, max_samples: int = 200,
                       size="single", height=None, output=None, formats=None):
    """Sample-sample correlation heatmap, samples ordered by average-linkage
    clustering (outliers and batches show up as blocks)."""
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    x = _expr(expr)
    if x.shape[1] > max_samples:
        x = x.iloc[:, np.linspace(0, x.shape[1] - 1, max_samples).astype(int)]
    v = x.to_numpy(float)
    if log:
        v = np.log10(np.clip(v, 0, None) + 1)
    c = pd.DataFrame(v, columns=x.columns).corr(method=method).to_numpy()
    d = np.clip(1 - c, 0, 2)
    np.fill_diagonal(d, 0)
    order = leaves_list(linkage(squareform(d, checks=False), "average"))
    with canvas(None, size, height or figsize(size)[0] / 0.0393701, 1.0) as (fig, ax):
        cm = c[np.ix_(order, order)]
        lo = float(np.nanmin(cm))
        im = ax.imshow(cm, cmap=DIVERGING if lo < 0 else "viridis",
                       vmin=-1 if lo < 0 else lo, vmax=1, interpolation="nearest")
        n = len(order)
        labels = [str(x.columns[i]) for i in order]
        if n <= 50:
            ax.set_xticks(range(n), labels, rotation=90, fontsize=4)
            ax.set_yticks(range(n), labels, fontsize=4)
        else:
            ax.set_xticks([])
            ax.set_yticks([])
        cb = fig.colorbar(im, ax=ax, shrink=0.7)
        cb.set_label(f"{method.capitalize()} r")
        cb.outline.set_linewidth(0.3)
        ax.set_title(f"{n} samples")
    return finish(fig, output, formats)
