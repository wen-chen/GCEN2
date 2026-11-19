"""Differential co-expression figures (``gcen2 diff``): edge categories,
correlation in condition A against B, the change against its significance,
and a network of the most changed edges."""

from __future__ import annotations

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.plot._common import canvas, figure, finish
from gcen2.plot._io import table as read_table
from gcen2.plot._style import OKABE_ITO, style

_CAT = {"A_only": OKABE_ITO[0], "B_only": OKABE_ITO[4], "shared_stable": "#999999",
        "shared_changed": OKABE_ITO[5], "gained": OKABE_ITO[4], "lost": OKABE_ITO[0]}


def _diff(x) -> pd.DataFrame:
    return read_table(x)


def diff_overview(diff, alpha: float = 0.05, size="double", height=None, output=None,
                  formats=None):
    """``gcen2 diff`` at a glance: edges per category, correlation in A
    against B (point density; edges of one condition only on the axes), and
    for the shared edges the change in correlation against -log10 p."""
    d = _diff(diff)
    ca = pd.to_numeric(d["cor_a"], errors="coerce").to_numpy(float)
    cb = pd.to_numeric(d["cor_b"], errors="coerce").to_numpy(float)
    with style():
        fig = figure(size, height or 58)
        a1, a2, a3 = fig.subplots(1, 3, width_ratios=[0.8, 1, 1])
        cnt = d["category"].value_counts()
        y = np.arange(len(cnt))
        a1.barh(y, cnt.to_numpy(), color=[_CAT.get(c, OKABE_ITO[3]) for c in cnt.index],
                height=0.7)
        a1.set_yticks(y, [c.replace("_", " ") for c in cnt.index])
        a1.invert_yaxis()
        a1.set_xlabel("Edges")
        for yi, v in zip(y, cnt.to_numpy()):
            a1.text(v, yi, f" {v:,}", va="center", fontsize=4.5)
        a1.set_xlim(0, cnt.max() * 1.35)
        a1.ticklabel_format(axis="x", style="sci", scilimits=(-3, 4), useMathText=True)
        both = np.isfinite(ca) & np.isfinite(cb)
        vals = np.r_[ca[both], cb[both]] if both.any() else np.array([-1.0, 1.0])
        lo, hi = float(np.nanmin(vals)), float(np.nanmax(vals))
        pad = 0.02 * (hi - lo or 1)
        lo, hi = lo - pad, hi + pad
        # square 2-D histogram on the same range for A and B
        from matplotlib.colors import LogNorm

        edges = np.linspace(lo, hi, 61)
        h, _, _ = np.histogram2d(ca[both], cb[both], bins=[edges, edges])
        h[h == 0] = np.nan
        hb = a2.pcolormesh(edges, edges, h.T, cmap="viridis", rasterized=True,
                           norm=LogNorm(vmin=1, vmax=max(np.nanmax(h) if both.any() else 1, 2)))
        a2.plot([lo, hi], [lo, hi], color="black", lw=0.4, ls=":")
        a2.set_xlim(lo, hi)
        a2.set_ylim(lo, hi)
        a2.set_aspect("equal")
        a2.set_xlabel("Correlation in A")
        a2.set_ylabel("Correlation in B")
        a2.set_title(f"{both.sum():,} shared edges")
        fig.colorbar(hb, ax=a2, shrink=0.6).set_label("Edges")
        p = pd.to_numeric(d.get("p_value"), errors="coerce").to_numpy(float) \
            if "p_value" in d.columns else np.full(len(d), np.nan)
        dr = cb - ca
        ok = both & np.isfinite(p)
        if ok.any():
            v = -np.log10(np.maximum(p[ok], 1e-300))
            sig = p[ok] < alpha
            s1 = a3.scatter(dr[ok][~sig], v[~sig], s=1, color="#BBBBBB", edgecolors="none")
            s2 = a3.scatter(dr[ok][sig], v[sig], s=2, color=OKABE_ITO[5], edgecolors="none",
                            label=f"p < {alpha} ({sig.sum():,})")
            for s in (s1, s2):
                s.set_rasterized(True)
            a3.axhline(-np.log10(alpha), color="#999999", lw=0.4, ls="--")
            a3.legend(loc="upper left", markerscale=3)
        a3.set_xlabel("Δ correlation (B − A)")
        a3.set_ylabel("−log10 p")
    return finish(fig, output, formats)


def diff_network(diff, top: int = 60, alpha: float = 0.05, seed: int = 0, size="single",
                 height=None, output=None, formats=None):
    """The ``top`` edges with the largest ``|change|`` among significant shared
    edges and one-condition edges: edge colour = gained in B (blue) or lost
    (orange), width = ``|change|``."""
    try:
        import networkx as nx
    except ImportError as exc:
        raise GCENInputError("diff_network needs networkx.") from exc
    d = _diff(diff).copy()
    d["ca"] = pd.to_numeric(d["cor_a"], errors="coerce").fillna(0)
    d["cb"] = pd.to_numeric(d["cor_b"], errors="coerce").fillna(0)
    d["dr"] = d["cb"].abs() - d["ca"].abs()
    if "p_value" in d.columns:
        pv = pd.to_numeric(d["p_value"], errors="coerce")
        d = d[(pv < alpha) | pv.isna()]
    d = d.reindex(d["dr"].abs().sort_values(ascending=False).index).head(top)
    if d.empty:
        raise GCENInputError("no differential edges.")
    g = nx.Graph()
    for _, r in d.iterrows():
        g.add_edge(str(r["gene_a"]), str(r["gene_b"]), dr=float(r["dr"]))
    from gcen2.plot._common import pack_components, repel_labels

    pos = pack_components(g, seed=seed, weight=None)
    m = max(abs(x) for x in d["dr"]) or 1
    with canvas(None, size, height or 85, 1.0) as (fig, ax):
        for u, v, a in g.edges(data=True):
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]],
                    color=OKABE_ITO[4] if a["dr"] > 0 else OKABE_ITO[0],
                    lw=0.3 + 1.6 * abs(a["dr"]) / m, zorder=1)
        deg = dict(g.degree())
        ns = list(g.nodes)
        ax.scatter([pos[n][0] for n in ns], [pos[n][1] for n in ns],
                   s=[8 + 6 * deg[n] for n in ns], color="white", edgecolors="black",
                   linewidths=0.4, zorder=2)
        ax.set_aspect("equal")
        ax.autoscale_view()
        top = sorted(ns, key=lambda x: -deg[x])[:12]
        repel_labels(ax, np.array([pos[n] for n in top]), top, point_px=3)
        from matplotlib.lines import Line2D

        ax.legend([Line2D([], [], color=OKABE_ITO[4]), Line2D([], [], color=OKABE_ITO[0])],
                  ["stronger in B", "stronger in A"], loc="upper center",
                  bbox_to_anchor=(0.5, 0.0), ncol=2)
        ax.axis("off")
        ax.set_title(f"{g.number_of_edges()} most changed edges")
    return finish(fig, output, formats)
