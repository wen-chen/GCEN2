"""Single-cell and spatial module figures (hdWGCNA): eigengenes on an
embedding or in space, by group (dot plot, violins, radar), differential
module eigengenes, the cell-type overview, consensus sets, module
preservation and projection."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.plot import _io
from gcen2.plot._common import canvas, figure, finish
from gcen2.plot._style import (
    DIVERGING,
    OKABE_ITO,
    figsize,
    is_unassigned,
    module_color,
    module_label,
    style,
)


def _modules_of(me: pd.DataFrame, modules=None) -> list[str]:
    cols = [c for c in me.columns if not is_unassigned(c) and c not in ("set", "sample")]
    if modules is None:
        return cols
    want = [("ME" + str(m)) if str(m).isdigit() else str(m) for m in modules]
    missing = [m for m in want if m not in cols]
    if missing:
        raise GCENInputError(f"modules not found: {', '.join(missing)}.")
    return want


def coordinates(coords, cells=None, basis: str = "X_umap") -> pd.DataFrame:
    """cells x 2 coordinates from a DataFrame (index = cells), an array (rows
    in the order of ``cells``), a TSV file (first column cells), or an
    AnnData / .h5ad (its ``obsm[basis]``)."""
    from gcen2._anndata import is_anndata

    if isinstance(coords, (str, Path)) and str(coords).endswith(".h5ad"):
        import anndata as ad

        coords = ad.read_h5ad(coords, backed="r")
    if is_anndata(coords):
        if basis not in coords.obsm:
            raise GCENInputError(f"no obsm['{basis}'] (have: {', '.join(coords.obsm)}).")
        xy = np.asarray(coords.obsm[basis])[:, :2]
        return pd.DataFrame(xy, index=coords.obs_names.map(str), columns=["x", "y"])
    if isinstance(coords, (str, Path)):
        t = _io.table(coords, index_col=0)
        t.index = t.index.map(str)
        return t.iloc[:, :2].set_axis(["x", "y"], axis=1).astype(float)
    if isinstance(coords, pd.DataFrame):
        t = coords.iloc[:, :2].copy()
        t.index = t.index.map(str)
        return t.set_axis(["x", "y"], axis=1).astype(float)
    xy = np.asarray(coords, dtype=float)[:, :2]
    if cells is None:
        raise GCENInputError("give the cell names of an array of coordinates.")
    return pd.DataFrame(xy, index=list(map(str, cells)), columns=["x", "y"])


def _panel_grid(n, ncols):
    ncols = max(1, min(ncols, n))
    return ncols, int(np.ceil(n / ncols))


# --------------------------------------------------------------------------
def me_embedding(eigengenes, coords, basis: str = "X_umap", modules=None, ncols: int = 4,
                 point_size: float | None = None, quantile: float = 0.99,
                 spatial: bool = False, flip_y: bool = True, size="double", height=None,
                 output=None, formats=None):
    """Eigengene of every cell on an embedding, one panel per module (hdWGCNA
    ModuleFeaturePlot): symmetric diverging colours around 0, clipped at the
    ``quantile`` of ``|eigengene|``, high ``|values|`` drawn on top. ``eigengenes``:
    cell_eigengenes.tsv (or cell_harmonized_eigengenes.tsv, spot_eigengenes.tsv,
    projected_eigengenes.tsv); ``coords``: an AnnData / .h5ad (``obsm[basis]``),
    a TSV or DataFrame of x, y by cell. ``spatial``: equal axes, y pointing
    down when ``flip_y`` (image coordinates)."""
    me = _io.eigengenes(eigengenes)
    me.index = me.index.map(str)
    xy = coordinates(coords, basis=basis).reindex(me.index)
    ok = xy.notna().all(axis=1).to_numpy()
    if ok.sum() == 0:
        raise GCENInputError("no cell of the eigengenes has coordinates.")
    mods = _modules_of(me, modules)
    ncols, nrows = _panel_grid(len(mods), ncols)
    n = int(ok.sum())
    if point_size is not None:
        s = point_size
    else:  # diameter ~ panel width / sqrt(points): spots tile, cells stay legible
        d_pt = figsize(size)[0] * 72 / ncols * 0.8 / np.sqrt(n) * (0.95 if spatial else 0.6)
        s = float(np.clip(d_pt ** 2, 0.2, 40))
    with style():
        w = figsize(size)[0] / 0.0393701
        fig = figure(size, height or min(170, nrows * (w / ncols) * (0.95 if spatial else 0.85)))
        axes = np.atleast_1d(fig.subplots(nrows, ncols, squeeze=False)).ravel()
        x, y = xy["x"].to_numpy()[ok], xy["y"].to_numpy()[ok]
        for ax, m in zip(axes, mods):
            v = me[m].to_numpy(float)[ok]
            lim = float(np.nanquantile(np.abs(v), quantile)) or 1.0
            order = np.argsort(np.abs(v))
            sc = ax.scatter(x[order], y[order], c=v[order], s=s, cmap=DIVERGING, vmin=-lim,
                            vmax=lim, edgecolors="none")
            if n > 2000:
                sc.set_rasterized(True)
            ax.set_title(module_label(m))
            ax.set_xticks([])
            ax.set_yticks([])
            for sp in ax.spines.values():
                sp.set_visible(False)
            if spatial:
                ax.set_aspect("equal")
                if flip_y:
                    ax.invert_yaxis()
            cb = fig.colorbar(sc, ax=ax, shrink=0.6, aspect=12, pad=0.01)
            cb.outline.set_linewidth(0.3)
            cb.ax.tick_params(labelsize=5)
            cb.ax.text(0.5, 1.03, "ME", transform=cb.ax.transAxes, ha="center", va="bottom",
                       fontsize=5)
        for ax in axes[len(mods):]:
            ax.axis("off")
    return finish(fig, output, formats)


def spatial_me(eigengenes, coords, basis: str = "spatial", **kwargs):
    """Spatial maps of spot eigengenes (:func:`me_embedding` with equal axes
    and image coordinates; ``coords`` an AnnData with ``obsm["spatial"]`` or
    x, y by spot)."""
    kwargs.setdefault("spatial", True)
    return me_embedding(eigengenes, coords, basis=basis, **kwargs)


# --------------------------------------------------------------------------
def _group_series(groups, index) -> pd.Series:
    if isinstance(groups, pd.Series):
        g = groups.copy()
        g.index = g.index.map(str)
        return g.reindex(index).astype(str)
    if isinstance(groups, (str, Path)):
        t = _io.table(groups, index_col=0)
        t.index = t.index.map(str)
        return t.iloc[:, 0].reindex(index).astype(str)
    g = np.asarray(groups)
    if g.size != len(index):
        raise GCENInputError("groups must give one label per cell.")
    return pd.Series(g.astype(str), index=index)


def groups_from(adata_or_table, column: str, cells) -> pd.Series:
    """Group labels of ``cells``: ``column`` of an AnnData / .h5ad obs or of a
    TSV (first column cells)."""
    from gcen2._anndata import is_anndata

    x = adata_or_table
    if isinstance(x, (str, Path)) and str(x).endswith(".h5ad"):
        import anndata as ad

        x = ad.read_h5ad(x, backed="r")
    obs = x.obs if is_anndata(x) else _io.table(x, index_col=0)
    if column not in obs.columns:
        raise GCENInputError(f"no column '{column}'.")
    s = obs[column].astype(str)
    s.index = s.index.map(str)
    return s.reindex(list(map(str, cells)))


def me_dotplot(eigengenes, groups, modules=None, scale: bool = True, size="single",
               height=None, output=None, formats=None):
    """Eigengenes by group (Seurat DotPlot of the eigengenes, as in the
    hdWGCNA tutorials): colour = mean eigengene of the group (z-scored
    across groups with ``scale``), dot size = percentage of the group's
    cells with an eigengene above 0. ``groups``: one label per cell (a
    Series indexed by cell)."""
    me = _io.eigengenes(eigengenes)
    me.index = me.index.map(str)
    g = _group_series(groups, me.index)
    keep = g.notna() & (g != "nan")
    me, g = me[keep.to_numpy()], g[keep]
    mods = _modules_of(me, modules)
    order = sorted(pd.unique(g))
    mean = me[mods].groupby(g.to_numpy()).mean().reindex(order)
    pct = (me[mods] > 0).groupby(g.to_numpy()).mean().reindex(order) * 100
    val = (mean - mean.mean()) / mean.std(ddof=1).replace(0, 1) if scale else mean
    with canvas(None, size, height or max(35, 8 + 5.5 * len(mods)), 0.8) as (fig, ax):
        xs, ys = np.meshgrid(np.arange(len(order)), np.arange(len(mods)))
        lim = float(np.nanmax(np.abs(val.to_numpy()))) or 1.0
        sc = ax.scatter(xs.ravel(), ys.ravel(), s=(pct.T.to_numpy().ravel() / 100) * 60 + 0.5,
                        c=val.T.to_numpy().ravel(), cmap=DIVERGING, vmin=-lim, vmax=lim,
                        edgecolors="black", linewidths=0.2)
        ax.set_xticks(range(len(order)), order, rotation=45, ha="right")
        ax.set_yticks(range(len(mods)), [module_label(m) for m in mods])
        ax.set_xlim(-0.6, len(order) - 0.4)
        ax.set_ylim(len(mods) - 0.5, -0.5)
        cb = fig.colorbar(sc, ax=ax, shrink=0.6, aspect=12, pad=0.02)
        cb.set_label("Scaled mean ME" if scale else "Mean ME")
        cb.outline.set_linewidth(0.3)
        for p in (25, 50, 100):
            ax.scatter([], [], s=(p / 100) * 60 + 0.5, color="white", edgecolors="black",
                       linewidths=0.3, label=f"{p}%")
        ax.legend(title="Cells > 0", loc="upper left", bbox_to_anchor=(1.18, 1.0),
                  labelspacing=0.8)
    return finish(fig, output, formats)


def me_violin(eigengenes, groups, modules=None, ncols: int = 3, size="double", height=None,
              output=None, formats=None):
    """Distribution of each module's eigengene by group (violins with the
    median marked), one panel per module."""
    me = _io.eigengenes(eigengenes)
    me.index = me.index.map(str)
    g = _group_series(groups, me.index)
    mods = _modules_of(me, modules)
    order = sorted(pd.unique(g[g != "nan"]))
    ncols, nrows = _panel_grid(len(mods), ncols)
    with style():
        fig = figure(size, height or min(170, 48 * nrows))
        axes = np.atleast_1d(fig.subplots(nrows, ncols, squeeze=False)).ravel()
        for ax, m in zip(axes, mods):
            data = [me.loc[(g == o).to_numpy(), m].dropna().to_numpy(float) for o in order]
            ok = [i for i, d in enumerate(data) if d.size > 1]
            parts = ax.violinplot([data[i] for i in ok], positions=ok, widths=0.8,
                                  showmedians=True, showextrema=False)
            for b in parts["bodies"]:
                b.set_facecolor(module_color(m))
                b.set_edgecolor("black")
                b.set_linewidth(0.3)
                b.set_alpha(0.8)
            parts["cmedians"].set_color("black")
            parts["cmedians"].set_linewidth(0.6)
            ax.axhline(0, color="#999999", lw=0.4, ls=":")
            ax.set_xticks(range(len(order)), order, rotation=45, ha="right")
            ax.set_title(module_label(m))
            ax.set_ylabel("ME")
        for ax in axes[len(mods):]:
            ax.axis("off")
    return finish(fig, output, formats)


def me_radar(eigengenes, groups, modules=None, ncols: int = 4, size="double", height=None,
             output=None, formats=None):
    """Each module's mean eigengene per group on a radar chart (hdWGCNA
    ModuleRadarPlot), scaled to 0-1 across the groups."""
    me = _io.eigengenes(eigengenes)
    me.index = me.index.map(str)
    g = _group_series(groups, me.index)
    mods = _modules_of(me, modules)
    order = sorted(pd.unique(g[g != "nan"]))
    if len(order) < 3:
        raise GCENInputError("a radar chart needs at least 3 groups.")
    mean = me[mods].groupby(g.to_numpy()).mean().reindex(order)
    rng = (mean.max() - mean.min()).replace(0, 1)
    sc = (mean - mean.min()) / rng
    ncols, nrows = _panel_grid(len(mods), ncols)
    ang = np.linspace(0, 2 * np.pi, len(order), endpoint=False)
    with style():
        w = figsize(size)[0] / 0.0393701
        fig = figure(size, height or min(170, nrows * w / ncols * 1.05))
        for i, m in enumerate(mods):
            ax = fig.add_subplot(nrows, ncols, i + 1, projection="polar")
            v = sc[m].to_numpy()
            col = module_color(m)
            ax.plot(np.append(ang, ang[0]), np.append(v, v[0]), color="black", lw=0.6)
            ax.fill(np.append(ang, ang[0]), np.append(v, v[0]), color=col, alpha=0.6)
            ax.set_xticks(ang, [_wrap(o) for o in order], fontsize=4.5)
            ax.tick_params(axis="x", pad=7)
            ax.set_yticks([0.5, 1.0], ["", ""])
            ax.set_ylim(0, 1.05)
            ax.grid(linewidth=0.3)
            ax.spines["polar"].set_linewidth(0.3)
            ax.set_title(module_label(m), pad=8)
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def _wrap(text: str, width: int = 12) -> str:
    """Break a long label at underscores, spaces or dashes."""
    import re

    if len(text) <= width:
        return text
    parts = re.split(r"(?<=[_\- ])", text)
    lines, cur = [], ""
    for p in parts:
        if cur and len(cur) + len(p) > width:
            lines.append(cur)
            cur = p
        else:
            cur += p
    lines.append(cur)
    return "\n".join(lines)


def _dmes(dmes) -> pd.DataFrame:
    d = _io.table(dmes)
    if "module" not in d.columns:
        d = d.reset_index().rename(columns={"index": "module"})
    if "group" not in d.columns:
        d["group"] = "group"
    return d


def dme_volcano(dmes, alpha: float = 0.05, ncols: int = 3, size="double", height=None,
                output=None, formats=None):
    """Differential module eigengenes (hdWGCNA PlotDMEsVolcano): log2 fold
    change against -log10 adjusted p-value, modules in their colours and
    named, one panel per group. ``dmes``: dmes.tsv (or find_dmes /
    FindDMEs output)."""
    d = _dmes(dmes)
    groups = list(pd.unique(d["group"]))
    ncols, nrows = _panel_grid(len(groups), ncols)
    floor = np.nextafter(0, 1)
    with style():
        fig = figure(size, height or min(170, 50 * nrows))
        axes = np.atleast_1d(fig.subplots(nrows, ncols, squeeze=False)).ravel()
        for ax, grp in zip(axes, groups):
            x = d[d["group"] == grp]
            fc = x["avg_log2FC"].to_numpy(float)
            raw = x["p_val_adj"].to_numpy(float)
            capped = raw <= floor  # p underflowed to 0: drawn at the cap, as triangles
            yv = -np.log10(np.maximum(raw, floor))
            sig = raw < alpha
            fin = np.isfinite(fc)
            lim = np.nanmax(np.abs(fc[fin])) * 1.15 if fin.any() else 1
            cols = np.array([module_color(m) for m in x["module"]], dtype=object)
            for mk, sel in (("o", fin & ~capped), ("^", fin & capped)):
                if sel.any():
                    ax.scatter(fc[sel], yv[sel], s=np.where(sig[sel], 22, 10), marker=mk,
                               c=list(cols[sel]), edgecolors="black", linewidths=0.3,
                               zorder=3)
            ax.axhline(-np.log10(alpha), color="#999999", lw=0.4, ls="--")
            if (fin & capped).any():
                from matplotlib.lines import Line2D

                ax.legend([Line2D([], [], marker="^", ls="", color="#888888", mec="black",
                                  mew=0.3, ms=4)], ["p below 1e-308 (capped)"], loc="best",
                          fontsize=4.5, handletextpad=0.2)
            ax.axvline(0, color="#999999", lw=0.4)
            ax.set_xlim(-lim, lim)
            ax.margins(y=0.25)
            ax.autoscale_view()
            ax.set_ylim(-0.04 * ax.get_ylim()[1], ax.get_ylim()[1])
            ax.set_xlabel("Average log2 fold change")
            ax.set_ylabel("−log10 adjusted p")
            ax.set_title(str(grp))
            from gcen2.plot._common import repel_labels

            repel_labels(ax, np.column_stack([fc[fin], yv[fin]]),
                         [module_label(m).split(" ")[0] for m in x["module"][fin]],
                         point_px=4, leader=False)
        for ax in axes[len(groups):]:
            ax.axis("off")
    return finish(fig, output, formats)


def dme_lollipop(dmes, modules=None, alpha: float = 0.05, ncols: int = 3, size="double",
                 height=None, output=None, formats=None):
    """Differential module eigengenes as a lollipop chart (hdWGCNA
    PlotDMEsLollipop): fold change per module from 0, dot size = module size
    (with ``modules``: modules.tsv), crosses for modules not significant at
    ``alpha``; one panel per group."""
    d = _dmes(dmes)
    sizes = None
    if modules is not None:
        sizes = _io.modules(modules).value_counts()
    groups = list(pd.unique(d["group"]))
    ncols, nrows = _panel_grid(len(groups), ncols)
    with style():
        n_mod = d.groupby("group").size().max()
        fig = figure(size, height or min(170, nrows * (16 + 7 * n_mod)))
        axes = np.atleast_1d(fig.subplots(nrows, ncols, squeeze=False)).ravel()
        for ax, grp in zip(axes, groups):
            x = d[d["group"] == grp].sort_values("avg_log2FC")
            y = np.arange(len(x))
            fc = x["avg_log2FC"].to_numpy(float)
            fc = np.where(np.isfinite(fc), fc, 0)
            sig = x["p_val_adj"].to_numpy(float) < alpha
            cols = [module_color(m) for m in x["module"]]
            ax.hlines(y, 0, fc, color=cols, lw=1.0)
            if sizes is not None:
                s = np.array([sizes.get(m, 0) for m in x["module"]], dtype=float)
                s = 12 + 60 * s / max(s.max(), 1)
            else:
                s = np.full(len(x), 30.0)
            ax.scatter(fc[sig], y[sig], s=s[sig], c=[c for c, k in zip(cols, sig) if k],
                       edgecolors="black", linewidths=0.3, zorder=3)
            ax.scatter(fc[~sig], y[~sig], s=s[~sig] * 0.6, marker="x", color="black",
                       linewidths=0.6, zorder=3)
            ax.axvline(0, color="black", lw=0.4)
            ax.set_yticks(y, [module_label(m) for m in x["module"]])
            ax.set_ylim(-0.7, len(x) - 0.3)
            if ax is axes[0]:
                from matplotlib.lines import Line2D

                hd = [Line2D([], [], marker="o", ls="", color="#888888", mec="black", mew=0.3,
                             ms=4, label=f"adj. p < {alpha:g}"),
                      Line2D([], [], marker="x", ls="", color="black", ms=4,
                             label="not significant")]
                if sizes is not None:
                    hd.append(Line2D([], [], ls="", label="dot size: genes in module"))
                fig.legend(handles=hd, loc="outside lower center", ncol=len(hd),
                           frameon=False, fontsize=5)
            ax.margins(x=0.12)
            ax.set_xlabel("Average log2 fold change")
            ax.set_title(str(grp))
        for ax in axes[len(groups):]:
            ax.axis("off")
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def celltype_overview(table, size="double", height=None, output=None, formats=None):
    """celltype-modules at a glance: per cell type the cells, metacells,
    genes, modules and genes in modules (celltypes.tsv), cell types from top
    to bottom in the order of the table (as in :func:`metacells` and the
    report's sections)."""
    t = _io.table(table)
    cols = [("n_cells", "Cells"), ("n_metacells", "Metacells"), ("n_genes", "Genes"),
            ("modules", "Modules"), ("genes_in_modules", "Genes in modules")]
    cols = [c for c in cols if c[0] in t.columns]
    with style():
        fig = figure(size, height or max(35, 8 + 4 * len(t)))
        axes = np.atleast_1d(fig.subplots(1, len(cols), sharey=True))
        y = np.arange(len(t))
        for i, (ax, (c, name)) in enumerate(zip(axes, cols)):
            v = t[c].to_numpy(float)
            ax.barh(y, v, color=OKABE_ITO[i % len(OKABE_ITO)], height=0.7)
            ax.set_title(name)
            ax.tick_params(axis="x", labelsize=5)
            for yi, vi in zip(y, v):
                ax.text(vi, yi, f" {vi:,.0f}", va="center", ha="left", fontsize=4.5)
            ax.set_xlim(0, v.max() * 1.3 if v.max() > 0 else 1)
        axes[0].set_yticks(y, t["cell_type"].astype(str))
        axes[0].invert_yaxis()  # shared y: first cell type on top in every panel
    return finish(fig, output, formats)


def set_eigengenes(table, modules=None, ncols: int = 4, size="double", height=None,
                   output=None, formats=None):
    """Consensus network: each module's eigengene in the metacells of every
    set (set_eigengenes.tsv), as boxes per set, to show the module is
    present in all sets."""
    t = _io.table(table)
    mods = _modules_of(t.drop(columns=["set", "sample"]), modules)
    sets = list(pd.unique(t["set"]))
    ncols, nrows = _panel_grid(len(mods), ncols)
    with style():
        fig = figure(size, height or min(170, 40 * nrows))
        axes = np.atleast_1d(fig.subplots(nrows, ncols, squeeze=False)).ravel()
        for ax, m in zip(axes, mods):
            data = [t.loc[t["set"] == s, m].to_numpy(float) for s in sets]
            bp = ax.boxplot(data, widths=0.6, patch_artist=True, showfliers=False)
            for b in bp["boxes"]:
                b.set_facecolor(module_color(m))
                b.set_linewidth(0.4)
            for k in ("medians", "whiskers", "caps"):
                for line in bp[k]:
                    line.set_linewidth(0.4)
                    line.set_color("black")
            ax.set_xticks(range(1, len(sets) + 1), sets, rotation=45, ha="right")
            ax.set_title(module_label(m))
            ax.set_ylabel("ME")
        for ax in axes[len(mods):]:
            ax.axis("off")
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def _pres(z):
    t = _io.table(z, index_col=0)
    t.index = t.index.map(str)
    return t


def preservation(z_table, observed=None, size="double", height=None, output=None,
                 formats=None):
    """Module preservation (hdWGCNA PlotModulePreservation, WGCNA's
    preservation figure): Zsummary.pres against module size with the
    thresholds 2 and 10, and the median rank (with ``observed``:
    preservation.tsv) against module size; modules in their colours and
    named, WGCNA's gold (random) module included as reference."""
    z = _pres(z_table)
    obs = _pres(observed) if observed is not None else None
    panels = [("Zsummary.pres", z, "Zsummary (preservation)")]
    if obs is not None and "medianRank.pres" in obs.columns:
        panels.append(("medianRank.pres", obs, "Median rank (preservation)"))
    with style():
        fig = figure(size if len(panels) > 1 else "single", height or 60)
        axes = np.atleast_1d(fig.subplots(1, len(panels)))
        for ax, (col, t, lab) in zip(axes, panels):
            x = t["moduleSize"].to_numpy(float)
            y = t[col].to_numpy(float)
            mods = list(t.index)
            ax.scatter(x, y, s=30, c=[module_color(m) for m in mods], edgecolors="black",
                       linewidths=0.4, zorder=3)
            for xi, yi, m in zip(x, y, mods):
                ax.annotate(module_label(m).split(" ")[0] if not is_unassigned(m) else "grey",
                            (xi, yi), xytext=(3, 2), textcoords="offset points", fontsize=4.5)
            ax.set_xscale("log")
            from gcen2.plot._common import plain_log_axis

            plain_log_axis(ax.xaxis, subs=(1.0, 2.0, 3.0, 5.0))
            ax.margins(x=0.12)
            ax.set_xlabel("Module size (genes)")
            ax.set_ylabel(lab)
            if col.startswith("Zsummary"):
                for zv, c, txt in ((2, "#0072B2", "Z = 2 (moderate)"),
                                   (10, "#009E73", "Z = 10 (strong)")):
                    ax.axhline(zv, color=c, lw=0.5, ls="--")
                    ax.text(1.0, zv, txt, transform=ax.get_yaxis_transform(), ha="right",
                            va="bottom", fontsize=4.5, color=c)
                ax.axhline(0, color="#999999", lw=0.3)
            else:
                ax.invert_yaxis()
    return finish(fig, output, formats)


def preservation_lollipop(z_table, size="single", height=None, output=None, formats=None):
    """Zsummary.pres of each module sorted (hdWGCNA
    PlotModulePreservationLollipop), with the thresholds 2 (moderate) and 10
    (strong)."""
    z = _pres(z_table)
    z = z.sort_values("Zsummary.pres")
    y = np.arange(len(z))
    v = z["Zsummary.pres"].to_numpy(float)
    cols = [module_color(m) for m in z.index]
    with canvas(None, size, height or max(30, 8 + 5 * len(z)), 0.7) as (fig, ax):
        ax.hlines(y, 0, v, color=cols, lw=1.0)
        ax.scatter(v, y, s=30, c=cols, edgecolors="black", linewidths=0.4, zorder=3)
        for zv, c, txt in ((2, "#0072B2", "Z = 2"), (10, "#009E73", "Z = 10")):
            ax.axvline(zv, color=c, lw=0.5, ls="--")
            ax.text(zv, 1.0, f" {txt}", transform=ax.get_xaxis_transform(), ha="left",
                    va="bottom", fontsize=4.5, color=c)
        ax.axvline(0, color="black", lw=0.4)
        ax.set_yticks(y, [module_label(m) if not is_unassigned(m) else "grey (unassigned)"
                          for m in z.index])
        ax.set_xlabel("Zsummary (preservation)")
    return finish(fig, output, formats)



def _members(members) -> pd.DataFrame:
    """metacell_members.tsv tables with a cell_type column: a file, a
    DataFrame, or a celltype-modules output directory (all its cell types)."""
    if isinstance(members, pd.DataFrame):
        t = members.copy()
    else:
        p = Path(members)
        if p.is_dir():
            tab = _io.table(p / "celltypes.tsv") if (p / "celltypes.tsv").is_file() else None
            dirs = ([(str(r["cell_type"]), p / str(r["dir"])) for _, r in tab.iterrows()]
                    if tab is not None else [(d.name, d) for d in sorted(p.iterdir())
                                             if d.is_dir()])
            parts = []
            for name, d in dirs:
                f = d / "metacell_members.tsv"
                if f.is_file():
                    parts.append(pd.read_csv(f, sep="\t").assign(cell_type=name))
            if not parts:
                raise GCENInputError(f"no metacell_members.tsv under {p}.")
            t = pd.concat(parts, ignore_index=True)
        else:
            t = _io.table(p)
    if not {"metacell", "cell"} <= set(t.columns):
        raise GCENInputError("metacell members need columns metacell and cell.")
    if "cell_type" not in t.columns:
        t["cell_type"] = "all"
    t["cell"] = t["cell"].astype(str)
    return t


def metacells(members, coords, basis: str = "X_umap", n_cells=None, size="double",
              height=None, output=None, formats=None):
    """Metacells on the cell embedding (hdWGCNA DimPlotMetacells): every cell
    in grey, each metacell at the mean position of its cells, coloured by
    cell type; and per cell type the number of metacells and the share of
    its cells that are in at least one. ``members``: metacell_members.tsv
    (or a celltype-modules output directory, all cell types); ``coords`` as
    in :func:`me_embedding`; ``n_cells``: cells per type (default: from
    celltypes.tsv of the directory, else not shown)."""
    t = _members(members)
    xy = coordinates(coords, basis=basis)
    xy.index = xy.index.map(str)
    t = t[t["cell"].isin(xy.index)]
    if t.empty:
        raise GCENInputError("no metacell member has coordinates.")
    pts = xy.loc[t["cell"], ["x", "y"]].to_numpy()
    t = t.assign(x=pts[:, 0], y=pts[:, 1])
    cen = t.groupby(["cell_type", "metacell"], sort=False)[["x", "y"]].mean().reset_index()
    types = list(pd.unique(t["cell_type"]))
    cols = {c: OKABE_ITO[i % len(OKABE_ITO)] for i, c in enumerate(types)}
    if n_cells is None and isinstance(members, (str, Path)) and Path(members).is_dir() \
            and (Path(members) / "celltypes.tsv").is_file():
        tab = _io.table(Path(members) / "celltypes.tsv")
        n_cells = dict(zip(tab["cell_type"].astype(str), tab["n_cells"]))
    used = t.groupby("cell_type", sort=False)["cell"].nunique()
    nmc = cen.groupby("cell_type", sort=False).size()
    with style():
        fig = figure(size, height or 80)
        a1, a2 = fig.subplots(1, 2, width_ratios=[2.2, 1])
        n = len(xy)
        a1.scatter(xy["x"], xy["y"], s=float(np.clip(3000 / max(n, 1), 0.3, 4)),
                   color="#D0D0D0", edgecolors="none", rasterized=n > 2000, label="cells")
        for c in types:
            m = cen[cen["cell_type"] == c]
            a1.scatter(m["x"], m["y"], s=7, color=cols[c], edgecolors="black",
                       linewidths=0.2, label=f"{c} ({len(m):,})", zorder=3)
        a1.set_xticks([])
        a1.set_yticks([])
        for sp in a1.spines.values():
            sp.set_visible(False)
        a1.set_xlabel(f"{basis.replace('X_', '').upper()} 1")
        a1.set_ylabel(f"{basis.replace('X_', '').upper()} 2")
        a1.set_title(f"{len(cen):,} metacells of {t['cell'].nunique():,} cells")
        a1.legend(loc="upper left", bbox_to_anchor=(0, -0.02), ncol=min(len(types) + 1, 4),
                  markerscale=1.5, title="metacells per cell type", frameon=False)
        y = np.arange(len(types))
        a2.barh(y, [nmc[c] for c in types], color=[cols[c] for c in types], height=0.7)
        for yi, c in zip(y, types):
            lab = f" {nmc[c]:,}"
            if n_cells is not None and c in n_cells and n_cells[c]:
                lab += f" ({used[c] / n_cells[c]:.0%} of cells)"
            a2.text(nmc[c], yi, lab, va="center", fontsize=4.5)
        a2.set_yticks(y, types)
        a2.invert_yaxis()
        a2.set_xlim(0, max(nmc) * 1.6)
        a2.set_xlabel("Metacells")
        a2.set_title("Metacells per cell type")
    return finish(fig, output, formats)


def moran(table, n_label: int = 8, size="double", height=None, output=None, formats=None):
    """Spatial autocorrelation of the genes (``spatial-build --moran``): the
    distribution of Moran's I (``filtered`` mode keeps the edges between
    genes with I > 0), I against mean expression with the most
    autocorrelated genes named, and, for ``weighted`` mode, I before and
    after the spatial smoothing."""
    t = _io.table(table)
    if not {"gene", "morans_i"} <= set(t.columns):
        raise GCENInputError("the Moran's I table needs columns gene and morans_i.")
    i = t["morans_i"].to_numpy(float)
    smooth = "morans_i_smoothed" in t.columns
    has_mean = "mean" in t.columns
    k = 1 + has_mean + smooth
    with style():
        fig = figure(size, height or 55)
        axes = np.atleast_1d(fig.subplots(1, k))
        a = axes[0]
        w = (i.max() - min(i.min(), 0)) / 50 or 0.01
        bins = np.arange(np.floor(min(i.min(), 0) / w), np.ceil(i.max() / w) + 1) * w  # 0 an edge
        a.hist(i[i <= 0], bins=bins, color="#BBBBBB", edgecolor="none", label="I ≤ 0")
        a.hist(i[i > 0], bins=bins, color=OKABE_ITO[4], edgecolor="none", label="I > 0")
        a.axvline(0, color="black", lw=0.5)
        a.set_xlabel("Moran's I")
        a.set_ylabel("Genes")
        a.set_title(f"{(i > 0).sum():,} of {len(i):,} genes with I > 0")
        a.legend(loc="upper right")
        if has_mean:
            a = axes[1]
            mu = t["mean"].to_numpy(float)
            pos = mu > 0
            a.scatter(mu[pos], i[pos], s=2, color=np.where(i[pos] > 0, OKABE_ITO[4], "#AAAAAA"),
                      edgecolors="none", rasterized=pos.sum() > 2000)
            a.set_xscale("log")
            from gcen2.plot._common import plain_log_axis

            plain_log_axis(a.xaxis)
            a.axhline(0, color="black", lw=0.5)
            a.set_xlabel("Mean expression")
            a.set_ylabel("Moran's I")
            top = t.loc[pos].nlargest(n_label, "morans_i")
            if len(top):
                from gcen2.plot._common import repel_labels

                repel_labels(a, top[["mean", "morans_i"]].to_numpy(float),
                             top["gene"].astype(str).tolist(), point_px=2)
            a.set_title("Most autocorrelated genes")
        if smooth:
            a = axes[-1]
            s2 = t["morans_i_smoothed"].to_numpy(float)
            a.scatter(i, s2, s=2, color=OKABE_ITO[2], edgecolors="none",
                      rasterized=len(i) > 2000)
            lo, hi = float(min(i.min(), s2.min())), float(max(i.max(), s2.max()))
            a.plot([lo, hi], [lo, hi], color="black", lw=0.4, ls=":")
            a.set_xlabel("Moran's I")
            a.set_ylabel("Moran's I after smoothing")
            a.set_title(f"Smoothing: median I {np.median(i):.2f} → {np.median(s2):.2f}")
    return finish(fig, output, formats)


__all__ = ["celltype_overview", "coordinates", "dme_lollipop", "dme_volcano", "groups_from",
           "me_dotplot", "me_embedding", "me_radar", "me_violin", "metacells", "moran",
           "preservation",
           "preservation_lollipop", "set_eigengenes", "spatial_me"]
