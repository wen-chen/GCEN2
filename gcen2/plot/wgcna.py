"""WGCNA / hdWGCNA module figures: soft-power choice, gene dendrogram with
module colours, eigengene network, module-trait heatmap, kME curves, module
membership vs gene significance, hub-gene networks, the module UMAP, the TOM
heatmap, the sample tree and module sizes."""

from __future__ import annotations

import json
import warnings

import numpy as np
import pandas as pd

from gcen2._matrix_io import require_file
from gcen2.exceptions import GCENInputError
from gcen2.plot import _io
from gcen2.plot._common import canvas, color_rows, draw_dendrogram, figure, finish, r_linkage
from gcen2.plot._style import (
    DIVERGING,
    GREY,
    figsize,
    is_unassigned,
    module_color,
    module_label,
    style,
)


# --------------------------------------------------------------------------
def module_sizes(modules, ax=None, size="single", height=None, include_unassigned=False,
                 output=None, formats=None):
    """Number of genes in each module (bars in the module colours, largest
    first). ``modules``: modules.tsv (or a gene -> module Series)."""
    mods = _io.modules(modules)
    counts = mods.value_counts()
    if not include_unassigned:
        counts = counts[[not is_unassigned(m) for m in counts.index]]
    counts = counts.sort_values(ascending=False)
    n = len(counts)
    h = height if height else (max(30, 6 + 4.2 * n) if n > 6 else 45)
    with canvas(ax, size, h, 0.6) as (fig, ax):
        y = np.arange(n)
        ax.barh(y, counts.to_numpy(), color=[module_color(m) for m in counts.index],
                edgecolor="black", linewidth=0.3, height=0.75)
        ax.set_yticks(y, [module_label(m) for m in counts.index])
        ax.invert_yaxis()
        ax.set_xlabel("Genes")
        for yi, v in zip(y, counts.to_numpy()):
            ax.text(v, yi, f" {v:,}", va="center", ha="left", fontsize=5)
        ax.margins(x=0.15)
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def soft_power(table, r2_cut: float | None = 0.8, size="single", height=None, output=None,
               formats=None):
    """Choice of the soft-thresholding power (WGCNA pickSoftThreshold /
    hdWGCNA PlotSoftPowers): signed scale-free fit R^2 (-sign(slope) x R^2)
    and mean connectivity for each power, the power used marked; one line per
    set for a consensus network. ``table``: soft_power.tsv."""
    t = _io.table(table)
    if t.empty:
        raise GCENInputError("the soft-power table is empty (the power was given).")
    sets = list(pd.unique(t["set"])) if "set" in t.columns else [None]
    with style():
        fig = figure(size, height or 45, layout=None)
        a1, a2 = fig.subplots(1, 2)
        from gcen2.plot._style import OKABE_ITO

        for i, s in enumerate(sets):
            d = t if s is None else t[t["set"] == s]
            c = "black" if s is None else OKABE_ITO[i % len(OKABE_ITO)]
            fit = -np.sign(d["slope"].to_numpy()) * d["r2"].to_numpy()
            for ax, y in ((a1, fit), (a2, d["mean_k"].to_numpy())):
                ax.plot(d["power"], y, color=c, lw=0.6, label=s if s is not None else None)
                ax.scatter(d["power"], y, s=4, color=c, zorder=3, linewidths=0)
                if "selected" in d.columns and d["selected"].astype(bool).any():
                    sel = d["selected"].astype(bool).to_numpy()
                    ax.scatter(d["power"][sel], y[sel], s=22, facecolors="none",
                               edgecolors="#D55E00", linewidths=0.8, zorder=4)
        if r2_cut is not None:
            a1.axhline(r2_cut, color="#D55E00", lw=0.5, ls="--")
        a1.set_xlabel("Soft power")
        a1.set_ylabel("Scale-free fit (signed R²)")
        a2.set_xlabel("Soft power")
        a2.set_ylabel("Mean connectivity")
        a2.set_yscale("log")
        sel = t[t.get("selected", pd.Series(False, index=t.index)).astype(bool)]
        from matplotlib.lines import Line2D

        handles = [Line2D([], [], marker="o", ls="", mfc="none", mec="#D55E00", ms=5,
                          label="power used")]
        if r2_cut is not None:
            handles.append(Line2D([], [], color="#D55E00", ls="--", lw=0.5,
                                  label=f"R² = {r2_cut:g}"))
        if sets != [None]:
            used = {str(r["set"]): r["power"] for _, r in sel.iterrows()}
            handles += [Line2D([], [], color=OKABE_ITO[i % len(OKABE_ITO)], lw=0.8,
                               marker="o", ms=2,
                               label=f"{s} (power {used[str(s)]:g})" if str(s) in used
                               else str(s)) for i, s in enumerate(sets)]
        elif len(sel):
            fig.suptitle(f"Soft power used: {sel['power'].iloc[0]:g}", fontsize=6)
        ncol = 2 if len(handles) > 2 else len(handles)
        rows = -(-len(handles) // ncol)
        fig.tight_layout(w_pad=1.0, rect=(0, 0.02 + 0.055 * rows, 1, 1))
        fig.legend(handles=handles, loc="lower center", ncol=ncol, frameon=False, fontsize=5)
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def _tree(tree):
    if isinstance(tree, dict):
        return tree
    return json.loads(require_file(tree).read_text(encoding="utf-8"))


def dendrogram(tree, show_unmerged: bool = True, size="double", height=None, output=None,
               formats=None):
    """Gene dendrogram (average linkage of 1 - TOM) with the module of every
    gene below it, after and before merging close modules (WGCNA
    plotDendroAndColors / hdWGCNA PlotDendrogram). ``tree``: dendrogram.json."""
    d = _tree(tree)
    z = r_linkage(d["merge"], d["height"])
    mods = np.asarray(d["modules"])
    unm = np.asarray(d.get("unmerged", d["modules"]))
    rows = [mods]
    labels = ["Modules"]
    if show_unmerged and not np.array_equal(mods, unm):
        rows = [unm, mods]
        labels = ["Unmerged", "Merged"]
    with style():
        fig = figure(size, height or 60, layout=None)
        gs = fig.add_gridspec(2, 1, height_ratios=[4, 0.45 * len(rows)], hspace=0.04)
        ax = fig.add_subplot(gs[0])
        order = draw_dendrogram(ax, z, lw=0.3 if len(mods) < 2000 else 0.15)
        ax.set_xticks([])
        ax.set_ylabel("Height (1 − TOM)")
        ax.spines["bottom"].set_visible(False)
        cax = fig.add_subplot(gs[1])
        colors = [[module_color(int(m)) if m >= 0 else GREY for m in r[order]] for r in rows]
        color_rows(cax, colors, labels)
        n_mod = len(set(mods[mods >= 0].tolist()))
        ax.set_title(f"{len(mods):,} genes, {n_mod} modules")
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def _me_frame(eigengenes, include_unassigned=False) -> pd.DataFrame:
    me = _io.eigengenes(eigengenes)
    if "set" in me.columns:
        me = me.drop(columns=[c for c in ("set", "sample") if c in me.columns])
    keep = [c for c in me.columns if include_unassigned or not is_unassigned(c)]
    return me[keep].astype(float)


def eigengene_network(eigengenes, traits=None, annotate: bool | None = None, size="single",
                      height=None, output=None, formats=None):
    """Eigengene network (WGCNA plotEigengeneNetworks / hdWGCNA
    ModuleCorrelogram): the eigengenes clustered by 1 - correlation (average
    linkage) above the heatmap of their correlations (values printed for up
    to 12 modules). ``traits`` (samples x numeric traits) are added as extra
    rows / columns. ``eigengenes``: eigengenes.tsv or cell_eigengenes.tsv."""
    me = _me_frame(eigengenes)
    if me.shape[1] < 2:
        raise GCENInputError("the eigengene network needs at least two modules.")
    data = me.copy()
    if traits is not None:
        t = _io.table(traits, index_col=0) if not isinstance(traits, pd.DataFrame) else traits
        t.index = t.index.map(str)
        data = data.join(t.select_dtypes("number"), how="left")
    c = data.corr().to_numpy()
    from scipy.cluster.hierarchy import leaves_list, linkage
    from scipy.spatial.distance import squareform

    dist = np.clip(1 - c, 0, 2)
    np.fill_diagonal(dist, 0)
    z = linkage(squareform(dist, checks=False), "average")
    order = leaves_list(z)
    names = list(data.columns)
    n = len(names)
    annotate = n <= 12 if annotate is None else annotate
    with style():
        w = figsize(size)[0]
        fig = figure(size, height or min(170, w / 0.0393701 * 1.05), layout=None)
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 3.2], width_ratios=[1, 0.04],
                              hspace=0.02, wspace=0.04)
        ax_t = fig.add_subplot(gs[0, 0])
        from scipy.cluster.hierarchy import dendrogram as sdend

        dd = sdend(z, no_plot=True)
        for xs, ys in zip(dd["icoord"], dd["dcoord"]):
            ax_t.plot(np.asarray(xs) / 10, ys, color="black", lw=0.5)
        ax_t.set_xlim(0, n)
        ax_t.set_xticks([])
        ax_t.set_ylabel("1 − r")
        ax_t.spines["bottom"].set_visible(False)
        ax = fig.add_subplot(gs[1, 0])
        cm = c[np.ix_(order, order)]
        im = ax.imshow(cm, cmap=DIVERGING, vmin=-1, vmax=1, aspect="auto",
                       extent=(0, n, n, 0), interpolation="nearest")
        lab = [module_label(names[i]) if names[i] in me.columns else names[i] for i in order]
        ax.set_xticks(np.arange(n) + 0.5, lab, rotation=90)
        ax.set_yticks(np.arange(n) + 0.5, lab)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        if annotate:
            for a in range(n):
                for b in range(n):
                    v = cm[a, b]
                    ax.text(b + 0.5, a + 0.5, f"{v:.2f}", ha="center", va="center", fontsize=5,
                            color="white" if abs(v) > 0.6 else "black")
        # module colour strips along the axes
        for k, i in enumerate(order):
            col = module_color(names[i]) if names[i] in me.columns else "white"
            ax.add_patch(_rect(k, -0.35, 1, 0.3, col, ax))
        ax.set_ylim(n, -0.4)
        cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 1]))
        cb.set_label("Pearson r")
        cb.outline.set_linewidth(0.4)
    return finish(fig, output, formats)


def _color_chips(ax, mods, fig):
    """A small square in each module's colour just left of its y label."""
    from matplotlib.transforms import ScaledTranslation

    tr = ax.get_yaxis_transform() + ScaledTranslation(-4.5 / 72, 0, fig.dpi_scale_trans)
    ax.scatter([0] * len(mods), range(len(mods)), marker="s", s=16,
               c=[module_color(m) for m in mods], edgecolors="black", linewidths=0.3,
               transform=tr, clip_on=False, zorder=5)


def _rect(x, y, w, h, color, ax):
    from matplotlib.patches import Rectangle

    return Rectangle((x, y), w, h, facecolor=color, edgecolor="black", linewidth=0.2,
                     clip_on=False, transform=ax.transData)


# --------------------------------------------------------------------------
def _stars(p):
    if not np.isfinite(p):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def module_trait(trait_cor, annotate: str = "p", use_fdr: bool | None = None, group=None,
                 size="single", height=None, output=None, formats=None):
    """Module-trait correlations (WGCNA labeledHeatmap / hdWGCNA
    PlotModuleTraitCorrelation): modules x traits coloured by r, each cell
    printed with r and the p-value (``annotate="p"``) or significance stars
    (``"stars"``) or nothing (``"none"``). ``trait_cor``: the table of
    ``module -e --trait-cor`` (module, trait, cor, p, n) or of celltype-modules
    (group, trait, module, cor, p, fdr, n; one panel per group, or only
    ``group``). ``use_fdr``: stars / values from the FDR (default: when
    present)."""
    t = _io.table(trait_cor)
    if annotate not in ("p", "stars", "none"):
        raise GCENInputError("annotate must be p, stars or none.")
    t = t[[not is_unassigned(m) for m in t["module"]]]
    use_fdr = ("fdr" in t.columns) if use_fdr is None else use_fdr
    pcol = "fdr" if use_fdr and "fdr" in t.columns else "p"
    groups = [None]
    if "group" in t.columns:
        groups = [group] if group is not None else list(pd.unique(t["group"]))
    mods = list(pd.unique(t["module"]))
    traits = list(pd.unique(t["trait"]))
    ncol = min(len(groups), 3)
    nrow = int(np.ceil(len(groups) / ncol))
    cell_w = 9 if annotate == "p" else 6
    w_mm = min(183, 28 + ncol * (len(traits) * cell_w + 6) + 14)
    h_mm = min(170, nrow * (len(mods) * (6.5 if annotate == "p" else 5) + 20) + 4)
    with style():
        fig = figure(w_mm if size == "single" else size, height or h_mm)
        from matplotlib.patches import Rectangle

        axes = np.atleast_1d(fig.subplots(nrow, ncol, squeeze=False, sharey=True)).ravel()
        im = None
        for n_ax, (ax, g) in enumerate(zip(axes, groups)):
            d = t if g is None else t[t["group"] == g]
            r = d.pivot(index="module", columns="trait", values="cor").reindex(index=mods,
                                                                              columns=traits)
            p = d.pivot(index="module", columns="trait", values=pcol).reindex(index=mods,
                                                                             columns=traits)
            im = ax.imshow(r.to_numpy(float), cmap=DIVERGING, vmin=-1, vmax=1, aspect="auto",
                           interpolation="nearest")
            for i in range(len(mods)):
                for j in range(len(traits)):
                    v, pv = r.iat[i, j], p.iat[i, j]
                    if not np.isfinite(v):
                        ax.add_patch(Rectangle((j - 0.5, i - 0.5), 1, 1, fc="#EEEEEE",
                                                   ec="none", lw=0))
                        ax.text(j, i, "NA", ha="center", va="center", fontsize=4.5,
                                color="#777777")
                        continue
                    if annotate == "p":
                        txt = f"{v:.2f}\n({pv:.0e})" if pv < 0.001 else f"{v:.2f}\n({pv:.2f})"
                    elif annotate == "stars":
                        txt = _stars(pv)
                    else:
                        txt = ""
                    ax.text(j, i, txt, ha="center", va="center", fontsize=4.5,
                            color="white" if abs(v) > 0.6 else "black", linespacing=0.9)
            ax.set_xticks(range(len(traits)), traits, rotation=45, ha="right")
            ax.set_yticks(range(len(mods)), [module_label(m) for m in mods])
            ax.tick_params(length=0)
            for sp in ax.spines.values():
                sp.set_visible(False)
            if n_ax % ncol == 0:
                _color_chips(ax, mods, fig)
                ax.tick_params(axis="y", pad=8)
            if g is not None:
                ax.set_title(str(g))
        for ax in axes[len(groups):]:
            ax.axis("off")
        cb = fig.colorbar(im, ax=axes[:len(groups)].tolist(), shrink=0.7, aspect=18)
        cb.set_label("Correlation")
        cb.outline.set_linewidth(0.4)
        if annotate == "stars":
            fig.supxlabel(f"* {pcol} < 0.05, ** < 0.01, *** < 0.001", fontsize=5)
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def _kme_columns(k: pd.DataFrame) -> dict:
    """module name -> kME column of a kme table (GCEN2 kME<k> / kME_unassigned,
    hdWGCNA kME_<colour>)."""
    out = {}
    for c in k.columns:
        if c.startswith("p."):
            continue
        if c.startswith("kME_"):
            name = c[4:]
        elif c.startswith("kME"):
            name = "ME" + c[3:]
        else:
            continue
        if not is_unassigned(name):
            out[name] = c
    return out


def _gene_module(k: pd.DataFrame) -> pd.Series:
    m = k["module"].astype(str)
    return m.where(~m.str.fullmatch(r"\d+"), "ME" + m)


def kme(kme_table, modules=None, n_label: int = 5, ncols: int = 4, size="double",
        height=None, output=None, formats=None):
    """kME of the genes of each module, sorted (hdWGCNA PlotKMEs), the top
    ``n_label`` hub genes named. ``kme_table``: kme.tsv or cell_kme.tsv (or
    hdWGCNA's modules table with kME_<colour> columns)."""
    k = _io.kme(kme_table)
    cols = _kme_columns(k)
    gm = _gene_module(k)
    mods = [m for m in cols if (gm == m).any()] if modules is None else list(modules)
    if not mods:
        raise GCENInputError("no module with kME values.")
    ncols = min(ncols, len(mods))
    nrows = int(np.ceil(len(mods) / ncols))
    with style():
        fig = figure(size, height or min(170, 30 * nrows + 8), layout=None)
        axes = np.atleast_1d(fig.subplots(nrows, ncols, squeeze=False)).ravel()
        for ax, m in zip(axes, mods):
            v = k.loc[gm == m, cols[m]].astype(float).sort_values(ascending=False)
            x = np.arange(1, len(v) + 1)
            col = module_color(m)
            ax.fill_between(x, 0, v.to_numpy(), color=col, alpha=0.35, lw=0)
            ax.plot(x, v.to_numpy(), color="black" if col.upper() == "#FFFF00" else col, lw=0.8)
            if n_label:
                ax.text(0.98, 0.97, "\n".join(map(str, v.index[:n_label])),
                        transform=ax.transAxes, ha="right", va="top", fontsize=4.5,
                        linespacing=1.1)
            ax.set_title(module_label(m))
            ax.set_xlabel("Gene rank")
            ax.set_ylabel("kME")
            ax.set_xlim(0, len(v) * 1.02)
        for ax in axes[len(mods):]:
            ax.axis("off")
        fig.tight_layout(h_pad=0.8, w_pad=0.6)
    return finish(fig, output, formats)


def mm_gs(kme_table, gs_table, module, trait, ax=None, size="single", height=None,
          output=None, formats=None):
    """Module membership (``|kME|``) against gene significance (``|GS|`` for a trait)
    for the genes of one module, with their correlation and p-value (WGCNA
    verboseScatterplot). ``gs_table``: the file of ``module -e --gs``."""
    from gcen2.tom_modules import cor_p_student

    k = _io.kme(kme_table)
    g = _io.table(gs_table, index_col=0)
    g.index = g.index.map(str)
    cols = _kme_columns(k)
    m = str(module)
    if m.isdigit():
        m = "ME" + m
    if m not in cols:
        raise GCENInputError(f"no kME for module {module}.")
    gcol = f"GS.{trait}"
    if gcol not in g.columns:
        raise GCENInputError(f"no gene significance for trait {trait}.")
    genes = k.index[_gene_module(k) == m]
    x = k.loc[genes, cols[m]].abs().to_numpy(float)
    y = g.reindex(genes)[gcol].abs().to_numpy(float)
    ok = np.isfinite(x) & np.isfinite(y)
    r = float(np.corrcoef(x[ok], y[ok])[0, 1]) if ok.sum() > 2 else np.nan
    p = float(cor_p_student(np.array([r]), ok.sum())[0]) if np.isfinite(r) else np.nan
    with canvas(ax, size, height or 60, 0.8) as (fig, ax):
        col = module_color(m)
        ax.scatter(x[ok], y[ok], s=4, color=col, edgecolors="black", linewidths=0.15)
        ax.set_xlabel(f"Module membership in {module_label(m)}")
        ax.set_ylabel(f"Gene significance for {trait}")
        ax.set_title(f"r = {r:.2f}, p = {p:.1e}")
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def _hub_order(k: pd.DataFrame, m: str) -> list[str]:
    cols = _kme_columns(k)
    v = k.loc[_gene_module(k) == m, cols[m]].astype(float)
    return list(v.sort_values(ascending=False).index)


def hub_network(tom_hubs, kme_table, module, n_inner: int = 10, n_outer: int = 15,
                edge_prop: float = 0.4, ax=None, size="single", height=None, output=None,
                formats=None):
    """Hub genes of one module (hdWGCNA ModuleNetworkPlot): the ``n_inner``
    genes of highest kME on an inner circle, the next ``n_outer`` on an
    outer one; the strongest ``edge_prop`` of their TOM connections drawn,
    darker and wider for higher TOM. ``tom_hubs``: tom_hubs.tsv
    (``--write-tom`` / ``--tom-hubs``); ``kme_table``: kme.tsv."""
    tom = _io.table(tom_hubs, index_col=0)
    tom.index = tom.index.map(str)
    k = _io.kme(kme_table)
    m = str(module)
    if m.isdigit():
        m = "ME" + m
    genes = [g for g in _hub_order(k, m) if g in tom.columns][:n_inner + n_outer]
    if len(genes) < 3:
        raise GCENInputError(f"too few hub genes of {module} in the TOM table.")
    sub = tom.loc[genes, genes].to_numpy(float)
    n_in = min(n_inner, len(genes))
    pos = {}
    for i, g in enumerate(genes[:n_in]):
        a = 2 * np.pi * i / n_in + np.pi / 2
        pos[g] = (0.58 * np.cos(a), 0.58 * np.sin(a))
    outer = genes[n_in:]
    for i, g in enumerate(outer):
        a = 2 * np.pi * i / max(len(outer), 1) + np.pi / 2 + np.pi / max(len(outer), 1)
        pos[g] = (np.cos(a), np.sin(a))
    col = module_color(m)
    with canvas(ax, size, height or 80, 1.0) as (fig, ax):
        iu = np.triu_indices(len(genes), 1)
        w = sub[iu]
        if w.size and np.nanmax(w) > 0:
            keep = np.argsort(-w)[:max(1, int(round(w.size * edge_prop)))]
            lo, hi = w[keep].min(), w[keep].max()
            for e in keep[::-1]:
                v = (w[e] - lo) / (hi - lo) if hi > lo else 1.0
                ga, gb = genes[iu[0][e]], genes[iu[1][e]]
                ax.plot([pos[ga][0], pos[gb][0]], [pos[ga][1], pos[gb][1]], color="#404040",
                        lw=0.15 + 0.9 * v, alpha=0.1 + 0.6 * v ** 2, zorder=1)
        xy = np.array([pos[g] for g in genes])
        ax.scatter(xy[:n_in, 0], xy[:n_in, 1], s=36, color=col, edgecolors="black",
                   linewidths=0.4, zorder=2)
        ax.scatter(xy[n_in:, 0], xy[n_in:, 1], s=16, color=col, edgecolors="black",
                   linewidths=0.3, zorder=2)
        ax.set_aspect("equal")
        ax.set_xlim(-1.3, 1.3)
        ax.set_ylim(-1.3, 1.3)
        ax.axis("off")
        # names outside each node, along the ray from the centre
        for i, (g, (x, y)) in enumerate(zip(genes, xy)):
            a = np.arctan2(y, x)
            r = 0.09 if i < n_in else 0.07
            c = np.cos(a)
            ax.text(x + r * c, y + r * np.sin(a), g, fontsize=4.5 if i < n_in else 4,
                    ha="left" if c > 0.3 else "right" if c < -0.3 else "center",
                    va="center", zorder=4, fontweight="bold" if i < n_in else "normal",
                    bbox={"boxstyle": "round,pad=0.1", "fc": "white", "ec": "none", "alpha": 0.8})
        ax.set_xlim(-1.45, 1.45)
        ax.set_ylim(-1.2, 1.2)
        ax.set_title(f"{module_label(m)}: top {len(genes)} genes by kME")
    return finish(fig, output, formats)


def hub_network_all(tom_hubs, kme_table, n_hubs: int = 6, n_other: int = 3,
                    edge_prop: float = 0.5, seed: int = 0, size="single", height=None,
                    output=None, formats=None):
    """Hub genes of all modules in one network (hdWGCNA HubGeneNetworkPlot):
    the ``n_hubs`` top-kME genes and ``n_other`` other genes of each module,
    the strongest ``edge_prop`` of their TOM connections, a force-directed
    layout (networkx); hub genes named."""
    try:
        import networkx as nx
    except ImportError as exc:
        raise GCENInputError("hub_network_all needs networkx (pip install "
                             "gcen2[visualization]).") from exc
    tom = _io.table(tom_hubs, index_col=0)
    tom.index = tom.index.map(str)
    k = _io.kme(kme_table)
    rng = np.random.default_rng(seed)
    gm = _gene_module(k)
    nodes, hubs, mod_of = [], set(), {}
    for m in _kme_columns(k):
        order = [g for g in _hub_order(k, m) if g in tom.index]
        h = [g for g in order if g in tom.columns][:n_hubs]
        rest = [g for g in order if g not in h]
        others = list(rng.choice(rest, size=min(n_other, len(rest)), replace=False)) if rest \
            else []
        for g in h + others:
            nodes.append(g)
            mod_of[g] = m
        hubs.update(h)
    del gm
    edges = []
    for a in nodes:
        for b in nodes:
            if a < b and b in tom.columns:
                edges.append((a, b, float(tom.at[a, b])))
            elif a < b and a in tom.columns:
                edges.append((a, b, float(tom.at[b, a])))
    edges.sort(key=lambda e: -e[2])
    edges = edges[:max(1, int(len(edges) * edge_prop))]
    graph = nx.Graph()
    graph.add_nodes_from(nodes)
    graph.add_weighted_edges_from(edges)
    pos = nx.spring_layout(graph, weight="weight", seed=seed, k=1.2 / np.sqrt(len(nodes)))
    wmax = max(e[2] for e in edges) if edges else 1
    with canvas(None, size, height or 85, 1.0) as (fig, ax):
        for a, b, w in sorted(edges, key=lambda e: e[2]):
            v = w / wmax
            ax.plot([pos[a][0], pos[b][0]], [pos[a][1], pos[b][1]], color="#555555",
                    lw=0.15 + 0.8 * v, alpha=0.08 + 0.5 * v, zorder=1)
        for g in nodes:
            col = module_color(mod_of[g])
            hub = g in hubs
            ax.scatter(*pos[g], s=40 if hub else 12, color=col, edgecolors="black",
                       linewidths=0.3, zorder=2)
        ax.axis("off")
        ax.set_aspect("equal")
        from gcen2.plot._common import repel_labels

        hl = [g for g in nodes if g in hubs]
        repel_labels(ax, np.array([pos[g] for g in hl]), hl, fontsize=4.5, point_px=3.5)
        mods = list(dict.fromkeys(mod_of[g] for g in nodes))
        from matplotlib.lines import Line2D

        ax.legend([Line2D([], [], marker="o", ls="", color=module_color(m), markersize=4,
                          markeredgecolor="black", markeredgewidth=0.3) for m in mods],
                  [module_label(m) for m in mods], loc="upper left",
                  bbox_to_anchor=(1.0, 1.0))
    return finish(fig, output, formats)


def module_umap(tom_hubs, kme_table, n_hubs: int = 10, n_label: int = 3,
                n_neighbors: int = 15, min_dist: float = 0.1, seed: int = 0,
                size="single", height=None, output=None, formats=None):
    """Genes of all modules in two dimensions (hdWGCNA RunModuleUMAP +
    ModuleUMAPPlot): UMAP of each gene's TOM with the ``n_hubs`` top-kME genes
    of every module; points in the module colours, sized by kME, the top
    ``n_label`` hubs of each module named. Needs umap-learn (else a PCA of
    the same values, with a warning)."""
    tom = _io.table(tom_hubs, index_col=0)
    tom.index = tom.index.map(str)
    k = _io.kme(kme_table)
    cols = _kme_columns(k)
    gm = _gene_module(k).reindex(tom.index)
    feats, labels = [], []
    for m in cols:
        h = [g for g in _hub_order(k, m) if g in tom.columns]
        feats.extend(h[:n_hubs])
        labels.extend(h[:n_label])
    x = tom[feats].to_numpy(float)
    try:
        import umap

        kw = {"n_neighbors": n_neighbors, "min_dist": min_dist, "metric": "cosine",
              "random_state": seed, "n_jobs": 1}
        try:
            model = umap.UMAP(compatibility_layout=True, **kw)
        except TypeError:  # umap-learn before compatibility_layout
            model = umap.UMAP(**kw)
        emb = model.fit_transform(x)
        axis = "UMAP"
    except ImportError:
        warnings.warn("umap-learn is not installed: a PCA instead (pip install umap-learn).",
                      stacklevel=2)
        xc = x - x.mean(axis=0)
        u, s, _ = np.linalg.svd(xc, full_matrices=False)
        emb = u[:, :2] * s[:2]
        axis = "PC"
    kv = np.array([float(k.at[g, cols[gm[g]]]) if gm[g] in cols else np.nan
                   for g in tom.index])
    kv = np.nan_to_num(kv, nan=0.0)
    with canvas(None, size, height or 80, 1.0) as (fig, ax):
        order = np.argsort(kv)
        sc = ax.scatter(emb[order, 0], emb[order, 1], s=1 + 10 * np.clip(kv[order], 0, 1) ** 2,
                        c=[module_color(gm.iat[i]) for i in order], edgecolors="none")
        if len(order) > 5000:
            sc.set_rasterized(True)
        pos = {g: i for i, g in enumerate(tom.index)}
        _stack_labels(ax, emb, pos, labels, k)
        ax.set_xlabel(f"{axis} 1")
        ax.set_ylabel(f"{axis} 2")
        ax.set_xticks([])
        ax.set_yticks([])
        from matplotlib.lines import Line2D

        ms = list(cols)
        ax.legend([Line2D([], [], marker="o", ls="", color=module_color(m), markersize=4)
                   for m in ms], [module_label(m) for m in ms], loc="upper left",
                  bbox_to_anchor=(1.0, 1.0))
    return finish(fig, output, formats)


def _stack_labels(ax, emb, pos, labels, k):
    """Hub-gene names of each module stacked in a column beside the module's
    hubs, with thin leader lines (no overlapping names)."""
    span = np.ptp(emb, axis=0)
    gm = _gene_module(k)
    by_mod = {}
    for g in labels:
        by_mod.setdefault(gm.get(g, "grey"), []).append(g)
    for gs in by_mod.values():
        pts = np.array([emb[pos[g]] for g in gs])
        cx, cy = pts.mean(axis=0)
        side = 1 if cx < emb[:, 0].mean() + span[0] * 0.25 else -1
        tx = cx + side * 0.12 * span[0]
        step = 0.045 * span[1]
        for i, g in enumerate(gs):
            ty = cy + (i - (len(gs) - 1) / 2) * step
            ax.plot([emb[pos[g], 0], tx], [emb[pos[g], 1], ty], color="#606060", lw=0.3)
            ax.text(tx, ty, f" {g} " if side > 0 else f"{g} ", fontsize=4.5, va="center",
                    ha="left" if side > 0 else "right")


# --------------------------------------------------------------------------
def tom_heatmap(expr, modules, power: float, signed: bool = True, n_genes: int = 400,
                seed: int = 0, size="single", height=None, output=None, formats=None):
    """TOM of a random sample of ``n_genes`` genes (WGCNA TOMplot):
    (1 - TOM)^7 ordered by their average-linkage dendrogram, module colours
    along both sides. ``expr``: genes x samples (metacells.tsv, the module -e
    input); ``power`` and ``signed`` as the network was built."""
    from gcen2.core import require_extension

    if isinstance(expr, pd.DataFrame):
        x = expr
    else:
        from gcen2._matrix_io import read_matrix

        m, g, s = read_matrix(expr)
        x = pd.DataFrame(m, index=g, columns=s)
    mods = _io.modules(modules)
    rng = np.random.default_rng(seed)
    pick = rng.choice(x.shape[0], size=min(n_genes, x.shape[0]), replace=False)
    sub = x.iloc[np.sort(pick)]
    sub = sub[sub.std(axis=1) > 0]
    d, _ = require_extension().module.tom_dissimilarity(np.ascontiguousarray(sub.to_numpy(float)),
                                            "pearson", signed, float(power), [], 2)
    merge, height_ = require_extension().module.hclust_average(d)
    z = r_linkage(merge, height_)
    from scipy.cluster.hierarchy import leaves_list

    order = leaves_list(z)
    names = list(sub.index)
    lab = [mods.get(names[i], "grey") for i in order]
    with style():
        fig = figure(size, height or figsize(size)[0] / 0.0393701, layout=None)
        gs = fig.add_gridspec(2, 3, width_ratios=[0.04, 1, 0.03], height_ratios=[0.04, 1],
                              wspace=0.02, hspace=0.01)
        ax = fig.add_subplot(gs[1, 1])
        pm = d[np.ix_(order, order)] ** 7
        np.fill_diagonal(pm, np.nan)
        lo, hi = np.nanquantile(pm, [0.01, 0.99])
        im = ax.imshow(pm, cmap="YlOrRd_r", vmin=lo, vmax=hi, aspect="auto",
                       interpolation="nearest")
        im.set_rasterized(True)
        ax.set_xticks([])
        ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        cols = [module_color(m) for m in lab]
        top = fig.add_subplot(gs[0, 1])
        color_rows(top, [cols], [""])
        left = fig.add_subplot(gs[1, 0])
        from matplotlib.colors import to_rgb

        left.imshow(np.array([[to_rgb(c)] for c in cols]), aspect="auto",
                    interpolation="nearest")
        left.axis("off")
        top.set_title(f"TOM of {len(names)} random genes")
        cb = fig.colorbar(im, cax=fig.add_subplot(gs[1, 2]))
        cb.set_label(r"$(1 - \mathrm{TOM})^{7}$ (dark = close)")
        cb.outline.set_linewidth(0.3)
    return finish(fig, output, formats)


# --------------------------------------------------------------------------
def sample_tree(expr, traits=None, cut_height: float | None = None, size="double",
                height=None, output=None, formats=None):
    """Clustering of the samples (average linkage, Euclidean distance) to spot
    outliers, with trait values as colour bars below (WGCNA
    plotDendroAndColors of the samples); ``cut_height`` drawn as a line.
    ``expr``: genes x samples."""
    from scipy.cluster.hierarchy import linkage
    from scipy.spatial.distance import pdist

    if isinstance(expr, pd.DataFrame):
        x = expr
    else:
        from gcen2._matrix_io import read_matrix

        m, g, s = read_matrix(expr)
        x = pd.DataFrame(m, index=g, columns=s)
    z = linkage(pdist(x.T.to_numpy(float)), "average")
    rows, labels = [], []
    if traits is not None:
        import matplotlib as mpl

        t = _io.table(traits, index_col=0) if not isinstance(traits, pd.DataFrame) else traits
        t.index = t.index.map(str)
        t = t.reindex(x.columns.map(str))
        cmap = mpl.colormaps["viridis"]
        for c in t.columns:
            v = pd.to_numeric(t[c], errors="coerce")
            if v.notna().sum() == 0:
                codes = pd.Categorical(t[c].astype(str))
                v = pd.Series(codes.codes, index=t.index).astype(float)
            lo, hi = np.nanmin(v), np.nanmax(v)
            norm = (v - lo) / (hi - lo) if hi > lo else v * 0
            rows.append([mpl.colors.to_hex(cmap(a)) if np.isfinite(a) else "#FFFFFF"
                         for a in norm])
            labels.append(str(c))
    with style():
        fig = figure(size, height or (55 + 3 * len(rows)), layout=None)
        if rows:
            gs = fig.add_gridspec(2, 1, height_ratios=[4, 0.4 * len(rows)], hspace=0.04)
            ax = fig.add_subplot(gs[0])
        else:
            ax = fig.add_subplot(111)
        order = draw_dendrogram(ax, z, lw=0.4)
        n = x.shape[1]
        if n <= 120:
            ax.set_xticks(np.arange(n) + 0.5, [str(x.columns[i]) for i in order], rotation=90,
                          fontsize=4)
        else:
            ax.set_xticks([])
        ax.set_ylabel("Height (Euclidean)")
        if cut_height is not None:
            ax.axhline(cut_height, color="#D55E00", lw=0.6, ls="--")
        ax.set_title(f"Sample clustering ({n} samples)")
        if rows:
            ax.set_xticks([])
            cax = fig.add_subplot(gs[1])
            color_rows(cax, [[r[i] for i in order] for r in rows], labels)
    return finish(fig, output, formats)
