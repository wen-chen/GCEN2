"""Co-expression network figures: degree distribution with the scale-free
fit, edge weights, a network view, Jaccard-module sizes, cell-type networks
and a network against its shuffled (degree-preserving) nulls."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from gcen2._matrix_io import require_file
from gcen2.exceptions import GCENInputError
from gcen2.plot import _io
from gcen2.plot._common import canvas, figure, finish
from gcen2.plot._style import OKABE_ITO, module_color, module_label, style


def _degrees(network) -> np.ndarray:
    """Node degrees (sorted, so the figure does not depend on the order the
    source lists them) of a network file (C++ reader), a Network or an edge
    DataFrame, or the degrees themselves."""
    return np.sort(_degrees_any(network))


def _degrees_any(network) -> np.ndarray:
    if isinstance(network, (np.ndarray, list, pd.Series)):
        return np.asarray(network, dtype=float)
    if isinstance(network, pd.DataFrame):
        a, b = network.columns[:2]
        return pd.concat([network[a], network[b]]).value_counts().to_numpy(float)
    if hasattr(network, "to_networkx"):  # a Network: degrees as networkx G.degree
        return np.array([d for _, d in network.to_networkx().degree()], dtype=float)
    p = require_file(network)
    from gcen2.network import summary_from_file

    s = summary_from_file(p)
    if hasattr(s, "degrees"):
        return np.asarray(s.degrees, dtype=float)
    return _degrees(s)


def scale_free_fit(k: np.ndarray, n_breaks: int = 10):
    """WGCNA scaleFreeFitIndex: degrees cut into ``n_breaks`` equal bins; the
    log10 frequency regressed on log10 mean degree. Returns (log10 k,
    log10 p(k), R^2, slope)."""
    k = np.asarray(k, dtype=float)
    k = k[k > 0]
    edges = np.linspace(k.min(), k.max(), n_breaks + 1)
    idx = np.clip(np.digitize(k, edges[1:-1], right=True), 0, n_breaks - 1)
    dk, pk = [], []
    for b in range(n_breaks):
        sel = idx == b
        mid = (edges[b] + edges[b + 1]) / 2
        dk.append(k[sel].mean() if sel.any() else mid)
        pk.append(sel.mean())
    dk, pk = np.array(dk), np.array(pk)
    ok = pk > 0
    x, y = np.log10(dk[ok]), np.log10(pk[ok] + 1e-9)
    if ok.sum() < 3:
        return x, y, np.nan, np.nan
    slope, icpt = np.polyfit(x, y, 1)
    r2 = 1 - ((y - (slope * x + icpt)) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    return x, y, float(r2), float(slope)


def degree_distribution(network, size="double", height=None, output=None, formats=None):
    """Degree distribution of a network: histogram of node degrees and the
    scale-free plot (log10 p(k) against log10 k with the fitted line, R^2 and
    slope; WGCNA scaleFreePlot). ``network``: an edge-list file (read by the
    C++ core), a Network, or the degrees."""
    k = _degrees(network)
    if k.size == 0:
        raise GCENInputError("the network has no edges.")
    x, y, r2, slope = scale_free_fit(k)
    with style():
        fig = figure(size, height or 55)
        a1, a2 = fig.subplots(1, 2)
        a1.hist(k, bins=min(60, max(10, int(np.sqrt(k.size)))), color=OKABE_ITO[4],
                edgecolor="white", linewidth=0.3)
        a1.set_xlabel("Degree")
        a1.set_ylabel("Genes")
        a1.set_title(f"{k.size:,} genes, mean degree {k.mean():.1f}")
        a2.scatter(x, y, s=12, color="black", zorder=3)
        if np.isfinite(slope):
            xs = np.array([x.min(), x.max()])
            a2.plot(xs, slope * xs + (y.mean() - slope * x.mean()), color=OKABE_ITO[5], lw=0.8)
        a2.set_xlabel("log10 degree k")
        a2.set_ylabel("log10 p(k)")
        a2.set_title(f"scale-free R² = {r2:.2f}, slope = {slope:.2f}")
    return finish(fig, output, formats)


def _edges(network, usecols=None) -> pd.DataFrame:
    if isinstance(network, pd.DataFrame):
        return network
    if hasattr(network, "edges"):
        return network.edges
    p = require_file(network)
    head = pd.read_csv(p, sep="\t", nrows=0)
    if not set(map(str, head.columns)) >= {"gene_a", "gene_b"} and not \
            str(head.columns[0]).startswith("gene"):
        t = pd.read_csv(p, sep="\t", header=None, comment="#")
        t.columns = ["gene_a", "gene_b", "correlation", "p_value"][:t.shape[1]] + \
            [f"c{i}" for i in range(4, t.shape[1])]
        return t
    return pd.read_csv(p, sep="\t", usecols=usecols)


def edge_weights(network, column: str = "correlation", size="single", height=None,
                 output=None, formats=None):
    """Distribution of the edge weights (correlations, or another column such
    as the TOM or a fused score)."""
    t = _edges(network)
    if column not in t.columns:
        raise GCENInputError(f"no column '{column}' (have: {', '.join(map(str, t.columns))}).")
    v = pd.to_numeric(t[column], errors="coerce").dropna().to_numpy(float)
    with canvas(None, size, height or 50, 0.7) as (fig, ax):
        name = column.replace("_", " ")
        if (v < 0).any() and (v > 0).any():
            # signed weights: |weight| on one axis, positive and negative stacked
            a = np.abs(v)
            bins = np.linspace(a.min(), a.max(), 60)
            ax.hist([a[v > 0], a[v < 0]], bins=bins, stacked=True,
                    color=[OKABE_ITO[5], OKABE_ITO[4]], edgecolor="none",
                    label=[f"positive ({(v > 0).sum():,})", f"negative ({(v < 0).sum():,})"])
            ax.set_xlabel(f"|{name}|")
            ax.legend(loc="upper left")
        else:
            ax.hist(v, bins=60, color=OKABE_ITO[2], edgecolor="none")
            ax.set_xlabel(name.capitalize())
        ax.set_ylabel("Edges")
        ax.set_title(f"{v.size:,} edges")
    return finish(fig, output, formats)


def network_view(network, genes=None, modules=None, max_nodes: int = 150,
                 weight: str = "correlation", n_label: int = 12, k_strongest: int = 5,
                 seed: int = 0, size="single", height=None, output=None, formats=None):
    """A drawing of (part of) a network: the given ``genes`` (or the
    ``max_nodes`` best-connected genes) with, for each gene, its
    ``k_strongest`` edges (so dense modules show their structure rather than
    a solid ball), force-directed layout (networkx), node size by degree,
    colours by module (``modules``: modules.tsv), the ``n_label``
    highest-degree genes named."""
    try:
        import networkx as nx
    except ImportError as exc:
        raise GCENInputError("network_view needs networkx (pip install "
                             "gcen2[visualization]).") from exc
    t = _edges(network)
    a, b = t.columns[:2]
    if genes is not None:
        keep = set(map(str, genes))
        t = t[t[a].astype(str).isin(keep) & t[b].astype(str).isin(keep)]
    else:
        deg = pd.concat([t[a], t[b]]).value_counts()
        keep = set(map(str, deg.index[:max_nodes]))
        t = t[t[a].astype(str).isin(keep) & t[b].astype(str).isin(keep)]
    if t.empty:
        raise GCENInputError("no edges among the selected genes.")
    t = t.assign(_w=t[weight].abs().to_numpy(float) if weight in t.columns else 1.0)
    if k_strongest:
        long = pd.concat([t.assign(_n=t[a].astype(str)), t.assign(_n=t[b].astype(str))])
        top = long.sort_values("_w", ascending=False).groupby("_n").head(k_strongest)
        t = top.drop(columns="_n").drop_duplicates(subset=[a, b])
    g = nx.Graph()
    w = t["_w"].to_numpy(float)
    for (u, v), x in zip(zip(t[a].astype(str), t[b].astype(str)), w):
        g.add_edge(u, v, weight=x)
    from gcen2.plot._common import pack_components, repel_labels

    pos = pack_components(g, seed=seed)
    mods = _io.modules(modules) if modules is not None else None
    deg = dict(g.degree())
    with canvas(None, size, height or 85, 1.0) as (fig, ax):
        from matplotlib.collections import LineCollection

        segs = [(pos[u], pos[v]) for u, v in g.edges()]
        ws = np.array([d["weight"] for _, _, d in g.edges(data=True)])
        wn = (ws - ws.min()) / (np.ptp(ws) or 1)
        ax.add_collection(LineCollection(segs, colors="#606060", linewidths=0.2 + 0.6 * wn,
                                         alpha=0.35, zorder=1))
        nodes = list(g.nodes())
        cols = [module_color(mods.get(n, "grey")) if mods is not None else OKABE_ITO[4]
                for n in nodes]
        dv = np.array([deg[n] for n in nodes], dtype=float)
        ax.scatter([pos[n][0] for n in nodes], [pos[n][1] for n in nodes],
                   s=6 + 40 * dv / dv.max(), c=cols, edgecolors="black", linewidths=0.25,
                   zorder=2)
        ax.set_aspect("equal")
        ax.axis("off")
        ax.autoscale_view()
        top = sorted(nodes, key=lambda x: -deg[x])[:n_label]
        repel_labels(ax, np.array([pos[n] for n in top]), top, point_px=3)
        ax.set_title(f"{g.number_of_nodes()} genes, {g.number_of_edges()} edges")
        if mods is not None:
            from matplotlib.lines import Line2D

            present = list(dict.fromkeys(mods.get(n, "grey") for n in nodes))[:12]
            ax.legend([Line2D([], [], marker="o", ls="", color=module_color(m), markersize=4,
                              markeredgecolor="black", markeredgewidth=0.3) for m in present],
                      [module_label(m) for m in present], loc="upper left",
                      bbox_to_anchor=(1.0, 1.0))
    return finish(fig, output, formats)


def module_size_hist(modules, size="single", height=None, output=None, formats=None):
    """Sizes of many modules (e.g. Jaccard modules of ``gcen2 module -i``):
    histogram on a log scale, with the number of modules and genes."""
    mods = _io.modules(modules)
    sz = mods.value_counts().to_numpy(float)
    with canvas(None, size, height or 50, 0.7) as (fig, ax):
        bins = np.logspace(np.log10(max(sz.min(), 1)), np.log10(sz.max() + 1), 30)
        ax.hist(sz, bins=bins, color=OKABE_ITO[0], edgecolor="white", linewidth=0.3)
        ax.set_xscale("log")
        ax.set_yscale("log")  # a few large modules stay visible next to many small ones
        from gcen2.plot._common import plain_log_axis

        plain_log_axis(ax.xaxis)
        plain_log_axis(ax.yaxis, subs=(1.0, 2.0, 5.0))
        ax.set_xlabel("Module size (genes)")
        ax.set_ylabel("Modules")
        ax.set_title(f"{sz.size:,} modules, {int(sz.sum()):,} genes")
    return finish(fig, output, formats)


def celltype_networks(directory, size="double", height=None, output=None, formats=None):
    """Networks of ``gcen2 celltype``: edges and genes per cell type, and the
    overlap of their edge sets (Jaccard index) between cell types."""
    d = Path(directory)
    tab = _io.table(d / "celltypes.tsv")
    names, sets = [], []
    for _, r in tab.iterrows():
        f = d / f"{r['dir']}.network" if "dir" in tab.columns else None
        if f is None or not f.is_file():
            cands = sorted(d.glob(f"{str(r['cell_type']).replace(' ', '_')}*.network"))
            f = cands[0] if cands else None
        if f is None:
            continue
        e = _edges(f, usecols=[0, 1])
        x, y = e.iloc[:, 0].astype(str).to_numpy(), e.iloc[:, 1].astype(str).to_numpy()
        lo, hi = np.minimum(x, y), np.maximum(x, y)
        sets.append(set(zip(lo, hi)))
        names.append(str(r["cell_type"]))
    if not names:
        raise GCENInputError(f"no cell-type networks in {d}.")
    n = len(names)
    jac = np.full((n, n), np.nan)
    for i in range(n):
        for j in range(i + 1, n):
            u = len(sets[i] | sets[j])
            jac[i, j] = jac[j, i] = len(sets[i] & sets[j]) / u if u else 0
    with style():
        fig = figure(size, height or max(50, 14 + 9 * n))
        a1, a2 = fig.subplots(1, 2, width_ratios=[1, 1])
        y = np.arange(n)
        a1.barh(y, [len(s) for s in sets], color=OKABE_ITO[4], height=0.7)
        a1.set_yticks(y, names)
        a1.invert_yaxis()
        a1.set_xlabel("Edges")
        off = jac[~np.eye(n, dtype=bool)]
        vmax = float(np.nanmax(off)) if n > 1 and np.isfinite(off).any() else 1.0
        im = a2.imshow(jac, cmap="viridis", vmin=0, vmax=max(vmax, 1e-6), aspect="auto")
        a2.set_xticks(y, names, rotation=45, ha="right")
        a2.set_yticks(y, names)
        for i in range(n):
            for j in range(n):
                if i != j:
                    a2.text(j, i, f"{jac[i, j]:.3f}", ha="center", va="center", fontsize=4.5,
                            color="white" if jac[i, j] < vmax * 0.6 else "black")
                else:
                    a2.text(j, i, "–", ha="center", va="center", fontsize=6, color="#999999")
        a1.set_title("Edges per cell type")
        a2.set_title("Shared edges between cell types")
        for yi, s in zip(y, sets):
            a1.text(len(s), yi, f" {len(s):,}", va="center", fontsize=4.5)
        a1.set_xlim(0, max(len(s) for s in sets) * 1.25)
        cb = fig.colorbar(im, ax=a2, shrink=0.7)
        cb.set_label("Edge overlap (Jaccard)")
        cb.outline.set_linewidth(0.3)
    return finish(fig, output, formats)


def null_distribution(network, shuffled, metric: str = "clustering_coefficient",
                      size="single", height=None, output=None, formats=None):
    """A network statistic (default the clustering coefficient) of the real
    network against its degree-preserving shuffles (``gcen2 shuffle``
    outputs): the shuffles' values as a histogram, the real value as a line,
    with the empirical p-value and z-score."""
    from gcen2.network import stats_from_file

    real = float(getattr(stats_from_file(network), metric))
    null = np.array([float(getattr(stats_from_file(p), metric)) for p in shuffled])
    if null.size == 0:
        raise GCENInputError("give at least one shuffled network.")
    p = (1 + (null >= real).sum()) / (1 + null.size)
    z = (real - null.mean()) / null.std(ddof=1) if null.size > 1 and null.std() > 0 else np.nan
    with canvas(None, size, height or 45, 0.6) as (fig, ax):
        if null.size >= 20:
            ax.hist(null, bins=min(30, null.size // 2), color="#BBBBBB", edgecolor="white",
                    linewidth=0.3, label=f"shuffled (n = {null.size})")
            ax.set_ylabel("Networks")
        else:  # few shuffles: show each one
            jit = np.random.default_rng(0).uniform(-0.15, 0.15, null.size)
            ax.scatter(null, jit, s=14, color="#999999", edgecolors="black", linewidths=0.3,
                       label=f"shuffled (n = {null.size})", zorder=3)
            ax.set_ylim(-0.6, 0.6)
            ax.set_yticks([])
            ax.spines["left"].set_visible(False)
        ax.axvline(real, color=OKABE_ITO[5], lw=1.0, label="observed")
        lo, hi = min(null.min(), real), max(null.max(), real)
        pad = 0.1 * (hi - lo or abs(hi) or 1)
        ax.set_xlim(lo - pad, hi + pad)
        ax.set_xlabel(metric.replace("_", " "))
        title = f"p = {p:.3g}" + (f", z = {z:.1f}" if np.isfinite(z) else "")
        if null.size < 20:
            title += f" (only {null.size} shuffle{'s' if null.size > 1 else ''})"
        ax.set_title(title)
        ax.legend(loc="upper left")
    return finish(fig, output, formats)


def _edge_keys(network) -> pd.Index:
    """Undirected edge keys ("a<TAB>b", a <= b) of a network, without duplicates."""
    e = _edges(network, usecols=[0, 1])
    x = e.iloc[:, 0].astype(str).to_numpy()
    y = e.iloc[:, 1].astype(str).to_numpy()
    lo, hi = np.minimum(x, y), np.maximum(x, y)
    return pd.Index(pd.Series(lo, dtype=object) + "\t" + hi).unique()


def _set_name(x, i: int) -> str:
    return Path(x).stem if isinstance(x, (str, Path)) else f"network {i + 1}"


def merge_overlap(networks, merged=None, names=None, top: int = 15, size="single",
                  height=None, output=None, formats=None):
    """Edges shared by the input networks of ``gcen2 merge``: a Venn diagram
    for 2 or 3 networks (areas not to scale), an UpSet-style bar chart for
    more; with ``merged`` (the merge output), how many edges of each
    region pass the cutoff on the summed correlation."""
    networks = list(networks)
    if len(networks) < 2:
        raise GCENInputError("merge_overlap needs at least 2 networks.")
    names = [str(n) for n in names] if names else [_set_name(x, i)
                                                   for i, x in enumerate(networks)]
    if len(names) != len(networks):
        raise GCENInputError(f"{len(names)} names for {len(networks)} networks.")
    parts = [pd.Series(1 << i, index=_edge_keys(x)) for i, x in enumerate(networks)]
    mask = pd.concat(parts).groupby(level=0).sum()
    counts = {int(k): int(v) for k, v in mask.value_counts().items()}
    extra = kept = None
    if merged is not None:
        mk = _edge_keys(merged)
        kept = {int(k): int(v) for k, v in mask[mask.index.isin(mk)].value_counts().items()}
    cols = [OKABE_ITO[i % len(OKABE_ITO)] for i in range(len(names))]
    total = len(mask)
    title = f"{total:,} distinct edges in {len(names)} networks"
    if kept is not None:
        title += f"; {sum(kept.values()):,} in the merged network"
    with style():
        if len(names) <= 3:
            fig = figure(size, height or (62 if len(names) == 3 else 45))
            ax = fig.add_subplot(111)
            if kept is not None:
                extra = {m: f"({kept.get(m, 0):,} kept)" for m in counts}
            from gcen2.plot._common import draw_venn

            draw_venn(ax, names, counts, cols, extra=extra)
            ax.set_title(title)
        else:
            fig = figure(size, height or 70)
            from gcen2.plot._common import draw_upset

            a1, _ = draw_upset(fig, names, counts, cols, top=top, extra=kept)
            a1.set_title(title)
    return finish(fig, output, formats)


def extract_overview(network, genes, top: int = 20, size="single", height=None,
                     output=None, formats=None):
    """The query genes of ``gcen2 extract``: how many are in the network (a
    Venn diagram of the query genes and the network's genes), and the
    degree of the best-connected ones, split into edges to other query
    genes and to the rest of the network. ``network``: the full network or
    the extracted subnetwork (it holds every edge of a query gene)."""
    if isinstance(genes, (str, Path)):
        from gcen2.io import read_gene_list

        genes = read_gene_list(genes)
    q = pd.Index([str(g) for g in genes]).unique()
    e = _edges(network, usecols=[0, 1])
    a = e.iloc[:, 0].astype(str).to_numpy()
    b = e.iloc[:, 1].astype(str).to_numpy()
    net_genes = pd.Index(np.unique(np.r_[a, b]))
    found = q[q.isin(net_genes)]
    ina, inb = pd.Index(a).isin(q), pd.Index(b).isin(q)
    inner = pd.Series(np.r_[a[ina & inb], b[ina & inb]]).value_counts()
    outer = pd.Series(np.r_[a[ina & ~inb], b[inb & ~ina]]).value_counts()
    deg = pd.DataFrame({"query": inner, "other": outer}).reindex(found).fillna(0).astype(int)
    deg = deg.assign(total=deg.sum(axis=1)).sort_values("total", ascending=False).head(top)
    counts = {1: len(q) - len(found), 2: len(net_genes) - len(found), 3: len(found)}
    with style():
        fig = figure(size, height or max(55, 12 + 2.6 * len(deg)))
        a1, a2 = fig.subplots(1, 2, width_ratios=[1.3, 1])
        from gcen2.plot._common import draw_venn

        draw_venn(a1, ["query genes", "network genes"], counts, [OKABE_ITO[5], OKABE_ITO[4]],
                  fontsize=5)
        a1.set_title(f"{len(found):,} of {len(q):,} query genes found")
        y = np.arange(len(deg))
        a2.barh(y, deg["query"], color=OKABE_ITO[5], height=0.7, label="to query genes")
        a2.barh(y, deg["other"], left=deg["query"], color=OKABE_ITO[4], height=0.7,
                label="to other genes")
        a2.set_yticks(y, deg.index, fontsize=4.5)
        a2.invert_yaxis()
        a2.set_xlabel("Edges")
        a2.set_title(f"Top {len(deg)} query genes by degree")
        a2.legend(loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=2)
    return finish(fig, output, formats)
