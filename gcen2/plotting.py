"""Visualization module for GCEN2."""

from __future__ import annotations

from pathlib import Path
from typing import Optional, Union

from gcen2.exceptions import GCENInputError
from gcen2.models import DiffCoExResult, EnrichmentResult, ExpressionMatrix, ModuleSet, Network


def _require_matplotlib():
    from gcen2.plot._common import _require_matplotlib

    _require_matplotlib()
    import matplotlib.pyplot as plt

    return plt


def _require_networkx():
    try:
        import networkx as nx

        return nx
    except ImportError as exc:
        raise GCENInputError(
            "networkx is required for network plotting. "
            "Install with: pip install gcen2[visualization]"
        ) from exc


def network_interactive(
    network: Network,
    output: Union[str, Path] = "network.html",
    color_by: str = "module",
    modules: Optional[ModuleSet] = None,
    annotation: Optional[dict] = None,
) -> None:
    """Generate an interactive HTML network visualization using pyvis.

    Parameters
    ----------
    network:
        Input network.
    output:
        Output HTML file path.
    color_by:
        'module' (requires modules) or 'degree'.
    modules:
        ModuleSet for coloring nodes by module membership.
    annotation:
        Optional dict mapping gene → list of top GO terms (for hover tooltips).
    """
    try:
        from pyvis.network import Network as PyvisNetwork
    except ImportError as exc:
        raise GCENInputError(
            "pyvis is required for interactive visualization. "
            "Install with: pip install gcen2[visualization]"
        ) from exc

    G = network.to_networkx()
    n_nodes = G.number_of_nodes()
    physics = n_nodes > 500

    # in_line: the JS inside the HTML (one self-contained file; "local" writes lib/ into
    # the current directory)
    pv = PyvisNetwork(height="800px", width="100%", bgcolor="#ffffff", font_color="black",
                      cdn_resources="in_line")
    pv.force_atlas_2based() if not physics else pv.barnes_hut()

    gene_module: dict[str, int] = {}
    if modules is not None:
        for i, mod in enumerate(modules.modules):
            for gene in mod:
                gene_module[gene] = i

    import matplotlib
    import matplotlib.colors as mcolors

    degrees = dict(G.degree())
    max_deg = max(degrees.values()) if degrees else 1
    colormap = matplotlib.colormaps["tab20"]  # cm.get_cmap is gone since matplotlib 3.9

    for node in G.nodes():
        deg = degrees.get(node, 1)
        size = 5 + 20 * (deg / max_deg) ** 0.5
        if color_by == "module" and node in gene_module:
            color = mcolors.to_hex(colormap(gene_module[node] % 20))
        else:
            color = mcolors.to_hex(matplotlib.colormaps["viridis"](deg / max_deg))

        top_go = ""
        if annotation and node in annotation:
            top_go = ", ".join(annotation[node][:3])
        title = f"<b>{node}</b><br>Degree: {deg}"
        if node in gene_module:
            title += f"<br>Module: {gene_module[node]}"
        if top_go:
            title += f"<br>GO: {top_go}"

        pv.add_node(node, label=node, title=title, size=size, color=color)

    for a, b, data in G.edges(data=True):
        cor = abs(data.get("correlation", 1.0))
        pv.add_edge(a, b, width=1 + 3 * cor, title=f"r={data.get('correlation', ''):.3f}")

    # pyvis writes with the platform's default encoding (cp1252 on Windows,
    # which cannot hold every gene name or GO term): write its HTML as UTF-8
    html = pv.generate_html(name=str(output), notebook=False)
    with open(output, "w", encoding="utf-8", newline="\n") as fh:
        fh.write(html)


def network_static(
    network: Network,
    output: Optional[Union[str, Path]] = None,
    modules: Optional[ModuleSet] = None,
    figsize: tuple[int, int] = (12, 10),
):
    """Generate a static publication-quality network figure via networkx + matplotlib."""
    plt = _require_matplotlib()
    nx = _require_networkx()

    G = network.to_networkx()
    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(G, seed=42)

    gene_module: dict[str, int] = {}
    if modules is not None:
        for i, mod in enumerate(modules.modules):
            for gene in mod:
                gene_module[gene] = i

    node_colors = [gene_module.get(n, -1) for n in G.nodes()]
    nx.draw_networkx(G, pos=pos, ax=ax, node_color=node_colors, cmap="tab20",
                     node_size=50, font_size=6, edge_color="gray", alpha=0.8)
    ax.axis("off")

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def dotplot(
    enrichment: EnrichmentResult,
    top_n: int = 20,
    output: Optional[Union[str, Path]] = None,
    title: str = "",
    figsize: tuple[int, int] = (8, 8),
):
    """Dot plot (bubble plot) for enrichment results."""
    plt = _require_matplotlib()

    df = enrichment.results.head(top_n)
    if df.empty:
        raise GCENInputError("No enrichment results to plot.")

    fig, ax = plt.subplots(figsize=figsize)
    scatter = ax.scatter(
        df["enrichment"],
        range(len(df)),
        s=df["study_count"] * 20,
        c=df["p_value"],
        cmap="Reds_r",
        alpha=0.8,
    )
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["name"] if "name" in df.columns else df["term_id"])
    ax.set_xlabel("Enrichment Ratio")
    ax.set_title(title or f"{enrichment.ontology_type} Enrichment")
    plt.colorbar(scatter, label="p-value")
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def barplot(
    enrichment: EnrichmentResult,
    top_n: int = 20,
    output: Optional[Union[str, Path]] = None,
    title: str = "",
    figsize: tuple[int, int] = (8, 8),
):
    """Horizontal bar chart for enrichment results."""
    plt = _require_matplotlib()
    import numpy as np

    df = enrichment.results.head(top_n)
    if df.empty:
        raise GCENInputError("No enrichment results to plot.")

    fig, ax = plt.subplots(figsize=figsize)
    ax.barh(
        range(len(df)),
        -np.log10(df["p_value"]),
        color="steelblue",
        alpha=0.8,
    )
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(df["name"] if "name" in df.columns else df["term_id"])
    ax.set_xlabel("-log10(p-value)")
    ax.set_title(title or f"{enrichment.ontology_type} Enrichment")
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def expression_boxplot(expr: ExpressionMatrix, output=None, figsize=(12, 6)):
    """Box plots comparing expression distributions across samples."""
    plt = _require_matplotlib()

    df = expr.data if hasattr(expr, "data") else expr
    fig, ax = plt.subplots(figsize=figsize)
    df.T.boxplot(ax=ax, rot=90)
    ax.set_ylabel("Expression")
    ax.set_title("Expression Distribution per Sample")
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def expression_violin(expr: ExpressionMatrix, output=None, figsize=(12, 6)):
    """Violin plots of expression distributions."""
    plt = _require_matplotlib()
    try:
        import seaborn as sns
    except ImportError as exc:
        raise GCENInputError(
            "seaborn is required for violin plots. "
            "Install with: pip install gcen2[visualization]"
        ) from exc


    df = expr.data if hasattr(expr, "data") else expr
    long_df = df.T.melt(var_name="gene", value_name="expression")
    long_df.rename(columns={"gene": "sample"}, inplace=True)

    fig, ax = plt.subplots(figsize=figsize)
    sns.violinplot(data=long_df, x="sample", y="expression", ax=ax)
    ax.set_xticklabels(ax.get_xticklabels(), rotation=90)
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def heatmap(expr: ExpressionMatrix, modules: ModuleSet, output=None, figsize=(12, 12)):
    """Heatmap (clustermap) of gene expression within modules."""
    _require_matplotlib()
    try:
        import seaborn as sns
    except ImportError as exc:
        raise GCENInputError(
            "seaborn is required for heatmaps. "
            "Install with: pip install gcen2[visualization]"
        ) from exc

    df = expr.data if hasattr(expr, "data") else expr
    module_genes = [g for mod in modules.modules for g in mod if g in df.index]
    fig = sns.clustermap(df.loc[module_genes], figsize=figsize)

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def degree_distribution(network: Network, output=None, figsize=(8, 6)):
    """Histogram of node degree distribution."""
    _require_matplotlib()
    _require_networkx()

    G = network.to_networkx()
    return degree_histogram([d for _, d in G.degree()], output=output, figsize=figsize)


def degree_histogram(degrees, output=None, figsize=(8, 6)):
    """:func:`degree_distribution` from the node degrees themselves (any
    order), e.g. those of :func:`gcen2.network.summary_from_file`."""
    plt = _require_matplotlib()

    fig, ax = plt.subplots(figsize=figsize)
    ax.hist(degrees, bins=50, color="steelblue", alpha=0.8, edgecolor="white")
    ax.set_xlabel("Degree")
    ax.set_ylabel("Count")
    ax.set_title("Degree Distribution")
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def _annotated_heatmap(ax, values, labels, rows, cols, title):
    """RdBu heatmap of correlations (-1..1) with a text label per cell."""
    import numpy as np

    im = ax.imshow(values, cmap="RdBu_r", vmin=-1, vmax=1, aspect="auto")
    ax.set_xticks(range(len(cols)))
    ax.set_xticklabels(cols, rotation=90)
    ax.set_yticks(range(len(rows)))
    ax.set_yticklabels(rows)
    if labels is not None:
        size = max(5, min(9, 180 // max(len(rows), len(cols), 1)))
        for i in range(values.shape[0]):
            for j in range(values.shape[1]):
                v = values[i, j]
                if np.isfinite(v):
                    ax.text(j, i, labels[i][j], ha="center", va="center", fontsize=size,
                            color="white" if abs(v) > 0.6 else "black")
    ax.set_title(title)
    return im


def eigengene_heatmap(eigengenes, output=None, figsize=None):
    """Correlations between module eigengenes (WGCNA
    ``plotEigengeneNetworks``' heatmap): ``eigengenes`` is a samples x
    eigengenes DataFrame (``gcen2 module -e --me``; ``ME_unassigned`` is
    left out). Eigengenes are ordered by average-linkage clustering on
    1 - r; values are printed for up to 25 modules."""
    import numpy as np

    plt = _require_matplotlib()
    me = eigengenes.drop(columns=[c for c in eigengenes.columns if c == "ME_unassigned"])
    if me.shape[1] == 0:
        raise GCENInputError("no module eigengenes to plot.")
    r = np.corrcoef(me.to_numpy(dtype=float), rowvar=False).reshape(me.shape[1], me.shape[1])
    order = np.arange(me.shape[1])
    if me.shape[1] > 2:
        from scipy.cluster.hierarchy import leaves_list, linkage
        from scipy.spatial.distance import squareform

        d = np.clip(1 - r, 0, 2)
        np.fill_diagonal(d, 0)
        order = leaves_list(linkage(squareform(d, checks=False), "average"))
    r = r[np.ix_(order, order)]
    names = [str(me.columns[i]) for i in order]
    n = len(names)
    labels = [[f"{v:.2f}" for v in row] for row in r] if n <= 25 else None
    figsize = figsize or (max(4, 0.45 * n + 2), max(3.5, 0.45 * n + 1.5))
    fig, ax = plt.subplots(figsize=figsize)
    im = _annotated_heatmap(ax, r, labels, names, names, "Eigengene correlation")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def module_trait_heatmap(trait_cor, output=None, figsize=None):
    """Module-trait relationships (WGCNA ``labeledHeatmap`` of the tutorial):
    ``trait_cor`` is the long table of ``gcen2 module -e --trait-cor``
    (module, trait, cor, p); each cell shows r and its p-value."""
    import numpy as np

    plt = _require_matplotlib()
    t = trait_cor[trait_cor["module"] != "ME_unassigned"]
    if t.empty:
        raise GCENInputError("no module-trait correlations to plot.")
    mods = list(dict.fromkeys(t["module"]))
    traits = list(dict.fromkeys(t["trait"]))
    r = t.pivot(index="module", columns="trait", values="cor").reindex(index=mods,
                                                                      columns=traits)
    p = t.pivot(index="module", columns="trait", values="p").reindex(index=mods, columns=traits)
    vals = r.to_numpy(dtype=float)
    labels = [[f"{vals[i, j]:.2f}\n({p.iat[i, j]:.1e})" if np.isfinite(vals[i, j]) else ""
               for j in range(len(traits))] for i in range(len(mods))]
    figsize = figsize or (max(4, 1.1 * len(traits) + 2.5), max(3, 0.5 * len(mods) + 1.5))
    fig, ax = plt.subplots(figsize=figsize)
    im = _annotated_heatmap(ax, vals, labels, [str(m) for m in mods], [str(c) for c in traits],
                            "Module-trait relationships")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    plt.tight_layout()
    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def spatial_eigengenes(coords, eigengenes, output=None, ncols: int = 4, point_size=None,
                       flip_y: bool = True):
    """Spatial maps of spot-level module eigengenes (hdWGCNA
    ``ModuleFeaturePlot``): ``coords`` (spots x 2: x, y, e.g. ``obsm
    ["spatial"]``), ``eigengenes`` (spots x modules DataFrame, e.g. ``gcen2
    spatial-modules`` spot_eigengenes.tsv, rows in the order of
    ``coords``). Each panel: one module, colour scale symmetric around 0
    (clipped at the 1st / 99th percentile of ``|ME|``). ``flip_y``: image
    coordinates (y down)."""
    import numpy as np

    plt = _require_matplotlib()
    xy = np.asarray(coords, dtype=float)
    me = eigengenes.drop(columns=[c for c in eigengenes.columns if c == "ME_unassigned"])
    if me.shape[1] == 0:
        raise GCENInputError("no module eigengenes to plot.")
    if xy.ndim != 2 or xy.shape[1] < 2 or len(xy) != len(me):
        raise GCENInputError("coords must be spots x 2, one row per eigengene row.")
    n = me.shape[1]
    ncols = max(1, min(ncols, n))
    nrows = -(-n // ncols)
    fig, axes = plt.subplots(nrows, ncols, figsize=(3.2 * ncols, 3.2 * nrows), squeeze=False)
    size = point_size or max(1.0, 25000 / max(len(xy), 1))
    for i, col in enumerate(me.columns):
        ax = axes[i // ncols][i % ncols]
        v = me[col].to_numpy(dtype=float)
        lim = float(np.nanpercentile(np.abs(v), 99)) or 1.0
        sc = ax.scatter(xy[:, 0], xy[:, 1], c=v, s=size, cmap="RdBu_r", vmin=-lim, vmax=lim,
                        linewidths=0)
        ax.set_title(str(col))
        ax.set_aspect("equal")
        if flip_y:
            ax.invert_yaxis()
        ax.set_xticks([])
        ax.set_yticks([])
        fig.colorbar(sc, ax=ax, fraction=0.046, pad=0.02)
    for j in range(n, nrows * ncols):
        axes[j // ncols][j % ncols].axis("off")
    plt.tight_layout()
    if output:
        fig.savefig(output, dpi=200, bbox_inches="tight")
    return fig


def diff_network(diff: DiffCoExResult, output=None, figsize=(10, 10)):
    """Visualize differential co-expression network with colored edges."""
    plt = _require_matplotlib()
    nx = _require_networkx()

    G = nx.Graph()
    color_map = {
        "A_only": "red",
        "B_only": "blue",
        "shared_changed": "orange",
        "shared_stable": "gray",
    }
    edge_colors = []
    for _, row in diff.results.iterrows():
        G.add_edge(row["gene_a"], row["gene_b"])
        edge_colors.append(color_map.get(row.get("category", "shared_stable"), "gray"))

    fig, ax = plt.subplots(figsize=figsize)
    pos = nx.spring_layout(G, seed=42)
    nx.draw_networkx(G, pos=pos, ax=ax, edge_color=edge_colors,
                     node_size=30, font_size=5)
    ax.axis("off")
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig


def trajectory_dynamics(traj, top: int = 10, output=None, figsize=(10, 7)):
    """Edges per window (top) and the correlation of the ``top`` edges that
    change most (largest ``|cor_change|``) along pseudotime (bottom).

    ``traj`` is a :class:`gcen2.trajectory.TrajectoryNetworks`.
    """
    plt = _require_matplotlib()
    import numpy as np

    table = traj.window_table()
    centers = table["center"].to_numpy()
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize, sharex=True,
                                   gridspec_kw={"height_ratios": [1, 2]})
    width = (np.min(np.diff(centers)) * 0.8) if len(centers) > 1 else 0.1
    ax1.bar(centers, table["edge_count"], width=width, color="steelblue", alpha=0.8)
    ax1.set_ylabel("Edges")
    ax1.set_title("Co-expression along pseudotime")

    dyn = traj.dynamics
    cor_cols = [f"cor_{w.index}" for w in traj.windows]
    if len(dyn) and top > 0:
        sel = dyn.loc[dyn["cor_change"].abs().sort_values(ascending=False).index[:top]]
        for _, row in sel.iterrows():
            ax2.plot(centers, row[cor_cols].to_numpy(dtype=float), marker="o", ms=3,
                     label=f"{row['gene_a']}–{row['gene_b']} ({row['pattern']})")
        ax2.legend(fontsize=7, loc="center left", bbox_to_anchor=(1.01, 0.5))
    ax2.axhline(0, color="gray", lw=0.5)
    ax2.set_xlabel("Pseudotime (window median)")
    ax2.set_ylabel("Correlation")
    plt.tight_layout()

    if output:
        fig.savefig(output, dpi=300, bbox_inches="tight")
    return fig
