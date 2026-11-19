"""``gcen2 plot``: one command per figure of :mod:`gcen2.plot`.

Every command is declared in :data:`_COMMANDS` (or by :func:`_plot` for the
few that need more than one call): its input files and options, mapped by
name onto the keyword arguments of the plot function, plus the options all
figures share (``-o``, ``--size``, ``--format``, ``--palette``). The plot
functions are imported when a command runs, which keeps ``import gcen2.cli``
free of matplotlib and pandas.
"""

from __future__ import annotations

import contextlib
import inspect
from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import plot_app

_SIZE_HELP = "Figure width: single (89 mm), onehalf (120 mm), double (183 mm) or mm."
_FMT_HELP = "Output format(s): pdf, png, svg, tiff, eps; repeatable [default: pdf, png, svg]."
_KW = inspect.Parameter.KEYWORD_ONLY


def _path(arg: str, *flags: str, help: str, required: bool = True,
          multiple: bool = False) -> inspect.Parameter:
    """An input file option, passed to the plot function as ``arg`` (left
    out when not given, so the function's default applies)."""
    ann = Optional[list[Path]] if multiple else (Path if required else Optional[Path])
    return inspect.Parameter(arg, _KW, annotation=ann,
                             default=typer.Option(... if required else None, *flags, help=help))


def _value(arg: str, annotation, default, *flags: str, help: str) -> inspect.Parameter:
    """Any other option, passed to the plot function as ``arg``."""
    return inspect.Parameter(arg, _KW, annotation=annotation,
                             default=typer.Option(default, *flags, help=help))


def _shared(size: str) -> list[inspect.Parameter]:
    return [
        _value("output", Path, ..., "-o", "--output", help="Output file stem (or file)."),
        _value("size", str, size, "--size", help=_SIZE_HELP),
        _value("formats", Optional[list[str]], None, "--format", help=_FMT_HELP),
        _value("palette", Optional[str], None, "--palette", help="wgcna or okabe-ito."),
    ]


def _draw(fn, output: Path, **kwargs) -> None:
    from gcen2.plot._style import get_defaults

    fn(output=output, **kwargs)
    typer.echo(f"Figure written: {output} ({', '.join(get_defaults()['formats'])})")


def _plot(name: str, target: str, doc: str, *params: inspect.Parameter, size: str = "single",
          run=None) -> None:
    """Register ``gcen2 plot <name>``: ``target`` is ``"<module>.<function>"``
    of :mod:`gcen2.plot`; ``run(fn, output, size, args)`` replaces the plain
    call ``fn(output=..., size=..., **args)`` when given."""
    optional_paths = {p.name for p in params if p.default.default is None
                      and p.annotation in (Optional[Path], Optional[list[Path]])}

    def cmd(**kw):
        import importlib

        from gcen2.plot import set_defaults

        module, fn_name = target.rsplit(".", 1)
        fn = getattr(importlib.import_module(f"gcen2.plot.{module}"), fn_name)
        # every command sets both, so one command's choice never leaks into
        # the next in the same process
        set_defaults(palette=kw.pop("palette") or "wgcna",
                     formats=kw.pop("formats") or ("pdf", "png", "svg"))
        sz = kw.pop("size")
        with contextlib.suppress(ValueError):
            sz = float(sz)
        output = kw.pop("output")
        args = {k: v for k, v in kw.items() if not (k in optional_paths and v is None)}
        if run is None:
            _draw(fn, output, size=sz, **args)
        else:
            run(fn, output, sz, args)

    cmd.__signature__ = inspect.Signature([*params, *_shared(size)])
    cmd.__doc__ = doc
    cmd.__name__ = f"plot_{name.replace('-', '_')}"
    plot_app.command(name)(cmd)


# -- figures that need more than one call -------------------------------------

def _hub_networks(fn, output, size, args) -> None:
    """One module, or every module of the kME table (one file each)."""
    if args["module"] is not None:
        _draw(fn, output, size=size, **args)
        return
    from gcen2.plot._io import kme
    from gcen2.plot.wgcna import _kme_columns

    for m in _kme_columns(kme(args["kme_table"])):
        _draw(fn, output.with_name(f"{output.name}_{m}"), size=size, **{**args, "module": m})


def _by_group(fn, output, size, args) -> None:
    """Eigengenes by group: the groups of the eigengenes' cells first."""
    from gcen2.plot._io import eigengenes
    from gcen2.plot.single_cell import groups_from

    cells = eigengenes(args["eigengenes"]).index
    args["groups"] = groups_from(args["groups"], args.pop("groupby"), cells)
    _draw(fn, output, size=size, **args)


def _preservation(fn, output, size, args) -> None:
    if args.pop("lollipop"):
        from gcen2.plot.single_cell import preservation_lollipop

        args.pop("observed", None)
        _draw(preservation_lollipop, output, size="single" if size == "double" else size,
              **args)
    else:
        _draw(fn, output, size=size, **args)


# -- the commands ----------------------------------------------------------------

_INPUT = ("-i", "--input")
_NET = _path("network", *_INPUT, help="Network edge-list file.")
_TOM = _path("tom_hubs", "-t", "--tom", help="tom_hubs.tsv (--write-tom / --tom-hubs).")
_KME = _path("kme_table", "-k", "--kme", help="kme.tsv.")
_MODULE_LIST = _value("modules", Optional[list[str]], None, "-m", "--module",
                      help="Only this module; repeatable.")
_TOP = "--top"

# WGCNA / hdWGCNA modules
_plot("module-sizes", "wgcna.module_sizes", "Genes per module (bar chart in the module colours).",
      _path("modules", "-m", "--modules", help="modules.tsv (module, gene)."),
      _value("include_unassigned", bool, False, "--unassigned",
             help="Also the unassigned genes."))
_plot("soft-power", "wgcna.soft_power",
      "Scale-free fit and mean connectivity of each soft power (PlotSoftPowers).",
      _path("table", *_INPUT, help="soft_power.tsv."),
      _value("r2_cut", Optional[float], 0.8, "--r2",
             help="R^2 threshold line (hdWGCNA 0.8, WGCNA 0.85)."))
_plot("dendrogram", "wgcna.dendrogram",
      "Gene dendrogram with module colours (PlotDendrogram, plotDendroAndColors).",
      _path("tree", *_INPUT, help="dendrogram.json."),
      _value("show_unmerged", bool, True, "--unmerged/--no-unmerged",
             help="Also the modules before merging."), size="double")
_plot("eigengene-network", "wgcna.eigengene_network",
      "Eigengene dendrogram and correlation heatmap (plotEigengeneNetworks, ModuleCorrelogram).",
      _path("eigengenes", *_INPUT, help="eigengenes.tsv / cell_eigengenes.tsv."),
      _path("traits", "--traits", help="Samples x numeric traits TSV.", required=False))
_plot("module-trait", "wgcna.module_trait",
      "Module-trait correlation heatmap (labeledHeatmap, PlotModuleTraitCorrelation).",
      _path("trait_cor", *_INPUT,
            help="Trait correlations: module -e --trait-cor or trait_cor.tsv."),
      _value("annotate", str, "p", "--annotate", help="p (r and p), stars or none."),
      _value("group", Optional[str], None, "--group",
             help="Only this group (celltype trait_cor)."))
_plot("kme", "wgcna.kme", "Sorted kME of each module's genes (PlotKMEs).",
      _path("kme_table", *_INPUT, help="kme.tsv or cell_kme.tsv."),
      _value("n_label", int, 5, "--n-label", help="Hub genes named per module."), size="double")
_plot("mm-gs", "wgcna.mm_gs", "Module membership vs gene significance (verboseScatterplot).",
      _path("kme_table", *_INPUT, help="kme.tsv."),
      _path("gs_table", "--gs", help="Gene significance (module -e --gs)."),
      _value("module", str, ..., "-m", "--module", help="Module (e.g. 0 or ME0)."),
      _value("trait", str, ..., "--trait", help="Trait name."))
_plot("hub-network", "wgcna.hub_network", "Hub-gene network of a module (ModuleNetworkPlot).",
      _TOM, _KME,
      _value("module", Optional[str], None, "-m", "--module",
             help="One module (e.g. 0 or ME0); default: all modules, one file each "
                  "(<output>_ME<k>)."),
      _value("n_inner", int, 10, "--n-inner", help="Genes on the inner circle."),
      _value("n_outer", int, 15, "--n-outer", help="Genes on the outer circle."),
      run=_hub_networks)
_plot("hub-network-all", "wgcna.hub_network_all",
      "Hub genes of all modules in one network (HubGeneNetworkPlot).",
      _path("tom_hubs", "-t", "--tom", help="tom_hubs.tsv."), _KME,
      _value("n_hubs", int, 6, "--n-hubs", help="Hub genes per module."),
      _value("n_other", int, 3, "--n-other", help="Other genes per module."),
      _value("edge_prop", float, 0.5, "--edge-prop", help="Fraction of the strongest edges drawn."),
      _value("seed", int, 0, "--seed", help="Layout seed."))
_plot("module-umap", "wgcna.module_umap",
      "UMAP of the module genes on their TOM with the hub genes (RunModuleUMAP).",
      _path("tom_hubs", "-t", "--tom", help="tom_hubs.tsv."), _KME,
      _value("n_hubs", int, 10, "--n-hubs", help="Hub genes per module used as UMAP features."),
      _value("n_label", int, 3, "--n-label", help="Hub genes named per module."),
      _value("seed", int, 0, "--seed", help="UMAP seed."))
_plot("tom-heatmap", "wgcna.tom_heatmap", "TOM heatmap of a random sample of genes (TOMplot).",
      _path("expr", "-e", "--expr", help="Genes x samples (metacells.tsv, module -e input)."),
      _path("modules", "-m", "--modules", help="modules.tsv."),
      _value("power", float, ..., "--power", help="Soft power of the network."),
      _value("signed", bool, True, "--signed/--unsigned", help="Signed network (hdWGCNA)."),
      _value("n_genes", int, 400, "--n-genes", help="Random genes shown."),
      _value("seed", int, 0, "--seed", help="Seed of the gene sample."))
_plot("sample-tree", "wgcna.sample_tree", "Sample clustering tree with trait colours (outlier "
                                          "check).",
      _path("expr", "-e", "--expr", help="Genes x samples expression TSV."),
      _path("traits", "--traits", help="Samples x traits TSV (colour bars).", required=False),
      _value("cut_height", Optional[float], None, "--cut-height", help="Draw this cut height."),
      size="double")

# single cell and spatial
_plot("me-embedding", "single_cell.me_embedding",
      "Eigengenes of every cell on an embedding or in space (ModuleFeaturePlot).",
      _path("eigengenes", *_INPUT, help="cell_eigengenes.tsv, cell_harmonized_eigengenes.tsv, "
                                       "spot_eigengenes.tsv or projected_eigengenes.tsv."),
      _path("coords", "-c", "--coords", help=".h5ad (obsm[--basis]) or TSV of x, y by cell."),
      _value("basis", str, "X_umap", "--basis",
             help="obsm key of an .h5ad (X_umap, spatial, ...)."),
      _value("spatial", bool, False, "--spatial", help="Spatial map: equal axes, image y axis."),
      _MODULE_LIST, _value("ncols", int, 4, "--ncols", help="Panels per row."), size="double")
for _name, _fn, _doc, _size in (
        ("me-dotplot", "me_dotplot", "Eigengenes by group: dot plot (DotPlot of MEs).", "single"),
        ("me-violin", "me_violin", "Eigengenes by group: violins.", "double"),
        ("me-radar", "me_radar", "Eigengenes by group: radar chart (ModuleRadarPlot).", "double")):
    _plot(_name, f"single_cell.{_fn}", _doc,
          _path("eigengenes", *_INPUT, help="cell_eigengenes.tsv (or harmonized)."),
          _path("groups", "-g", "--groups", help=".h5ad or TSV (first column cells)."),
          _value("groupby", str, ..., "--groupby", help="Column with the groups."),
          _MODULE_LIST, size=_size, run=_by_group)
_plot("dme-volcano", "single_cell.dme_volcano",
      "Differential module eigengenes: volcano plot (PlotDMEsVolcano).",
      _path("dmes", *_INPUT, help="dmes.tsv."),
      _value("alpha", float, 0.05, "--alpha", help="Significance threshold (adjusted p)."),
      size="double")
_plot("dme-lollipop", "single_cell.dme_lollipop",
      "Differential module eigengenes: lollipop chart (PlotDMEsLollipop).",
      _path("dmes", *_INPUT, help="dmes.tsv."),
      _path("modules", "-m", "--modules", help="modules.tsv (dot size).", required=False),
      _value("alpha", float, 0.05, "--alpha", help="Significance threshold (adjusted p)."),
      size="double")
_plot("celltype-overview", "single_cell.celltype_overview",
      "Cells, metacells, genes and modules of every cell type.",
      _path("table", *_INPUT, help="celltypes.tsv of celltype-modules."), size="double")
_plot("set-eigengenes", "single_cell.set_eigengenes", "Consensus modules: eigengenes in every set.",
      _path("table", *_INPUT, help="set_eigengenes.tsv (consensus network)."), size="double")
_plot("preservation", "single_cell.preservation",
      "Module preservation: Zsummary and median rank against module size "
      "(PlotModulePreservation), or a lollipop chart (PlotModulePreservationLollipop).",
      _path("z_table", *_INPUT, help="preservation_Z.tsv."),
      _path("observed", "--observed", help="preservation.tsv (adds the median rank panel).",
            required=False),
      _value("lollipop", bool, False, "--lollipop", help="Sorted Zsummary lollipop chart instead."),
      size="double", run=_preservation)
_plot("metacells", "single_cell.metacells",
      "Metacells on the cell embedding (DimPlotMetacells), per cell type.",
      _path("members", *_INPUT, help="metacell_members.tsv, or a celltype-modules output "
                                    "directory."),
      _path("coords", "-c", "--coords", help="AnnData (.h5ad) or TSV of x, y by cell."),
      _value("basis", str, "X_umap", "--basis", help="obsm key of the embedding."), size="double")
_plot("moran", "single_cell.moran",
      "Moran's I of the genes (spatial-build --moran): distribution, top genes, smoothing.",
      _path("table", *_INPUT, help="spatial-build --moran output."), size="double")

# networks
_plot("degree-distribution", "network.degree_distribution",
      "Degree histogram and scale-free plot of a network.", _NET)
_plot("edge-weights", "network.edge_weights", "Distribution of the edge weights.", _NET,
      _value("column", str, "correlation", "--column", help="Weight column."))
_plot("network-view", "network.network_view",
      "Drawing of a network (best-connected genes, or --genes), coloured by module.", _NET,
      _path("modules", "-m", "--modules", help="modules.tsv / module.txt.", required=False),
      _value("max_nodes", int, 150, "--max-nodes", help="Genes drawn."),
      _value("n_label", int, 12, "--n-label", help="Genes named."))
_plot("module-size-hist", "network.module_size_hist",
      "Histogram of module sizes (many modules, e.g. gcen2 module -i).",
      _path("modules", "-m", "--modules", help="Module file (module, gene)."))
_plot("celltype-networks", "network.celltype_networks",
      "Edges and edge overlap of gcen2 celltype networks.",
      _path("directory", *_INPUT, help="Output directory of gcen2 celltype."), size="double")
_plot("null-distribution", "network.null_distribution",
      "A network statistic against degree-preserving shuffles (gcen2 shuffle).", _NET,
      _path("shuffled", "-s", "--shuffled", help="Shuffled network; repeatable.",
            required=False, multiple=True),
      _value("metric", str, "clustering_coefficient", "--metric",
             help="clustering_coefficient, transitivity, ..."))
_plot("merge-overlap", "network.merge_overlap",
      "Edges shared by the inputs of gcen2 merge (Venn for 2-3 networks, UpSet for more).",
      _path("networks", *_INPUT, help="Input network of gcen2 merge; repeat for each.",
            multiple=True),
      _path("merged", "-m", "--merged", help="gcen2 merge output (edges kept per region).",
            required=False),
      _value("names", Optional[list[str]], None, "--name",
             help="Network names, one per -i (default: file names)."))
_plot("extract-overview", "network.extract_overview",
      "Query genes of gcen2 extract: found in the network, and their degree.",
      _path("network", *_INPUT, help="Network, or the gcen2 extract output."),
      _path("genes", "-g", "--genes", help="Gene list given to gcen2 extract."),
      _value("top", int, 20, _TOP, help="Genes shown."))

# expression, enrichment, function prediction, other analyses
_plot("expression-qc", "qc.expression_distribution",
      "Per-sample expression distributions, before and after normalization.",
      _path("expr", "-e", "--expr", help="Expression TSV (genes x samples)."),
      _path("normalized", "--normalized", help="Normalized expression (gcen2 normalize).",
            required=False), size="double")
_plot("mean-variance", "qc.mean_variance",
      "Mean against coefficient of variation; genes kept by gcen2 filter in colour.",
      _path("expr", "-e", "--expr", help="Expression TSV."),
      _path("kept", "--kept", help="Filtered expression (gcen2 filter output).", required=False))
_plot("sample-correlation", "qc.sample_correlation",
      "Sample-sample correlation heatmap (clustered).",
      _path("expr", "-e", "--expr", help="Expression TSV."))
_plot("enrich-bar", "enrichment.enrich_bar", "Top enriched terms of one result (bars).",
      _path("table", *_INPUT, help="gcen2 enrich output / one annotate file."),
      _value("top", int, 15, _TOP, help="Terms shown."))
_plot("enrich-dot", "enrichment.enrich_dot",
      "Top terms of every module (dot plot; EnrichrDotPlot).",
      _path("results", *_INPUT, help="gcen2 annotate output directory."),
      _value("top", int, 3, _TOP, help="Terms per module."), size="double")
_plot("enrich-heatmap", "enrichment.enrich_heatmap", "Modules x top terms heatmap of -log10 p.",
      _path("results", *_INPUT, help="gcen2 annotate output directory."),
      _value("top", int, 3, _TOP, help="Terms per module."), size="double")
_plot("rwr-rank", "function.rwr_rank", "Top genes by random-walk score (seeds marked).",
      _path("table", *_INPUT, help="gcen2 rwr / rwr-multiplex output."),
      _path("seeds", "-g", "--genes", help="Seed gene list.", required=False),
      _value("top", int, 30, _TOP, help="Genes shown."))
_plot("rwr-network", "function.rwr_subnetwork", "The seeds and the top candidates in the network.",
      _path("network", "-n", "--network", help="Network edge-list file."),
      _path("table", *_INPUT, help="gcen2 rwr output."),
      _path("seeds", "-g", "--genes", help="Seed gene list."),
      _value("top", int, 20, _TOP, help="Candidates shown."))
_plot("diff-overview", "diff.diff_overview",
      "Differential co-expression at a glance (categories, A vs B, change vs p).",
      _path("diff", *_INPUT, help="gcen2 diff output."), size="double")
_plot("diff-network", "diff.diff_network", "Network of the most changed edges.",
      _path("diff", *_INPUT, help="gcen2 diff output."),
      _value("top", int, 60, _TOP, help="Edges drawn."))
_plot("trajectory", "trajectory.trajectory_overview",
      "Networks along pseudotime (windows, edge patterns, dynamic genes).",
      _path("directory", *_INPUT, help="gcen2 trajectory output directory."), size="double")
_plot("fusion", "multiomics.fusion_overview", "A fused multi-layer network (gcen2 fuse).",
      _path("fused", *_INPUT, help="gcen2 fuse output."), size="double")
