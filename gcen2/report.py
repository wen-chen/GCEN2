"""Automated HTML/PDF research report generation."""

from __future__ import annotations

import io
from pathlib import Path
from typing import TYPE_CHECKING, Optional, Union

from gcen2._matrix_io import require_file
from gcen2.exceptions import GCENInputError, GCENIOError
from gcen2.models import EnrichmentResult, ExpressionMatrix, ModuleSet, Network, RWRResult

if TYPE_CHECKING:
    from gcen2.network import NetworkSummary

_DEFAULT_TEMPLATE = Path(__file__).parent / "templates" / "report_default.html.jinja2"


def generate(
    network: Optional[Union[Network, NetworkSummary]] = None,
    modules: Optional[ModuleSet] = None,
    enrichment: Optional[dict[str, EnrichmentResult]] = None,
    rwr: Optional[RWRResult] = None,
    expr: Optional[ExpressionMatrix] = None,
    output: Union[str, Path] = "report.html",
    template: Optional[Union[str, Path]] = None,
    wgcna: Optional[Union[dict, list]] = None,
    figures_dir: Optional[Union[str, Path]] = None,
    celltypes: Optional[dict] = None,
    moran=None,
    merge: Optional[dict] = None,
    extract: Optional[dict] = None,
) -> None:
    """Generate a comprehensive HTML research report.

    Parameters
    ----------
    network:
        Co-expression network, or its :class:`gcen2.network.NetworkSummary`
        (:func:`gcen2.network.summary_from_file`), which gives the same report
        without the edge table.
    modules:
        Identified gene modules.
    enrichment:
        Enrichment results dict (gene/module → EnrichmentResult).
    rwr:
        RWR prediction results.
    expr:
        Expression matrix (for QC section).
    output:
        Output file path. Use .pdf extension for PDF output (requires weasyprint).
    template:
        Custom Jinja2 template path. Uses the built-in default if None.
    wgcna:
        WGCNA-style module results: a section from :func:`wgcna_section`
        (``gcen2 module -e`` outputs) or a list of them
        (:func:`celltype_module_sections`), each shown with its module table
        and hub genes and the figures its files allow (gcen2.plot: module
        sizes, soft power, dendrogram, eigengene network, module-trait, kME,
        hub network, differential eigengenes, consensus sets).
    figures_dir:
        Also save every report figure there as PDF, PNG and SVG (named by
        section and figure), ready for a manuscript.
    celltypes:
        :func:`celltype_overview` of a ``celltype-modules`` directory: cell
        types at a glance and, with an embedding, the metacells on it.
    moran:
        Moran's I table of ``spatial-build --moran`` (path or DataFrame).
    merge:
        ``{"networks": [...], "merged": path or None, "names": [...] or
        None}``: the inputs (and output) of ``gcen2 merge``.
    extract:
        ``{"network": path, "genes": path or list}``: the network (or the
        extracted subnetwork) and query genes of ``gcen2 extract``.

    Notes
    -----
    Figures are embedded as vector SVG (so a PDF report keeps them vector),
    each with a caption.
    """
    try:
        from jinja2 import Environment, FileSystemLoader, select_autoescape
    except ImportError as exc:
        raise GCENInputError(
            "Jinja2 is required for report generation. "
            "Install with: pip install gcen2[report]"
        ) from exc

    output = Path(output)
    template_path = Path(template) if template else _DEFAULT_TEMPLATE
    if not template_path.is_file():
        raise GCENIOError(f"Report template not found: {template_path}")

    figs = _Figures(figures_dir)
    context = _build_context(network, modules, enrichment, rwr, expr, figs)
    context["wgcna"] = _wgcna_context(wgcna, figs)
    _extra_figures(context, figs, celltypes, moran, merge, extract)

    env = Environment(
        loader=FileSystemLoader(str(template_path.parent)),
        autoescape=select_autoescape(["html"]),
    )
    tmpl = env.get_template(template_path.name)
    html_content = tmpl.render(**context)

    if output.suffix.lower() == ".pdf":
        _render_pdf(html_content, output)
    else:
        try:
            with open(output, "w", encoding="utf-8", newline="\n") as fh:
                fh.write(html_content)
        except Exception as exc:
            raise GCENIOError(f"Failed to write report to {output}: {exc}") from exc


def _build_context(network, modules, enrichment, rwr, expr, figs=None) -> dict:
    import gcen2

    figs = figs or _Figures(None)
    context: dict = {
        "gcen2_version": gcen2.__version__,
        "network": None,
        "network_stats": None,
        "modules": None,
        "enrichment": None,
        "rwr": None,
        "expr_qc": None,
        "figures": {},
    }

    if network is not None:
        from gcen2.network import NetworkSummary

        summary = isinstance(network, NetworkSummary)
        context["network"] = network
        context["network_stats"] = network.stats if summary else network.stats()
        figs.add(context["figures"], "degree_distribution", "network_degree",
                 lambda p: p.degree_distribution(network.degrees if summary else network),
                 "Degree distribution of the co-expression network. Left: genes by degree. "
                 "Right: WGCNA's scale-free plot, log10 frequency against log10 degree in ten "
                 "bins, with the fitted line, its R² and slope (R² near 1 = scale-free).")

    if modules is not None:
        context["modules"] = modules
        figs.add(context["figures"], "module_sizes", "network_module_sizes",
                 lambda p: p.module_size_hist(_module_series(modules)),
                 "Sizes of the network modules (Jaccard clustering), log scale.")

    if enrichment is not None:
        context["enrichment"] = enrichment
        first_key = next(iter(enrichment))
        res = enrichment[first_key]
        table = getattr(res, "results", res)
        figs.add(context["figures"], "dotplot", "enrichment",
                 lambda p: p.enrich_bar(table),
                 f"Top enriched terms of {first_key} ({getattr(res, 'ontology_type', '')}): "
                 "bars of −log10 p, coloured by GO namespace, with the number of genes.")

    if rwr is not None:
        context["rwr"] = rwr
        figs.add(context["figures"], "rwr", "rwr_rank",
                 lambda p: p.rwr_rank(rwr.ranked_genes, seeds=rwr.seed_genes),
                 "Genes ranked by random walk with restart from the seed genes (diamonds).")

    if expr is not None:
        context["expr_qc"] = {
            "gene_count": len(expr.gene_names),
            "sample_count": len(expr.sample_names),
            "is_normalized": expr.is_normalized,
            "normalization_method": expr.normalization_method,
        }
        figs.add(context["figures"], "expr_qc", "expression_qc",
                 lambda p: p.expression_distribution(expr),
                 "Expression of every sample, log10(x + 1), as boxes (median, quartiles, "
                 "whiskers at 1.5 IQR).")

    return context


def celltype_overview(directory, coords=None, basis: str = "X_umap") -> dict:
    """The cell-type overview of a ``gcen2 celltype-modules`` directory for
    :func:`generate` (``celltypes``): cells, metacells, genes and modules per
    cell type and, with ``coords`` (an AnnData / .h5ad with ``obsm[basis]``,
    or a TSV of x, y by cell), the metacells on that embedding."""
    d = Path(directory)
    require_file(d / "celltypes.tsv")
    return {"directory": d, "coords": coords, "basis": basis}


def _extra_figures(context, figs, celltypes, moran, merge, extract) -> None:
    f = context["figures"]
    context["celltype_figures"] = []
    if celltypes is not None:
        d = Path(celltypes["directory"])
        for key, draw, cap in (
                ("celltype_overview", lambda p: p.celltype_overview(d / "celltypes.tsv"),
                 "Cell types at a glance: cells, metacells, genes, modules and genes in "
                 "modules of each cell type."),
                ("metacells", lambda p: p.metacells(d, celltypes["coords"],
                                                    basis=celltypes.get("basis", "X_umap")),
                 "Metacells on the cell embedding: every cell in grey, each metacell at the "
                 "mean position of its cells, coloured by cell type; right, metacells per cell "
                 "type and the share of its cells in at least one metacell.")):
            if key == "metacells" and celltypes.get("coords") is None:
                continue
            fig = figs.make(f"celltypes_{key}", draw, cap)
            if fig is not None:
                context["celltype_figures"].append(fig)
    if merge is not None:
        figs.add(f, "merge_overlap", "merge_overlap",
                 lambda p: p.merge_overlap(merge["networks"], merged=merge.get("merged"),
                                           names=merge.get("names")),
                 "Edges shared by the merged networks (Venn diagram, areas not to scale; "
                 "UpSet bars for more than three networks)"
                 + ("; in brackets, the edges of each region that pass the cutoff on the "
                    "summed correlation and are in the merged network." if merge.get("merged")
                    is not None else "."))
    if extract is not None:
        figs.add(f, "extract_overview", "extract_overview",
                 lambda p: p.extract_overview(extract["network"], extract["genes"]),
                 "Query genes of the extracted subnetwork: left, query genes found among the "
                 "network's genes; right, the best-connected query genes, with their edges to "
                 "other query genes and to the rest of the network.")
    if moran is not None:
        figs.add(f, "moran", "moran",
                 lambda p: p.moran(moran),
                 "Spatial autocorrelation of the genes: distribution of Moran's I (the filtered "
                 "spatial network keeps genes with I > 0), I against mean expression with the "
                 "most autocorrelated genes named, and, for the weighted network, I before and "
                 "after the spatial smoothing.")


def _module_series(modules):
    import pandas as pd

    rows = [(g, f"ME{k}") for k, mod in enumerate(modules.modules) for g in mod]
    return pd.Series(dict(rows))


class _Figures:
    """Report figures: drawn with gcen2.plot, embedded as SVG, optionally also
    written to a directory (PDF, PNG, SVG)."""

    def __init__(self, directory):
        self.dir = Path(directory) if directory is not None else None
        if self.dir is not None:
            self.dir.mkdir(parents=True, exist_ok=True)

    def make(self, name: str, draw, caption: str) -> Optional[dict]:
        import re
        import warnings

        try:
            import matplotlib  # noqa: F401
        except ImportError:
            return None
        from gcen2 import plot

        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                fig = draw(plot)
        except (GCENInputError, GCENIOError, ValueError, KeyError) as exc:
            warnings.warn(f"report figure {name} skipped: {exc}", stacklevel=3)
            return None
        svg = _fig_to_svg(fig)
        if self.dir is not None:
            safe = re.sub(r"[^A-Za-z0-9._-]+", "_", name).strip("_")
            plot.save_figure(fig, self.dir / safe, formats=("pdf", "png", "svg"))
        import matplotlib.pyplot as plt

        plt.close(fig)
        return {"svg": svg, "caption": caption, "name": name}

    def add(self, target: dict, key: str, name: str, draw, caption: str) -> None:
        f = self.make(name, draw, caption)
        if f is not None:
            target[key] = f


def _read_tsv(x, **kw):
    import pandas as pd

    if x is None or isinstance(x, pd.DataFrame):
        return x
    return pd.read_csv(require_file(x), sep="\t", **kw)


def wgcna_section(modules, eigengenes=None, kme=None, hubs=None, trait_cor=None,
                  title: str = "WGCNA-style modules", n_hubs: int = 10,
                  note: Optional[str] = None, soft_power=None, dendrogram=None,
                  tom_hubs=None, dmes=None, cell_trait_cor=None, set_eigengenes=None) -> dict:
    """A report section of ``gcen2 module -e`` results (paths or DataFrames):
    ``modules`` (module, gene), and optionally ``eigengenes`` (``--me``),
    ``kme`` (``--kme``), ``hubs`` (``--hubs``; else the top ``n_hubs`` genes
    by kME when ``kme`` is given) and ``trait_cor`` (``--trait-cor``). The
    files for its figures: ``soft_power`` (soft_power.tsv), ``dendrogram``
    (dendrogram.json), ``tom_hubs`` (tom_hubs.tsv), and of celltype-modules
    ``dmes`` (dmes.tsv), ``cell_trait_cor`` (trait_cor.tsv) and
    ``set_eigengenes`` (set_eigengenes.tsv)."""
    mods = _read_tsv(modules)
    if mods is None or not {"module", "gene"} <= set(mods.columns):
        raise GCENInputError("modules must have the columns module and gene.")
    me = _read_tsv(eigengenes, index_col=0)
    kme_df = _read_tsv(kme, index_col=0, dtype={"module": str}, keep_default_na=False)
    hub_df = _read_tsv(hubs)
    tc = _read_tsv(trait_cor)
    hub_lists: dict = {}
    if hub_df is not None:
        for m, grp in hub_df.sort_values(["module", "rank"]).groupby("module"):
            hub_lists[str(m)] = list(grp["gene"].astype(str))[:n_hubs]
    elif kme_df is not None:
        member = kme_df["module"].astype(str)
        for m in member.unique():
            col = f"kME{m}"
            if m == "NA" or col not in kme_df.columns:
                continue
            v = kme_df.loc[member == m, col].astype(float).sort_values(ascending=False)
            hub_lists[str(m)] = [str(g) for g in v.index[:n_hubs]]
    rows = []
    for m, grp in mods.groupby("module", sort=True):
        genes = sorted(grp["gene"].astype(str))
        rows.append({"module": m, "size": len(genes),
                     "hubs": hub_lists.get(str(m)), "genes": genes[:10],
                     "more": len(genes) > 10})
    return {"title": title, "note": note, "rows": rows, "n_modules": len(rows),
            "n_genes": int(len(mods)), "eigengenes": me, "trait_cor": tc,
            "files": {"modules": mods, "kme": kme_df, "soft_power": soft_power,
                      "dendrogram": dendrogram, "tom_hubs": tom_hubs, "dmes": dmes,
                      "cell_trait_cor": cell_trait_cor, "set_eigengenes": set_eigengenes}}


def celltype_module_sections(directory, n_hubs: int = 10) -> list:
    """One :func:`wgcna_section` per cell type of a ``gcen2
    celltype-modules`` output directory (its ``celltypes.tsv``)."""
    d = Path(directory)
    table = _read_tsv(d / "celltypes.tsv")
    out = []
    for row in table.itertuples():
        sub = d / str(row.dir)

        def f(name, sub=sub):
            return sub / name if (sub / name).is_file() else None

        hubs = sub / "hub_genes.tsv"
        power = row.power if isinstance(row.power, str) else f"{row.power:g}"
        note = (f"{row.n_cells} cells, {row.n_metacells} metacells, {row.n_genes} genes, "
                f"soft power {power}")
        out.append(wgcna_section(sub / "modules.tsv", eigengenes=f("eigengenes.tsv"),
                                 kme=f("kme.tsv"), hubs=hubs if hubs.is_file() else None,
                                 title=str(row.cell_type), n_hubs=n_hubs, note=note,
                                 soft_power=f("soft_power.tsv"),
                                 dendrogram=f("dendrogram.json"), tom_hubs=f("tom_hubs.tsv"),
                                 dmes=f("dmes.tsv"), cell_trait_cor=f("trait_cor.tsv"),
                                 set_eigengenes=f("set_eigengenes.tsv")))
    return out


_CAPTIONS = {
    "sizes": "Genes in each module, in the module colours (ME<k> is module k of modules.tsv; "
             "the colour names follow WGCNA's labels2colors).",
    "soft_power": "Choice of the soft-thresholding power: scale-free fit (signed R²) and mean "
                  "connectivity for each power tried; the circle marks the power used, the "
                  "dashed line the R² threshold.",
    "dendrogram": "Gene dendrogram (average linkage of 1 − TOM) with the module of every gene "
                  "below it, before and after merging modules with similar eigengenes; grey: "
                  "unassigned genes.",
    "eigengene_network": "Eigengene network: the module eigengenes clustered by 1 − correlation "
                         "(top) and their correlations (heatmap).",
    "module_trait": "Correlation of each module eigengene with each trait; cells give r and the "
                    "p-value.",
    "kme": "Module membership (kME, correlation of each gene with its module's eigengene), "
           "genes sorted; the hub genes are named.",
    "hub_network": "Hub genes of all modules: the top-kME genes of each module (large, named) and "
                   "a few other module genes, with their strongest TOM connections.",
    "dmes": "Differential module eigengenes: average log2 fold change and −log10 adjusted "
            "p-value of each module in each group against the other cells.",
    "cell_trait": "Correlation of the cell eigengenes with cell traits (r and p-value), for all "
                  "cells and per group.",
    "sets": "Consensus network: each module's eigengene in the metacells of every set "
            "(boxes: median and quartiles).",
}


def _section_figures(sec, figs) -> list:
    files = sec.get("files") or {}
    title = sec.get("title", "wgcna")
    out = []

    def add(key, draw):
        f = figs.make(f"{title}_{key}", draw, _CAPTIONS[key])
        if f is not None:
            out.append(f)

    mods = files.get("modules")
    if mods is not None:
        add("sizes", lambda p: p.module_sizes(mods.astype({"module": str})
                                              .assign(module=lambda t: "ME" + t["module"])
                                              .set_index("gene")["module"]))
    if files.get("soft_power") is not None:
        add("soft_power", lambda p: p.soft_power(files["soft_power"]))
    if files.get("dendrogram") is not None:
        add("dendrogram", lambda p: p.dendrogram(files["dendrogram"]))
    me = sec.get("eigengenes")
    if me is not None and sum(c != "ME_unassigned" for c in me.columns) >= 2:
        add("eigengene_network", lambda p: p.eigengene_network(me))
    tc = sec.get("trait_cor")
    if tc is not None and len(tc):
        add("module_trait", lambda p: p.module_trait(tc))
    kme = files.get("kme")
    if kme is not None:
        add("kme", lambda p: p.kme(kme.assign(module=kme["module"].astype(str))))
    if files.get("tom_hubs") is not None and kme is not None:
        add("hub_network", lambda p: p.hub_network_all(files["tom_hubs"], kme))
    if files.get("dmes") is not None:
        add("dmes", lambda p: p.dme_volcano(files["dmes"]))
    if files.get("cell_trait_cor") is not None:
        add("cell_trait", lambda p: p.module_trait(files["cell_trait_cor"], annotate="stars"))
    if files.get("set_eigengenes") is not None:
        add("sets", lambda p: p.set_eigengenes(files["set_eigengenes"]))
    return out


def _wgcna_context(wgcna, figs=None) -> Optional[list]:
    if wgcna is None:
        return None
    figs = figs or _Figures(None)
    sections = [wgcna] if isinstance(wgcna, dict) else list(wgcna)
    out = []
    for sec in sections:
        sec = dict(sec)
        sec["figures"] = _section_figures(sec, figs)
        tc = sec.get("trait_cor")
        sec["trait_top"] = (None if tc is None else
                            tc[tc["module"] != "ME_unassigned"].sort_values("p").head(10)
                            .to_dict("records"))
        out.append(sec)
    return out


def _fig_to_svg(fig) -> str:
    """The figure as an inline SVG element (text kept as text)."""
    import matplotlib as mpl

    from gcen2.plot._style import rc_params

    buf = io.StringIO()
    with mpl.rc_context(rc_params()):
        fig.savefig(buf, format="svg", bbox_inches="tight", pad_inches=0.02,
                    metadata={"Date": None})
    text = buf.getvalue()
    return text[text.find("<svg"):]


def _render_pdf(html: str, output: Path) -> None:
    try:
        import weasyprint
    except ImportError as exc:
        raise GCENInputError(
            "weasyprint is required for PDF export. "
            "Install with: pip install weasyprint"
        ) from exc
    except OSError as exc:  # installed, but its system libraries (Pango) are missing
        raise GCENInputError(
            f"weasyprint cannot load its system libraries (Pango): {exc}. See "
            "https://doc.courtbouillon.org/weasyprint/stable/first_steps.html"
        ) from exc
    try:
        weasyprint.HTML(string=html).write_pdf(str(output))
    except Exception as exc:
        raise GCENIOError(f"Failed to render PDF {output}: {exc}") from exc
