"""Single-cell commands: celltype (one network per cell type) and
celltype-modules (the hdWGCNA workflow per cell type)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import THREADS_HELP, app, fail, fmt_power, read_h5ad


@app.command("celltype")
def cmd_celltype(
    input: Path = typer.Option(..., "-i", "--input", help="AnnData (.h5ad) file."),
    groupby: str = typer.Option(..., "--groupby",
                                help="adata.obs column with the cell types."),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    groups: Optional[list[str]] = typer.Option(
        None, "--group", help="Only this cell type; repeatable (default: all)."),
    strategy: str = typer.Option("metacell", "--strategy",
                                 help="Samples per cell type: metacell, or pseudo_replicates (deprecated: "
                                      "random groups rank gene pairs at chance level)."),
    k: int = typer.Option(25, "-k", "--k", help="Cells per metacell."),
    max_shared: Optional[int] = typer.Option(
        None, "--max-shared", help="Max cells shared by two metacells (default min(10, k-1))."),
    max_metacells: int = typer.Option(1000, "--max-metacells",
                                      help="Max metacells per cell type."),
    min_cells: int = typer.Option(100, "--min-cells",
                                  help="Skip cell types with fewer cells."),
    use_rep: Optional[str] = typer.Option(
        None, "--use-rep", help="adata.obsm key for the neighbours (default X_pca, else PCA)."),
    n_pcs: int = typer.Option(30, "--n-pcs", help="Components used for the neighbours."),
    replicates: int = typer.Option(10, "--replicates",
                                   help="pseudo_replicates: random groups per cell type."),
    layer: Optional[str] = typer.Option(None, "--layer", help="adata.layers entry instead of X."),
    normalize: Optional[bool] = typer.Option(
        None, "--normalize/--no-normalize",
        help="Log-normalise metacells (default: when the data are integer counts)."),
    hvg: bool = typer.Option(False, "--hvg",
                             help="Only genes flagged in adata.var['highly_variable']."),
    method: str = typer.Option("spearman", "-m", "--method",
                               help="pearson, spearman, pearson_nonzero or "
                                    "pearson_pairwise_complete."),
    cor: float = typer.Option(0.5, "-c", "--cor", help="Minimum correlation cutoff."),
    pval: float = typer.Option(0.001, "-p", "--pval", help="P-value cutoff."),
    fdr: bool = typer.Option(True, "--fdr/--no-fdr", "-f",
                             help="Benjamini-Hochberg FDR over all gene pairs, -p applying to it "
                                  "(default on; off before 2026-11-18)."),
    signed: bool = typer.Option(False, "--signed", help="Keep positive correlations only."),
    log: Optional[str] = typer.Option(None, "-l", "--log",
                                      help="log(x+1) transform first: log, log2 or log10."),
    diff: bool = typer.Option(False, "--diff",
                              help="Also compare every pair of cell-type networks (Fisher's z on "
                                   "the networks' effective sample sizes)."),
    diff_pval: float = typer.Option(0.05, "--diff-pval",
                                    help="P-value cutoff of the differential test."),
    seed: int = typer.Option(42, "--seed", help="Random seed."),
    effective: str = typer.Option(
        "overlap", "--effective-samples",
        help="Metacells: sample size of the p-values. overlap: per gene pair, from the cells "
             "the metacells share; cells: n_cells // k for every pair (before 2026-11-18)."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
) -> None:
    """Cell-type-specific co-expression networks from single-cell data."""
    from gcen2.exceptions import GCENInputError
    from gcen2.single_cell import build_by_celltype, write_celltype_networks

    adata = read_h5ad(input)
    if groupby not in adata.obs.columns:
        fail(f"'{groupby}' is not a column of adata.obs.")
    if groups:
        labels = adata.obs[groupby].astype(str)
        missing = sorted(set(groups) - set(labels))
        if missing:
            fail(f"no cells with {groupby} = {', '.join(missing)}.")
        adata = adata[labels.isin(groups).to_numpy()]
    networks = build_by_celltype(
        adata, groupby, use_highly_variable=hvg, min_cells=min_cells, strategy=strategy, k=k,
        max_shared=max_shared, max_metacells=max_metacells, use_rep=use_rep, n_pcs=n_pcs,
        layer=layer, normalize=normalize, pseudo_replicates=replicates, seed=seed,
        effective=effective,
        method=method, cor=cor, pval=pval, fdr=fdr, signed=signed, log=log, threads=threads)
    if not networks:
        raise GCENInputError("no cell type had enough cells (see --min-cells).")
    table = write_celltype_networks(networks, output, diff=diff, pval=diff_pval, threads=threads)
    typer.echo(f"{len(table)} cell-type networks ({strategy}) → {output}/")
    for row in table.itertuples():
        typer.echo(f"  {row.cell_type}: {row.n_cells} cells, {row.n_samples} samples, "
                   f"{row.node_count} nodes, {row.edge_count} edges")
    if diff:
        import pandas as pd

        summary = pd.read_csv(output / "diff" / "summary.tsv", sep="\t")
        typer.echo(f"Differential co-expression ({len(summary)} pairs) → {output}/diff/")
        for row in summary.itertuples():
            typer.echo(f"  {row.cell_type_a} vs {row.cell_type_b}: "
                       f"{row.shared_changed} changed, {row.A_only} / {row.B_only} specific")


@app.command("celltype-modules")
def cmd_celltype_modules(
    input: Path = typer.Option(..., "-i", "--input", help="AnnData (.h5ad) file."),
    groupby: str = typer.Option(..., "--groupby",
                                help="adata.obs column with the cell types."),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    groups: Optional[list[str]] = typer.Option(
        None, "--group", help="Only this cell type; repeatable (default: all)."),
    fraction: float = typer.Option(
        0.05, "--fraction", help="Genes with a value > 0 in at least this fraction of the "
                                 "type's cells (hdWGCNA gene_select = 'fraction')."),
    hvg: bool = typer.Option(False, "--hvg",
                             help="Only genes flagged in adata.var['highly_variable']."),
    k: int = typer.Option(25, "-k", "--k", help="Cells per metacell."),
    max_shared: Optional[int] = typer.Option(
        None, "--max-shared", help="Max cells shared by two metacells (default min(10, k-1))."),
    max_metacells: int = typer.Option(1000, "--max-metacells",
                                      help="Max metacells per cell type."),
    min_cells: int = typer.Option(100, "--min-cells",
                                  help="Skip cell types with fewer cells."),
    use_rep: Optional[str] = typer.Option(
        None, "--use-rep", help="adata.obsm key for the neighbours (default X_pca, else PCA)."),
    n_pcs: int = typer.Option(30, "--n-pcs", help="Components used for the neighbours."),
    layer: Optional[str] = typer.Option(None, "--layer", help="adata.layers entry instead of X."),
    normalize: Optional[bool] = typer.Option(
        None, "--normalize/--no-normalize",
        help="Log-normalise metacells (default: when the data are integer counts)."),
    preset: str = typer.Option(
        "hdwgcna", "--preset",
        help="hdwgcna (hdWGCNA ConstructNetwork defaults, bicor power) or wgcna "
             "(blockwiseModules defaults)."),
    power: str = typer.Option("auto", "--power",
                              help="Soft-threshold power, or 'auto' (the preset's rule)."),
    min_size: Optional[int] = typer.Option(None, "--min-size",
                                           help="Minimum module size (preset default)."),
    deep_split: Optional[int] = typer.Option(None, "--deep-split",
                                             help="Dynamic tree cut deepSplit, 0-4."),
    merge_cut: Optional[float] = typer.Option(
        None, "--merge-cut", help="Merge modules whose eigengenes have 1 - r below this."),
    n_hubs: int = typer.Option(10, "--n-hubs", help="Hub genes per module."),
    write_metacells: bool = typer.Option(
        False, "--write-metacells",
        help="Also write each type's metacells.tsv (genes x metacells; gcen2 module -e input)."),
    write_tom: bool = typer.Option(
        False, "--write-tom",
        help="Also write each type's tom_hubs.tsv (TOM of the module genes with each module's "
             "top-kME genes; for the hub-network and module-UMAP figures)."),
    seed: int = typer.Option(42, "--seed", help="Random seed (metacell centres)."),
    cell_level: bool = typer.Option(
        True, "--cell-eigengenes/--no-cell-eigengenes",
        help="Eigengenes and kME of every cell (hdWGCNA ModuleEigengenes / "
             "ModuleConnectivity): cell_eigengenes.tsv, cell_kme.tsv."),
    reassign: bool = typer.Option(
        True, "--reassign/--no-reassign",
        help="Move genes with a negative cell kME for their module (hdWGCNA "
             "ModuleConnectivity); needs the cell eigengenes."),
    backed: bool = typer.Option(
        False, "--backed",
        help="Leave the expression on disk and read one cell type's cells at a time "
             "(less memory for an atlas; X must be stored as CSR)."),
    consensus_by: Optional[str] = typer.Option(
        None, "--consensus-by",
        help="adata.obs column (sample, study, ...): a consensus network over its values "
             "(hdWGCNA ConstructNetwork(consensus = TRUE)); metacells within each value."),
    consensus_quantile: Optional[float] = typer.Option(
        None, "--consensus-quantile",
        help="With --consensus-by: quantile of the sets' TOMs (0 = minimum; preset default: "
             "hdwgcna 0.3, wgcna 0)."),
    harmonize_by: Optional[str] = typer.Option(
        None, "--harmonize-by",
        help="adata.obs column (sample, batch, ...): also harmonized cell eigengenes "
             "(hdWGCNA ModuleEigengenes(group.by.vars = ...), Harmony): "
             "cell_harmonized_eigengenes.tsv; the cell kME then uses them."),
    npcs: int = typer.Option(50, "--npcs",
                             help="With --harmonize-by: principal components per module."),
    dme_by: Optional[str] = typer.Option(
        None, "--dme-by",
        help="adata.obs column: differential module eigengenes of each value against the "
             "type's other cells (hdWGCNA FindAllDMEs): dmes.tsv."),
    traits: Optional[list[str]] = typer.Option(
        None, "--trait",
        help="adata.obs column to correlate with the cell eigengenes (hdWGCNA "
             "ModuleTraitCorrelation); repeatable: trait_cor.tsv."),
    trait_groupby: Optional[str] = typer.Option(
        None, "--trait-groupby", help="With --trait: also per value of this adata.obs column."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
) -> None:
    """WGCNA-style modules per cell type (the hdWGCNA workflow): genes
    expressed in enough cells, metacells, TOM modules, eigengenes, kME and
    hub genes in OUTPUT/<cell type>/. With --consensus-by, one consensus
    network over the samples of each type (sets.tsv, set_eigengenes.tsv);
    with --harmonize-by, harmonized cell eigengenes; --dme-by and --trait
    test and correlate the (harmonized) cell eigengenes."""
    from gcen2.exceptions import GCENInputError
    from gcen2.single_cell import iter_modules_by_celltype, write_celltype_modules

    if consensus_quantile is not None and consensus_by is None:
        fail("--consensus-quantile needs --consensus-by.")
    if not cell_level and (harmonize_by is not None or dme_by is not None or traits):
        fail("--harmonize-by, --dme-by and --trait need the cell eigengenes.")
    if trait_groupby is not None and not traits:
        fail("--trait-groupby needs --trait.")
    adata = read_h5ad(input, backed=backed)
    given = dict(min_size=min_size, deep_split=deep_split, merge_cut=merge_cut,
                 consensus_quantile=consensus_quantile)
    # one cell type at a time: written, then released
    res = iter_modules_by_celltype(
        adata, groupby, groups=groups, fraction=fraction, use_highly_variable=hvg, k=k,
        max_shared=max_shared, max_metacells=max_metacells, min_cells=min_cells,
        use_rep=use_rep, n_pcs=n_pcs, layer=layer, normalize=normalize, seed=seed,
        threads=threads, preset=preset, cell_level=cell_level, reassign=reassign,
        consensus_by=consensus_by, harmonize_by=harmonize_by, npcs=npcs, dme_by=dme_by,
        traits=traits, trait_groupby=trait_groupby, power=power,
        **{key: v for key, v in given.items() if v is not None})
    table = write_celltype_modules(res, output, n_hubs=n_hubs, write_metacells=write_metacells,
                                   write_tom=write_tom, threads=threads)
    if table.empty:
        raise GCENInputError("no cell type had enough cells (see --min-cells).")
    typer.echo(f"Modules of {len(table)} cell types ({preset}) → {output}/")
    for row in table.itertuples():
        typer.echo(f"  {row.cell_type}: {row.n_cells} cells, {row.n_metacells} metacells, "
                   f"{row.n_genes} genes, power {fmt_power(row.power)}, {row.modules} modules "
                   f"({row.genes_in_modules} genes; {row.reassigned} reassigned)")
