"""WGCNA-style module commands: ``module -e`` (:func:`wgcna_modules`),
consensus-modules, module-preservation and project-modules."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import (
    THREADS_HELP,
    app,
    check_groups,
    echo_tom,
    fail,
    h5ad_only,
    module_options,
    read_h5ad,
    read_module_table,
)


def _expression_input(expr: Path, layer, hvg, groupby, group):
    """(matrix, genes, samples) of a genes x samples TSV or an .h5ad."""
    if expr.suffix == ".h5ad":
        from gcen2._anndata import adata_cell_names, adata_matrix

        check_groups(groupby, group)
        adata = read_h5ad(expr)
        matrix, genes = adata_matrix(adata, layer=layer, use_highly_variable=hvg,
                                     groupby=groupby, groups=group)
        return matrix, genes, adata_cell_names(adata, groupby, group)
    h5ad_only(layer=layer, hvg=hvg, groupby=groupby, group=group)
    from gcen2._matrix_io import read_matrix

    return read_matrix(expr)


def _read_traits(path: Path):
    import pandas as pd

    from gcen2._matrix_io import require_file

    return pd.read_csv(require_file(path), sep="\t", index_col=0)


def wgcna_modules(expr: Path, output: Path, *, method, signed, power, min_size, deep_split,
                  merge_cut, log, threads, preset, me_out, kme_out, hubs_out, n_hubs, traits,
                  trait_out, layer, hvg, groupby, group, soft_power_out, dendrogram_out,
                  tom_hubs_out, gs_out) -> None:
    """``gcen2 module -e``: WGCNA modules of an expression matrix and the
    files asked for."""
    from gcen2.tom_modules import HDWGCNA, gene_significance, tom_modules, write_module_files

    opts = module_options(preset, HDWGCNA, signed=signed, min_size=min_size,
                          deep_split=deep_split, merge_cut=merge_cut)
    if gs_out is not None and traits is None:
        fail("--gs needs --traits.")
    if (traits is None) != (trait_out is None) and gs_out is None:
        fail("--traits and --trait-cor go together.")
    trait_df = None if traits is None else _read_traits(traits)
    matrix, genes, samples = _expression_input(expr, layer, hvg, groupby, group)
    res = tom_modules(matrix, gene_names=genes, sample_names=samples, method=method,
                      power=power, log=log, threads=threads, **opts)
    del matrix
    echo_tom({"power": res.power, "estimated": res.power_estimated, "fits": res.fits,
              "rule": res.power_rule})
    write_module_files(
        res, modules=output, eigengenes=me_out, kme=kme_out, hubs=hubs_out, n_hubs=n_hubs,
        trait_table=None if trait_out is None else res.trait_cor(trait_df),
        trait_cor=trait_out, soft_power=soft_power_out, dendrogram=dendrogram_out,
        tom_hubs=tom_hubs_out, gs=gs_out,
        gs_table=None if gs_out is None else gene_significance(res, trait_df),
        threads=threads)
    n_genes = int((res.labels > 0).sum())
    typer.echo(f"{int(res.labels.max(initial=0))} modules ({n_genes} of {len(res.genes)} genes) "
               f"written to {output}")
    if me_out is not None:
        var = ", ".join(f"ME{k - 1 if k else '_unassigned'} {v:.2f}"
                        for k, v in zip(res.me_labels, res.var_explained))
        typer.echo(f"Eigengenes written to {me_out} (variance explained: {var})")
    for path, what in ((kme_out, "kME"), (hubs_out, "Hub genes"),
                       (trait_out, "Eigengene-trait correlations"),
                       (soft_power_out, "Soft-power fits"), (dendrogram_out, "Gene tree"),
                       (tom_hubs_out, "Hub TOM"), (gs_out, "Gene significance")):
        if path is not None:
            typer.echo(f"{what} written to {path}")


@app.command("consensus-modules")
def cmd_consensus_modules(
    expr: list[Path] = typer.Option(
        ..., "-e", "--expr", help="Expression matrix of one set (genes x samples TSV); "
                                  "repeat for each set (same genes)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    set_names: Optional[list[str]] = typer.Option(
        None, "--set-name", help="Name of each set, in the order of -e (default: file names)."),
    preset: str = typer.Option(
        "wgcna", "--preset",
        help="wgcna (blockwiseConsensusModules defaults) or hdwgcna (hdWGCNA "
             "ConstructNetwork(consensus = TRUE): signed, consensus quantile 0.3, deepSplit 4, "
             "no PAM stage, min size 50, merge cut 0.2, each set's power from bicor with "
             "R^2 >= 0.8 and power > 3)."),
    signed: bool = typer.Option(False, "--signed", help="Signed networks (wgcna preset)."),
    power: str = typer.Option(
        "auto", "--power", help="Soft-threshold power: 'auto' (each set's own), one number, or "
                                "one per set (comma-separated)."),
    consensus_quantile: Optional[float] = typer.Option(
        None, "--consensus-quantile",
        help="Quantile of the sets' TOMs (0 = minimum) [default: 0; hdwgcna 0.3]."),
    min_size: Optional[int] = typer.Option(None, "--min-size", help="Minimum module size."),
    deep_split: Optional[int] = typer.Option(None, "--deep-split",
                                             help="Dynamic tree cut deepSplit, 0-4."),
    merge_cut: Optional[float] = typer.Option(
        None, "--merge-cut", help="Merge modules whose consensus eigengene dissimilarity is "
                                  "below this."),
    log: Optional[str] = typer.Option(None, "-l", "--log",
                                      help="log(x+1) transform first: log, log2 or log10."),
    n_hubs: int = typer.Option(10, "--n-hubs", help="Hub genes per module."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
) -> None:
    """Consensus modules of several expression sets (WGCNA
    blockwiseConsensusModules): modules.tsv, eigengenes.tsv (all samples),
    kme.tsv, hub_genes.tsv, sets.tsv (powers, calibration) and
    set_eigengenes.tsv in OUTPUT/."""
    import numpy as np

    from gcen2._matrix_io import read_matrix
    from gcen2.consensus import HDWGCNA_CONSENSUS, consensus_tom_modules
    from gcen2.exceptions import GCENInputError
    from gcen2.network import _LOG_FUNCS
    from gcen2.single_cell.modules import _write_sets
    from gcen2.tom_modules import write_module_files

    opts = module_options(preset, HDWGCNA_CONSENSUS, signed=signed, min_size=min_size,
                          deep_split=deep_split, merge_cut=merge_cut,
                          consensus_quantile=consensus_quantile)
    pw = [x.strip() for x in power.split(",")]
    opts["power"] = pw[0] if len(pw) == 1 else pw
    if log is not None and log not in _LOG_FUNCS:
        raise GCENInputError(f"Unknown log base '{log}'. Choose from: {', '.join(_LOG_FUNCS)}")
    mats, samples, genes = [], [], None
    for path in expr:
        m, g, smp = read_matrix(path)
        if log is not None:
            m = _LOG_FUNCS[log](np.asarray(m, dtype=np.float64) + 1.0)
        if genes is None:
            genes = list(g)
        elif list(g) != genes:
            raise GCENInputError(f"{path} does not have the genes of {expr[0]} in the "
                                 "same order.")
        mats.append(m)
        samples.append(list(smp))
    res = consensus_tom_modules(mats, gene_names=genes,
                                set_names=list(set_names) if set_names else [p.stem for p in expr],
                                sample_names=samples, threads=threads, **opts)
    output.mkdir(parents=True, exist_ok=True)
    write_module_files(res, modules=output / "modules.tsv",
                       eigengenes=output / "eigengenes.tsv", kme=output / "kme.tsv",
                       hubs=output / "hub_genes.tsv", n_hubs=n_hubs,
                       soft_power=output / "soft_power.tsv",
                       dendrogram=output / "dendrogram.json")
    _write_sets(res, output)
    md = res.metadata
    pw_txt = ", ".join(f"{s} {md['set_powers'][s]:g}" for s in md["sets"])
    typer.echo(f"{len(md['sets'])} sets (powers: {pw_txt}): "
               f"{int(res.labels.max(initial=0))} consensus modules "
               f"({int((res.labels > 0).sum())} of {len(res.genes)} genes) → {output}/")


@app.command("module-preservation")
def cmd_module_preservation(
    ref: Path = typer.Option(
        ..., "-r", "--ref", help="Reference expression the modules were built from (genes x "
                                 "samples TSV, e.g. metacells.tsv of celltype-modules "
                                 "--write-metacells)."),
    modules: Path = typer.Option(
        ..., "-m", "--modules", help="Reference modules: GCEN2 modules.tsv or a table with gene "
                                     "and module columns; genes not in it are unassigned."),
    query: Path = typer.Option(
        ..., "-q", "--query", help="Test expression: genes x samples TSV, or an .h5ad with counts "
                                   "(metacells are then built from it, as hdWGCNA's query)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    n_permutations: int = typer.Option(500, "--permutations",
                                       help="Permutations (hdWGCNA default 500)."),
    seed: int = typer.Option(12345, "--seed", help="Random seed (R's set.seed; hdWGCNA 12345)."),
    network_type: str = typer.Option("unsigned", "--network-type",
                                     help="unsigned (WGCNA's and hdWGCNA's), signed or "
                                          "'signed hybrid'."),
    groupby: Optional[str] = typer.Option(None, "--groupby",
                                          help=".h5ad query: obs column to select cells."),
    groups: Optional[list[str]] = typer.Option(
        None, "--group", help=".h5ad query: only this --groupby value; repeatable."),
    k: int = typer.Option(25, "-k", "--k", help=".h5ad query: cells per metacell."),
    max_shared: Optional[int] = typer.Option(None, "--max-shared",
                                             help=".h5ad query: max cells shared by metacells."),
    max_metacells: int = typer.Option(1000, "--max-metacells", help=".h5ad query."),
    use_rep: Optional[str] = typer.Option(None, "--use-rep",
                                          help=".h5ad query: obsm key for the neighbours."),
    n_pcs: int = typer.Option(30, "--n-pcs", help=".h5ad query: components for neighbours."),
    metacell_seed: int = typer.Option(42, "--metacell-seed", help=".h5ad query."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
) -> None:
    """Preservation of reference modules in another data set (WGCNA
    modulePreservation, as hdWGCNA's ModulePreservation): preservation.tsv
    (observed statistics) and preservation_Z.tsv (Z scores, Zsummary.pres)."""
    import pandas as pd

    from gcen2._matrix_io import read_matrix
    from gcen2.preservation import module_preservation

    mods = read_module_table(modules)
    m, g, s = read_matrix(ref)
    ref_df = pd.DataFrame(m.T, index=s, columns=g)
    if query.suffix == ".h5ad":
        from gcen2.single_cell import metacells

        adata = read_h5ad(query)
        check_groups(groupby, groups)
        mc = metacells(adata, groupby=groupby, k=k, max_shared=max_shared,
                       max_metacells=max_metacells, use_rep=use_rep, n_pcs=n_pcs,
                       seed=metacell_seed, threads=threads, groups=groups)
        keep = [x for x in ref_df.columns if x in set(map(str, mc.var_names))]
        x = mc[:, keep].X
        test_df = pd.DataFrame(x.toarray() if hasattr(x, "toarray") else x,
                               index=mc.obs_names, columns=keep)
        test_df = test_df.loc[:, test_df.std(axis=0) > 0]
    else:
        if groupby is not None or groups:
            fail("--groupby and --group need an .h5ad query.")
        m, g, s = read_matrix(query)
        test_df = pd.DataFrame(m.T, index=s, columns=g)
    obs, z = module_preservation(ref_df, test_df, mods, n_permutations=n_permutations,
                                 seed=seed, network_type=network_type, threads=threads)
    output.mkdir(parents=True, exist_ok=True)
    obs.to_csv(output / "preservation.tsv", sep="\t", float_format="%.10g", na_rep="NA", lineterminator="\n")
    z.to_csv(output / "preservation_Z.tsv", sep="\t", float_format="%.10g", na_rep="NA", lineterminator="\n")
    typer.echo(f"{len(z)} modules (gold and grey included), {n_permutations} permutations "
               f"→ {output}/")
    for name, row in z.iterrows():
        typer.echo(f"  {name}: {int(row['moduleSize'])} genes, Zsummary.pres "
                   f"{row['Zsummary.pres']:.2f}")


@app.command("project-modules")
def cmd_project_modules(
    input: Path = typer.Option(..., "-i", "--input", help="Query AnnData (.h5ad), counts."),
    modules: Path = typer.Option(
        ..., "-m", "--modules", help="Reference modules: GCEN2 modules.tsv (module, gene) or a "
                                     "table with gene and module columns (hdWGCNA)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    groupby: Optional[str] = typer.Option(None, "--groupby",
                                          help="adata.obs column to select query cells."),
    groups: Optional[list[str]] = typer.Option(
        None, "--group", help="Only cells with this --groupby value; repeatable."),
    min_overlap: float = typer.Option(
        0.5, "--min-overlap", help="Project modules with at least this fraction of their genes "
                                   "in the query (hdWGCNA overlap_proportion)."),
    harmonize_by: Optional[str] = typer.Option(
        None, "--harmonize-by", help="adata.obs column: also harmonized eigengenes (Harmony)."),
    layer: Optional[str] = typer.Option(None, "--layer", help="adata.layers entry instead of X."),
    npcs: int = typer.Option(50, "--npcs", help="With --harmonize-by: components per module."),
    seed: int = typer.Option(0, "--seed", help="Random seed (Harmony)."),
) -> None:
    """Eigengenes of reference modules on the cells of another data set
    (hdWGCNA ProjectModules): projected_eigengenes.tsv, with --harmonize-by
    projected_harmonized_eigengenes.tsv, and projected_modules.tsv."""
    from gcen2.exceptions import GCENInputError
    from gcen2.single_cell import project_modules

    check_groups(groupby, groups)
    mods = read_module_table(modules)
    adata = read_h5ad(input)
    cells = None
    if groupby is not None:
        if groupby not in adata.obs.columns:
            raise GCENInputError(f"'{groupby}' not found in adata.obs.")
        cells = adata.obs[groupby].astype(str).isin(groups).to_numpy()
    me, hme, table = project_modules(adata, mods, min_overlap=min_overlap, layer=layer,
                                     cells=cells, harmonize_by=harmonize_by, npcs=npcs,
                                     seed=seed)
    output.mkdir(parents=True, exist_ok=True)
    me.rename_axis("cell").to_csv(output / "projected_eigengenes.tsv", sep="\t",
                                  float_format="%.10g", lineterminator="\n")
    if hme is not None:
        hme.rename_axis("cell").to_csv(output / "projected_harmonized_eigengenes.tsv",
                                       sep="\t", float_format="%.10g", lineterminator="\n")
    table.to_csv(output / "projected_modules.tsv", sep="\t", index=False, lineterminator="\n")
    typer.echo(f"{me.shape[1]} of {len(table)} modules projected onto {me.shape[0]} cells "
               f"→ {output}/")
