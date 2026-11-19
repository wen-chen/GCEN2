"""Spatial commands: spatial-modules (the hdWGCNA spatial workflow) and
spatial-build (networks filtered by, or smoothed with, the spatial
neighbourhood, or per region)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import THREADS_HELP, app, fail, read_h5ad


@app.command("spatial-modules")
def cmd_spatial_modules(
    input: Path = typer.Option(..., "-i", "--input",
                               help="Spatial AnnData (.h5ad) with counts: Visium (obs array_row / "
                                    "array_col) or any platform with obsm['spatial']."),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    groupby: Optional[str] = typer.Option(
        None, "--groupby", help="adata.obs column (regions / clusters); metaspots are built "
                                "within each group (default: all spots together)."),
    fraction: float = typer.Option(
        0.05, "--fraction", help="Genes with counts in at least round(fraction * spots) spots "
                                 "(hdWGCNA gene_select = 'fraction')."),
    hvg: bool = typer.Option(False, "--hvg",
                             help="Only genes flagged in adata.var['highly_variable']."),
    genes_file: Optional[Path] = typer.Option(
        None, "--genes", help="Use these genes (one per line) instead of --fraction / --hvg."),
    min_spots: int = typer.Option(50, "--min-spots", help="Skip groups with fewer spots."),
    mode: str = typer.Option("sum", "--mode", help="Metaspot aggregation: sum or average."),
    metaspot_method: str = typer.Option(
        "auto", "--metaspots",
        help="grid (Visium, hdWGCNA MetaspotsByGroups), knn (each centre with its k - 1 "
             "nearest spots in obsm[--spatial-key]; any platform) or auto (grid when "
             "obs has a Visium hexagonal grid in --row-key / --col-key, else knn)."),
    row_key: str = typer.Option("array_row", "--row-key", help="obs column: Visium grid row."),
    col_key: str = typer.Option("array_col", "--col-key",
                                help="obs column: Visium grid column."),
    spatial_key: str = typer.Option("spatial", "--spatial-key",
                                    help="obsm entry with the coordinates (knn metaspots, "
                                         "--plot)."),
    k: int = typer.Option(
        7, "-k", "--k-neighbors",
        help="knn metaspots: spots per metaspot (the centre included); 7 matches a Visium "
             "hexagon, about 100 for 10 um Slide-seq beads."),
    max_shared: Optional[int] = typer.Option(
        None, "--max-shared",
        help="knn metaspots: at most this many spots shared by two metaspots "
             "(default round(2k/7), 2 for k = 7)."),
    seed: int = typer.Option(42, "--seed", help="knn metaspots: seed for the centre order."),
    layer: Optional[str] = typer.Option(None, "--layer", help="adata.layers entry instead of X."),
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
    spot_level: bool = typer.Option(
        True, "--spot-eigengenes/--no-spot-eigengenes",
        help="Eigengenes of every spot (hdWGCNA ModuleEigengenes), spot kME and group "
             "correlations."),
    reassign: bool = typer.Option(
        True, "--reassign/--no-reassign",
        help="Move genes with a negative spot kME for their module (hdWGCNA "
             "ModuleConnectivity); needs the spot eigengenes."),
    plot: bool = typer.Option(False, "--plot",
                              help="Also draw spatial_eigengenes.png (needs obsm[--spatial-key])."),
    write_metaspots: bool = typer.Option(
        False, "--write-metaspots",
        help="Also write metaspots.tsv (genes x metaspots; gcen2 module -e input)."),
    write_tom: bool = typer.Option(
        False, "--write-tom",
        help="Also write tom_hubs.tsv (for the hub-network and module-UMAP figures)."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
) -> None:
    """WGCNA-style modules of a spatial sample (the hdWGCNA spatial workflow):
    metaspots within groups (Visium grid, or nearest spots for any
    platform), TOM modules, eigengenes, kME, hub genes and spot-level
    eigengenes in OUTPUT/."""
    from gcen2.spatial import spatial_modules, write_spatial_modules

    adata = read_h5ad(input)
    genes = None
    if genes_file is not None:
        if not genes_file.is_file():
            fail(f"genes file not found: {genes_file}")
        genes = [g.strip() for g in genes_file.read_text(encoding="utf-8").splitlines()
                 if g.strip()]
    if plot and spatial_key not in adata.obsm:
        fail(f"--plot needs adata.obsm['{spatial_key}'].")
    given = dict(min_size=min_size, deep_split=deep_split, merge_cut=merge_cut)
    res = spatial_modules(
        adata, groupby, fraction=fraction, use_highly_variable=hvg, genes=genes,
        min_spots=min_spots, mode=mode, row_key=row_key, col_key=col_key, layer=layer,
        threads=threads, preset=preset, spot_level=spot_level, reassign=reassign,
        method=metaspot_method, spatial_key=spatial_key, k=k, max_shared=max_shared, seed=seed,
        power=power, **{key: v for key, v in given.items() if v is not None})
    summary = write_spatial_modules(
        res, output, groups=adata.obs[groupby].astype(str).to_numpy() if groupby else None,
        n_hubs=n_hubs, write_metaspots=write_metaspots, write_tom=write_tom, threads=threads)
    me = res.metadata.get("spot_eigengenes")
    if plot and me is not None and me.shape[1]:
        import matplotlib.pyplot as plt
        import numpy as np
        import pandas as pd

        from gcen2.plot.single_cell import spatial_me

        me = me.rename(columns={c: f"ME{int(c[2:]) - 1}" for c in me.columns})
        xy = pd.DataFrame(np.asarray(adata.obsm[spatial_key])[:, :2], index=me.index)
        plt.close(spatial_me(me, xy, output=output / "spatial_eigengenes",
                             formats=("pdf", "png", "svg")))
    typer.echo(f"{summary['n_spots']} spots, {summary['n_metaspots']} {summary['metaspots']} "
               f"metaspots, {summary['n_genes']} genes, power {summary['power']:g}: "
               f"{summary['modules']} modules ({summary['genes_in_modules']} genes; "
               f"{summary['reassigned']} with a negative spot kME reassigned) → {output}/")


@app.command("spatial-build")
def cmd_spatial_build(
    input: Path = typer.Option(..., "-i", "--input",
                               help="Input AnnData file (.h5ad) with obsm['spatial']."),
    output: Path = typer.Option(..., "-o", "--output",
                                help="Output .network file (filtered/weighted) or directory "
                                     "(regions)."),
    mode: str = typer.Option("filtered", "-m", "--mode",
                             help="Analysis mode: filtered | weighted | regions."),
    n_neighbors: Optional[int] = typer.Option(
        6, "--n-neighbors", help="KNN neighbours. Use --radius to switch to radius mode."),
    radius: Optional[float] = typer.Option(None, "--radius",
                                           help="Radius cutoff (overrides --n-neighbors)."),
    bandwidth: Optional[float] = typer.Option(
        None, "--bandwidth", help="Gaussian bandwidth for weighted mode (auto if omitted)."),
    region_key: Optional[str] = typer.Option(
        None, "--region-key", help="obs column for region labels (regions mode)."),
    n_regions: int = typer.Option(5, "--n-regions",
                                  help="Number of auto clusters when --region-key is omitted."),
    cor_cutoff: float = typer.Option(0.6, "-c", "--cor-cutoff",
                                     help="Minimum absolute correlation."),
    p_cutoff: float = typer.Option(0.05, "-p", "--p-cutoff", help="Maximum p-value."),
    method: str = typer.Option("spearman", "--method", help="Correlation: spearman | pearson."),
    layer: Optional[str] = typer.Option(None, "--layer", help="adata.layers entry instead of X."),
    fraction: float = typer.Option(
        0.05, "--fraction", help="Genes with a value > 0 in at least round(fraction * spots) "
                                 "spots (0: all genes)."),
    normalize: Optional[bool] = typer.Option(
        None, "--normalize/--no-normalize",
        help="Log-normalize (counts per 10,000; default: when the data are integer counts)."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
    kernel: str = typer.Option("auto", "--kernel",
                               help="Correlation kernel: auto, pairwise or gemm (see build)."),
    moran: Optional[Path] = typer.Option(
        None, "--moran",
        help="Also write each gene's Moran's I (gene, morans_i, pval, fdr, mean, detected; "
             "weighted mode: also morans_i_smoothed) to this TSV (gcen2 plot moran)."),
    moran_fdr: float = typer.Option(
        0.05, "--moran-fdr",
        help="filtered: keep edges between genes with Moran's I > 0 and FDR <= this "
             "(one-sided test; 1: I > 0 only)."),
) -> None:
    """Spatial transcriptomics co-expression network analysis."""
    from gcen2 import spatial

    if mode not in ("filtered", "weighted", "regions"):
        fail(f"Unknown mode '{mode}'. Choose: filtered | weighted | regions")
    adata = read_h5ad(input)
    typer.echo(f"Loaded {adata.n_obs} spots × {adata.n_vars} genes from {input}")
    # radius overrides n_neighbors
    neighbours = ({"n_neighbors": None, "radius": radius} if radius is not None
                  else {"n_neighbors": n_neighbors})
    common = dict(neighbours, layer=layer, fraction=fraction, normalize=normalize,
                  method=method, threads=threads, kernel=kernel)
    if mode == "regions":
        # build_spatial_regions passes its keywords on to gcen2.network.build
        output.mkdir(parents=True, exist_ok=True)
        nets = spatial.build_spatial_regions(
            adata, cor=cor_cutoff, pval=p_cutoff, region_key=region_key,
            **({} if region_key is not None else {"n_regions": n_regions}), **common)
        for name, net in nets.items():
            net.to_edgelist(output / f"{name}.network")
        typer.echo(f"Regions network: {len(nets)} regions → {output}/")
        for name, net in sorted(nets.items()):
            typer.echo(f"  {name}: {net.node_count} nodes, {net.edge_count} edges")
    else:
        if mode == "filtered":
            net = spatial.build_spatial_filtered(adata, cor_cutoff=cor_cutoff,
                                                 p_cutoff=p_cutoff, moran_fdr=moran_fdr,
                                                 **common)
        else:
            net = spatial.build_spatial_weighted(
                adata, cor_cutoff=cor_cutoff, p_cutoff=p_cutoff,
                **({} if bandwidth is None else {"bandwidth": bandwidth}), **common)
        net.to_edgelist(output)
        typer.echo(f"{mode.capitalize()} network: {net.node_count} nodes, {net.edge_count} "
                   f"edges → {output}")
    if moran is not None:
        mi = spatial.spatial_autocorrelation(
            adata, layer=layer, fraction=fraction, normalize=normalize,
            n_neighbors=neighbours["n_neighbors"], radius=radius,
            smoothed=mode == "weighted", bandwidth=bandwidth)
        mi.to_csv(moran, sep="\t", index=False, float_format="%.6g", lineterminator="\n")
        sig = int(((mi["morans_i"] > 0) & (mi["fdr"] <= moran_fdr)).sum())
        typer.echo(f"Moran's I of {len(mi)} genes ({int((mi['morans_i'] > 0).sum())} "
                   f"with I > 0, {sig} with FDR <= {moran_fdr:g}) → {moran}")
