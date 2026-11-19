"""The report command: an HTML or PDF report with vector figures."""

from __future__ import annotations

import os
from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import app, fail


@app.command("report")
def cmd_report(
    network: Optional[Path] = typer.Option(None, "-n", "--network", help="Input network TSV."),
    modules: Optional[Path] = typer.Option(None, "-m", "--modules", help="Module TSV file."),
    output: Path = typer.Option("report.html", "-o", "--output", help="Output report file."),
    me: Optional[Path] = typer.Option(
        None, "--me", help="Eigengenes of -m (gcen2 module -e --me): WGCNA module section."),
    kme: Optional[Path] = typer.Option(None, "--kme", help="kME file of -m (module -e --kme)."),
    hubs: Optional[Path] = typer.Option(None, "--hubs",
                                        help="Hub genes of -m (module -e --hubs)."),
    trait_cor: Optional[Path] = typer.Option(
        None, "--trait-cor", help="Eigengene-trait correlations of -m (module -e --trait-cor)."),
    celltype_modules: Optional[Path] = typer.Option(
        None, "--celltype-modules",
        help="Output directory of gcen2 celltype-modules: one module section per cell type."),
    n_hubs: int = typer.Option(10, "--n-hubs", help="Hub genes shown per module."),
    soft_power: Optional[Path] = typer.Option(
        None, "--soft-power", help="Soft-power fits of -m (module -e --soft-power)."),
    dendrogram: Optional[Path] = typer.Option(
        None, "--dendrogram", help="Gene tree of -m (module -e --dendrogram)."),
    tom_hubs: Optional[Path] = typer.Option(
        None, "--tom-hubs", help="Hub TOM of -m (module -e --tom-hubs; hub-network figure)."),
    figures_dir: Optional[Path] = typer.Option(
        None, "--figures-dir",
        help="Also save every report figure there as PDF, PNG and SVG (for a manuscript)."),
    coords: Optional[Path] = typer.Option(
        None, "--coords",
        help="With --celltype-modules: cells' embedding (.h5ad with obsm[--basis], or TSV of "
             "x, y by cell) for the metacell figure."),
    basis: str = typer.Option("X_umap", "--basis", help="obsm key of --coords."),
    moran: Optional[Path] = typer.Option(
        None, "--moran", help="Moran's I table of spatial-build --moran."),
    merge_inputs: Optional[list[Path]] = typer.Option(
        None, "--merge-input", help="Input network of gcen2 merge; repeat for each."),
    merged: Optional[Path] = typer.Option(
        None, "--merged", help="With --merge-input: the gcen2 merge output."),
    extract_genes: Optional[Path] = typer.Option(
        None, "--extract-genes",
        help="Gene list of gcen2 extract; the network is -n (full network or subnetwork)."),
) -> None:
    """Generate a comprehensive HTML research report (or PDF with a .pdf
    output) with vector figures and captions. With --me / --kme / --hubs /
    --trait-cor / --soft-power / --dendrogram / --tom-hubs, -m is shown as
    WGCNA-style modules with their figures (module sizes, soft power,
    dendrogram, eigengene network, module-trait, kME, hub network)."""
    # figures only go into the file: no GUI backend, whose display probe
    # costs ~1.5 s over SSH X forwarding (an explicit MPLBACKEND wins)
    os.environ.setdefault("MPLBACKEND", "Agg")
    from gcen2 import report
    from gcen2._matrix_io import require_file

    wgcna_files = (me, kme, hubs, trait_cor, soft_power, dendrogram, tom_hubs)
    extras = any(v is not None for v in wgcna_files)
    if extras and modules is None:
        fail("--me, --kme, --hubs, --trait-cor, --soft-power, --dendrogram and --tom-hubs "
             "describe -m: give -m too.")
    if coords is not None and celltype_modules is None:
        fail("--coords is for the metacells of --celltype-modules.")
    if merged is not None and not merge_inputs:
        fail("--merged needs the networks merged (--merge-input).")
    if extract_genes is not None and network is None:
        fail("--extract-genes needs the network (-n).")
    for p in [moran, merged, extract_genes, *(merge_inputs or [])]:
        if p is not None:
            require_file(p)

    net_obj = None
    if network:
        # statistics and degrees read in C++ (no edge table); same report as
        # with read_network(network)
        from gcen2.network import summary_from_file

        net_obj = summary_from_file(network)
    wgcna, mod_obj, celltypes = [], None, None
    if modules and extras:
        wgcna.append(report.wgcna_section(modules, eigengenes=me, kme=kme, hubs=hubs,
                                          trait_cor=trait_cor, n_hubs=n_hubs,
                                          soft_power=soft_power, dendrogram=dendrogram,
                                          tom_hubs=tom_hubs))
    elif modules:
        import pandas as pd

        from gcen2.models import ModuleSet

        mod_df = pd.read_csv(require_file(modules), sep="\t")
        mod_obj = ModuleSet(modules=mod_df.groupby("module")["gene"].apply(set).tolist())
    if celltype_modules is not None:
        wgcna.extend(report.celltype_module_sections(celltype_modules, n_hubs=n_hubs))
        celltypes = report.celltype_overview(celltype_modules, coords=coords, basis=basis)
    report.generate(
        network=net_obj, modules=mod_obj, output=output, wgcna=wgcna or None,
        figures_dir=figures_dir, celltypes=celltypes, moran=moran,
        merge=({"networks": list(merge_inputs), "merged": merged} if merge_inputs else None),
        extract=({"network": network, "genes": extract_genes} if extract_genes else None))
    typer.echo(f"Report written to {output}")
    if figures_dir is not None:
        n = len(list(Path(figures_dir).glob("*.pdf")))
        typer.echo(f"{n} figures written to {figures_dir}/ (PDF, PNG, SVG)")
