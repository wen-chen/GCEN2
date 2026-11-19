"""The trajectory command: co-expression networks along pseudotime."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import THREADS_HELP, app, fail, read_h5ad


@app.command("trajectory")
def cmd_trajectory(
    input: Path = typer.Option(..., "-i", "--input",
                               help="AnnData (.h5ad) or expression TSV (genes x samples)."),
    pseudotime: str = typer.Option(
        ..., "--pseudotime",
        help="adata.obs column (.h5ad input), or a TSV file of sample / pseudotime pairs."),
    output: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    window: float = typer.Option(0.2, "-w", "--window",
                                 help="Cells per window (>= 3), or a fraction (< 1)."),
    step: Optional[float] = typer.Option(None, "-s", "--step",
                                         help="Shift between windows, count or fraction "
                                              "(default: half a window)."),
    metacell: int = typer.Option(1, "--metacell",
                                 help="Average this many consecutive cells into a metacell."),
    layer: Optional[str] = typer.Option(None, "--layer", help="AnnData layer instead of X."),
    method: str = typer.Option("spearman", "-m", "--method",
                               help="Correlation method: pearson or spearman."),
    cor: float = typer.Option(0.5, "-c", "--cor", help="Minimum correlation cutoff."),
    pval: float = typer.Option(0.001, "-p", "--pval", help="P-value cutoff."),
    fdr: bool = typer.Option(False, "-f", "--fdr", help="Apply FDR correction."),
    signed: bool = typer.Option(False, "--signed", help="Keep positive correlations only."),
    log: Optional[str] = typer.Option(None, "-l", "--log",
                                      help="log(x+1) transform first: log, log2 or log10."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
    plot: bool = typer.Option(False, "--plot", help="Also write trajectory.png."),
) -> None:
    """Co-expression networks in sliding windows along pseudotime."""
    import gcen2

    def as_count(v):
        return int(v) if v >= 1 else float(v)

    if input.suffix == ".h5ad":
        data = read_h5ad(input)
        pt = pseudotime
        if Path(pseudotime).is_file():
            pt = _read_pseudotime(Path(pseudotime)).reindex(data.obs_names.astype(str))
    else:
        data = gcen2.read_tsv(input).data
        pt = _read_pseudotime(Path(pseudotime))
    traj = gcen2.trajectory.build_trajectory(
        data, pt, window_size=as_count(window), step=None if step is None else as_count(step),
        metacell_size=metacell, layer=layer, method=method, cor=cor, pval=pval, signed=signed,
        fdr=fdr, log=log, threads=threads)
    traj.write(output)
    if plot:
        import matplotlib.pyplot as plt

        from gcen2.plot.trajectory import trajectory_overview

        plt.close(trajectory_overview(output, output=output / "trajectory",
                                      formats=("pdf", "png", "svg")))
    typer.echo(f"{len(traj.windows)} windows of {traj.params['window_size']} cells "
               f"(step {traj.params['step']}) → {output}/")
    for w in traj.windows:
        typer.echo(f"  window {w.index} [{w.start:.4g}, {w.end:.4g}]: "
                   f"{w.network.node_count} nodes, {w.network.edge_count} edges")
    counts = traj.dynamics["pattern"].value_counts()
    typer.echo("Edge patterns: " + ", ".join(
        f"{p} {int(counts.get(p, 0))}" for p in gcen2.trajectory.PATTERNS))


def _read_pseudotime(path: Path):
    """Sample / pseudotime pairs (first two columns; a header line is skipped
    when its second field is not a number)."""
    import pandas as pd

    from gcen2._matrix_io import require_file

    df = pd.read_csv(require_file(path), sep="\t", header=None, dtype=str, comment="#")
    if df.shape[1] < 2:
        fail(f"{path} needs two tab-separated columns: sample, pseudotime.")
    values = pd.to_numeric(df[1], errors="coerce")
    if len(df) and pd.isna(values.iloc[0]) and df[1].iloc[0].strip().lower() not in ("nan", ""):
        df, values = df.iloc[1:], values.iloc[1:]
    return pd.Series(values.to_numpy(dtype=float), index=df[0].astype(str).to_numpy())
