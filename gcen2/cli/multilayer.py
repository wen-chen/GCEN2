"""Multi-omics commands: fuse (weighted fusion of network layers) and
rwr-multiplex (random walk on a multiplex network)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import app, fail


def _layers(inputs: list[Path], names: Optional[list[str]]) -> dict:
    """``{name: path}`` for repeated -i / -n (and optional --name) options;
    the files are read (in C++ when possible) by gcen2.multilayer."""
    from gcen2._matrix_io import require_file

    if names and len(names) != len(inputs):
        fail(f"{len(names)} --name for {len(inputs)} input files.")
    names = names or [p.stem for p in inputs]
    if len(set(names)) != len(names):
        fail(f"layer names must be unique (got {names}); use --name.")
    return {n: require_file(p) for n, p in zip(names, inputs)}


@app.command("fuse")
def cmd_fuse(
    inputs: list[Path] = typer.Option(..., "-i", "--input",
                                      help="Layer edge-list files (co-expression, PPI, "
                                           "regulatory, ...); repeat for each layer."),
    output: Path = typer.Option(..., "-o", "--output", help="Output fused network TSV."),
    names: Optional[list[str]] = typer.Option(None, "--name",
                                              help="Layer names, one per -i (default: "
                                                   "file names)."),
    weights: Optional[list[float]] = typer.Option(None, "-w", "--weight",
                                                  help="Layer weights, one per -i "
                                                       "(default 1)."),
    normalize: str = typer.Option("max", "--normalize",
                                  help="Per-layer weight scaling: max, rank or none."),
    min_layers: int = typer.Option(1, "--min-layers",
                                   help="Keep edges found in at least this many layers."),
    min_score: float = typer.Option(0.0, "--min-score", help="Minimum fused score."),
    threads: int = typer.Option(4, "-t", "--threads", help="Threads for writing the table."),
) -> None:
    """Weighted fusion of heterogeneous networks (multi-omics integration)."""
    from gcen2.multilayer import fuse_to_file

    layers = _layers(inputs, names)
    if weights and len(weights) != len(inputs):
        fail(f"{len(weights)} --weight for {len(inputs)} -i files.")
    # the table is written by the C++ extension from gene codes
    r = fuse_to_file(layers, output, weights=dict(zip(layers, weights)) if weights else None,
                     normalize=normalize, min_layers=min_layers, min_score=min_score,
                     threads=threads)
    typer.echo(f"Fused network: {r['node_count']} nodes, {r['edge_count']} edges → {output}")
    typer.echo("Edges by number of layers: "
               + ", ".join(f"{k}: {v}" for k, v in r["layer_counts"].items()))


@app.command("rwr-multiplex")
def cmd_rwr_multiplex(
    inputs: list[Path] = typer.Option(..., "-n", "--network",
                                      help="Layer edge-list files; repeat for each layer."),
    genes: Path = typer.Option(..., "-g", "--genes", help="Seed genes list."),
    output: Path = typer.Option(..., "-o", "--output", help="Output TSV file."),
    names: Optional[list[str]] = typer.Option(None, "--name",
                                              help="Layer names, one per -n."),
    gamma: float = typer.Option(0.5, "--gamma", help="Restart probability."),
    delta: float = typer.Option(0.5, "--delta",
                                help="Weight of the jumps between layers (their probability "
                                     "with --transition layer)."),
    weighted: bool = typer.Option(False, "--weighted",
                                  help="Use layer edge weights for transitions."),
    aggregate: str = typer.Option("geometric", "--aggregate",
                                  help="Per-gene score: geometric, arithmetic or sum."),
    prob_cutoff: float = typer.Option(0.0, "--prob-cutoff", help="Minimum score."),
    transition: str = typer.Option(
        "supra", "--transition",
        help="supra: normalise the whole supra-adjacency matrix (RandomWalkRestartMH); "
             "layer: each layer first, every gene jumps with probability delta."),
) -> None:
    """Random walk with restart on a multiplex (multi-layer) network (RWR-M)."""
    from gcen2.io import read_gene_list
    from gcen2.multilayer import rwr_multiplex

    layers = _layers(inputs, names)
    res = rwr_multiplex(layers, read_gene_list(genes), gamma=gamma, delta=delta,
                        weighted=weighted, aggregate=aggregate, prob_cutoff=prob_cutoff,
                        transition=transition)
    res.ranked_genes.to_csv(output, sep="\t", index=False, lineterminator="\n")
    typer.echo(f"Multiplex RWR ({len(layers)} layers, {res.steps} steps): "
               f"{len(res.ranked_genes)} genes → {output}")
