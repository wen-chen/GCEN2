"""The GCEN v1 commands: normalize, filter, build, module, annotate, enrich,
rwr, merge, extract, shuffle, stats, diff. The network commands read and
write their files in C++ where possible (no pandas)."""

from __future__ import annotations

from pathlib import Path
from typing import Optional

import typer

from gcen2.cli._common import THREADS_HELP, app, check_groups, echo_tom, fail, h5ad_only


@app.command("normalize")
def cmd_normalize(
    input: Path = typer.Option(..., "-i", "--input", help="Input expression TSV file."),
    output: Path = typer.Option(..., "-o", "--output", help="Output normalized TSV file."),
    method: str = typer.Option("tmm", "-m", "--method",
                               help="Normalization method: upqt, median, deseq, tmm, hkg."),
    log: Optional[str] = typer.Option(None, "-l", "--log",
                                      help="Log transform: log, log2, log10."),
) -> None:
    """Normalize gene expression data."""
    import gcen2

    norm = gcen2.preprocessing.normalize(gcen2.read_tsv(input), method=method)
    if log:
        norm = gcen2.preprocessing.log_transform(norm, base=log)
    gcen2.write_tsv(norm.data, output)
    typer.echo(f"Normalized expression written to {output}")


@app.command("filter")
def cmd_filter(
    input: Path = typer.Option(..., "-i", "--input", help="Input expression TSV file."),
    output: Path = typer.Option(..., "-o", "--output", help="Output filtered TSV file."),
    percent_mean: float = typer.Option(1.0, "-p", "--percent-mean",
                                       help="Retain top fraction by mean expression."),
    percent_sd: float = typer.Option(1.0, "--percent-sd",
                                     help="Retain top fraction by standard deviation."),
    cutoff_mean: float = typer.Option(0.0, "--cutoff-mean", help="Minimum mean expression."),
    cutoff_sd: float = typer.Option(0.0, "--cutoff-sd", help="Minimum standard deviation."),
) -> None:
    """Filter genes by expression statistics."""
    import gcen2

    filtered = gcen2.preprocessing.filter_genes(
        gcen2.read_tsv(input), cutoff_mean=cutoff_mean, cutoff_sd=cutoff_sd,
        percent_mean=percent_mean, percent_sd=percent_sd)
    gcen2.write_tsv(filtered.data, output)
    typer.echo(f"Filtered expression written to {output}")


@app.command("build")
def cmd_build(
    input: Path = typer.Option(..., "-i", "--input",
                               help="Input expression TSV file, or AnnData (.h5ad)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output network TSV file."),
    method: str = typer.Option("spearman", "-m", "--method",
                                help="Correlation method: pearson, spearman, pearson_nonzero "
                                     "(only cells where both genes are non-zero) or "
                                     "pearson_pairwise_complete."),
    cor: float = typer.Option(0.1, "-c", "--cor", help="Minimum correlation cutoff."),
    pval: float = typer.Option(0.001, "-p", "--pval", help="P-value cutoff."),
    fdr: bool = typer.Option(False, "-f", "--fdr", help="Apply FDR correction."),
    signed: bool = typer.Option(
        False, "-s", "--signed",
        help="Signed network: keep positive correlations only (as GCEN v1 -s)."),
    all_edges: bool = typer.Option(False, "-a", "--all",
                                   help="Output all edges without any cutoff (-c/-p ignored)."),
    log: Optional[str] = typer.Option(None, "-l", "--log",
                                      help="log(x+1) transform first: log, log2 or log10."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
    max_pairs: Optional[int] = typer.Option(
        None, "--max-pairs",
        help="Gene pairs computed per block; edges are written block by block "
             "(default 500000)."),
    layer: Optional[str] = typer.Option(None, "--layer",
                                        help=".h5ad: adata.layers entry instead of X."),
    hvg: bool = typer.Option(False, "--hvg",
                             help=".h5ad: only genes flagged in adata.var['highly_variable']."),
    groupby: Optional[str] = typer.Option(None, "--groupby",
                                          help=".h5ad: adata.obs column to select cells by."),
    group: Optional[list[str]] = typer.Option(
        None, "--group", help=".h5ad: keep cells whose --groupby value is this; repeatable."),
    tom: bool = typer.Option(
        False, "--tom",
        help="Add WGCNA's topological overlap as a 'tom' column (pearson/spearman; "
             "-s: signed adjacency). Holds a genes x genes float matrix in memory."),
    power: str = typer.Option(
        "auto", "--power",
        help="--tom: soft-threshold power, or 'auto' (WGCNA pickSoftThreshold)."),
    tom_min: float = typer.Option(0.0, "--tom-min",
                                  help="--tom: drop edges with TOM below this."),
    kernel: str = typer.Option(
        "auto", "--kernel",
        help="pearson/spearman kernel: pairwise (as GCEN v1), gemm (matrix product, ~10x "
             "faster with hundreds of samples; r equal to rounding) or auto (gemm from 64 "
             "samples)."),
) -> None:
    """Build a gene co-expression network (streamed to the output file)."""
    opts = dict(method=method, cor=cor, pval=pval, fdr=fdr, signed=signed, log=log,
                all_edges=all_edges, threads=threads, tom=tom, power=power, tom_min=tom_min,
                kernel=kernel)
    if input.suffix == ".h5ad":
        from gcen2._anndata import adata_matrix
        from gcen2._build_file import MAX_PAIRS
        from gcen2.cli._common import read_h5ad
        from gcen2.network import build_to_file

        check_groups(groupby, group)
        matrix, genes = adata_matrix(read_h5ad(input), layer=layer, use_highly_variable=hvg,
                                     groupby=groupby, groups=group)
        info = build_to_file(matrix, output, gene_names=genes,
                             max_pairs_per_block=max_pairs or MAX_PAIRS, **opts)
        typer.echo(f"{matrix.shape[1]} cells x {matrix.shape[0]} genes")
    else:
        h5ad_only(layer=layer, hvg=hvg, groupby=groupby, group=group)
        # plain matrices are built in C++ without numpy or pandas
        from gcen2._build_file import build_from_file

        info = build_from_file(input, output, max_pairs_per_block=max_pairs, **opts)
    echo_tom(info.get("tom"))
    typer.echo(f"Network with {info['node_count']} nodes and {info['edge_count']} edges "
               f"written to {output}")


@app.command("module")
def cmd_module(
    input: Optional[Path] = typer.Option(None, "-i", "--input", help="Input network TSV file."),
    output: Path = typer.Option(..., "-o", "--output", help="Output module TSV file."),
    similarity: float = typer.Option(0.5, "-s", "--similarity",
                                     help="Jaccard similarity cutoff."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
    expr: Optional[Path] = typer.Option(
        None, "-e", "--expr",
        help="Expression matrix instead of -i: WGCNA modules (TOM, average-linkage "
             "clustering, dynamic tree cut, eigengene merging; as blockwiseModules)."),
    method: str = typer.Option("pearson", "-m", "--method",
                               help="-e: correlation for the TOM: pearson or spearman."),
    signed: bool = typer.Option(False, "--signed",
                                help="-e: signed network ((1 + r) / 2)^power, else |r|^power."),
    power: str = typer.Option("auto", "--power",
                              help="-e: soft-threshold power, or 'auto' (pickSoftThreshold)."),
    min_size: Optional[int] = typer.Option(
        None, "--min-size", help="-e: minimum module size [default: 20; hdwgcna 50]."),
    deep_split: Optional[int] = typer.Option(
        None, "--deep-split", help="-e: dynamic tree cut deepSplit, 0-4 [default: 2; hdwgcna 4]."),
    merge_cut: Optional[float] = typer.Option(
        None, "--merge-cut", help="-e: merge modules whose eigengenes have 1 - r below this "
                                  "[default: 0.15; hdwgcna 0.2]."),
    preset: str = typer.Option(
        "wgcna", "--preset",
        help="-e: wgcna (blockwiseModules defaults) or hdwgcna (hdWGCNA ConstructNetwork "
             "defaults: signed, deepSplit 4, no PAM stage, minKMEtoStay 0.2, min size 50, "
             "merge cut 0.2, power from bicor with R^2 >= 0.8 and power > 3)."),
    log: Optional[str] = typer.Option(None, "-l", "--log",
                                      help="-e: log(x+1) transform first: log, log2 or log10."),
    me_out: Optional[Path] = typer.Option(
        None, "--me", help="-e: write the module eigengenes (samples x ME<k> for module k "
                           "of -o, and ME_unassigned) to this TSV file."),
    kme_out: Optional[Path] = typer.Option(
        None, "--kme", help="-e: write every gene's kME (correlation with each eigengene) "
                            "and its p-values to this TSV file."),
    hubs_out: Optional[Path] = typer.Option(
        None, "--hubs", help="-e: write the hub genes (highest kME) of each module to this "
                             "TSV file."),
    n_hubs: int = typer.Option(10, "--n-hubs", help="-e: hub genes per module."),
    traits: Optional[Path] = typer.Option(
        None, "--traits", help="-e: sample traits TSV (samples as rows, numeric columns); "
                               "needs --trait-cor."),
    trait_out: Optional[Path] = typer.Option(
        None, "--trait-cor", help="-e: write the eigengene-trait correlations "
                                  "(module, trait, cor, p, n) to this TSV file."),
    layer: Optional[str] = typer.Option(None, "--layer",
                                        help="-e .h5ad: adata.layers entry instead of X."),
    hvg: bool = typer.Option(False, "--hvg",
                             help="-e .h5ad: only genes in adata.var['highly_variable']."),
    groupby: Optional[str] = typer.Option(None, "--groupby",
                                          help="-e .h5ad: adata.obs column to select cells by."),
    group: Optional[list[str]] = typer.Option(
        None, "--group", help="-e .h5ad: keep cells whose --groupby value is this; repeatable."),
    soft_power_out: Optional[Path] = typer.Option(
        None, "--soft-power", help="-e: write the scale-free fit of every power tried to this "
                                   "TSV file (for gcen2 plot soft-power)."),
    dendrogram_out: Optional[Path] = typer.Option(
        None, "--dendrogram", help="-e: write the gene tree and modules to this JSON file (for "
                                   "gcen2 plot dendrogram)."),
    tom_hubs_out: Optional[Path] = typer.Option(
        None, "--tom-hubs", help="-e: write the TOM of the module genes with each module's 25 "
                                 "top-kME genes to this TSV file (hub networks, module UMAP)."),
    gs_out: Optional[Path] = typer.Option(
        None, "--gs", help="-e: write the gene significance (correlation of every gene with "
                           "every trait, and p) to this TSV file; needs --traits."),
) -> None:
    """Identify gene modules from a co-expression network (-i, Jaccard + DFS)
    or from an expression matrix or .h5ad (-e, WGCNA TOM + hierarchical
    clustering)."""
    if (input is None) == (expr is None):
        fail("give one of -i (network) or -e (expression matrix).")
    if expr is not None:
        from gcen2.cli.modules import wgcna_modules

        wgcna_modules(expr, output, method=method, signed=signed, power=power,
                      min_size=min_size, deep_split=deep_split, merge_cut=merge_cut, log=log,
                      threads=threads, preset=preset, me_out=me_out, kme_out=kme_out,
                      hubs_out=hubs_out, n_hubs=n_hubs, traits=traits, trait_out=trait_out,
                      layer=layer, hvg=hvg, groupby=groupby, group=group,
                      soft_power_out=soft_power_out, dendrogram_out=dendrogram_out,
                      tom_hubs_out=tom_hubs_out, gs_out=gs_out)
        return
    wgcna_only = (me_out, kme_out, hubs_out, traits, trait_out, layer, groupby, min_size,
                  deep_split, merge_cut, soft_power_out, dendrogram_out, tom_hubs_out, gs_out)
    if any(v is not None for v in wgcna_only) or hvg or group or signed or preset != "wgcna":
        fail("the WGCNA options (--signed, --min-size, --me, ...) and the .h5ad options "
             "need -e.")
    from gcen2._matrix_io import write_modules
    from gcen2.module import modules_from_file

    modules = modules_from_file(input, similarity=similarity, threads=threads)
    write_modules(output, [sorted(m) for m in modules])
    typer.echo(f"{len(modules)} modules written to {output}")


@app.command("annotate")
def cmd_annotate(
    go_obo: Optional[Path] = typer.Option(None, "-g", "--go", "--go-obo", help="GO OBO file."),
    kegg: Optional[Path] = typer.Option(
        None, "-k", "--kegg", help="K2ko file (KEGG); ignored if -g is given."),
    assoc: Path = typer.Option(..., "-a", "--assoc", help="Gene-term association file."),
    network: Optional[Path] = typer.Option(None, "-n", "--network", help="Input network file."),
    module_file: Optional[Path] = typer.Option(
        None, "-m", "--module", help="Module file (ignored if -n is given)."),
    output_dir: Path = typer.Option(..., "-o", "--output", help="Output directory."),
    pval: float = typer.Option(0.05, "-p", "--pval", help="P-value cutoff."),
    threads: int = typer.Option(2, "-t", "--thread", "--threads", help=THREADS_HELP),
) -> None:
    """GO/KEGG functional annotation of a network or of modules (GCEN v1 rules).

    Network: genes without annotation are annotated from their annotated
    neighbours and written to <gene>.go / <gene>.kegg. Modules: written to
    module_<i>.go / .kegg plus module_<i>.txt with the module genes. Only
    genes/modules with significant terms produce files. Tables end with
    p_value (one-sided hypergeometric) and fdr (Benjamini-Hochberg over every
    term annotated in the background); -p filters on p_value.
    """
    # pandas/scipy-free with the compiled extension (C++ network reader)
    from gcen2 import annotation
    from gcen2.core import extension
    from gcen2.io import read_modules

    if go_obo is None and kegg is None:
        raise typer.BadParameter("one of -g/--go or -k/--kegg is required")
    if network is not None:
        if extension() is not None:
            target = annotation.NetworkAdjacency.from_file(network)
        else:
            from gcen2.io import read_network

            target = read_network(network, genes_only=True)
    elif module_file is not None:
        target = read_modules(module_file)
    else:
        raise typer.BadParameter("one of -n/--network or -m/--module is required")
    written = annotation.annotate_to_dir(
        target, output_dir, assoc, go=go_obo, kegg=None if go_obo is not None else kegg,
        pval=pval, threads=threads)
    typer.echo(f"{written} annotation files written to {output_dir}")


@app.command("enrich")
def cmd_enrich(
    query_genes: Path = typer.Option(..., "-e", "--query", help="Query gene list."),
    background: Path = typer.Option(..., "-b", "--background", help="Background gene list."),
    go_obo: Path = typer.Option(..., "-g", "--go-obo", help="GO OBO file."),
    assoc: Path = typer.Option(..., "-a", "--assoc", help="Gene-term association file."),
    output: Path = typer.Option(..., "-o", "--output", help="Output TSV file."),
    pval: float = typer.Option(0.05, "-p", "--pval", help="P-value cutoff."),
) -> None:
    """Standalone enrichment analysis for an arbitrary gene list.

    Columns end with p_value (one-sided hypergeometric) and fdr
    (Benjamini-Hochberg over every term annotated in the background);
    -p filters on p_value.
    """
    from gcen2 import annotation
    from gcen2.io import read_gene_list

    ontology = annotation.load_go(go_obo)
    result = annotation.enrich(read_gene_list(query_genes), read_gene_list(background),
                               ontology, annotation.load_association(assoc, ontology),
                               pval=pval)
    result.results.to_csv(output, sep="\t", index=False, lineterminator="\n")
    typer.echo(f"Enrichment results written to {output}")


@app.command("rwr")
def cmd_rwr(
    network: Path = typer.Option(..., "-n", "--network", help="Input network TSV file."),
    genes: Path = typer.Option(..., "-g", "--genes", help="Seed genes list (one per line; "
                                                          "GCEN v1 seed files accepted)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output TSV file."),
    gamma: float = typer.Option(0.5, "--gamma", help="Restart probability."),
    prob_cutoff: float = typer.Option(0.01, "--prob-cutoff", help="Minimum probability."),
) -> None:
    """Gene function prediction using Random Walk with Restart."""
    # pandas/scipy-free path (the C++ extension reads the network file)
    from gcen2._matrix_io import tsv_field
    from gcen2.io import read_gene_list
    from gcen2.rwr import predict_from_file

    ranked, prob = predict_from_file(network, read_gene_list(genes), gamma=gamma,
                                     prob_cutoff=prob_cutoff)
    with open(output, "w", encoding="utf-8", newline="") as fh:  # as pandas.to_csv
        fh.write("gene\tprobability\n")
        fh.write("".join(f"{tsv_field(g)}\t{p!r}\n" for g, p in zip(ranked, prob.tolist())))
    typer.echo(f"RWR predictions written to {output}")


@app.command("merge")
def cmd_merge(
    inputs: list[Path] = typer.Option(..., "-i", "--input",
                                      help="Input network files (repeat for each network)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output merged network TSV file."),
    cor_cutoff: float = typer.Option(0.5, "-c", "--cor-cutoff",
                                     help="Minimum absolute aggregated correlation."),
) -> None:
    """Merge two or more co-expression networks."""
    from gcen2.network import merge_files

    r = merge_files(inputs, output, cor_cutoff=cor_cutoff)
    typer.echo(f"Merged network: {r['node_count']} nodes, {r['edge_count']} edges → {output}")


@app.command("extract")
def cmd_extract(
    input: Path = typer.Option(..., "-i", "--input", help="Input network TSV file."),
    genes: Path = typer.Option(..., "-g", "--genes",
                               help="Gene list file (one gene per line)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output subnetwork TSV file."),
) -> None:
    """Extract a subnetwork for a list of genes."""
    from gcen2.io import read_gene_list
    from gcen2.network import extract_file

    r = extract_file(input, read_gene_list(genes), output)
    typer.echo(f"Subnetwork: {r['node_count']} nodes, {r['edge_count']} edges → {output}")


@app.command("shuffle")
def cmd_shuffle(
    input: Path = typer.Option(..., "-i", "--input", help="Input network TSV file."),
    output: Path = typer.Option(..., "-o", "--output", help="Output shuffled network TSV file."),
    swap_factor: int = typer.Option(100, "--swap-factor",
                                    help="Number of swaps = swap_factor × edge_count."),
    seed: int = typer.Option(42, "--seed", help="Random seed."),
) -> None:
    """Degree-preserving edge shuffle (null model generation)."""
    # read, swapped (numpy's random stream reproduced in C++) and written by
    # the extension; same file as shuffle(read_network(input)).to_edgelist
    from gcen2.network import shuffle_file

    r = shuffle_file(input, output, swap_factor=swap_factor, seed=seed)
    typer.echo(f"Shuffled network: {r['node_count']} nodes, {r['edge_count']} edges "
               f"→ {output}")


@app.command("stats")
def cmd_stats(
    input: Path = typer.Option(..., "-i", "--input", help="Input network TSV file."),
    threads: int = typer.Option(2, "-t", "--threads", help=THREADS_HELP),
) -> None:
    """Print network statistics."""
    # read and counted in C++ (no pandas/scipy); same values as
    # gcen2.network.stats(read_network(input))
    from gcen2.network import stats_from_file

    s = stats_from_file(input, threads=threads)
    for label, value in (
            ("Nodes", s.node_count), ("Edges", s.edge_count), ("Density", f"{s.density:.6f}"),
            ("Average degree", f"{s.average_degree:.2f}"), ("Max degree", s.max_degree),
            ("Connected components", s.connected_components),
            ("Giant component size", s.giant_component_size),
            ("Heterogeneity", f"{s.heterogeneity:.4f}"),
            ("Clustering coefficient", f"{s.clustering_coefficient:.4f}"),
            ("Transitivity", f"{s.transitivity:.4f}")):
        typer.echo(f"{label + ':':<24}{value}")


@app.command("diff")
def cmd_diff(
    network_a: Path = typer.Option(..., "--network-a", "-A",
                                   help="Network from condition A."),
    network_b: Path = typer.Option(..., "--network-b", "-B",
                                   help="Network from condition B."),
    output: Path = typer.Option(..., "-o", "--output",
                                help="Output differential results TSV file."),
    pval: float = typer.Option(0.05, "-p", "--pval", help="P-value cutoff."),
    n_a: int = typer.Option(30, "--n-a", help="Number of samples in condition A."),
    n_b: int = typer.Option(30, "--n-b", help="Number of samples in condition B."),
    method: str = typer.Option(
        "pearson", "-m", "--method",
        help="Correlation of the networks: pearson, or spearman (variance 1.06 / (n - 3))."),
    fdr: bool = typer.Option(False, "-f", "--fdr",
                             help="Call changed edges by FDR (Benjamini-Hochberg) <= -p."),
    threads: int = typer.Option(2, "-t", "--threads",
                                help="Threads for reading and writing the files."),
) -> None:
    """Differential co-expression analysis between two conditions (Fisher's z test)."""
    # edges read, joined and written in C++; same file as
    # differential.compare(read_network(A), read_network(B)).results.to_csv
    from gcen2.differential import compare_files

    counts = compare_files(network_a, network_b, output, pval=pval, n_a=n_a, n_b=n_b,
                           fdr=fdr, correlation=method, threads=threads)
    typer.echo(f"Differential analysis complete → {output}")
    for cat, n in sorted(counts.items()):
        typer.echo(f"  {cat}: {n}")
