"""Single-cell analysis functions: pseudo-bulk aggregation and cell-type networks."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd

from gcen2._anndata import (  # noqa: F401  (is_anndata re-exported)
    adata_matrix,
    is_anndata,
    require_anndata,
)
from gcen2.exceptions import GCENInputError
from gcen2.models import Network


def pseudobulk(
    adata: Any,
    groupby: str,
    layer: str | None = None,
    method: str = "sum",
    min_cells: int = 10,
) -> pd.DataFrame:
    """Aggregate single-cell data into pseudo-bulk samples.

    Parameters
    ----------
    adata:
        AnnData object with cells as rows and genes as columns.
    groupby:
        Key in ``adata.obs`` to group cells by (e.g., cell type or cluster).
    layer:
        Which layer to use for aggregation. Uses ``adata.X`` if None.
    method:
        Aggregation method: 'sum' or 'mean'.
    min_cells:
        Minimum cells per group. Groups below this threshold are excluded
        with a warning.

    Returns
    -------
    pd.DataFrame with genes as rows and pseudo-bulk samples as columns.
    """
    require_anndata()

    if groupby not in adata.obs.columns:
        raise GCENInputError(f"groupby key '{groupby}' not found in adata.obs.")

    import scipy.sparse as sp

    X = adata.layers[layer] if layer else adata.X

    groups = adata.obs[groupby].unique()
    result_cols = {}

    for group in groups:
        mask = (adata.obs[groupby] == group).to_numpy()
        n_cells = mask.sum()
        if n_cells < min_cells:
            warnings.warn(
                f"Group '{group}' has only {n_cells} cells (< min_cells={min_cells}). Skipping.",
                stacklevel=2,
            )
            continue

        X_group = X[mask, :]
        if sp.issparse(X_group):
            X_group = X_group.toarray()

        if method == "sum":
            agg = X_group.sum(axis=0)
        elif method == "mean":
            agg = X_group.mean(axis=0)
        else:
            raise GCENInputError(f"Unknown aggregation method '{method}'. Choose 'sum' or 'mean'.")

        result_cols[str(group)] = np.asarray(agg).flatten()

    if not result_cols:
        raise GCENInputError("No groups passed the min_cells filter. Cannot build pseudo-bulk.")

    return pd.DataFrame(result_cols, index=adata.var_names)


def build_by_celltype(
    adata: Any,
    celltype_key: str,
    use_highly_variable: bool = False,
    min_cells: int = 10,
    strategy: str = "metacell",
    k: int = 25,
    max_shared: int | None = None,
    max_metacells: int = 1000,
    use_rep: str | None = None,
    n_pcs: int = 30,
    layer: str | None = None,
    normalize: bool | None = None,
    pseudo_replicates: int = 10,
    seed: int = 42,
    max_genes_per_chunk: int | None = None,
    effective: str = "overlap",
    fdr: bool = True,
    **build_kwargs,
) -> dict[str, Network]:
    """Build cell-type-specific co-expression networks.

    Parameters
    ----------
    adata:
        AnnData object.
    celltype_key:
        Key in ``adata.obs`` identifying cell types.
    use_highly_variable:
        If True, subset to highly variable genes before building networks.
    min_cells:
        Minimum cells per cell type (for metacells, at least ``k``).
    strategy:
        How cells become correlation samples:

        * ``metacell`` (default): kNN metacells within each cell type (see
          :func:`gcen2.single_cell.metacells`; ``k``, ``max_shared``,
          ``max_metacells``, ``use_rep``, ``n_pcs``), each a mean of ``k``
          similar cells, so the samples follow the cell states.
        * ``pseudo_replicates`` (deprecated, kept for compatibility with
          GCEN v1; warns): cells randomly split into ``pseudo_replicates``
          groups whose sums are the samples. Random groups have the same
          expected composition, so their sums estimate the cell-level
          correlations from a few samples: on PBMCs, 10 groups ranked gene
          pairs at chance level, and their edges were about the number
          expected by chance (``benchmarks/hdwgcna_pbmc``).
    layer:
        Expression from ``adata.layers[layer]`` instead of ``X``.
    normalize:
        Metacells only: log-normalise each metacell (see
        :func:`gcen2.single_cell.metacells`; default: for integer counts).
    seed:
        Random seed (metacell centres or the random split).
    effective:
        Metacells only: the sample size of the p-values. ``overlap``
        (default): per gene pair, from the metacells' shared cells (see
        :func:`metacell_effective_samples`); ``cells``: ``n_cells // k``
        for every pair (GCEN2 before 2026-11-18).
    fdr:
        Benjamini-Hochberg FDR over all gene pairs, the ``pval`` cutoff
        applying to it (default; ``gcen2 celltype --no-fdr`` turns it off).
        Unlike :func:`gcen2.network.build`, whose default is off: with
        calibrated p-values, ``pval`` alone keeps about ``pval`` x 30-40
        million gene pairs of a cell type as chance edges, which dominate the
        networks and modules of small cell types (GCEN2 before 2026-11-18:
        off).
    max_genes_per_chunk:
        Not used (kept for compatibility): the correlation kernel works
        through the gene pairs in blocks itself, and splitting the genes
        would have adjusted the FDR within each block instead of over all
        gene pairs.
    **build_kwargs:
        Additional keyword arguments passed to ``gcen2.network.build()``.

    Returns
    -------
    dict mapping cell type name → Network. ``metadata`` holds
    ``cell_type``, ``n_cells``, ``n_samples``, ``strategy`` (and ``k`` and
    ``p_value_samples`` for metacells).

    Metacells overlap and each averages ``k`` cells, so they are far fewer
    independent samples than their number: p-values (and FDR) of metacell
    networks are computed on an effective sample size, and the ``pval`` /
    ``fdr`` cutoffs apply to those. Counted on all metacells they would be
    far too small (on simulated data, thousands of spurious edges between
    module genes and unrelated genes). The effective size depends on the
    gene: a gene that varies smoothly over overlapping metacells has fewer
    independent values than a noisy one (:func:`metacell_effective_samples`;
    ``metadata["gene_effective_samples"]``, and ``p_value_samples`` their
    median).
    """
    from gcen2 import network as net_module

    require_anndata()
    del max_genes_per_chunk  # see the docstring
    build_kwargs["fdr"] = fdr
    if effective not in ("overlap", "cells"):
        raise GCENInputError(f"effective must be 'overlap' or 'cells', got '{effective}'.")
    if strategy not in ("metacell", "pseudo_replicates"):
        raise GCENInputError(
            f"strategy must be 'metacell' or 'pseudo_replicates', got '{strategy}'.")
    if strategy == "pseudo_replicates":
        warnings.warn(
            "strategy='pseudo_replicates' is deprecated and will be removed: sums of random "
            "groups of cells rank gene pairs at chance level; use the default metacells.",
            FutureWarning, stacklevel=2)
    if use_highly_variable and "highly_variable" not in adata.var.columns:
        raise GCENInputError(
            "use_highly_variable=True but 'highly_variable' not in adata.var. "
            "Run sc.pp.highly_variable_genes(adata) first."
        )
    if celltype_key not in adata.obs.columns:
        raise GCENInputError(f"celltype_key '{celltype_key}' not found in adata.obs.")

    samples = {}  # cell type -> (genes x samples DataFrame, n_cells)
    labels = adata.obs[celltype_key].astype(str)
    if strategy == "metacell":
        from gcen2.single_cell.metacell import metacells

        counts = labels.value_counts()
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")  # small types reported below
            meta = metacells(adata, groupby=celltype_key, k=k, max_shared=max_shared,
                             max_metacells=max_metacells, min_cells=max(min_cells, k),
                             use_rep=use_rep, n_pcs=n_pcs, layer=layer,
                             normalize=normalize, seed=seed,
                             threads=build_kwargs.get("threads", 1)) \
                if (counts >= max(min_cells, k)).any() else None
        for ct in labels.unique():
            if counts[ct] < max(min_cells, k):
                warnings.warn(f"Cell type '{ct}' has {counts[ct]} cells "
                              f"(< {max(min_cells, k)}). Skipping.", stacklevel=2)
                continue
            m, genes = adata_matrix(meta, use_highly_variable=use_highly_variable,
                                    groupby=celltype_key, groups=ct)
            if m.shape[1] < 10:
                warnings.warn(f"Cell type '{ct}': only {m.shape[1]} metacells; "
                              "correlations will be unreliable.", stacklevel=2)
            rows = np.flatnonzero((meta.obs[celltype_key].astype(str) == ct).to_numpy())
            samples[ct] = (pd.DataFrame(m, index=genes), int(counts[ct]),
                           meta.uns["metacell_members"][rows])
    else:
        import scipy.sparse as sp

        rng = np.random.default_rng(seed)
        X_all = adata.layers[layer] if layer is not None else adata.X
        gene_mask = (adata.var["highly_variable"].to_numpy(dtype=bool)
                     if use_highly_variable else np.ones(adata.n_vars, dtype=bool))
        genes = np.asarray(adata.var_names, dtype=str)[gene_mask]
        for ct in labels.unique():
            cells = np.flatnonzero((labels == ct).to_numpy())
            if len(cells) < min_cells:
                warnings.warn(f"Cell type '{ct}' has {len(cells)} cells (< {min_cells}). "
                              "Skipping.", stacklevel=2)
                continue
            cols = []
            for indices in _partition_cells(len(cells), pseudo_replicates, rng):
                X = X_all[cells[indices]][:, gene_mask]
                if sp.issparse(X):
                    X = X.toarray()
                cols.append(np.asarray(X.sum(axis=0), dtype=np.float64).ravel())
            samples[ct] = (pd.DataFrame(np.column_stack(cols), index=genes,
                                        columns=[f"rep_{i}" for i in range(len(cols))]),
                           len(cells), None)

    networks: dict[str, Network] = {}
    for ct, (expr, n_cells, members) in samples.items():
        try:
            ct_network = net_module.build(expr, **build_kwargs)
        except GCENInputError as exc:
            warnings.warn(f"Failed to build network for cell type '{ct}': {exc}",
                          stacklevel=2)
            continue
        ct_network.metadata.update(cell_type=str(ct), n_cells=n_cells,
                                   n_samples=expr.shape[1], strategy=strategy)
        if strategy == "metacell":
            ct_network.metadata["k"] = k
            n_eff = effective_samples(ct_network)
            if n_eff < 4:
                warnings.warn(f"Cell type '{ct}': {n_eff} effective samples "
                              f"({n_cells} cells / k = {k}); too few. Skipping.",
                              stacklevel=2)
                continue
            if effective == "overlap":
                n_eff = metacell_effective_samples(
                    expr, members, k, method=build_kwargs.get("method", "spearman"),
                    log=build_kwargs.get("log"))
            ct_network = _effective_pvalues(ct_network, n_eff, len(expr), build_kwargs)
            ct_network.metadata["effective"] = effective
        networks[str(ct)] = ct_network

    return networks


def build_dropout_aware(
    expr: Any,
    method: str = "pearson",
    zero_handling: str = "both_nonzero",
    **build_kwargs,
) -> Network:
    """Build a co-expression network with dropout-aware correlation.

    Parameters
    ----------
    expr:
        Expression matrix (genes × samples DataFrame) or AnnData (any input
        of :func:`gcen2.network.build`; ``layer`` / ``use_highly_variable``
        go in ``build_kwargs``).
    method:
        Correlation method; the two zero-masking modes support ``pearson``
        only (the ranks would change with every pair's cell subset), so
        log-transform the data first (e.g. ``log="log2"``).
    zero_handling:
        'both_nonzero': per gene pair, Pearson over the cells where both
        genes are non-zero, zeros being treated as missing values (technical
        dropouts); p-value on that number of cells
        (``method="pearson_nonzero"`` of :func:`gcen2.network.build`, C++).
        Suited to zero-inflated data; when zeros mostly reflect low
        expression (UMI counts), plain Pearson on metacells is better.
        'pairwise_complete': ignore only the cells where both genes are zero
        (``pearson_pairwise_complete``). Warns: on sparse data the cells
        left with one gene at zero make unrelated genes look anti-correlated.
        'impute': apply MAGIC imputation before computing correlation.
        Warns: imputed cells are not independent samples, so the p-values
        (and FDR), which count every cell, are not valid. Diffusion also
        makes unrelated genes correlated: on PBMCs, 12-27% of the pairs
        with a permuted null gene became edges (``benchmarks/hdwgcna_pbmc``).
        Use it to rank pairs at most; metacells for networks.
    """
    from gcen2 import network as net_module

    masked = {"both_nonzero": "pearson_nonzero",
              "pairwise_complete": "pearson_pairwise_complete"}
    if zero_handling in masked:
        if method not in ("pearson", masked[zero_handling]):
            raise GCENInputError(
                f"zero_handling='{zero_handling}' supports method='pearson' only "
                "(ranks would change with every gene pair's cell subset).")
        if zero_handling == "pairwise_complete":
            warnings.warn(
                "pairwise_complete keeps the cells where one gene is 0 and the other is "
                "not, which biases correlations towards negative values on sparse data; "
                "consider zero_handling='both_nonzero' or metacells.", stacklevel=2)
        net = net_module.build(expr, method=masked[zero_handling], **build_kwargs)
        net.metadata["zero_handling"] = zero_handling
        return net
    if zero_handling == "impute":
        try:
            import magic
        except ImportError as exc:
            raise GCENInputError(
                "MAGIC imputation requires the 'magic-impute' package. "
                "Install with: pip install magic-impute"
            ) from exc
        if is_anndata(expr):
            m, genes = adata_matrix(expr, layer=build_kwargs.pop("layer", None),
                                    use_highly_variable=build_kwargs.pop(
                                        "use_highly_variable", False))
            expr = pd.DataFrame(m, index=genes)
        warnings.warn(
            "zero_handling='impute': MAGIC-imputed cells are not independent samples, so the "
            "p-values and FDR are not valid, and imputation correlates unrelated genes; "
            "build networks from metacells (build_by_celltype) instead.", stacklevel=2)
        expr = magic.MAGIC().fit_transform(expr.T).T
        net = net_module.build(expr, method=method, **build_kwargs)
        net.metadata["zero_handling"] = zero_handling
        return net
    raise GCENInputError(
        f"Unknown zero_handling '{zero_handling}'. "
        "Choose 'both_nonzero', 'pairwise_complete' or 'impute'.")


def _effective_pvalues(net: Network, n_eff, n_genes: int, build_kwargs) -> Network:
    """Recompute the p-values (and FDR) of ``net`` with ``n_eff`` samples
    (a number, or a Series of every gene's effective samples: a pair gets the
    smaller of its two genes') and re-apply the cutoffs of ``build_kwargs``
    (build's defaults otherwise).

    The kernel filtered on p-values from more samples (n_eff never exceeds
    the number of metacells), which are smaller, so its edges are a superset
    of the result: every edge with an effective p-value <= pval is present
    and the BH ranks of those are exact.
    """
    import inspect

    from scipy import stats

    from gcen2 import network as net_module

    defaults = {k: v.default for k, v in inspect.signature(net_module.build).parameters.items()}
    opt = {k: build_kwargs.get(k, defaults[k])
           for k in ("cor", "pval", "signed", "fdr", "all_edges")}
    edges = net.edges
    r = edges["correlation"].to_numpy(dtype=np.float64)
    meta = dict(net.metadata)
    if isinstance(n_eff, pd.Series):
        ne = n_eff.to_numpy(dtype=np.float64)
        idx = n_eff.index
        pair = np.minimum(ne[idx.get_indexer(edges["gene_a"].astype(str))],
                          ne[idx.get_indexer(edges["gene_b"].astype(str))])
        df = pair - 2
        meta.update(p_value_samples=int(np.median(ne)), gene_effective_samples=n_eff)
    else:
        df = n_eff - 2
        meta.update(p_value_samples=n_eff)
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.abs(r) * np.sqrt(df / np.maximum(1.0 - r * r, 0.0))
    p = np.where(np.abs(r) >= 1.0, 0.0, 2.0 * stats.t.sf(t, df))
    edges = edges.assign(p_value=p)
    q = None
    if opt["fdr"]:
        n_tests = n_genes * (n_genes - 1) // 2
        q = np.full(len(p), np.nan)
        # exact for the edges with p <= pval (all of them are present)
        ok = p <= opt["pval"] if not opt["all_edges"] else np.ones(len(p), dtype=bool)
        q[ok] = net_module._benjamini_hochberg(p[ok], n_tests)
        edges = edges.assign(fdr=q)
    if not opt["all_edges"]:
        mask = net_module._final_mask(r, p, None if q is None else np.nan_to_num(q, nan=2.0),
                                      opt["cor"], opt["pval"], opt["signed"])
        edges = edges[mask].reset_index(drop=True)
    return Network(edges=edges, is_signed=net.is_signed, method=net.method, metadata=meta)


def metacell_effective_samples(expr: Any, members: Any, k: int, method: str = "pearson",
                               log: str | None = None) -> pd.Series:
    """Every gene's effective number of independent metacells.

    Metacells that share cells are not independent: a gene's metacell values
    are correlated wherever metacells overlap. If gene x is unrelated to gene
    y and to the choice of metacells, its metacell values have the
    covariance of the overlap, ``C[a, b] = |cells(a) & cells(b)| / k``, and
    the correlation of x with y has variance ``s(y) / M``, where
    ``s(y) = y'Cy / y'y`` (y centred) and M is the number of metacells. The
    effective sample size of the pair is then ``M / s(y)``. A gene that
    varies smoothly across overlapping metacells (cell-state genes) has a
    large s and few effective samples; a noisy one has s near 1. With even,
    complete coverage of the cells, a perfectly smooth gene gets
    ``n_cells / k``, the former fixed value.

    ``expr``: genes x metacells (as given to :func:`gcen2.network.build`,
    which applies ``log``; Spearman networks: ranks); ``members``: the
    member cells of every metacell (M x k, ``metacells().uns``). Returns a
    Series of ``M / s``, between 3 and M, indexed by gene; a pair is tested
    on the smaller of its two genes' values. On PBMCs these agree with the
    variance of correlations of permuted null genes (median 351 against 372
    for the half of the genes with the smallest s, 110 against 102 for the
    top 1%; ``benchmarks/hdwgcna_pbmc``).
    """
    import scipy.sparse as sp

    from gcen2 import network as net_module

    x, genes = net_module._prepare_matrix(expr, "pearson" if method.startswith("pearson")
                                          else method, log)
    if method == "spearman":
        from scipy.stats import rankdata

        x = rankdata(x, axis=1)
    members = np.asarray(members)
    n_meta = x.shape[1]
    if members.shape[0] != n_meta:
        raise GCENInputError(f"{members.shape[0]} member lists for {n_meta} metacells.")
    cells, inv = np.unique(members, return_inverse=True)
    o = sp.csr_matrix((np.ones(members.size), (np.repeat(np.arange(n_meta), members.shape[1]),
                                               inv.ravel())), shape=(n_meta, len(cells)))
    c = (o @ o.T).toarray() / k
    y = x - x.mean(axis=1, keepdims=True)
    norm2 = np.einsum("gm,gm->g", y, y)
    num = np.einsum("gm,gm->g", y, y @ c)
    s = np.divide(num, norm2, out=np.ones_like(norm2), where=norm2 > 0)
    n_eff = np.clip(n_meta / np.maximum(s, 1e-12), 3.0, n_meta)
    return pd.Series(n_eff, index=pd.Index(genes, dtype=object))


def effective_samples(network: Network) -> int:
    """One number of samples for a :func:`build_by_celltype` network, for
    the check for too few samples and Fisher's z test between cell types
    when the network has no per-gene sizes (:func:`write_celltype_networks`): ``n_cells // k`` for metacells (a gene that
    varies smoothly across evenly overlapping metacells has about that many
    independent values; the network's own p-values use every gene's value,
    :func:`metacell_effective_samples`); pseudo-replicates are disjoint and
    count fully."""
    md = network.metadata
    n = int(md.get("n_samples", 0))
    if md.get("strategy") == "metacell" and md.get("k"):
        n = min(n, int(md.get("n_cells", 0)) // int(md["k"]))
    return n


def write_celltype_networks(
    networks: dict[str, Network],
    output_dir: Any,
    diff: bool = False,
    pval: float = 0.05,
    threads: int = 2,
) -> pd.DataFrame:
    """Write :func:`build_by_celltype` networks to ``output_dir``.

    Writes ``<cell type>.network`` per cell type (names made file-safe),
    ``celltypes.tsv`` (``cell_type, file, n_cells, n_samples,
    effective_samples, node_count, edge_count``) and, with ``diff``,
    ``diff/<A>__vs__<B>.tsv`` for every pair of cell types
    (:func:`gcen2.differential.compare_files`, Fisher's z test, Spearman's
    variance when both networks are Spearman networks) plus
    ``diff/summary.tsv`` (category counts per pair; ``n_a``, ``n_b`` the
    sample sizes, for per-gene sizes their median);
    ``threads`` for reading and writing them.

    The sample sizes of the test are those of the networks' own p-values:
    for metacell networks with per-gene effective sizes
    (``metadata["gene_effective_samples"]``, ``effective="overlap"``) every
    edge is tested on the smaller of its two genes' sizes in each cell type
    (edges with fewer than 4 in either: no test, ``shared_stable``); otherwise
    :func:`effective_samples` per network, and pairs of cell types with fewer
    than 4 are skipped with a warning.
    Returns the ``celltypes.tsv`` table (``effective_samples``: the
    networks' ``p_value_samples``, the median gene's for per-gene sizes).
    """
    import re
    from itertools import combinations
    from pathlib import Path

    from gcen2.differential import CATEGORIES, compare_files

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows, files, used = [], {}, set()
    for ct, net in networks.items():
        stem = re.sub(r"[^A-Za-z0-9._-]+", "_", str(ct)).strip("._") or "celltype"
        name, i = stem, 1
        while name.lower() in used:
            i += 1
            name = f"{stem}_{i}"
        used.add(name.lower())
        path = out / f"{name}.network"
        net.to_edgelist(path)
        files[ct] = path
        md = net.metadata
        rows.append({"cell_type": ct, "file": path.name, "n_cells": md.get("n_cells"),
                     "n_samples": md.get("n_samples"),
                     "effective_samples": int(md.get("p_value_samples",
                                                     effective_samples(net))),
                     "node_count": net.node_count, "edge_count": net.edge_count})
    table = pd.DataFrame(rows, columns=["cell_type", "file", "n_cells", "n_samples",
                                        "effective_samples", "node_count", "edge_count"])
    table.to_csv(out / "celltypes.tsv", sep="\t", index=False, lineterminator="\n")
    if diff:
        ddir = out / "diff"
        ddir.mkdir(exist_ok=True)
        summary = []
        for a, b in combinations(networks, 2):
            ga = networks[a].metadata.get("gene_effective_samples")
            gb = networks[b].metadata.get("gene_effective_samples")
            if ga is not None and gb is not None:  # per edge (the networks' own p-values)
                n_a, n_b = ga, gb
                na, nb = (int(networks[x].metadata["p_value_samples"]) for x in (a, b))
            else:
                n_a = na = effective_samples(networks[a])
                n_b = nb = effective_samples(networks[b])
                if na < 4 or nb < 4:
                    warnings.warn(f"Skipping diff {a} vs {b}: fewer than 4 effective samples "
                                  f"({na}, {nb}).", stacklevel=2)
                    continue
            path = ddir / f"{files[a].stem}__vs__{files[b].stem}.tsv"
            spearman = networks[a].method == networks[b].method == "spearman"
            counts = compare_files(files[a], files[b], path, pval=pval, n_a=n_a, n_b=n_b,
                                   correlation="spearman" if spearman else "pearson",
                                   threads=threads)
            summary.append({"cell_type_a": a, "cell_type_b": b, "file": path.name,
                            "n_a": na, "n_b": nb, **{c: counts.get(c, 0) for c in CATEGORIES}})
        pd.DataFrame(summary, columns=["cell_type_a", "cell_type_b", "file", "n_a", "n_b",
                                       *CATEGORIES]).to_csv(ddir / "summary.tsv", sep="\t",
                                                            index=False, lineterminator="\n")
    return table


def _partition_cells(n_cells: int, k: int, rng: np.random.Generator) -> list[np.ndarray]:
    indices = rng.permutation(n_cells)
    return [indices[i::k] for i in range(k)]
