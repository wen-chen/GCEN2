"""WGCNA-style modules from spatial data, as hdWGCNA builds them:
metaspots (each centre spot summed with its neighbours, within groups such
as regions or clusters), log-normalized (``NormalizeMetacells``), then
:func:`gcen2.tom_modules.tom_modules` on all metaspots (by default with
hdWGCNA's ConstructNetwork settings) and module eigengenes of the original
spots (``ModuleEigengenes``).

Metaspots come from the Visium grid (hdWGCNA ``MetaspotsByGroups``) or, for
platforms without one (Slide-seq, Visium HD, Stereo-seq, imaging-based
single-cell data), from the k nearest neighbours in physical coordinates
(``method="knn"``)."""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np

from gcen2.exceptions import GCENInputError


def _grid(adata: Any, row_key: str, col_key: str) -> np.ndarray:
    for key in (row_key, col_key):
        if key not in adata.obs.columns:
            raise GCENInputError(
                f"'{key}' not found in adata.obs: metaspots need the Visium array "
                "coordinates (array_row / array_col of tissue_positions).")
    rc = adata.obs[[row_key, col_key]].to_numpy(dtype=np.float64)
    if not np.isfinite(rc).all() or not np.array_equal(rc, np.round(rc)):
        raise GCENInputError(f"{row_key} / {col_key} must be integer grid coordinates.")
    return rc.astype(np.int64)


def _group_metaspots(rc: np.ndarray):
    """hdWGCNA ConstructMetaspots on the grid coordinates (rows x 2: row,
    col) of one group's spots: (centre index, member indices) per metaspot,
    members = [centre] + neighbours (a neighbour can appear twice, as in
    hdWGCNA, when the grid has no diagonal neighbours). None if hdWGCNA's
    aggregation fails (fewer than two centre columns)."""
    row, col = rc[:, 0], rc[:, 1]
    bounds = np.arange(col.min(), col.max() + 1)
    if len(bounds) >= 4:
        bounds = bounds[3::4]
    elif len(bounds) == 3:
        bounds = bounds[1:2]
    else:
        bounds = bounds[:1]
    if len(bounds) <= 1:
        return None
    even = bool(np.all(bounds % 2 == 0))
    rows = np.arange(row.min(), row.max() + 1)
    rows = rows[(rows % 2 == 0) == even]
    # smallest non-zero distance between two spots of the group (squared;
    # integer coordinates, so the equality test below is exact)
    from scipy.spatial import cKDTree

    uniq = np.unique(rc, axis=0)
    if len(uniq) < 2:
        return None
    d, nb = cKDTree(uniq).query(uniq, k=2)
    i = int(np.argmin(d[:, 1]))
    min_d2 = int(((uniq[i] - uniq[nb[i, 1]]) ** 2).sum())
    where = {}
    for i, (r, c) in enumerate(zip(row.tolist(), col.tolist())):
        where.setdefault((r, c), []).append(i)
    ucols = np.unique(col)
    col_pos = {int(c): i for i, c in enumerate(ucols)}
    by_row: dict = {}
    for i, r in enumerate(row.tolist()):
        by_row.setdefault(r, []).append(i)
    out = []
    for x in bounds.tolist():
        for y in rows.tolist():
            for centre in where.get((y, x), ()):
                d2 = (row - y) ** 2 + (col - x) ** 2
                near = np.flatnonzero(d2 == min_d2).tolist()
                ix = col_pos[x]  # unique_cols[c(ix - 2, ix + 2)] (1-based in R)
                side = {int(ucols[k]) for k in (ix - 2, ix + 2) if 0 <= k < len(ucols)}
                other = [i for i in by_row.get(y, ()) if int(col[i]) in side]
                members = near + other
                if len(members) >= 2:
                    out.append((centre, [centre] + members))
    return out


METASPOT_METHODS = ("auto", "grid", "knn")


def _has_grid(adata: Any, row_key: str, col_key: str) -> bool:
    if row_key not in adata.obs.columns or col_key not in adata.obs.columns:
        return False
    try:
        rc = adata.obs[[row_key, col_key]].to_numpy(dtype=np.float64)
    except (TypeError, ValueError):
        return False
    if not (np.isfinite(rc).all() and np.array_equal(rc, np.round(rc))):
        return False
    # a Visium hexagonal grid: row + col has one parity (a square grid, such
    # as Visium HD bins, has both)
    return len(np.unique(rc.sum(axis=1).astype(np.int64) % 2)) == 1


def _coords(adata: Any, spatial_key: str) -> np.ndarray:
    if spatial_key not in adata.obsm:
        raise GCENInputError(
            f"'{spatial_key}' not found in adata.obsm: knn metaspots need the spot / cell "
            "coordinates (and grid metaspots need obs array_row / array_col).")
    xy = np.asarray(adata.obsm[spatial_key], dtype=np.float64)
    if xy.ndim != 2 or xy.shape[1] not in (2, 3) or len(xy) != adata.n_obs:
        raise GCENInputError(f"adata.obsm['{spatial_key}'] must be n_obs x 2 (or 3).")
    if not np.isfinite(xy).all():
        raise GCENInputError(f"adata.obsm['{spatial_key}'] has missing coordinates.")
    return xy


def _knn_metaspots(xy: np.ndarray, k: int, max_shared: int,
                   rng: np.random.Generator, threads: int = 1):
    """(centre index, member indices) per metaspot of one group: each
    centre with its k - 1 nearest spots, centres tried in random order and
    kept when they share at most ``max_shared`` spots with every metaspot
    kept so far (:func:`gcen2.single_cell.metacell._select`)."""
    from gcen2._knn import knn_indices
    from gcen2.single_cell.metacell import _select

    n = len(xy)
    # itself first, equidistant spots by index: the same on every platform
    nn = knn_indices(xy, k, workers=threads)
    chosen = _select(nn, max_shared, n, rng)
    return [(int(c), nn[c].tolist()) for c in sorted(chosen)]


def metaspots(adata: Any, groupby: str | None = None, row_key: str = "array_row",
              col_key: str = "array_col", min_spots: int = 50, layer: str | None = None,
              mode: str = "sum", method: str = "auto", spatial_key: str = "spatial",
              k: int = 7, max_shared: int | None = None, seed: int = 42, threads: int = 1):
    """Metaspots within each group of ``groupby`` (sorted; all spots if
    None) with at least ``min_spots`` spots, each centre aggregated
    (``mode`` "sum" or "average") with its neighbours. Needs counts
    (``layer`` or ``X``).

    ``method="grid"`` is hdWGCNA ``MetaspotsByGroups``: every 4th grid
    column and the rows of matching parity give the centre spots, and each
    centre is aggregated with its nearest grid neighbours and the spots two
    columns away in its row (the six neighbours of a Visium hexagon). Needs
    the Visium array coordinates in ``obs`` (``row_key``, ``col_key``).

    ``method="knn"`` needs only coordinates (``obsm[spatial_key]``, any unit)
    and so works without a grid: each centre with its ``k - 1`` nearest
    spots of the same group (``k=7``: a centre and its ring, the size of a
    Visium hexagon); centres are tried in random order (``seed``) and kept
    when they share at most ``max_shared`` spots with every metaspot kept so
    far (default ``round(2 * k / 7)``: 2 for k = 7, the overlap of
    neighbouring grid metaspots). Units much smaller than a Visium spot
    (Slide-seq beads, Visium HD bins, single cells) need a larger ``k`` for
    a comparable number of counts per metaspot: about 100 Slide-seq beads
    (10 um) cover the tissue a Visium metaspot captures.

    ``method="auto"`` uses the grid when ``obs`` has integer ``row_key`` /
    ``col_key`` columns on a hexagonal (Visium) grid, i.e. row + col of one
    parity, else knn (also for the square grid of Visium HD bins).

    Returns an AnnData (metaspots x genes; obs: the centre's obs and
    ``n_spots``; index: the centre barcode) with the member spot indices
    in ``uns["metaspot_members"]`` and the method in
    ``uns["metaspot_params"]``.
    """
    import anndata as ad
    import scipy.sparse as sp

    from gcen2.single_cell.metacell import _expression

    if mode not in ("sum", "average"):
        raise GCENInputError("mode must be 'sum' or 'average'.")
    if method not in METASPOT_METHODS:
        raise GCENInputError(f"method must be one of {', '.join(METASPOT_METHODS)}.")
    if groupby is not None and groupby not in adata.obs.columns:
        raise GCENInputError(f"groupby '{groupby}' not found in adata.obs.")
    if method == "auto":
        method = "grid" if _has_grid(adata, row_key, col_key) else "knn"
    if method == "grid":
        rc = _grid(adata, row_key, col_key)
    else:
        if k < 2:
            raise GCENInputError(f"k must be >= 2, got {k}.")
        if max_shared is None:
            max_shared = round(2 * k / 7)
        if not 0 <= max_shared < k:
            raise GCENInputError(f"max_shared must be in [0, k), got {max_shared}.")
        xy = _coords(adata, spatial_key)
        rng = np.random.default_rng(seed)
    X = _expression(adata, layer)
    labels = (np.full(adata.n_obs, "all", dtype=object) if groupby is None
              else adata.obs[groupby].astype(str).to_numpy())
    rows, cols, vals, centres, members_all = [], [], [], [], []
    for g in sorted(set(labels)):
        idx = np.flatnonzero(labels == g)
        if len(idx) < max(min_spots, k if method == "knn" else 0):
            warnings.warn(f"Group '{g}' has {len(idx)} spots (< {min_spots}). Skipping.",
                          stacklevel=2)
            continue
        res = (_group_metaspots(rc[idx]) if method == "grid"
               else _knn_metaspots(xy[idx], k, max_shared, rng, threads))
        if not res:
            warnings.warn(f"Metaspot aggregation failed for group '{g}'.", stacklevel=2)
            continue
        for centre, members in res:
            row = len(centres)
            w = 1.0 / len(members) if mode == "average" else 1.0
            rows.extend([row] * len(members))
            cols.extend(idx[members].tolist())
            vals.extend([w] * len(members))
            centres.append(idx[centre])
            members_all.append(idx[members])
    if not centres:
        raise GCENInputError("no group gave metaspots.")
    S = sp.csr_matrix((vals, (rows, cols)), shape=(len(centres), adata.n_obs))
    S.sum_duplicates()
    M = S @ (X.astype(np.float64) if sp.issparse(X) else np.asarray(X, dtype=np.float64))
    obs = adata.obs.iloc[centres].copy()
    obs["n_spots"] = [len(m) for m in members_all]
    out = ad.AnnData(X=sp.csr_matrix(M) if sp.issparse(M) else M, obs=obs,
                     var=adata.var.copy())
    out.uns["metaspot_members"] = members_all
    out.uns["metaspot_params"] = ({"method": "grid"} if method == "grid" else
                                  {"method": "knn", "k": k, "max_shared": max_shared,
                                   "seed": seed, "spatial_key": spatial_key})
    return out


def log_normalize(X, columns=None, scale: float = 1e4) -> np.ndarray:
    """Seurat LogNormalize of the rows of X (log1p of counts per ``scale``
    counts of the row, all columns), dense, of the given ``columns`` only."""
    import scipy.sparse as sp

    lib = np.asarray(X.sum(axis=1, dtype=np.float64), dtype=np.float64).ravel()
    lib[lib == 0] = 1.0
    if columns is not None:
        X = X[:, columns]
    Y = X.astype(np.float64).toarray() if sp.issparse(X) else np.array(X, dtype=np.float64)
    # x * scale / library rounds once (x * scale is exact for counts), so equal
    # ratios give equal values and the Spearman ranks of the result do not
    # depend on the platform's log1p (CPUs differ in its last bit)
    Y *= scale
    Y /= lib[:, None]
    return np.log1p(Y)


def spot_eigengenes(data: np.ndarray, labels: np.ndarray, levels) -> np.ndarray:
    """hdWGCNA ``ModuleEigengenes`` on spots: for each module the genes of
    ``data`` (spots x genes, log-normalized) are scaled (Seurat ScaleData:
    z-scores, clipped above at 10) and the eigengene is the score on the
    first principal component of that matrix (Seurat RunPCA: no
    re-centring), signed to correlate positively with the module's
    average (unscaled) expression. Returns spots x modules."""
    out = np.empty((data.shape[0], len(levels)))
    for j, k in enumerate(levels):
        out[:, j] = module_eigengene_scores(data[:, labels == k])
    return out


def module_pc_scores(x: np.ndarray, npcs: int = 50) -> np.ndarray:
    """Seurat ``RunPCA`` of one module as hdWGCNA's ``ModuleEigengenes``
    runs it: the rows of ``x`` (rows x the module's genes, log-normalized)
    scaled as :func:`module_eigengene_scores` does, then the scores (u * d)
    on the first ``min(npcs, genes - 1)`` principal components, PC1 signed
    as the eigengene. Rows x components."""
    from scipy.linalg import eigh

    x = np.asfortranarray(x, dtype=np.float64)
    sd = x.std(axis=0, ddof=1)
    sd[sd == 0] = 1.0
    z = np.minimum((x - x.mean(axis=0)) / sd, 10.0)
    g = z.shape[1]
    k = max(1, min(npcs, g - 1))
    if g <= z.shape[0]:
        _, v = eigh(z.T @ z, subset_by_index=[g - k, g - 1])
        pcs = z @ v[:, ::-1]
    else:
        gram = z @ z.T
        m = gram.shape[0]
        w, u = eigh(gram, subset_by_index=[m - k, m - 1])
        pcs = u[:, ::-1] * np.sqrt(np.maximum(w[::-1], 0.0))
    c = np.corrcoef(x.mean(axis=1), pcs[:, 0])[0, 1]
    if np.isfinite(c) and c < 0:
        pcs[:, 0] = -pcs[:, 0]
    return pcs


def module_eigengene_scores(x: np.ndarray) -> np.ndarray:
    """The eigengene of one module on the rows of ``x`` (rows x the module's
    genes, log-normalized): Seurat ScaleData (z-scores with the n - 1 SD,
    clipped above at 10), the first principal component score without
    re-centring, signed to correlate positively with the rows' average.
    Computed in column-major order whatever the input layout (the BLAS
    products round differently by layout), so the result does not depend
    on how the caller built ``x``. The first column of :func:`module_pc_scores`."""
    return module_pc_scores(x, npcs=1)[:, 0]


def spatial_modules(
    adata: Any,
    groupby: str | None = None,
    fraction: float = 0.05,
    use_highly_variable: bool = False,
    genes=None,
    min_spots: int = 50,
    mode: str = "sum",
    row_key: str = "array_row",
    col_key: str = "array_col",
    layer: str | None = None,
    threads: int = 2,
    preset: str = "hdwgcna",
    spot_level: bool = True,
    reassign: bool = True,
    method: str = "auto",
    spatial_key: str = "spatial",
    k: int = 7,
    max_shared: int | None = None,
    seed: int = 42,
    **tom_kwargs,
):
    """WGCNA-style modules of a spatial sample as hdWGCNA builds them.

    Genes: those with counts in at least ``round(fraction * n)`` of the n
    spots (:func:`gcen2.single_cell.modules.fraction_genes`; and in
    ``adata.var["highly_variable"]`` with ``use_highly_variable``), or the list ``genes``.
    Metaspots within each
    ``groupby`` group (:func:`metaspots`: the Visium grid, or the ``k``
    nearest spots in ``obsm[spatial_key]`` for data without a grid;
    ``method``, ``max_shared`` and ``seed`` as there), log-normalized; then
    :func:`gcen2.tom_modules.tom_modules` on the genes x metaspots matrix of
    all groups with the ``preset`` options (``"hdwgcna"``:
    :data:`gcen2.tom_modules.HDWGCNA`; ``"wgcna"``: blockwiseModules
    defaults) updated by ``tom_kwargs``. Eigengenes and kME of the result
    are those of the metaspots; with ``spot_level`` the eigengene of every
    module on every spot (:func:`spot_eigengenes`, hdWGCNA
    ``ModuleEigengenes``) is in ``metadata["spot_eigengenes"]`` (spots x
    ME<k> DataFrame), for spatial maps, and the Pearson correlation of every
    gene with them on the spots (hdWGCNA ``ModuleConnectivity``) in
    ``metadata["spot_kme"]`` (genes x kME<k>); with ``reassign`` (and at
    least two modules), as ``ModuleConnectivity``, genes with a negative
    spot kME for their module move to the module of their largest spot kME
    (:func:`reassign_negative_kme`; ``metadata["reassigned"]`` lists the
    genes considered). Eigengenes are not recomputed afterwards (nor are
    they in hdWGCNA).

    ``metadata`` also holds ``n_spots``, ``n_metaspots``, ``n_genes``,
    ``groups`` and ``metaspots`` (genes x metaspots DataFrame, the
    ``gcen2 module -e`` input).
    """
    import pandas as pd

    from gcen2._anndata import require_anndata
    from gcen2.single_cell.metacell import _expression, _is_counts
    from gcen2.single_cell.modules import cell_eigengenes, fraction_genes
    from gcen2.tom_modules import HDWGCNA, tom_modules

    require_anndata()
    if preset not in ("hdwgcna", "wgcna"):
        raise GCENInputError("preset must be 'hdwgcna' or 'wgcna'.")
    if not 0 <= fraction <= 1:
        raise GCENInputError("fraction must be between 0 and 1.")
    X = _expression(adata, layer)
    if not _is_counts(X):
        raise GCENInputError("spatial_modules needs raw counts (adata.X or a layer).")
    var = np.asarray(adata.var_names, dtype=str)
    if genes is not None:
        genes = list(dict.fromkeys(str(g) for g in genes))
        missing = sorted(set(genes) - set(var))
        if missing:
            raise GCENInputError(f"{len(missing)} genes not in adata.var_names "
                                 f"(e.g. {', '.join(missing[:5])}).")
        pos = {g: i for i, g in enumerate(var)}
        gene_idx = np.array([pos[g] for g in genes], dtype=np.int64)
    else:
        keep = fraction_genes(X, fraction)
        if use_highly_variable:
            if "highly_variable" not in adata.var.columns:
                raise GCENInputError("use_highly_variable=True but 'highly_variable' not "
                                     "in adata.var.")
            keep &= adata.var["highly_variable"].to_numpy(dtype=bool)
        gene_idx = np.flatnonzero(keep)
    if len(gene_idx) < 3:
        raise GCENInputError("fewer than 3 genes selected.")
    opts = dict(HDWGCNA) if preset == "hdwgcna" else {}
    opts.update(tom_kwargs)

    ms = metaspots(adata, groupby=groupby, row_key=row_key, col_key=col_key,
                   min_spots=min_spots, layer=layer, mode=mode, method=method,
                   spatial_key=spatial_key, k=k, max_shared=max_shared, seed=seed,
                   threads=threads)
    data = log_normalize(ms.X, gene_idx)
    names = [str(s) for s in ms.obs_names]
    gene_names = [str(var[i]) for i in gene_idx]
    m = np.ascontiguousarray(data.T)
    res = tom_modules(m, gene_names=gene_names, sample_names=names, method="pearson",
                      threads=threads, **opts)
    md = res.metadata
    md.update(n_spots=int(adata.n_obs), n_metaspots=len(names), n_genes=len(gene_names),
              metaspot_params=dict(ms.uns["metaspot_params"]),
              groups=(sorted(set(ms.obs[groupby].astype(str))) if groupby else ["all"]),
              metaspots=pd.DataFrame(m, index=gene_names, columns=names))
    if spot_level:
        levels = [k for k in res.me_labels if k > 0]
        good = np.flatnonzero(res._good)
        # in row blocks on the sparse counts (the single-cell code): no dense
        # spots x genes matrix, which for Visium HD bins is several GB
        me, kme = cell_eigengenes(X, gene_idx[good], res.labels[good], levels)
        md["spot_eigengenes"] = pd.DataFrame(
            me, index=[str(s) for s in adata.obs_names], columns=[f"ME{k}" for k in levels])
        good_genes = [res.genes[i] for i in good]
        md["spot_kme"] = pd.DataFrame(kme, index=good_genes,
                                      columns=[f"kME{k}" for k in levels])
        md["reassigned"] = []
        if reassign and len(levels) >= 2:
            lab = res.labels[good]
            new, moved = reassign_negative_kme(lab, kme, levels)
            res.labels[good] = new
            md["reassigned"] = [good_genes[i] for i in moved]
    return res


def reassign_negative_kme(labels: np.ndarray, kme: np.ndarray, levels):
    """hdWGCNA ``ReassignModules`` (run by ``ModuleConnectivity``): a gene
    of module k (> 0) whose kME for k is negative moves to the module with
    its largest kME (the first on ties), or to 0 (unassigned) if all its
    kMEs are negative. ``kme``: genes x ``levels``. Returns the new labels
    and the indices of the genes considered (kME < 0 for their module)."""
    col = {k: j for j, k in enumerate(levels)}
    new = labels.copy()
    moved = []
    for i, k in enumerate(labels):
        if k > 0 and kme[i, col[int(k)]] < 0:
            moved.append(i)
            j = int(np.argmax(kme[i]))
            new[i] = levels[j] if kme[i, j] >= 0 else 0
    return new, moved


def group_trait_cor(spot_me, groups):
    """Correlation of every spot-level eigengene with each group indicator
    (spot in the group: 1, else 0), as :meth:`TomModules.trait_cor`: a long
    table module, trait, cor, p, n."""
    import pandas as pd

    from gcen2.tom_modules import cor_pairwise

    g = pd.Series(np.asarray(groups, dtype=str), index=spot_me.index)
    levels = sorted(g.unique())
    ind = np.array([(g == lv).to_numpy(dtype=np.float64) for lv in levels])
    r, p, n = cor_pairwise(spot_me.to_numpy(dtype=np.float64).T, ind)
    k, m = r.shape
    return pd.DataFrame({"module": np.repeat(list(spot_me.columns), m),
                         "trait": np.tile(levels, k), "cor": r.ravel(), "p": p.ravel(),
                         "n": n.ravel()})


def write_spatial_modules(res, output_dir: Any, groups=None, n_hubs: int = 10,
                          write_metaspots: bool = False, write_tom: bool = False,
                          threads: int = 2):
    """Write :func:`spatial_modules` results to ``output_dir``: as ``gcen2
    module -e`` (:func:`gcen2.tom_modules.write_module_files`)
    ``modules.tsv``, ``eigengenes.tsv`` (metaspots), ``kme.tsv`` and
    ``hub_genes.tsv``; ``spot_eigengenes.tsv`` (spots x ME0..ME<K-1>) and
    ``spot_kme.tsv`` (gene, module, kME0..: hdWGCNA's kME), if computed; with ``groups`` (one
    label per spot, e.g. the ``groupby``
    column) ``group_cor.tsv``, their correlations (:func:`group_trait_cor`,
    the ``--trait-cor`` format of ``gcen2 report``); with
    ``write_metaspots`` ``metaspots.tsv`` (genes x metaspots); and
    ``summary.tsv``. Returns the summary as a dict."""
    from pathlib import Path

    import pandas as pd

    from gcen2.tom_modules import write_module_files

    d = Path(output_dir)
    d.mkdir(parents=True, exist_ok=True)
    write_module_files(res, modules=d / "modules.tsv", eigengenes=d / "eigengenes.tsv",
                       kme=d / "kme.tsv", hubs=d / "hub_genes.tsv", n_hubs=n_hubs,
                       soft_power=d / "soft_power.tsv", dendrogram=d / "dendrogram.json",
                       tom_hubs=d / "tom_hubs.tsv" if write_tom else None, threads=threads)
    md = res.metadata
    me = md.get("spot_eigengenes")
    if me is not None:
        me = me.rename(columns={c: f"ME{int(c[2:]) - 1}" for c in me.columns})
        me.rename_axis("spot").to_csv(d / "spot_eigengenes.tsv", sep="\t",
                                      float_format="%.10g", lineterminator="\n")
        sk = md.get("spot_kme")
        if sk is not None:
            lab = dict(zip(res.genes, res.labels))
            sk = sk.rename(columns={c: f"kME{int(c[3:]) - 1}" for c in sk.columns})
            sk.insert(0, "module", [str(lab[g] - 1) if lab[g] else "NA" for g in sk.index])
            sk.rename_axis("gene").to_csv(d / "spot_kme.tsv", sep="\t", float_format="%.10g",
                                          na_rep="NA", lineterminator="\n")
        if groups is not None and me.shape[1]:
            group_trait_cor(me, groups).to_csv(d / "group_cor.tsv", sep="\t", index=False,
                                               float_format="%.10g", na_rep="NA", lineterminator="\n")
    if write_metaspots:
        md["metaspots"].to_csv(d / "metaspots.tsv", sep="\t", float_format="%.10g", lineterminator="\n")
    summary = {"n_spots": md.get("n_spots"), "n_metaspots": md.get("n_metaspots"),
               "n_genes": md.get("n_genes"), "power": res.power,
               "modules": int(res.labels.max(initial=0)),
               "genes_in_modules": int((res.labels > 0).sum()),
               "reassigned": len(md.get("reassigned", [])),
               "groups": ",".join(md.get("groups", [])),
               "metaspots": md.get("metaspot_params", {}).get("method", "grid")}
    pd.DataFrame([summary]).to_csv(d / "summary.tsv", sep="\t", index=False, lineterminator="\n")
    return summary
