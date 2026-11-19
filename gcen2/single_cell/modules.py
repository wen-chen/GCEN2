"""WGCNA-style modules per cell type from single-cell data (the hdWGCNA
workflow): genes expressed in enough of the type's cells, kNN metacells
within the type, then :func:`gcen2.tom_modules.tom_modules` on the metacells
(by default with hdWGCNA's ConstructNetwork settings, ``HDWGCNA``)."""

from __future__ import annotations

import re
import warnings
from pathlib import Path
from typing import Any

import numpy as np

from gcen2.exceptions import GCENInputError


def expressed_fraction(adata: Any, cells: np.ndarray, layer: str | None = None) -> np.ndarray:
    """Fraction of the given cells (boolean mask) with a value > 0, per gene
    (hdWGCNA ``SetupForWGCNA(gene_select = "fraction")`` per group)."""
    X = adata.layers[layer] if layer is not None else adata.X
    X = X[cells]
    n = X.shape[0]
    if hasattr(X, "toarray"):  # scipy sparse
        pos = np.asarray((X > 0).sum(axis=0)).ravel()
    else:
        pos = (np.asarray(X) > 0).sum(axis=0)
    return pos / max(n, 1)


def fraction_genes(X, fraction: float) -> np.ndarray:
    """Genes (columns) with a value > 0 in at least ``round(fraction * n)``
    of the n rows of X: hdWGCNA ``SetupForWGCNA(gene_select = "fraction")``
    (R's round, half to even)."""
    from gcen2.tom_modules import _r_round

    return _positive_counts(X) >= _r_round(fraction * X.shape[0])


def _positive_counts(X) -> np.ndarray:
    """Per column, the number of rows with a value > 0. Sparse input is
    counted in blocks of its stored values (no copy of the matrix)."""
    import scipy.sparse as sp

    if not sp.issparse(X):
        return (np.asarray(X) > 0).sum(axis=0)
    if X.format != "csr":
        X = sp.csr_matrix(X)
    if not X.has_canonical_format:
        X = X.copy()
        X.sum_duplicates()
    pos = np.zeros(X.shape[1], dtype=np.int64)
    step = 1 << 24
    for a in range(0, X.nnz, step):
        idx = X.indices[a:a + step]
        pos += np.bincount(idx[X.data[a:a + step] > 0], minlength=X.shape[1])
    return pos


def modules_by_celltype(*args, **kwargs) -> dict:
    """WGCNA-style modules of each cell type, as hdWGCNA builds them: a dict
    cell type -> :class:`gcen2.tom_modules.TomModules` of
    :func:`iter_modules_by_celltype` (same arguments). On an atlas, iterate
    over :func:`iter_modules_by_celltype` instead and keep only what is
    needed: every result holds its genes x metacells matrix."""
    return dict(iter_modules_by_celltype(*args, **kwargs))


def iter_modules_by_celltype(
    adata: Any,
    celltype_key: str,
    groups=None,
    fraction: float = 0.05,
    use_highly_variable: bool = False,
    k: int = 25,
    max_shared: int | None = None,
    max_metacells: int = 1000,
    min_cells: int = 100,
    use_rep: str | None = None,
    n_pcs: int = 30,
    layer: str | None = None,
    normalize: bool | None = None,
    seed: int = 42,
    threads: int = 2,
    preset: str = "hdwgcna",
    cell_level: bool = True,
    reassign: bool = True,
    consensus_by: str | None = None,
    harmonize_by: str | None = None,
    npcs: int = 50,
    harmony_kwargs: dict | None = None,
    dme_by: str | None = None,
    traits=None,
    trait_groupby: str | None = None,
    **tom_kwargs,
):
    """WGCNA-style modules of each cell type, as hdWGCNA builds them; yields
    ``(cell type, TomModules)`` one cell type at a time.

    For each cell type (``groups``: only these) with at least
    ``max(min_cells, k)`` cells: the genes with a value > 0 in at least
    ``round(fraction * n)`` of its n cells (:func:`fraction_genes`; and in
    ``adata.var["highly_variable"]`` with ``use_highly_variable``), metacells of its cells
    (:func:`gcen2.single_cell.metacells`: ``k``, ``max_shared``,
    ``max_metacells``, ``use_rep``, ``n_pcs``, ``layer``, ``normalize``,
    ``seed``; computed on all genes, so the neighbour space does not depend
    on the gene selection), and :func:`gcen2.tom_modules.tom_modules` on the
    genes x metacells matrix with the ``preset`` options (``"hdwgcna"``:
    :data:`gcen2.tom_modules.HDWGCNA`; ``"wgcna"``: blockwiseModules
    defaults) updated by ``tom_kwargs`` (``power``, ``min_size``, ...).

    Metacell members are chosen for all cell types first (cheap; one random
    generator, so the result does not depend on how the output is consumed),
    their profiles only when the cell type is reached, so memory holds one
    cell type's metacells at a time.

    Each :class:`gcen2.tom_modules.TomModules` has a ``metadata`` with ``cell_type``,
    ``n_cells``, ``n_metacells``,
    ``n_genes`` and the metacell matrix (``"metacells"``, genes x metacells
    DataFrame). The result's own eigengenes and kME are those of the
    metacells. With ``cell_level`` (default), as hdWGCNA's
    ``ModuleEigengenes`` and ``ModuleConnectivity``, the eigengenes and kME
    on the type's single cells (:func:`cell_eigengenes`) are in
    ``metadata["cell_eigengenes"]`` (cells x ME<k>) and
    ``metadata["cell_kme"]`` (genes x kME<k>); with ``reassign`` (default;
    at least two modules) genes with a negative cell kME for their module
    move to the module of their largest one (hdWGCNA ``ReassignModules``;
    ``metadata["reassigned"]`` lists them). Eigengenes are not recomputed
    afterwards (nor are they in hdWGCNA).

    ``consensus_by`` (an ``adata.obs`` column, e.g. sample or study): hdWGCNA's
    consensus network (``MetacellsByGroups(group.by = c(type, consensus_by))``,
    ``SetMultiExpr``, ``TestSoftPowersConsensus``, ``ConstructNetwork(consensus
    = TRUE)``): metacells within each cell type and value of the column, and
    :func:`gcen2.consensus.consensus_tom_modules` over the values (the sets,
    in sorted order, the first being the calibration reference; values with
    no metacell of the type left out) with the ``preset`` options (``"hdwgcna"``:
    :data:`gcen2.consensus.HDWGCNA_CONSENSUS`). ``metadata`` then also has
    ``sets``, ``set_powers`` and ``set_eigengenes`` (the metacells of each set).

    ``harmonize_by`` (an ``adata.obs`` column; needs ``cell_level``):
    hdWGCNA's harmonized eigengenes (``ModuleEigengenes(group.by.vars =
    ...)``; :func:`cell_eigengenes` with ``npcs``, ``seed`` and
    ``harmony_kwargs``) in ``metadata["cell_harmonized_eigengenes"]``; the
    cell kME and the reassignment then use them, as ``ModuleConnectivity``
    does.

    ``dme_by`` (an ``adata.obs`` column; needs ``cell_level``): hdWGCNA's
    ``FindAllDMEs`` within each type (each value against the type's other
    cells; :func:`gcen2.single_cell.find_all_dmes`) in ``metadata["dmes"]``.
    ``traits`` (``adata.obs`` columns; needs ``cell_level``): hdWGCNA's
    ``ModuleTraitCorrelation`` on the type's cells, for all of them and per
    value of ``trait_groupby`` (:func:`gcen2.single_cell.module_trait_cor`)
    in ``metadata["trait_cor"]``. Both use the harmonized eigengenes when
    there are any, as hdWGCNA does.
    Cell types that are too small or yield no TOM are skipped with a warning.
    """
    import pandas as pd

    from gcen2._anndata import require_anndata
    from gcen2.consensus import HDWGCNA_CONSENSUS, consensus_tom_modules
    from gcen2.single_cell.metacell import (
        _is_counts,
        _metacell_plan,
        _metacell_profiles,
        _rows,
    )
    from gcen2.tom_modules import HDWGCNA, tom_modules

    require_anndata()
    traits = [traits] if isinstance(traits, str) else list(traits or [])
    _check_celltype_options(adata, celltype_key, preset, fraction, use_highly_variable,
                            cell_level, consensus_by, harmonize_by, dme_by, traits,
                            trait_groupby)
    presets = HDWGCNA if consensus_by is None else HDWGCNA_CONSENSUS
    opts = dict(presets) if preset == "hdwgcna" else {}
    opts.update(tom_kwargs)

    labels = adata.obs[celltype_key].astype(str)
    need = max(min_cells, k)
    keep = _cell_types(labels, groups, need, celltype_key)
    if not keep:
        return
    # no AnnData view here: every .X of a view copies the whole matrix
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        plan = _metacell_plan(adata,
                              celltype_key if consensus_by is None
                              else [celltype_key, consensus_by],
                              k, max_shared, max_metacells, need,
                              use_rep, n_pcs, layer, seed, threads, keep)
    X = plan["X"]
    if normalize is None:
        normalize = _is_counts(X)
    meta_type = np.array([g[0] for g in plan["groups"]], dtype=object)
    if consensus_by is not None:
        meta_set = np.array([g[1] for g in plan["groups"]], dtype=object)
        # sorted (C locale): the order in which hdWGCNA's ConstructNetwork takes
        # the powers of TestSoftPowersConsensus (dplyr::group_split)
        set_order = sorted(str(v) for v in pd.unique(adata.obs[consensus_by].astype(str)))
    batch_all = (None if harmonize_by is None
                 else adata.obs[harmonize_by].astype(str).to_numpy())
    all_labels = labels.to_numpy()
    var = np.asarray(adata.var_names, dtype=str)
    cell_names = np.asarray(adata.obs_names, dtype=str)
    hv = (adata.var["highly_variable"].to_numpy(dtype=bool) if use_highly_variable
          else np.ones(adata.n_vars, dtype=bool))
    for t in keep:
        cells = np.flatnonzero(all_labels == t)
        # the type's rows, read once (from disk for a backed AnnData) and used
        # for the gene selection, its metacells and the cell level
        Xc = _rows(X, cells)
        gene_mask = fraction_genes(Xc, fraction) & hv
        rows = np.flatnonzero(meta_type == t)
        if not rows.size:
            continue
        local = np.searchsorted(cells, plan["members"][rows])  # members are of type t
        prof = _metacell_profiles(Xc, local, "mean", normalize)
        names = np.array([f"metacell_{i}" for i in rows])
        m = np.ascontiguousarray(prof[:, gene_mask].T)
        del prof
        genes = var[gene_mask].tolist()
        try:
            if consensus_by is None:
                res = tom_modules(m, gene_names=genes, sample_names=names, method="pearson",
                                  threads=threads, **opts)
            else:
                sets = [v for v in set_order if (meta_set[rows] == v).any()]
                cols = [np.flatnonzero(meta_set[rows] == v) for v in sets]
                res = consensus_tom_modules([m[:, c] for c in cols], gene_names=genes,
                                            set_names=sets,
                                            sample_names=[names[c].tolist() for c in cols],
                                            threads=threads, **opts)
                order = np.concatenate(cols)  # metacells as in res.samples
                m, names = m[:, order], names[order]
        except GCENInputError as exc:
            warnings.warn(f"Cell type '{t}': {exc} Skipping.", stacklevel=2)
            continue
        members = _members_table(plan["members"][rows], rows, cell_names)
        res.metadata.update(cell_type=t, n_cells=int(cells.size), n_metacells=m.shape[1],
                            n_genes=len(genes),
                            metacells=pd.DataFrame(m, index=genes, columns=names, copy=False),
                            metacell_members=members)
        del m
        res.metadata["reassigned"] = []
        if cell_level:
            _cell_level(res, Xc, np.flatnonzero(gene_mask), normalize, reassign,
                        cell_names[cells],
                        batches=None if batch_all is None else batch_all[cells],
                        npcs=npcs, seed=seed, harmony_kwargs=harmony_kwargs)
            _cell_stats(res, adata.obs.iloc[cells], dme_by, traits, trait_groupby)
        del Xc
        yield t, res


def _check_celltype_options(adata, celltype_key, preset, fraction, use_highly_variable,
                            cell_level, consensus_by, harmonize_by, dme_by, traits,
                            trait_groupby) -> None:
    """The argument checks of :func:`iter_modules_by_celltype`."""
    for key in (consensus_by, harmonize_by, dme_by, trait_groupby, *traits):
        if key is not None and key not in adata.obs.columns:
            raise GCENInputError(f"'{key}' not found in adata.obs.")
    if (dme_by is not None or traits) and not cell_level:
        raise GCENInputError("dme_by and traits need cell_level=True.")
    if trait_groupby is not None and not traits:
        raise GCENInputError("trait_groupby needs traits.")
    if consensus_by == celltype_key:
        raise GCENInputError("consensus_by must differ from celltype_key.")
    if harmonize_by is not None and not cell_level:
        raise GCENInputError("harmonize_by needs cell_level=True.")
    if celltype_key not in adata.obs.columns:
        raise GCENInputError(f"celltype_key '{celltype_key}' not found in adata.obs.")
    if preset not in ("hdwgcna", "wgcna"):
        raise GCENInputError("preset must be 'hdwgcna' or 'wgcna'.")
    if not 0 <= fraction <= 1:
        raise GCENInputError("fraction must be between 0 and 1.")
    if use_highly_variable and "highly_variable" not in adata.var.columns:
        raise GCENInputError("use_highly_variable=True but 'highly_variable' not in "
                             "adata.var. Run sc.pp.highly_variable_genes(adata) first.")


def _cell_types(labels, groups, need: int, celltype_key: str) -> list:
    """The cell types to build (``groups``, else all, in order of first
    appearance) with at least ``need`` cells; the others are skipped with a
    warning."""
    import pandas as pd

    if groups is not None:
        groups = [str(g) for g in ([groups] if isinstance(groups, str) else groups)]
        missing = sorted(set(groups) - set(labels))
        if missing:
            raise GCENInputError(f"no cells with {celltype_key} = {', '.join(missing)}.")
    counts = labels.value_counts()
    keep = []
    for t in pd.unique(labels):
        if groups is not None and t not in groups:
            continue
        if counts[t] < need:
            warnings.warn(f"Cell type '{t}' has {counts[t]} cells (< {need}). Skipping.",
                          stacklevel=3)
        else:
            keep.append(t)
    return keep


def _members_table(mem: np.ndarray, rows: np.ndarray, cell_names: np.ndarray):
    """metacell, cell, center (the first member) of the metacells ``rows``."""
    import pandas as pd

    return pd.DataFrame({
        "metacell": np.repeat(np.array([f"metacell_{i}" for i in rows]), mem.shape[1]),
        "cell": cell_names[mem.ravel()],
        "center": np.tile(np.arange(mem.shape[1]) == 0, len(rows))})


def _cell_stats(res, obs, dme_by, traits, trait_groupby) -> None:
    """FindAllDMEs / ModuleTraitCorrelation of one type's cell eigengenes."""
    from gcen2.single_cell.dme import find_all_dmes, module_trait_cor

    md = res.metadata
    me = md.get("cell_harmonized_eigengenes", md.get("cell_eigengenes"))
    if me is None or me.shape[1] == 0:
        return
    if dme_by is not None:
        md["dmes"] = find_all_dmes(me, obs[dme_by].astype(str).to_numpy())
    if traits:
        t = obs[traits].copy()
        t.index = me.index
        g = None if trait_groupby is None else obs[trait_groupby].astype(str).to_numpy()
        md["trait_cor"] = module_trait_cor(me, t, groups=g)


def _cell_level(res, Xc, gene_idx, normalize: bool, reassign: bool, cell_names,
                batches=None, npcs: int = 50, seed: int = 0, harmony_kwargs=None) -> None:
    """hdWGCNA ModuleEigengenes / ModuleConnectivity on the cells of one type:
    ``metadata["cell_eigengenes"]`` (cells x ME<k>), ``["cell_kme"]``
    (genes x kME<k>), with ``batches`` ``["cell_harmonized_eigengenes"]``
    (the kME then against them) and, with ``reassign``, hdWGCNA's
    ReassignModules of genes with a negative kME for their module (labels
    changed in place)."""
    import pandas as pd

    from gcen2.spatial.modules import reassign_negative_kme

    levels = [k for k in res.me_labels if k > 0]
    good = np.flatnonzero(res._good)
    if not levels or not good.size:
        return
    out = cell_eigengenes(Xc, gene_idx[good], res.labels[good], levels,
                          normalize=normalize, batches=batches, npcs=npcs, seed=seed,
                          harmony_kwargs=harmony_kwargs)
    me, kme = out[0], out[1]
    md = res.metadata
    md["cell_eigengenes"] = pd.DataFrame(me, index=cell_names,
                                         columns=[f"ME{k}" for k in levels])
    if batches is not None:
        md["cell_harmonized_eigengenes"] = pd.DataFrame(out[2], index=cell_names,
                                                        columns=[f"ME{k}" for k in levels])
    good_genes = [res.genes[i] for i in good]
    md["cell_kme"] = pd.DataFrame(kme, index=good_genes, columns=[f"kME{k}" for k in levels])
    if reassign and len(levels) >= 2:
        new, moved = reassign_negative_kme(res.labels[good], kme, levels)
        res.labels[good] = new
        md["reassigned"] = [good_genes[i] for i in moved]


def cell_eigengenes(X, gene_idx, labels, levels, normalize: bool = True,
                    scale: float = 1e4, batches=None, npcs: int = 50, seed: int = 0,
                    harmony_kwargs: dict | None = None):
    """Module eigengenes and kME on single cells, as hdWGCNA's
    ``ModuleEigengenes`` and ``ModuleConnectivity`` compute them.

    ``X``: cells x all genes (counts, sparse or dense); ``gene_idx``: the
    columns of the module genes; ``labels``: their module (0 = unassigned);
    ``levels``: the modules (> 0). With ``normalize``, Seurat LogNormalize
    first (library size over all genes); else ``X`` is used as it is.
    Eigengenes: per module, ScaleData (clipped at 10) and the first
    principal component score, signed by the average expression
    (:func:`gcen2.spatial.module_eigengene_scores`, the spatial code). kME:
    Pearson correlation of every gene in ``gene_idx`` with every eigengene,
    computed on the sparse matrix (as hdWGCNA's ``corSparse``), so a large
    cell type needs no dense cells x genes matrix. Returns (cells x levels,
    genes x levels).

    With ``batches`` (one label per cell, or a list of such vectors), also
    hdWGCNA's harmonized eigengenes (``ModuleEigengenes(group.by.vars =
    ...)``): per module, Harmony (:func:`gcen2.harmony.harmony_integrate`,
    ``seed`` and ``harmony_kwargs``) on the module's first
    ``min(npcs, genes - 1)`` principal components (Seurat ``RunPCA``), and the
    first corrected component, signed as the eigengene. The kME is then
    computed against the harmonized eigengenes, as hdWGCNA's
    ``ModuleConnectivity`` does when they exist. Returns (eigengenes, kME,
    harmonized eigengenes)."""
    import scipy.sparse as sp

    from gcen2.spatial.modules import module_eigengene_scores, module_pc_scores

    if batches is not None:
        from gcen2.harmony import harmony_integrate

    X = sp.csr_matrix(X)
    gene_idx = np.asarray(gene_idx)
    labels = np.asarray(labels)
    n = X.shape[0]
    if normalize:
        # x * scale / library, one rounding: as spatial.modules.log_normalize
        lib = np.asarray(X.sum(axis=1, dtype=np.float64), dtype=np.float64).ravel()
        lib[lib == 0] = 1.0
        f = lib
    else:
        f = None
    # blocks of rows with about 1e7 stored values: memory holds one block (and
    # one module's dense cells x genes block), not a copy of X
    per_row = max(1.0, X.nnz / max(n, 1))
    step = max(256, int(1e7 / per_row))
    blocks = [(r, min(r + step, n)) for r in range(0, n, step)]

    def dense(r0, r1, cols):
        """Rows r0:r1, columns cols: x * scale / library, log1p (Seurat
        LogNormalize), float64."""
        d = X[r0:r1][:, cols].toarray().astype(np.float64)
        if f is not None:
            d *= scale
            d /= f[r0:r1, None]
            np.log1p(d, out=d)
        return d

    me = np.empty((n, len(levels)))
    hme = None if batches is None else np.empty((n, len(levels)))
    for j, k in enumerate(levels):
        cols = gene_idx[labels == k]
        x = np.empty((n, cols.size), order="F")  # the layout of module_eigengene_scores
        for r0, r1 in blocks:
            x[r0:r1] = dense(r0, r1, cols)
        me[:, j] = module_eigengene_scores(x)
        if batches is not None:
            pcs = module_pc_scores(x, npcs)
            if pcs.shape[1] < 2:  # RunHarmony needs two dimensions
                raise GCENInputError(f"module {k} has too few genes for Harmony.")
            hme[:, j] = harmony_integrate(pcs, batches, seed=seed,
                                          **(harmony_kwargs or {}))[:, 0]
        del x
    # kME: Pearson of every gene with every eigengene, accumulated over the
    # row blocks: sum x (y - mean y), sum x, sum x^2
    yc = (me if hme is None else hme)
    yc = yc - yc.mean(axis=0)
    p = gene_idx.size
    cov = np.zeros((p, len(levels)))
    s1 = np.zeros(p)
    s2 = np.zeros(p)
    for r0, r1 in blocks:
        sub = X[r0:r1][:, gene_idx]
        d = sub.data.astype(np.float64)
        if f is not None:
            d *= scale
            d /= np.repeat(f[r0:r1], np.diff(sub.indptr))
            np.log1p(d, out=d)
        blk = sp.csr_matrix((d, sub.indices, sub.indptr), shape=sub.shape)
        cov += np.asarray(blk.T @ yc[r0:r1])
        s1 += np.bincount(sub.indices, weights=d, minlength=p)
        s2 += np.bincount(sub.indices, weights=d * d, minlength=p)
    ss = s2 - s1 * s1 / n
    with np.errstate(divide="ignore", invalid="ignore"):
        r = cov / (np.sqrt(np.maximum(ss, 0.0))[:, None] * np.linalg.norm(yc, axis=0)[None, :])
    r[~np.isfinite(r)] = np.nan
    if hme is not None:
        return me, np.clip(r, -1.0, 1.0), hme
    return me, np.clip(r, -1.0, 1.0)


def _write_cell_level(res, d: Path) -> None:
    """cell_eigengenes.tsv (cells x ME0..ME<K-1>, the module numbers of
    modules.tsv) and cell_kme.tsv (gene, module, kME0..), if computed."""
    md = res.metadata
    me = md.get("cell_eigengenes")
    if me is None:
        return
    me = me.rename(columns={c: f"ME{int(c[2:]) - 1}" for c in me.columns})
    me.rename_axis("cell").to_csv(d / "cell_eigengenes.tsv", sep="\t", float_format="%.10g", lineterminator="\n")
    hme = md.get("cell_harmonized_eigengenes")
    if hme is not None:
        hme = hme.rename(columns={c: f"ME{int(c[2:]) - 1}" for c in hme.columns})
        hme.rename_axis("cell").to_csv(d / "cell_harmonized_eigengenes.tsv", sep="\t",
                                       float_format="%.10g", lineterminator="\n")
    kme = md["cell_kme"]
    lab = dict(zip(res.genes, res.labels))
    kme = kme.rename(columns={c: f"kME{int(c[3:]) - 1}" for c in kme.columns})
    kme.insert(0, "module", [str(lab[g] - 1) if lab[g] else "NA" for g in kme.index])
    kme.rename_axis("gene").to_csv(d / "cell_kme.tsv", sep="\t", float_format="%.10g",
                                   na_rep="NA", lineterminator="\n")


def _write_cell_stats(res, d: Path) -> None:
    """dmes.tsv and trait_cor.tsv (modules numbered as modules.tsv), if computed."""
    md = res.metadata

    def rename(name: str) -> str:
        return f"ME{int(name[2:]) - 1}" if name.startswith("ME") and name[2:].isdigit() else name

    if md.get("dmes") is not None:
        df = md["dmes"].assign(module=lambda x: x["module"].map(rename))
        df[["group", "module", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj"]].to_csv(
            d / "dmes.tsv", sep="\t", index=False, float_format="%.10g", na_rep="NA", lineterminator="\n")
    if md.get("trait_cor") is not None:
        df = md["trait_cor"].assign(module=lambda x: x["module"].map(rename))
        df.to_csv(d / "trait_cor.tsv", sep="\t", index=False, float_format="%.10g",
                  na_rep="NA", lineterminator="\n")


def _write_sets(res, d: Path) -> None:
    """sets.tsv (set, samples, power, calibration quantile and power) and
    set_eigengenes.tsv (set, metacell, ME0.. of each set's metacells) of a
    consensus network."""
    import pandas as pd

    md = res.metadata
    pd.DataFrame({
        "set": md["sets"],
        "samples": [len(md["set_samples"][s]) for s in md["sets"]],
        "power": [md["set_powers"][s] for s in md["sets"]],
        "calibration_quantile": [md["calibration_quantiles"][s] for s in md["sets"]],
        "calibration_power": [md["calibration_powers"][s] for s in md["sets"]],
    }).to_csv(d / "sets.tsv", sep="\t", index=False, float_format="%.10g", lineterminator="\n")
    frames = []
    for s in md["sets"]:
        e = md["set_eigengenes"][s]
        e = e.rename(columns={c: f"ME{int(c[2:]) - 1}" if c != "ME0" else "ME_unassigned"
                              for c in e.columns})
        frames.append(e.rename_axis("sample").reset_index().assign(set=s))
    df = pd.concat(frames, ignore_index=True)
    df = df[["set", "sample"] + [c for c in df.columns if c not in ("set", "sample")]]
    df.to_csv(d / "set_eigengenes.tsv", sep="\t", index=False, float_format="%.10g", lineterminator="\n")


def _file_stem(name: str, used: set) -> str:
    stem = re.sub(r"[^A-Za-z0-9._-]+", "_", str(name)).strip("._") or "celltype"
    out, i = stem, 1
    while out.lower() in used:
        i += 1
        out = f"{stem}_{i}"
    used.add(out.lower())
    return out


def write_celltype_modules(results: dict, output_dir: Any, n_hubs: int = 10,
                           write_metacells: bool = False, write_tom: bool = False,
                           threads: int = 2):
    """Write :func:`modules_by_celltype` results (a dict, or the pairs of
    :func:`iter_modules_by_celltype`, written and released one at a time):
    per cell type a directory
    ``<cell type>/`` (name made file-safe) with ``modules.tsv``,
    ``eigengenes.tsv``, ``kme.tsv`` and ``hub_genes.tsv`` (as ``gcen2
    module -e`` writes them: :func:`gcen2.tom_modules.write_module_files`;
    eigengenes and kME of the metacells), ``cell_eigengenes.tsv`` and
    ``cell_kme.tsv`` (on the single cells, if computed; also
    ``cell_harmonized_eigengenes.tsv``, ``dmes.tsv``, ``trait_cor.tsv``,
    ``sets.tsv`` and ``set_eigengenes.tsv`` when computed)
    ``metacell_members.tsv`` (metacell, cell, center: the cells averaged in
    each metacell, its centre cell first; ``gcen2 plot metacells``)
    and, with ``write_metacells``, ``metacells.tsv`` (genes x metacells,
    the input of ``gcen2 module -e``); and ``celltypes.tsv`` (cell_type,
    dir, n_cells, n_metacells, n_genes, power, modules, genes_in_modules,
    reassigned).
    Returns that table."""
    import pandas as pd

    from gcen2.tom_modules import write_module_files

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    rows, used = [], set()
    for t, res in (results.items() if isinstance(results, dict) else results):
        d = out / _file_stem(t, used)
        d.mkdir(exist_ok=True)
        consensus = bool(res.metadata.get("consensus"))
        if write_tom and consensus:
            warnings.warn(f"{t}: no tom_hubs.tsv for a consensus network.", stacklevel=2)
        write_module_files(res, modules=d / "modules.tsv", eigengenes=d / "eigengenes.tsv",
                           kme=d / "kme.tsv", hubs=d / "hub_genes.tsv", n_hubs=n_hubs,
                           soft_power=d / "soft_power.tsv", dendrogram=d / "dendrogram.json",
                           tom_hubs=d / "tom_hubs.tsv" if write_tom and not consensus else None,
                           threads=threads)
        md = res.metadata
        if write_metacells:
            md["metacells"].to_csv(d / "metacells.tsv", sep="\t", float_format="%.10g", lineterminator="\n")
        if md.get("metacell_members") is not None:
            md["metacell_members"].to_csv(d / "metacell_members.tsv", sep="\t", index=False, lineterminator="\n")
        _write_cell_level(res, d)
        if md.get("consensus"):
            _write_sets(res, d)
        _write_cell_stats(res, d)
        rows.append({"cell_type": t, "dir": d.name, "n_cells": md.get("n_cells"),
                     "n_metacells": md.get("n_metacells"), "n_genes": md.get("n_genes"),
                     "power": (",".join(f"{p:g}" for p in md["set_powers"].values())
                               if md.get("consensus") else res.power),
                     "modules": int(res.labels.max(initial=0)),
                     "genes_in_modules": int((res.labels > 0).sum()),
                     "reassigned": len(md.get("reassigned", []))})
        del res, md
    table = pd.DataFrame(rows, columns=["cell_type", "dir", "n_cells", "n_metacells", "n_genes",
                                        "power", "modules", "genes_in_modules", "reassigned"])
    table.to_csv(out / "celltypes.tsv", sep="\t", index=False, lineterminator="\n")
    return table
