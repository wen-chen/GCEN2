"""Projecting modules onto another data set (hdWGCNA ``ProjectModules``).

:func:`project_modules`: the eigengenes of reference modules on the cells
of a query data set, as hdWGCNA computes them: modules with fewer than
``min_overlap`` (0.5) of their genes in the query are left out; the others
keep their genes found in the query; then ``ModuleEigengenes`` on the query
cells (Seurat LogNormalize, ScaleData of the module's genes, first principal
component signed by the average expression;
:func:`gcen2.single_cell.cell_eigengenes`), optionally harmonized
(``group.by.vars``).
"""

from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError


def _module_table(modules) -> pd.Series:
    """gene -> module (str) from a Series / dict, a DataFrame with gene and
    module columns, or a TomModules (modules 1..K)."""
    if hasattr(modules, "labels") and hasattr(modules, "genes"):
        s = pd.Series(modules.labels, index=modules.genes)
        s = s[s > 0].map(lambda k: f"ME{k}")
    elif isinstance(modules, pd.DataFrame):
        cols = {c.lower(): c for c in modules.columns}
        g = cols.get("gene") or cols.get("gene_name")
        m = cols.get("module")
        if g is None or m is None:
            raise GCENInputError("the modules table needs gene and module columns.")
        s = pd.Series(modules[m].to_numpy(), index=modules[g].astype(str).to_numpy())
    else:
        s = pd.Series(modules)
    s.index = s.index.map(str)
    s = s.astype(str)
    return s[~s.index.duplicated()]


def project_modules(adata: Any, modules, min_overlap: float = 0.5, layer: str | None = None,
                    normalize: bool | None = None, cells=None, batches=None,
                    harmonize_by: str | None = None, npcs: int = 50, seed: int = 0,
                    include_unassigned: bool = False, **harmony_kwargs):
    """Eigengenes of reference ``modules`` on the cells of ``adata`` (the
    query), as hdWGCNA's ``ProjectModules``.

    ``modules``: gene -> module (a Series or dict, a table with gene and
    module columns such as ``modules.tsv`` of GCEN2 or hdWGCNA's
    ``GetModules``, or a :class:`gcen2.tom_modules.TomModules`); modules
    named "grey" (hdWGCNA's unassigned genes) are skipped unless
    ``include_unassigned`` (hdWGCNA projects them too). ``cells``: only these
    query cells (a boolean mask or names). ``batches`` (one label per cell)
    or ``harmonize_by`` (an obs column): also the harmonized eigengenes.

    Returns ``(eigengenes, harmonized, table)``: cells x modules DataFrames
    (``harmonized`` None without batches) and a table module, genes,
    genes_found, projected.
    """
    from gcen2._anndata import require_anndata
    from gcen2.single_cell.metacell import _expression, _is_counts, _rows
    from gcen2.single_cell.modules import cell_eigengenes

    require_anndata()
    if not 0 <= min_overlap <= 1:
        raise GCENInputError("min_overlap must be between 0 and 1.")
    mods = _module_table(modules)
    if not include_unassigned:
        mods = mods[mods != "grey"]
    if mods.empty:
        raise GCENInputError("no modules to project.")
    var = {g: i for i, g in enumerate(map(str, adata.var_names))}
    rows = []
    order = list(pd.unique(mods.to_numpy()))
    for m in order:
        genes = mods.index[mods == m]
        found = [g for g in genes if g in var]
        rows.append({"module": m, "genes": len(genes), "genes_found": len(found),
                     "projected": len(found) / len(genes) >= min_overlap and len(found) >= 2})
    table = pd.DataFrame(rows)
    skipped = table.loc[~table["projected"], "module"].tolist()
    if skipped:
        warnings.warn(f"modules not projected (too few genes in the query): "
                      f"{', '.join(skipped)}", stacklevel=2)
    keep = table.loc[table["projected"], "module"].tolist()
    if not keep:
        raise GCENInputError("no module has enough genes in the query.")

    if cells is None:
        idx = np.arange(adata.n_obs)
    else:
        c = np.asarray(cells)
        if c.dtype == bool:
            idx = np.flatnonzero(c)
        else:
            pos = {n: i for i, n in enumerate(map(str, adata.obs_names))}
            missing = [str(x) for x in c if str(x) not in pos]
            if missing:
                raise GCENInputError(f"{len(missing)} cells not in adata (e.g. {missing[0]}).")
            idx = np.array([pos[str(x)] for x in c], dtype=np.int64)
    if idx.size < 3:
        raise GCENInputError("need at least 3 query cells.")
    if harmonize_by is not None:
        if harmonize_by not in adata.obs.columns:
            raise GCENInputError(f"'{harmonize_by}' not found in adata.obs.")
        batches = adata.obs[harmonize_by].astype(str).to_numpy()[idx]
    elif batches is not None:
        batches = np.asarray(batches).astype(str)
        if batches.size == adata.n_obs and idx.size != adata.n_obs:
            batches = batches[idx]
        if batches.size != idx.size:
            raise GCENInputError("batches must give one label per cell.")

    X = _expression(adata, layer)
    Xc = _rows(X, idx)
    if normalize is None:
        normalize = _is_counts(Xc)
    gene_idx, labels = [], []
    for k, m in enumerate(keep, start=1):
        for g in mods.index[mods == m]:
            if g in var:
                gene_idx.append(var[g])
                labels.append(k)
    out = cell_eigengenes(Xc, np.array(gene_idx), np.array(labels),
                          list(range(1, len(keep) + 1)), normalize=normalize,
                          batches=batches, npcs=npcs, seed=seed,
                          harmony_kwargs=harmony_kwargs or None)
    names = np.asarray(adata.obs_names, dtype=str)[idx]
    me = pd.DataFrame(out[0], index=names, columns=keep)
    hme = None if batches is None else pd.DataFrame(out[2], index=names, columns=keep)
    return me, hme, table
