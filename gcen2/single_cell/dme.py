"""Differential module eigengenes and module-trait correlations on single
cells, as hdWGCNA's ``FindDMEs`` / ``FindAllDMEs`` and
``ModuleTraitCorrelation`` compute them from the cell eigengenes
(harmonized ones when there are any).

``find_dmes``: Seurat ``FindMarkers`` on the eigengenes clipped at 0 (as
hdWGCNA does) with its defaults there (``test.use = "wilcox"``,
``logfc.threshold = 0``, ``min.pct = 0``, ``pseudocount.use = 0``): R's
``wilcox.test`` p-value (exact for fewer than 50 cells per group without
ties, else the normal approximation with tie and continuity corrections),
``avg_log2FC`` = log2 of the ratio of the group means of expm1(eigengene),
``pct.1`` / ``pct.2`` the fractions of cells above 0 (3 decimals), sorted by
p-value and then by ``|pct.1 - pct.2|`` (descending), ``p_val_adj`` Bonferroni
over all modules.

``module_trait_cor``: Pearson (or Spearman) correlation of every trait with
every eigengene on the pairwise complete cells (``Hmisc::rcorr``), its
Student t p-value, and Benjamini-Hochberg FDR over the traits x modules
table, for all cells and for each group with at least 5 cells.
Categorical traits become their level number (1, 2, ...), as hdWGCNA's
``as.numeric(factor)``.
"""

from __future__ import annotations

import warnings

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError


def _wilcox_p(x: np.ndarray, y: np.ndarray) -> float:
    """R ``wilcox.test(x, y)$p.value`` (two-sided, correct = TRUE)."""
    from scipy.stats import mannwhitneyu, norm, rankdata

    nx, ny = x.size, y.size
    r = rankdata(np.concatenate([x, y]))
    stat = r[:nx].sum() - nx * (nx + 1) / 2
    ties = np.unique(r).size < r.size
    if nx < 50 and ny < 50 and not ties:
        return float(mannwhitneyu(x, y, alternative="two-sided", method="exact").pvalue)
    _, counts = np.unique(r, return_counts=True)
    n = nx + ny
    z = stat - nx * ny / 2
    sigma = np.sqrt((nx * ny / 12) * ((n + 1) - (counts ** 3 - counts).sum() / (n * (n - 1))))
    if sigma == 0:
        return float("nan")
    z = (z - np.sign(z) * 0.5) / sigma
    return float(2 * min(norm.cdf(z), norm.sf(z)))


def _eigengene_frame(eigengenes) -> pd.DataFrame:
    if not isinstance(eigengenes, pd.DataFrame):
        raise GCENInputError("eigengenes must be a cells x modules DataFrame.")
    me = eigengenes.drop(columns=[c for c in eigengenes.columns
                                  if str(c) in ("grey", "ME_unassigned")])
    me.index = me.index.map(str)
    return me


def find_dmes(eigengenes, cells1, cells2, only_pos: bool = False,
              logfc_threshold: float = 0.0, min_pct: float = 0.0,
              pseudocount: float = 0.0) -> pd.DataFrame:
    """hdWGCNA ``FindDMEs``: modules whose eigengene differs between the
    cells ``cells1`` and ``cells2`` (names in the index of ``eigengenes``,
    a cells x modules DataFrame). Returns a table (index: module) with
    p_val, avg_log2FC, pct.1, pct.2, p_val_adj and module."""
    me = _eigengene_frame(eigengenes)
    c1 = [str(c) for c in cells1]
    c2 = [str(c) for c in cells2]
    missing = sorted((set(c1) | set(c2)) - set(me.index))
    if missing:
        raise GCENInputError(f"{len(missing)} cells not in the eigengenes "
                             f"(e.g. {', '.join(missing[:3])}).")
    if set(c1) & set(c2):
        raise GCENInputError("cells1 and cells2 overlap.")
    if len(c1) < 3 or len(c2) < 3:
        raise GCENInputError("each group needs at least 3 cells.")
    x1 = np.clip(me.loc[c1].to_numpy(dtype=np.float64), 0, None)
    x2 = np.clip(me.loc[c2].to_numpy(dtype=np.float64), 0, None)
    with np.errstate(divide="ignore", invalid="ignore"):
        fc = (np.log2((np.expm1(x1).sum(axis=0) + pseudocount) / x1.shape[0])
              - np.log2((np.expm1(x2).sum(axis=0) + pseudocount) / x2.shape[0]))
    pct1 = np.round((x1 > 0).mean(axis=0), 3)
    pct2 = np.round((x2 > 0).mean(axis=0), 3)
    keep = np.maximum(pct1, pct2) >= min_pct
    keep &= (fc >= logfc_threshold) if only_pos else (np.abs(fc) >= logfc_threshold)
    rows = []
    for j in np.flatnonzero(keep):
        rows.append({"p_val": _wilcox_p(x1[:, j], x2[:, j]), "avg_log2FC": fc[j],
                     "pct.1": pct1[j], "pct.2": pct2[j], "module": str(me.columns[j])})
    df = pd.DataFrame(rows, columns=["p_val", "avg_log2FC", "pct.1", "pct.2", "module"])
    if only_pos:
        df = df[df["avg_log2FC"] > 0]
    df = df.assign(_d=-(df["pct.1"] - df["pct.2"]).abs()) \
        .sort_values(["p_val", "_d"], kind="mergesort").drop(columns="_d")
    df.insert(4, "p_val_adj", np.minimum(df["p_val"] * me.shape[1], 1.0))
    return df.set_index(df["module"].to_numpy())


def find_all_dmes(eigengenes, groups, **kwargs) -> pd.DataFrame:
    """hdWGCNA ``FindAllDMEs``: :func:`find_dmes` of the cells of each group
    against all other cells (``groups``: one label per row of
    ``eigengenes``, or a Series indexed by cell), groups in order of first
    appearance. One table with a ``group`` column; groups with fewer than 3
    cells (or fewer than 3 others) are skipped with a warning."""
    me = _eigengene_frame(eigengenes)
    g = (groups.reindex(me.index) if isinstance(groups, pd.Series)
         else pd.Series(np.asarray(groups), index=me.index))
    if g.isna().any() or len(g) != len(me):
        raise GCENInputError("groups must give a label for every cell.")
    g = g.astype(str)
    out = []
    for value in pd.unique(g):
        inside = g.index[g == value]
        rest = g.index[g != value]
        if len(inside) < 3 or len(rest) < 3:
            warnings.warn(f"group '{value}' has too few cells (or too few others); "
                          "skipped.", stacklevel=2)
            continue
        d = find_dmes(me, inside, rest, **kwargs)
        d["group"] = value
        d.index = [f"{value}.{m}" for m in d.index]
        out.append(d)
    if not out:
        return pd.DataFrame(columns=["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj",
                                     "module", "group"])
    return pd.concat(out)


def _trait_matrix(traits: pd.DataFrame) -> pd.DataFrame:
    out = {}
    for c in traits.columns:
        v = traits[c]
        if isinstance(v.dtype, pd.CategoricalDtype):
            codes = v.cat.codes.to_numpy().astype(np.float64)
            codes[codes < 0] = np.nan
            out[c] = codes + 1
        elif pd.api.types.is_bool_dtype(v):
            out[c] = v.astype(np.float64).to_numpy()
        elif pd.api.types.is_numeric_dtype(v):
            out[c] = v.to_numpy(dtype=np.float64)
        else:
            cat = pd.Categorical(v.astype("string").to_numpy(dtype=object))
            warnings.warn(f"trait '{c}' is not numeric: its levels "
                          f"{', '.join(map(str, cat.categories))} are used as 1, 2, ...",
                          stacklevel=3)
            codes = cat.codes.astype(np.float64)
            codes[codes < 0] = np.nan
            out[c] = codes + 1
    return pd.DataFrame(out, index=traits.index)


def _bh(p: np.ndarray) -> np.ndarray:
    """R p.adjust(p, "fdr") (NA kept, not counted)."""
    out = np.full(p.shape, np.nan)
    ok = np.isfinite(p)
    q = p[ok]
    n = q.size
    if n:
        order = np.argsort(-q, kind="mergesort")
        adj = np.minimum.accumulate(q[order] * n / np.arange(n, 0, -1))
        res = np.empty(n)
        res[order] = np.minimum(adj, 1.0)
        out[ok] = res
    return out


def _cor_table(t: np.ndarray, m: np.ndarray, method: str):
    from scipy.stats import rankdata
    from scipy.stats import t as tdist

    k, q = t.shape[1], m.shape[1]
    r = np.full((k, q), np.nan)
    n = np.zeros((k, q), dtype=np.int64)
    for i in range(k):
        for j in range(q):
            ok = np.isfinite(t[:, i]) & np.isfinite(m[:, j])
            n[i, j] = ok.sum()
            if n[i, j] < 2:
                continue
            a, b = t[ok, i], m[ok, j]
            if method == "spearman":
                a, b = rankdata(a), rankdata(b)
            a = a - a.mean()
            b = b - b.mean()
            den = np.sqrt((a @ a) * (b @ b))
            if den > 0:
                r[i, j] = min(1.0, max(-1.0, (a @ b) / den))
    df = n - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        stat = np.abs(r) * np.sqrt(df) / np.sqrt(1 - r * r)
        p = 2 * tdist.sf(stat, df)
    p[np.abs(r) == 1] = 0.0
    return r, p, n


def module_trait_cor(eigengenes, traits, groups=None, method: str = "pearson",
                     min_cells: int = 5) -> pd.DataFrame:
    """hdWGCNA ``ModuleTraitCorrelation``: correlation of every trait
    (columns of ``traits``, rows matched to the cells of ``eigengenes`` by
    name) with every module eigengene, its p-value and the FDR, for all
    cells (group ``all_cells``) and, with ``groups`` (one label per cell or
    a Series indexed by cell), for each group with at least ``min_cells``
    cells (groups in sorted order). Returns a long table group, trait,
    module, cor, p, fdr, n."""
    if method not in ("pearson", "spearman"):
        raise GCENInputError("method must be pearson or spearman.")
    me = _eigengene_frame(eigengenes)
    if not isinstance(traits, pd.DataFrame):
        raise GCENInputError("traits must be a DataFrame (cells x traits).")
    traits = traits.copy()
    traits.index = traits.index.map(str)
    missing = me.index.difference(traits.index)
    if len(missing):
        raise GCENInputError(f"{len(missing)} cells have no trait values "
                             f"(e.g. {', '.join(missing[:3])}).")
    tm = _trait_matrix(traits.loc[me.index])
    blocks = [("all_cells", np.ones(len(me), dtype=bool))]
    if groups is not None:
        g = (groups.reindex(me.index) if isinstance(groups, pd.Series)
             else pd.Series(np.asarray(groups), index=me.index)).astype(str)
        for value in sorted(pd.unique(g)):
            sel = (g == value).to_numpy()
            if sel.sum() >= min_cells:
                blocks.append((value, sel))
    rows = []
    mods = [str(c) for c in me.columns]
    for name, sel in blocks:
        r, p, n = _cor_table(tm.to_numpy()[sel], me.to_numpy(dtype=np.float64)[sel], method)
        fdr = _bh(p.ravel()).reshape(p.shape)
        for i, tr in enumerate(tm.columns):
            for j, mod in enumerate(mods):
                rows.append((name, str(tr), mod, r[i, j], p[i, j], fdr[i, j], int(n[i, j])))
    return pd.DataFrame(rows, columns=["group", "trait", "module", "cor", "p", "fdr", "n"])
