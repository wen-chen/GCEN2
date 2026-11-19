"""Module preservation (WGCNA ``modulePreservation``, as hdWGCNA's
``ModulePreservation`` calls it): how well the modules of a reference data
set hold together in a test data set.

:func:`module_preservation` follows WGCNA 1.74 for expression data with one
reference and one test set, without interpolation (``useInterpolation =
FALSE``), cluster coefficients or accuracy statistics (the test set has no
modules of its own):

- a "gold" module of ``min(max_gold_module_size, genes / 2)`` random genes
  is added (WGCNA's baseline), and modules (grey included) with more than
  ``max_module_size`` genes are subsampled;
- quality (in the reference) and preservation (in the test set) statistics
  of every module: propVarExplained, meanSignAwareKME, separability,
  meanSignAwareCorDat, meanAdj and meanMAR (adjacency ``|r|^6`` for "unsigned",
  ((1 + r) / 2)^12 for "signed"), and cor.kIM, cor.kME, cor.kMEall, cor.cor
  and cor.MAR between the sets;
- ``n_permutations`` permutations of the genes; Z = (observed - mean) / sd
  over the permutations; Zdensity.pres = median Z of propVarExplained,
  meanSignAwareKME, meanSignAwareCorDat and meanAdj; Zconnectivity.pres =
  median Z of cor.kIM, cor.kME and cor.cor; Zsummary.pres = their mean
  (and the median ranks of the observed statistics).

The random draws (the gold genes, the subsampling and the permutations) are
those of R's ``sample()`` after ``set.seed(seed)``, in WGCNA's order
(:mod:`gcen2._rrng`), so with the same data the statistics are WGCNA's (to
about 1e-8 relative: meanMAR subtracts nearly equal sums, which magnifies
last-bit differences of the correlations); the permutations run in parallel.
hdWGCNA uses ``n_permutations = 500``, ``seed = 12345`` and an unsigned
network (WGCNA's default).
"""

from __future__ import annotations

from concurrent.futures import ThreadPoolExecutor

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError

_QUAL = ["propVarExplained", "meanSignAwareKME", "separability", "meanSignAwareCorDat",
         "meanAdj", "meanClusterCoeff", "meanMAR"]
_INTER = ["cor.kIM", "cor.kME", "cor.kMEall", "cor.cor", "cor.clusterCoeff", "cor.MAR"]


def _std_cols(x: np.ndarray) -> np.ndarray:
    xc = x - x.mean(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        return xc / np.linalg.norm(xc, axis=0)


def _cor(x: np.ndarray, y: np.ndarray | None = None) -> np.ndarray:
    """Pearson correlation of the columns of x (with those of y)."""
    a = _std_cols(x)
    b = a if y is None else _std_cols(y)
    return np.clip(a.T @ b, -1.0, 1.0)


def _cor1(a: np.ndarray, b: np.ndarray) -> float:
    ok = np.isfinite(a) & np.isfinite(b)
    if ok.sum() < 2:
        return np.nan
    a, b = a[ok] - a[ok].mean(), b[ok] - b[ok].mean()
    den = np.sqrt((a @ a) * (b @ b))
    return float(np.clip((a @ b) / den, -1, 1)) if den > 0 else np.nan


def _eigengenes(dat: np.ndarray, colors: np.ndarray, levels) -> np.ndarray:
    """WGCNA moduleEigengenes (samples x levels) of dat (samples x genes)."""
    from gcen2.tom_modules import module_eigengene

    return np.column_stack([module_eigengene(dat[:, colors == lv].T) for lv in levels])


def _lower(m: np.ndarray) -> np.ndarray:
    return m[np.tril_indices(m.shape[0], -1)]


def _core(dat_ref, dat_ref_p, dat_test, colors, levels, n_type, permuted):
    """WGCNA .coreCalcForExpr (densityOnly = FALSE, no cluster coefficients):
    quality (levels x 7), intra (levels x 7), inter (levels x 6)."""
    n_mods = len(levels)
    sizes = np.array([(colors == lv).sum() for lv in levels])
    act = sizes > 1
    me1 = _eigengenes(dat_ref, colors, levels)
    kme1 = _cor(dat_ref, me1)
    if permuted:
        me2 = _eigengenes(dat_ref_p, colors, levels)
        kme2 = _cor(dat_ref_p, me2)
    else:
        me2, kme2 = me1, kme1
    me3 = _eigengenes(dat_test, colors, levels)
    kme3 = _cor(dat_test, me3)
    mod = [np.flatnonzero(colors == lv) for lv in levels]
    qual = np.full((n_mods, 7), np.nan)
    intra = np.full((n_mods, 7), np.nan)
    inter = np.full((n_mods, 6), np.nan)
    for j in range(n_mods):
        if not act[j]:
            continue
        g = mod[j]
        inter[j, 1] = abs(_cor1(kme1[g, j], kme3[g, j]))
        inter[j, 2] = abs(_cor1(kme1[:, j], kme3[:, j]))
        qual[j, 0] = np.nanmean(kme2[g, j] ** 2)
        intra[j, 0] = np.nanmean(kme3[g, j] ** 2)
        qual[j, 1] = abs(np.nanmean(np.abs(kme2[g, j])))
        intra[j, 1] = abs(np.nanmean(np.sign(kme1[g, j]) * kme3[g, j]))
    excl = np.array([lv in ("gold", "grey") for lv in levels])
    for k, me in enumerate((me2, me3)):
        c = _cor(me)
        if n_type == 0:  # as WGCNA (never: network types are 1..3)
            c = np.abs(c)
        np.fill_diagonal(c, 0.0)
        with np.errstate(invalid="ignore"):
            sep = 1 - np.max(np.where(np.isnan(c[:, ~excl]), -np.inf, c[:, ~excl]), axis=1) \
                if (~excl).any() else np.full(n_mods, np.inf)
        (qual if k == 0 else intra)[:, 2] = sep
    for j in range(n_mods):
        if not act[j]:
            continue
        g = mod[j]
        c1 = _cor(dat_ref[:, g])
        c2 = _cor(dat_ref_p[:, g]) if permuted else c1
        c3 = _cor(dat_test[:, g])
        s2 = np.abs(c2) if n_type == 1 else c2
        s3 = np.sign(c1) * c3
        qual[j, 3] = np.nanmean(_lower(s2))
        intra[j, 3] = np.nanmean(_lower(s3))
        inter[j, 3] = _cor1(_lower(c1), _lower(c3))
        if n_type == 1:
            a1, a2, a3 = c1 ** 6, c2 ** 6, c3 ** 6
        elif n_type == 2:
            a1, a2, a3 = ((1 + c1) / 2) ** 12, ((1 + c2) / 2) ** 12, ((1 + c3) / 2) ** 12
        else:
            a1, a2, a3 = c1 ** 6, c2 ** 6, c3 ** 6
            a1[c1 < 0] = 0
            a2[c2 < 0] = 0
            a3[c3 < 0] = 0

        def mar(a):
            den = a.sum(axis=0) - 1
            with np.errstate(divide="ignore", invalid="ignore"):
                m = ((a ** 2).sum(axis=0) - 1) / den
            m[den == 0] = np.nan
            return m

        m1, m2, m3 = mar(a1), mar(a2), mar(a3)
        qual[j, 6] = np.mean(m2)
        intra[j, 6] = np.mean(m3)
        inter[j, 0] = _cor1(np.nansum(a1, axis=0), np.nansum(a3, axis=0))
        inter[j, 5] = _cor1(m1, m3)
        qual[j, 4] = np.nanmean(_lower(a2))
        intra[j, 4] = np.nanmean(_lower(a3))
    return sizes, qual, intra, inter


def _draws(rng, colors: np.ndarray, gold_size: int, max_module_size: int, permuted: bool):
    """The random draws of one .modulePreservationInternal call, in WGCNA's
    order: gold genes (reference and test), the genes kept of large modules,
    and the permutation (0-based indices)."""
    n = colors.size
    if permuted:
        rng.sample(n, n)  # colorRefAcc = sample(colorRefAcc)
        rng.sample(n, n)  # colorTestAcc = sample(colorRefAcc)
    gold_r = rng.sample(n, gold_size) - 1
    gold_t = rng.sample(n, gold_size) - 1 if permuted else gold_r
    keep = np.ones(n, dtype=bool)
    for lv in sorted(set(colors.tolist())):
        inm = np.flatnonzero(colors == lv)
        if inm.size > max_module_size:
            sam = rng.sample(inm.size, max_module_size) - 1
            keep[inm] = False
            keep[inm[sam]] = True
    perm = None
    if keep.sum() < n:
        if permuted:
            perm = rng.sample(int(keep.sum()), int(keep.sum())) - 1
    elif permuted:
        perm = rng.sample(n, n) - 1
    return gold_r, gold_t, keep, perm


def _internal(dat_ref, dat_test, colors, gold_size, draws, permuted, n_type):
    gold_r, gold_t, keep, perm = draws
    n = colors.size
    gref, grefp, gtest = dat_ref[:, gold_r], dat_ref[:, gold_t], dat_test[:, gold_t]
    if keep.sum() < n:
        c = colors[keep]
        dref = dat_ref[:, keep]
        if permuted:
            dtest = dat_test[:, perm]  # WGCNA: datTest[, keepPerm] (not subset first)
            drefp = dref[:, perm]
        else:
            dtest = dat_test[:, keep]
            drefp = dref
    else:
        c = colors
        dref = dat_ref
        if permuted:
            drefp, dtest = dat_ref[:, perm], dat_test[:, perm]
        else:
            drefp, dtest = dat_ref, dat_test
    dref = np.hstack([dref, gref])
    drefp = np.hstack([drefp, grefp])
    dtest = np.hstack([dtest, gtest])
    c2 = np.concatenate([c.astype(object), np.array(["gold"] * gold_size, dtype=object)])
    levels = sorted(set(c2.tolist()))
    return levels, _core(dref, drefp, dtest, c2, levels, n_type, permuted)


def _ranks(x: np.ndarray) -> np.ndarray:
    """R rank(-x, na.last = "keep") per column (average ties)."""
    from scipy.stats import rankdata

    out = np.full(x.shape, np.nan)
    for j in range(x.shape[1]):
        ok = np.isfinite(x[:, j])
        out[ok, j] = rankdata(-x[ok, j])
    return out


def module_preservation(ref, test, modules, n_permutations: int = 500, seed: int = 12345,
                        network_type: str = "unsigned", max_gold_module_size: int = 1000,
                        max_module_size: int = 1000, threads: int = 2):
    """Preservation of the modules of ``ref`` in ``test`` (see the module
    docstring), as hdWGCNA's ``ModulePreservation``.

    ``ref`` / ``test``: samples (e.g. metacells) x genes DataFrames; the
    genes they share are used (in the order of ``ref``). ``modules``: the
    module of every gene of ``ref`` (a Series indexed by gene; genes absent
    or labelled "grey" (or 0) are unassigned). Returns ``(observed, Z)``:
    DataFrames indexed by module (with "gold" and "grey"), WGCNA's quality
    and preservation tables side by side (all-missing columns dropped, as
    WGCNA's ``discardInvalidOutput``), the Z table with ``Zsummary.pres``.
    """
    if network_type not in ("unsigned", "signed", "signed hybrid"):
        raise GCENInputError("network_type must be unsigned, signed or signed hybrid.")
    n_type = {"unsigned": 1, "signed": 2, "signed hybrid": 3}[network_type]
    if not isinstance(ref, pd.DataFrame) or not isinstance(test, pd.DataFrame):
        raise GCENInputError("ref and test must be samples x genes DataFrames.")
    if n_permutations < 2:
        raise GCENInputError("n_permutations must be at least 2.")
    in_test = set(test.columns)
    genes = [g for g in ref.columns if g in in_test]
    if len(genes) < 10:
        raise GCENInputError(f"ref and test share only {len(genes)} genes.")
    dref = ref[genes].to_numpy(dtype=np.float64)
    dtest = test[genes].to_numpy(dtype=np.float64)
    for name, d in (("ref", dref), ("test", dtest)):
        bad = ~np.isfinite(d).all(axis=0) | (d.std(axis=0) == 0)
        if bad.any():
            raise GCENInputError(f"{int(bad.sum())} genes are constant or not finite in "
                                 f"{name} (WGCNA's goodSamplesGenes would stop); remove them.")
    mods = pd.Series(modules)
    mods.index = mods.index.map(str)
    lab = mods.reindex(genes).astype(object)
    lab = lab.where(lab.notna(), "grey").astype(str)
    lab[lab.isin(["0", "grey", "nan"])] = "grey"
    colors = lab.to_numpy(dtype=object)
    if len(set(colors.tolist())) < 2:
        raise GCENInputError("give at least one module besides the unassigned genes.")
    if "gold" in set(colors.tolist()):
        raise GCENInputError("'gold' is reserved for WGCNA's random module.")
    n = colors.size
    gold_size = int(min(max_gold_module_size, n / 2))

    from gcen2._rrng import RRandom

    rng = RRandom(seed)
    obs_draws = _draws(rng, colors, gold_size, max_module_size, False)
    perm_draws = [_draws(rng, colors, gold_size, max_module_size, True)
                  for _ in range(n_permutations)]
    from gcen2.tom_modules import _blas_threads

    with _blas_threads(1):
        levels, (sizes, qual, intra, inter) = _internal(dref, dtest, colors, gold_size,
                                                        obs_draws, False, n_type)

        def one(dr):
            _, (_, q, i, e) = _internal(dref, dtest, colors, gold_size, dr, True, n_type)
            return np.hstack([q, i, e])

        with ThreadPoolExecutor(max_workers=max(1, int(threads))) as ex:
            perm = np.stack(list(ex.map(one, perm_draws)), axis=2)  # levels x 20 x perms

    reg_obs = np.hstack([qual, intra, inter])
    import warnings

    with np.errstate(invalid="ignore", divide="ignore"), warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)  # all-missing statistics
        z = (reg_obs - np.nanmean(perm, axis=2)) / np.nanstd(perm, axis=2, ddof=1)
    qnames = [f"{s}.qual" for s in _QUAL]
    pnames = [f"{s}.pres" for s in _QUAL]
    z_qual = pd.DataFrame(z[:, :7], index=levels, columns=[f"Z.{s}" for s in qnames])
    z_pres = pd.DataFrame(z[:, 7:], index=levels,
                          columns=[f"Z.{s}" for s in pnames + _INTER])

    # observed tables (separability apart), median ranks
    q = pd.DataFrame(qual, index=levels, columns=qnames).drop(columns="separability.qual")
    rq = _ranks(q.iloc[:, [0, 1, 2, 3]].to_numpy())
    obs_qual = pd.concat([pd.DataFrame({"moduleSize": sizes,
                                        "medianRank.qual": np.nanmedian(rq, axis=1)},
                                       index=levels), q], axis=1)
    pres = pd.concat([pd.DataFrame(intra, index=levels, columns=pnames)
                      .drop(columns="separability.pres"),
                      pd.DataFrame(inter, index=levels, columns=_INTER)], axis=1)
    rd = _ranks(pres.iloc[:, [0, 1, 2, 3]].to_numpy())
    rc = _ranks(pres.iloc[:, [6, 7, 9]].to_numpy())
    obs_pres = pd.concat([pd.DataFrame({
        "medianRank.pres": np.nanmedian(np.hstack([rd, rc]), axis=1),
        "medianRankDensity.pres": np.nanmedian(rd, axis=1),
        "medianRankConnectivity.pres": np.nanmedian(rc, axis=1)}, index=levels), pres], axis=1)

    zq = z_qual.drop(columns="Z.separability.qual")
    zq.insert(0, "Zsummary.qual", np.nanmedian(zq.iloc[:, :4].to_numpy(), axis=1))
    zp = z_pres.drop(columns="Z.separability.pres")
    zdens = np.nanmedian(zp.iloc[:, [0, 1, 2, 3]].to_numpy(), axis=1)
    zconn = np.nanmedian(zp.iloc[:, [6, 7, 9]].to_numpy(), axis=1)
    zp.insert(0, "Zconnectivity.pres", zconn)
    zp.insert(0, "Zdensity.pres", zdens)
    zp.insert(0, "Zsummary.pres", np.nanmean(np.column_stack([zdens, zconn]), axis=1))

    def tidy(df):
        return df.loc[:, df.notna().any(axis=0)]

    observed = tidy(pd.concat([obs_qual, obs_pres], axis=1))
    zt = tidy(pd.concat([pd.DataFrame({"moduleSize": sizes}, index=levels), zq, zp], axis=1))
    observed.index.name = zt.index.name = "module"
    return observed, zt
