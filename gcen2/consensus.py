"""Consensus modules of several expression sets (WGCNA
``blockwiseConsensusModules``, as hdWGCNA's ``ConstructNetwork(consensus =
TRUE)`` calls it).

:func:`consensus_tom_modules` follows WGCNA 1.74 on one block (all genes in
one TOM):

1. the TOM of each set, with its own soft power (``power="auto"``: each
   set's power as :func:`gcen2.tom_modules.tom_modules` picks it);
2. calibration ("single quantile"): each set's TOM is raised to the power
   that maps its ``calibration_quantile`` quantile onto the first set's
   (quantile type 8 of ``calibration_sample_factor / 0.05`` entries drawn by
   R's ``sample()`` after ``set.seed(calibration_seed)``, as WGCNA's
   ``consensusTOM``: :mod:`gcen2._rrng`);
3. the consensus TOM: per gene pair the ``consensus_quantile`` quantile over
   the sets (0: the minimum), by WGCNA's own quantile routine
   (``rowQuantileC``), so with its rounding;
4. average-linkage clustering of 1 - consensus TOM and the dynamic tree cut
   (as :func:`gcen2.tom_modules.tom_modules`);
5. the kME checks on the consensus kME (the ``consensus_quantile`` quantile
   over the sets of each gene's kME in each set, against that set's
   eigengene); WGCNA's consensus reassignment compares lower-tail
   chi-square p-values and never moves a gene, so there is none;
6. ``mergeCloseModules`` on the consensus eigengene dissimilarity (per pair
   of modules, the ``1 - consensus_quantile`` quantile over the sets of
   1 - cor of their eigengenes in the set).

The result is a :class:`gcen2.tom_modules.TomModules` whose samples are the
samples of all sets (eigengenes, kME and hub genes on all of them together);
``metadata`` has the sets' names, powers, calibration and eigengenes.
"""

from __future__ import annotations

import math

import numpy as np

from gcen2.core import require_extension
from gcen2.exceptions import GCENInputError
from gcen2.tom_modules import (
    HDWGCNA,
    TomModules,
    _blas_threads,
    _check_module_options,
    _eigengenes,
    _merge_modules,
    _pick_power,
    _relabel_by_size,
    cutree_hybrid,
    module_eigengenes,
)

# hdWGCNA 0.4.12 ConstructNetwork(consensus = TRUE): its defaults, with the
# power of each set by TestSoftPowersConsensus (pickSoftThreshold with
# bicor per set; the lowest power above 3 with R^2 >= 0.8)
HDWGCNA_CONSENSUS = dict(HDWGCNA, consensus_quantile=0.3)


def r_quantile7(x: np.ndarray, p: float, axis: int = 0) -> np.ndarray:
    """R ``quantile(x, p, type = 7, na.rm = TRUE)`` along ``axis`` (same
    arithmetic: (1 - h) * x[lo] + h * x[hi])."""
    x = np.sort(np.moveaxis(np.asarray(x, dtype=np.float64), axis, 0), axis=0)  # NaN last
    n = np.isfinite(x).sum(axis=0)
    out = np.full(x.shape[1:], np.nan)
    for m in np.unique(n):
        if m == 0:
            continue
        sel = n == m
        index = 1 + (m - 1) * p
        lo, hi = math.floor(index), math.ceil(index)
        qs = x[lo - 1][sel]
        xhi = x[hi - 1][sel]
        h = index - lo
        do = (index > lo) & (xhi != qs)
        qs = np.where(do, (1 - h) * qs + h * xhi, qs)
        out[sel] = qs
    return out


def r_quantile8(x: np.ndarray, p: float) -> float:
    """R ``quantile(x, p, type = 8)`` of a vector without missing values."""
    x = np.sort(np.asarray(x, dtype=np.float64))
    n = x.size
    a = b = 1.0 / 3.0
    fuzz = 4 * np.finfo(float).eps
    nppm = a + p * (n + 1 - a - b)
    j = math.floor(nppm + fuzz)
    h = nppm - j
    if abs(h) < fuzz:
        h = 0.0
    xx = np.concatenate([[x[0], x[0]], x, [x[-1], x[-1]]])
    qs = xx[j + 1]  # R x[j + 2L], 1-based
    if h == 1:
        return float(xx[j + 2])
    if 0 < h < 1 and xx[j + 1] != xx[j + 2]:
        return float((1 - h) * xx[j + 1] + h * xx[j + 2])
    return float(qs)


def _upper(m: np.ndarray, out: np.ndarray) -> None:
    """The upper triangle of the square matrix m, row by row (R's dist order
    of the lower triangle) into out."""
    n = m.shape[0]
    pos = 0
    for i in range(n - 1):
        k = n - 1 - i
        out[pos:pos + k] = m[i, i + 1:]
        pos += k


def _fill_diss(cons: np.ndarray, n: int) -> np.ndarray:
    """n x n matrix 1 - cons (cons in dist order) in both triangles, diagonal 0."""
    d = np.zeros((n, n))
    pos = 0
    for i in range(n - 1):
        k = n - 1 - i
        d[i, i + 1:] = 1.0 - cons[pos:pos + k]
        d[i + 1:, i] = d[i, i + 1:]
        pos += k
    return d


def _cor_cols(x: np.ndarray, me: np.ndarray) -> np.ndarray:
    """Pearson correlation of the rows of x (genes x samples) with the vector me."""
    xc = x - x.mean(axis=1, keepdims=True)
    mc = me - me.mean()
    with np.errstate(divide="ignore", invalid="ignore"):
        return (xc @ mc) / (np.linalg.norm(xc, axis=1) * np.linalg.norm(mc))


def consensus_post_process(xs, labels: np.ndarray, signed: bool, min_size: float,
                           consensus_quantile: float, min_core_kme: float = 0.5,
                           min_kme_to_stay: float = 0.2) -> np.ndarray:
    """blockwiseConsensusModules after cutreeDynamic: the consensus kME
    checks (``xs``: genes x samples of each set)."""
    labels = labels.copy()
    levels = [int(k) for k in np.unique(labels[labels != 0])]
    if not levels:
        return labels
    mes = [_eigengenes(x, labels, levels) for x in xs]
    for i, k in enumerate(levels):
        genes = np.flatnonzero(labels == k)
        kme = np.column_stack([_cor_cols(x[genes], me[i]) for x, me in zip(xs, mes)])
        if not signed:
            kme = np.abs(kme)
        cons = r_quantile7(kme, consensus_quantile, axis=1)
        if (cons > min_core_kme).sum() < min_size / 3:
            labels[genes] = 0
        elif (cons < min_kme_to_stay).any():
            labels[genes[cons < min_kme_to_stay]] = 0
            if (labels[genes] > 0).sum() < min_size:
                labels[genes] = 0
    return labels


def consensus_merge_close_modules(xs, labels: np.ndarray, cut_height: float,
                                  consensus_quantile: float) -> np.ndarray:
    """WGCNA mergeCloseModules on several sets (iterate, relabel by size):
    modules are merged when the ``1 - consensus_quantile`` quantile over the
    sets of 1 - cor of their eigengenes is below ``cut_height`` (average
    linkage)."""
    def dissimilarity(labels, levels):
        diss = [1.0 - np.corrcoef(_eigengenes(x, labels, levels)) for x in xs]
        return r_quantile7(np.array(diss), 1 - consensus_quantile, axis=0)

    return _merge_modules(labels, dissimilarity, cut_height)


def _set_matrices(sets, gene_names, set_names, sample_names):
    """Validated genes x samples float64 matrices with one gene list."""
    from gcen2.network import _prepare_matrix
    from gcen2.tom_modules import _sample_names

    if isinstance(sets, dict):
        if set_names is None:
            set_names = list(sets)
        sets = list(sets.values())
    sets = list(sets)
    if len(sets) < 1:
        raise GCENInputError("give at least one expression set.")
    if set_names is None:
        set_names = [f"set_{i + 1}" for i in range(len(sets))]
    set_names = [str(s) for s in set_names]
    if len(set_names) != len(sets) or len(set(set_names)) != len(sets):
        raise GCENInputError("set_names must give one distinct name per set.")
    mats, samples = [], []
    genes = None
    for i, expr in enumerate(sets):
        names = _sample_names(expr, sample_names[i] if sample_names is not None else None)
        m, g = _prepare_matrix(expr, "pearson", None, gene_names)
        if genes is None:
            genes = list(g)
        elif list(g) != genes:
            if sorted(g) != sorted(genes) or len(set(g)) != len(g):
                raise GCENInputError(f"set '{set_names[i]}' does not have the genes of "
                                     f"set '{set_names[0]}'.")
            pos = {name: j for j, name in enumerate(g)}
            m = np.ascontiguousarray(m[[pos[name] for name in genes]])
        if names is None:
            names = [f"{set_names[i]}_{j + 1}" for j in range(m.shape[1])]
        elif len(names) != m.shape[1]:
            raise GCENInputError(f"set '{set_names[i]}': one sample name per column needed.")
        mats.append(m)
        samples.append([str(s) for s in names])
    return mats, genes, set_names, samples


def _set_toms(xs, set_names, set_power, pick, signed: bool, threads: int):
    """The TOM of every set (upper triangle, float64) and how its power was
    found: ``{"power", "estimated", "fits", "rule"}``. A power of 0 is picked
    as pickSoftThreshold does, or by ``pick`` (method, R^2 cut, minimum
    power, powers; :func:`gcen2.tom_modules._pick_power`) when given."""
    from gcen2._build_file import _FIT_COLUMNS

    core = require_extension()
    n = xs[0].shape[0]
    toms, found = [], []
    for name, x, p in zip(set_names, xs, set_power):
        if x.shape[1] < 3:
            raise GCENInputError(f"set '{name}' has fewer than 3 samples.")
        rule = None
        if p == 0 and pick is not None:
            method, r2, min_power, powers = pick
            p, meta, rule = _pick_power(x, method, signed, r2, min_power, powers, threads)
        tom, info = core.module.tom_dissimilarity(x, "pearson", signed, p, [], threads,
                                                  similarity=True)
        if rule is None:
            meta = {"power": info["power"], "estimated": info["estimated"],
                    "fits": [dict(zip(_FIT_COLUMNS, f)) for f in info["fits"]]}
        t = np.empty(n * (n - 1) // 2)
        _upper(tom, t)
        del tom
        toms.append(t)
        found.append({"power": meta["power"], "estimated": meta["estimated"],
                      "fits": meta["fits"], "rule": rule})
    return toms, found


def _consensus_tom(toms, consensus_quantile, calibration, calibration_quantile,
                   sample_for_calibration, calibration_sample_factor, calibration_seed,
                   threads):
    """WGCNA consensusTOM: each set's TOM calibrated to the first one ("single
    quantile": raised to the power that matches their ``calibration_quantile``
    quantiles, on a random sample of the pairs as WGCNA draws it), then the
    ``consensus_quantile`` quantile over the sets (0: the minimum). The TOMs
    are modified in place. Returns (consensus, quantiles, powers)."""
    n_sets, n_pairs = len(toms), toms[0].size
    scale_quant = [1.0] * n_sets
    scale_power = [1.0] * n_sets
    if calibration == "single quantile" and n_sets > 1:
        idx = None
        if sample_for_calibration:
            from gcen2._rrng import RRandom

            qx = min(calibration_quantile, 1 - calibration_quantile)
            n_sc = int(min(calibration_sample_factor * 1 / qx, n_pairs))
            idx = RRandom(calibration_seed).sample(n_pairs, n_sc) - 1
        for s in range(n_sets):
            vals = toms[s] if idx is None else toms[s][idx]
            scale_quant[s] = r_quantile8(vals, calibration_quantile)
            if s > 0:
                scale_power[s] = math.log(scale_quant[0]) / math.log(scale_quant[s])
                np.power(toms[s], scale_power[s], out=toms[s])
    if n_sets == 1:
        cons = toms[0]
    elif consensus_quantile == 0:
        cons = toms[0]
        for t in toms[1:]:
            np.minimum(cons, t, out=cons)
    else:
        cons = require_extension().module.row_quantile(toms, consensus_quantile, threads)
    return cons, scale_quant, scale_power


def consensus_tom_modules(sets, gene_names=None, set_names=None, sample_names=None,
                          signed: bool = False, power="auto", consensus_quantile: float = 0.0,
                          calibration: str = "single quantile",
                          calibration_quantile: float = 0.95,
                          sample_for_calibration: bool = True,
                          calibration_sample_factor: float = 1000,
                          calibration_seed: int = 54321,
                          min_size: int = 20, deep_split: int = 2, merge_cut: float = 0.15,
                          detect_cut: float = 0.995, pam_stage: bool = True,
                          min_core_kme: float = 0.5, min_kme_to_stay: float = 0.2,
                          reassign_threshold: float | None = None,
                          power_method: str | None = None, power_r2: float | None = None,
                          min_power: float = 0, powers=None,
                          threads: int = 2) -> TomModules:
    """Consensus modules of several expression sets, as WGCNA
    ``blockwiseConsensusModules`` (see the module docstring).

    ``sets``: a list (or dict name -> set) of genes x samples matrices
    (DataFrames, AnnData, or arrays with ``gene_names``) with the same genes;
    ``set_names`` and ``sample_names`` (a list per set) name them. ``power``:
    "auto", one power for all sets, or one per set; ``power_method``,
    ``power_r2``, ``min_power`` and ``powers`` as in
    :func:`gcen2.tom_modules.tom_modules`, applied to each set.
    ``consensus_quantile`` (0: the minimum over the sets),
    ``calibration`` ("single quantile" or "none"), ``calibration_quantile``,
    ``sample_for_calibration``, ``calibration_sample_factor`` and the
    clustering options are blockwiseConsensusModules' (``minModuleSize``
    20, ``deepSplit`` 2, ``mergeCutHeight`` 0.15, ``minKMEtoStay`` 0.2);
    ``calibration_seed`` is ``consensusTOM``'s random seed (54321).
    ``reassign_threshold`` is accepted for the options of
    :func:`gcen2.tom_modules.tom_modules` and has no effect (see the module
    docstring). ``consensus_tom_modules(sets, **HDWGCNA_CONSENSUS)`` is
    hdWGCNA's ``ConstructNetwork(consensus = TRUE)`` with the powers of
    ``TestSoftPowersConsensus``.

    Memory: one genes x genes double matrix plus half of one per set.
    """
    from gcen2._build_file import _tom_power

    del reassign_threshold
    mats, genes, set_names, samples = _set_matrices(sets, gene_names, set_names, sample_names)
    n_sets = len(mats)
    if not 0 <= consensus_quantile <= 1:
        raise GCENInputError("consensus_quantile must be between 0 and 1.")
    if calibration not in ("single quantile", "none"):
        raise GCENInputError("calibration must be 'single quantile' or 'none'.")
    if not 0 < calibration_quantile < 1:
        raise GCENInputError("calibration_quantile must be between 0 and 1.")
    _check_module_options(power_method, min_size, merge_cut)
    if isinstance(power, (list, tuple, np.ndarray)):
        if len(power) != n_sets:
            raise GCENInputError("give one power per set (or one for all).")
        set_power = [_tom_power(True, p, "pearson") for p in power]
    else:
        set_power = [_tom_power(True, power, "pearson")] * n_sets
    custom_pick = (power_method not in (None, "pearson") or power_r2 is not None
                   or min_power or powers is not None)

    # goodGenesMS: genes finite and not constant in every set
    good = np.ones(len(genes), dtype=bool)
    for m in mats:
        good &= np.isfinite(m).all(axis=1) & (m.std(axis=1) > 0)
    if good.sum() < 3:
        raise GCENInputError("fewer than 3 genes with finite, non-constant expression in "
                             "every set.")
    xs = [np.ascontiguousarray(m[good]) for m in mats]
    del mats
    n = xs[0].shape[0]
    pick = None
    if custom_pick:
        pick = (power_method or "pearson", power_r2, min_power, powers)
    toms, sets_power = _set_toms(xs, set_names, set_power, pick, signed, threads)
    used_power = [float(m["power"]) for m in sets_power]
    cons, scale_quant, scale_power = _consensus_tom(
        toms, consensus_quantile, calibration, calibration_quantile, sample_for_calibration,
        calibration_sample_factor, calibration_seed, threads)
    del toms
    diss = _fill_diss(cons, n)
    del cons

    if min_size > n / 2:  # checkMinModuleSize
        min_size = n / 2
    merge, height = require_extension().module.hclust_average_inplace(diss)
    with _blas_threads(threads):
        lab = cutree_hybrid(merge, height, diss, min_cluster_size=min_size,
                            deep_split=deep_split, cut_height=detect_cut, pam_stage=pam_stage)
        del diss
        lab = consensus_post_process(xs, lab, signed, min_size, consensus_quantile,
                                     min_core_kme=min_core_kme, min_kme_to_stay=min_kme_to_stay)
        labels = np.zeros(len(genes), dtype=np.int64)
        unmerged = labels.copy()
        unmerged[good] = _relabel_by_size(lab)
        labels[good] = consensus_merge_close_modules(xs, lab, merge_cut, consensus_quantile)
        x_all = np.hstack(xs)
        me_labels, me, var = module_eigengenes(x_all, labels[good])
        set_me = {}
        for name, x, smp in zip(set_names, xs, samples):
            import pandas as pd

            lv, e, _ = module_eigengenes(x, labels[good])
            set_me[name] = pd.DataFrame(e.T, index=smp, columns=[f"ME{k}" for k in lv])
    first = sets_power[0]
    res = TomModules(genes=genes, labels=labels, unmerged=unmerged, power=used_power[0],
                     power_estimated=bool(first["estimated"]), fits=first["fits"],
                     power_rule=first["rule"],
                     samples=[s for smp in samples for s in smp], me_labels=me_labels,
                     eigengenes=me, var_explained=var, signed=signed, _x=x_all, _good=good)
    res.metadata["tree"] = {"merge": np.asarray(merge), "height": np.asarray(height)}
    res.metadata.update(
        consensus=True, sets=set_names, set_powers=dict(zip(set_names, used_power)),
        set_fits={s: m["fits"] for s, m in zip(set_names, sets_power)}, set_samples=dict(zip(set_names, samples)),
        consensus_quantile=consensus_quantile,
        calibration=calibration, calibration_quantiles=dict(zip(set_names, scale_quant)),
        calibration_powers=dict(zip(set_names, scale_power)), set_eigengenes=set_me)
    return res
