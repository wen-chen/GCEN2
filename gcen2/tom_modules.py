"""WGCNA-style modules: TOM, average-linkage clustering, dynamic tree cut.

:func:`tom_modules` follows WGCNA 1.74 ``blockwiseModules`` on one block
(every gene in one TOM, as ``maxBlockSize`` >= the number of genes):

1. 1 - TOM of all genes (the ``tom`` edge score of ``gcen2 build --tom``;
   ``power="auto"`` picks the power as ``pickSoftThreshold``);
2. ``fastcluster::hclust(method = "average")`` on it (C++, same merges and
   heights);
3. ``dynamicTreeCut::cutreeHybrid`` (1.63-1) with ``deepSplit``,
   ``minClusterSize``, ``cutHeight = 0.995`` and the PAM stage
   (:func:`cutree_hybrid`, a line-by-line port);
4. module eigengenes (first principal component of the scaled module
   genes, aligned with their average) and kME checks: a module with fewer
   than ``min_size / 3`` genes at kME > 0.5 is dropped, genes with kME < 0.3
   leave their module, and a module left with fewer than ``min_size`` genes
   is dropped;
5. genes whose kME is much higher for another module (Fisher p-value ratio
   < 1e-6) are reassigned to it (as WGCNA does, by the eigengene's column
   number);
6. ``mergeCloseModules``: modules whose eigengenes correlate above
   ``1 - merge_cut`` (average linkage on 1 - cor) are merged, repeatedly.

Eigengenes and kME use the Pearson correlation of the expression values
(WGCNA's ``corType = "pearson"``), also when the TOM is from Spearman
correlations. Genes with a constant or non-finite row are left out (label 0),
as WGCNA's ``goodSamplesGenes`` removes them; missing values are not imputed.
Modules are numbered 1, 2, ... by decreasing size, 0 = unassigned.

Downstream (as the WGCNA tutorial): :attr:`TomModules.eigengenes`
(``moduleEigengenes`` of the merged modules, ME0 for the unassigned genes
included, with ``varExplained``), :meth:`TomModules.kme` (``signedKME`` and
``corPvalueStudent``) and :meth:`TomModules.trait_cor` (``corAndPvalue`` of
the eigengenes with sample traits, pairwise complete observations).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from gcen2.core import require_extension
from gcen2.exceptions import GCENInputError

if TYPE_CHECKING:
    import pandas as pd

# hdWGCNA 0.4.12 ConstructNetwork defaults (blockwiseConsensusModules on one
# set: deepSplit 4, no PAM stage, minKMEtoStay 0.2; its reassignment compares
# lower-tail chi-square p-values and never moves a gene) with the power of
# TestSoftPowers (pickSoftThreshold with bicor, powers 1..10, 12..30; the
# lowest power above 3 with R^2 >= 0.8)
HDWGCNA = dict(signed=True, min_size=50, deep_split=4, merge_cut=0.2, pam_stage=False,
               min_kme_to_stay=0.2, reassign_threshold=0.0, power_method="bicor",
               power_r2=0.8, min_power=3,
               powers=tuple(range(1, 11)) + tuple(range(12, 31, 2)))

# dynamicTreeCut defaults (deepSplit 0..4)
_DEF_MCS = (0.64, 0.73, 0.82, 0.91, 0.95)


@dataclass
class TomModules:
    """Result of :func:`tom_modules`."""

    genes: list[str]
    labels: np.ndarray          # merged modules: 0 unassigned, 1..K by size
    unmerged: np.ndarray        # before mergeCloseModules (after the kME steps)
    power: float
    power_estimated: bool       # power="auto" reached scale-free fit R^2 > 0.85
    fits: list[dict] = field(default_factory=list)
    power_rule: str | None = None   # how power="auto" picked it, if not pickSoftThreshold's
    samples: list[str] = field(default_factory=list)
    me_labels: list[int] = field(default_factory=list)   # eigengene k: ME<k>
    eigengenes: np.ndarray | None = None                 # modules x samples
    var_explained: np.ndarray | None = None
    signed: bool = False
    metadata: dict = field(default_factory=dict)        # e.g. cell type, cells, metacells
    _x: np.ndarray | None = field(default=None, repr=False)     # good genes x samples
    _good: np.ndarray | None = field(default=None, repr=False)

    def modules(self) -> list[list[str]]:
        """Genes of modules 1..K (sorted names)."""
        out = []
        for k in range(1, int(self.labels.max(initial=0)) + 1):
            out.append(sorted(g for g, lab in zip(self.genes, self.labels) if lab == k))
        return out

    def hub_genes(self, n: int = 10):
        """The ``n`` genes of each module with the highest kME for it (``|kME|``
        for unsigned networks; as hdWGCNA ``GetHubGenes``): a table module
        (0-based, as :func:`write_module_files`), gene, kME, rank."""
        import pandas as pd

        kme, _ = self.kme()
        col = {k: i for i, k in enumerate(self.me_labels)}
        rows = []
        for k in range(1, int(self.labels.max(initial=0)) + 1):
            idx = np.flatnonzero(self.labels == k)
            v = kme[idx, col[k]]
            score = v if self.signed else np.abs(v)
            order = idx[np.lexsort((np.array(self.genes, dtype=object)[idx], -score))][:n]
            for r, g in enumerate(order, start=1):
                rows.append((k - 1, self.genes[g], kme[g, col[k]], r))
        return pd.DataFrame(rows, columns=["module", "gene", "kME", "rank"])

    def eigengene_frame(self):
        """Eigengenes as a samples x ME<k> DataFrame (WGCNA ``MEs``)."""
        import pandas as pd

        return pd.DataFrame(self.eigengenes.T, index=self.samples,
                            columns=[f"ME{k}" for k in self.me_labels])

    def kme(self) -> tuple[np.ndarray, np.ndarray]:
        """WGCNA ``signedKME`` (Pearson correlation of every gene with every
        eigengene; genes x eigengenes, NaN for the genes left out) and its
        ``corPvalueStudent`` p-values."""
        n = len(self.genes)
        kme = np.full((n, len(self.me_labels)), np.nan)
        if self.me_labels:
            kme[self._good] = _cor(self._x, self.eigengenes)
        return kme, cor_p_student(kme, len(self.samples))

    def trait_cor(self, traits) -> pd.DataFrame:
        """Eigengene-trait correlations as WGCNA ``corAndPvalue(MEs, traits,
        use = "p")``: a long table module, trait, cor, p, n. ``traits``: a
        DataFrame with samples as rows (matched by name; samples missing from
        it count as missing values) and numeric columns."""
        import pandas as pd

        if not isinstance(traits, pd.DataFrame):
            raise GCENInputError("traits must be a DataFrame with samples as rows.")
        traits = traits.copy()
        traits.index = traits.index.map(str)
        if not traits.index.is_unique:
            raise GCENInputError("trait sample names must be unique.")
        if not traits.index.isin(self.samples).any():
            raise GCENInputError("no trait sample matches an expression sample.")
        bad = [c for c in traits.columns if not pd.api.types.is_numeric_dtype(traits[c])]
        if bad:
            raise GCENInputError(f"traits must be numeric (not: {', '.join(map(str, bad))}).")
        t = traits.reindex(self.samples).to_numpy(dtype=np.float64).T
        r, p, nobs = cor_pairwise(self.eigengenes, t)
        k, m = r.shape
        return pd.DataFrame({
            "module": np.repeat([f"ME{x}" for x in self.me_labels], m),
            "trait": np.tile([str(c) for c in traits.columns], k),
            "cor": r.ravel(), "p": p.ravel(), "n": nobs.ravel()})


def _blas_threads(n: int):
    """Limit numpy's BLAS to n threads (threadpoolctl, a dependency; skipped
    if it is missing): its default of one thread per logical CPU makes the
    small products here many times slower on CPUs with many cores. scipy
    ships its own OpenBLAS: it is loaded first, as the limit applies to the
    libraries loaded when it is set."""
    try:
        from threadpoolctl import threadpool_limits
    except ImportError:
        import contextlib

        return contextlib.nullcontext()
    import scipy.linalg  # noqa: F401

    return threadpool_limits(limits=max(1, int(n)), user_api="blas")


def _r_round(x: float) -> int:
    return int(round(x))  # both round half to even


def _interpolate(data, index):
    i = _r_round(index)
    n = len(data)
    if i < 1:
        return data[0]
    if i >= n:
        return data[n - 1]
    r = index - i
    return data[i - 1] * (1 - r) + data[i] * r


def _core_size(branch_size: int, min_cluster_size: int) -> int:
    base = min_cluster_size / 2 + 1
    if base < branch_size:
        return int(base + math.sqrt(branch_size - base))
    return branch_size


def _core_scatter(dist, core) -> float:
    sub = np.asarray(dist[np.ix_(core, core)], dtype=np.float64)
    return float(np.mean(sub.sum(axis=0) / (len(core) - 1)))


def cutree_hybrid(merge, height, dist, min_cluster_size: int = 20, deep_split: int = 1,
                  cut_height: float | None = None, pam_stage: bool = True,
                  pam_respects_dendro: bool = True,
                  max_pam_dist: float | None = None) -> np.ndarray:
    """``dynamicTreeCut::cutreeHybrid(...)$labels`` (1.63-1; no external
    split functions, ``useMedoids = FALSE``, ``respectSmallClusters = TRUE``).

    ``merge`` / ``height``: an R-convention dendrogram (as
    ``_gcen_core.module.hclust_average`` returns); ``dist``: the n x n
    dissimilarities. ``max_pam_dist`` defaults to ``cut_height`` as given
    (``cutreeDynamic`` passes it before it is capped at the tree height).
    Returns labels 0 (unassigned), 1, 2, ... by decreasing cluster size.
    """
    merge = np.asarray(merge)
    height = np.asarray(height, dtype=np.float64)
    n_merge = len(height)
    n_points = n_merge + 1
    if max_pam_dist is None:
        max_pam_dist = cut_height
    ref_merge = max(1, _r_round(n_merge * 0.05))
    ref_height = height[ref_merge - 1]
    if cut_height is None:
        cut_height = 0.99 * (height.max() - ref_height) + ref_height
    elif cut_height > height.max():
        cut_height = float(height.max())
    if max_pam_dist is None:
        max_pam_dist = cut_height
    if int((height <= cut_height).sum()) < min_cluster_size:
        return np.zeros(n_points, dtype=np.int64)

    ds = deep_split + 1
    if not 1 <= ds <= len(_DEF_MCS):
        raise GCENInputError(f"deep_split must be 0..{len(_DEF_MCS) - 1}.")
    max_core_scatter = _interpolate(_DEF_MCS, ds)
    min_gap = _interpolate([(1 - x) * 3 / 4 for x in _DEF_MCS], ds)
    max_abs_core_scatter = ref_height + max_core_scatter * (cut_height - ref_height)
    min_abs_gap = min_gap * (cut_height - ref_height)
    min_abs_split_height = ref_height  # minSplitHeight 0

    # branches 1.. (index 0 unused), as in R
    is_basic = [True]
    is_top_basic = [True]
    fail_size = [False]
    size = [2]
    n_singletons = [2]
    attach_height = [math.nan]
    singletons: list[list[int]] = [[]]
    basic_clusters: list[list[int]] = [[]]
    n_branches = 0

    def new_branch():
        is_basic.append(True)
        is_top_basic.append(True)
        fail_size.append(False)
        size.append(2)
        n_singletons.append(2)
        attach_height.append(math.nan)
        singletons.append([])
        basic_clusters.append([])
        return len(is_basic) - 1

    ind_merge_to_branch = np.zeros(n_merge + 1, dtype=np.int64)  # by merge step 1..
    on_branch = np.zeros(n_points, dtype=np.int64)

    for m in range(1, n_merge + 1):
        h = height[m - 1]
        if not h <= cut_height:
            continue
        a, b = int(merge[m - 1, 0]), int(merge[m - 1, 1])
        if a < 0 and b < 0:
            n_branches = new_branch()
            singletons[n_branches] = [-a - 1, -b - 1]
            ind_merge_to_branch[m] = n_branches
        elif (a < 0) != (b < 0):
            clust = ind_merge_to_branch[max(a, b)]
            gene = -min(a, b) - 1
            if is_basic[clust]:
                singletons[clust].append(gene)
            else:
                on_branch[gene] = clust
            size[clust] += 1
            n_singletons[clust] += 1
            ind_merge_to_branch[m] = clust
        else:
            c1, c2 = ind_merge_to_branch[a], ind_merge_to_branch[b]
            small, large = (c1, c2) if size[c1] <= size[c2] else (c2, c1)

            def ave_dist(br):
                if not is_basic[br]:
                    return 0.0
                cs = _core_size(n_singletons[br], min_cluster_size)
                return _core_scatter(dist, singletons[br][:cs])

            sm_ave, lg_ave = ave_dist(small), ave_dist(large)
            sm = (size[small] < min_cluster_size, sm_ave > max_abs_core_scatter,
                  h - sm_ave < min_abs_gap, h < min_abs_split_height)
            if is_basic[small] and any(sm):
                do_merge = True
                smaller_fail_size = not (sm[1] or sm[2])
            else:
                lg = (size[large] < min_cluster_size, lg_ave > max_abs_core_scatter,
                      h - lg_ave < min_abs_gap, h < min_abs_split_height)
                if is_basic[large] and any(lg):
                    do_merge = True
                    smaller_fail_size = not (lg[1] or lg[2])
                    small, large = large, small
                else:
                    do_merge = False
            if do_merge:
                fail_size[small] = smaller_fail_size
                attach_height[small] = h
                is_top_basic[small] = False
                if is_basic[large]:
                    singletons[large].extend(singletons[small][:n_singletons[small]])
                    n_singletons[large] += n_singletons[small]
                else:
                    if not is_basic[small]:
                        raise RuntimeError("cutree_hybrid: merging two composite clusters")
                    on_branch[singletons[small]] = large
                size[large] += size[small]
                ind_merge_to_branch[m] = large
            else:
                if is_basic[large] and not is_basic[small]:
                    small, large = large, small
                if is_basic[large] or (pam_stage and pam_respects_dendro):
                    nb = new_branch()
                    n_branches = nb
                    attach_height[large] = attach_height[small] = h
                    add = [small] if is_basic[small] else list(basic_clusters[small])
                    add += [large] if is_basic[large] else list(basic_clusters[large])
                    is_basic[nb] = False
                    is_top_basic[nb] = False
                    basic_clusters[nb] = add
                    size[nb] = size[small] + size[large]
                    ind_merge_to_branch[m] = nb
                else:
                    add = [small] if is_basic[small] else list(basic_clusters[small])
                    basic_clusters[large].extend(add)
                    size[large] += size[small]
                    attach_height[small] = h
                    ind_merge_to_branch[m] = large

    # clusters and small labels
    is_cluster = [False] * (n_branches + 1)
    small_labels = np.zeros(n_points, dtype=np.int64)
    for c in range(1, n_branches + 1):
        if math.isnan(attach_height[c]):
            attach_height[c] = cut_height
        if is_top_basic[c]:
            cs = _core_size(n_singletons[c], min_cluster_size)
            scatter = _core_scatter(dist, singletons[c][:cs])
            is_cluster[c] = (size[c] >= min_cluster_size and scatter < max_abs_core_scatter
                             and attach_height[c] - scatter > min_abs_gap)
        if fail_size[c]:
            small_labels[singletons[c]] = c
    colors = np.zeros(n_points, dtype=np.int64)
    branch_labels = np.zeros(n_branches + 1, dtype=np.int64)
    color = 0
    for c in range(1, n_branches + 1):
        if is_cluster[c]:
            color += 1
            colors[singletons[c]] = color
            small_labels[singletons[c]] = 0
            branch_labels[c] = color
    n_proper = color

    if pam_stage and (colors == 0).any() and n_proper > 0:
        diam = np.zeros(n_proper + 1)
        for k in range(1, n_proper + 1):
            members = np.flatnonzero(colors == k)
            if len(members) > 1:
                sub = np.asarray(dist[np.ix_(members, members)], dtype=np.float64)
                diam[k] = np.max(sub.sum(axis=0) / (len(members) - 1))
        colors_x = colors.copy()

        def nearest(dist_rows_mean, use):
            """Mean distance per label of the objects `use` -> (label, dist)."""
            labs = colors_x[use]
            levels = np.unique(labs)
            means = np.array([dist_rows_mean[labs == lv].mean() for lv in levels])
            k = int(np.argmin(means))  # which.min: first minimum
            return int(levels[k]), float(means[k])

        def assignable(label, d):
            # R would stop on label 0 (ClusterDiam[0] is empty): leave unassigned
            return label > 0 and (d < diam[label] or d < max_pam_dist)

        if not pam_respects_dendro:
            raise GCENInputError("only pam_respects_dendro=True is supported.")
        for sclust in np.unique(small_labels[small_labels != 0]):
            in_cluster = np.flatnonzero(small_labels == sclust)
            on_br = np.unique(on_branch[in_cluster])
            if len(on_br) > 1:
                raise RuntimeError("cutree_hybrid: a small cluster lies on several branches")
            if on_br[0] > 0:
                labels_on = np.unique(branch_labels[basic_clusters[on_br[0]]])
                use = np.flatnonzero(np.isin(colors_x, labels_on))
                mean_dist = np.asarray(dist[np.ix_(in_cluster, use)], np.float64).mean(axis=0)
                label, d = nearest(mean_dist, use)
                colors[in_cluster] = label if assignable(label, d) else -1
        for obj in np.flatnonzero(colors == 0):
            br = on_branch[obj]
            if br <= 0:
                continue
            labels_on = np.unique(branch_labels[basic_clusters[br]])
            use = np.flatnonzero(np.isin(colors_x, labels_on))
            label, d = nearest(np.asarray(dist[use, obj], np.float64), use)
            if assignable(label, d):
                colors[obj] = label
    colors[colors < 0] = 0
    return _relabel_by_size(colors)


def _relabel_by_size(colors: np.ndarray) -> np.ndarray:
    """Labels 1.. by decreasing size (ties: smaller old label first), 0 stays 0."""
    out = np.zeros_like(colors)
    levels, counts = np.unique(colors[colors != 0], return_counts=True)
    order = np.lexsort((levels, -counts))
    for new, k in enumerate(order, start=1):
        out[colors == levels[k]] = new
    return out


def module_eigengene(x: np.ndarray) -> np.ndarray:
    """WGCNA moduleEigengenes for one module: genes x samples -> samples.

    The first right singular vector of the scaled genes, from the top
    eigenvector of the smaller Gram matrix only (same vector as the SVD up
    to sign and rounding), aligned with the average scaled expression.
    """
    from scipy.linalg import eigh

    z = (x - x.mean(axis=1, keepdims=True)) / x.std(axis=1, ddof=1, keepdims=True)
    g, n = z.shape
    if g <= n:
        _, u = eigh(z @ z.T, subset_by_index=[g - 1, g - 1])
        pc = z.T @ u[:, 0]
        pc /= np.linalg.norm(pc)
    else:
        _, v = eigh(z.T @ z, subset_by_index=[n - 1, n - 1])
        pc = v[:, 0]
    c = np.corrcoef(z.mean(axis=0), pc)[0, 1]
    return -pc if np.isfinite(c) and c < 0 else pc


def _cor(x: np.ndarray, y: np.ndarray) -> np.ndarray:
    """Pearson correlation of the rows of x (genes x samples) with the rows
    of y: x.shape[0] x y.shape[0]."""
    xc = x - x.mean(axis=1, keepdims=True)
    yc = y - y.mean(axis=1, keepdims=True)
    xc /= np.linalg.norm(xc, axis=1, keepdims=True)
    yc /= np.linalg.norm(yc, axis=1, keepdims=True)
    return np.clip(xc @ yc.T, -1.0, 1.0)


def _fisher_p(r: np.ndarray, n_samples: int) -> np.ndarray:
    """WGCNA corPvalueFisher (two-sided)."""
    from scipy.special import erfc

    with np.errstate(divide="ignore"):
        z = np.abs(0.5 * np.log((1 + r) / (1 - r)) * math.sqrt(n_samples - 3))
    return erfc(z / math.sqrt(2.0))


def cor_p_student(r: np.ndarray, n) -> np.ndarray:
    """WGCNA corPvalueStudent: two-sided Student t p-value of correlations
    from n observations (n broadcasts with r)."""
    from scipy.special import stdtr

    r = np.asarray(r, dtype=np.float64)
    df = np.asarray(n, dtype=np.float64) - 2
    with np.errstate(divide="ignore", invalid="ignore"):
        t = np.abs(r) * np.sqrt(df) / np.sqrt(1 - r * r)
        return 2 * stdtr(df, -t)


def cor_pairwise(x: np.ndarray, y: np.ndarray):
    """Pearson correlation of the rows of x with the rows of y on pairwise
    complete observations (R ``cor(use = "p")``), with WGCNA corAndPvalue's
    Student p-values and observation counts."""
    k, m = x.shape[0], y.shape[0]
    r = np.full((k, m), np.nan)
    nobs = np.zeros((k, m), dtype=np.int64)
    fx = np.isfinite(x)
    for j in range(m):
        for i in range(k):
            ok = fx[i] & np.isfinite(y[j])
            nobs[i, j] = ok.sum()
            if nobs[i, j] < 2:
                continue
            a = x[i, ok] - x[i, ok].mean()
            b = y[j, ok] - y[j, ok].mean()
            den = math.sqrt((a @ a) * (b @ b))
            if den > 0:
                r[i, j] = min(1.0, max(-1.0, (a @ b) / den))
    return r, cor_p_student(r, nobs), nobs


def module_eigengenes(x: np.ndarray, labels: np.ndarray):
    """WGCNA moduleEigengenes of every label present (0 included):
    (levels, eigengenes as levels x samples, varExplained)."""
    levels = [int(k) for k in np.unique(labels)]
    me = _eigengenes(x, labels, levels)
    var = np.array([np.mean(_cor(x[labels == k], me[i:i + 1]) ** 2)
                    for i, k in enumerate(levels)])
    return levels, me, var


def _eigengenes(x: np.ndarray, labels: np.ndarray, levels) -> np.ndarray:
    return np.array([module_eigengene(x[labels == k]) for k in levels])


def _cutree_height(merge: np.ndarray, height: np.ndarray, h: float) -> np.ndarray:
    """R cutree(tree, h = h) group of each leaf (0-based group ids)."""
    n = len(height) + 1
    parent = list(range(n))

    def find(i):
        while parent[i] != i:
            parent[i] = parent[parent[i]]
            i = parent[i]
        return i

    rep = {}
    for step, (a, b) in enumerate(merge, start=1):
        ra = find(-a - 1) if a < 0 else rep[a]
        rb = find(-b - 1) if b < 0 else rep[b]
        if height[step - 1] <= h:
            parent[rb] = ra
            rep[step] = ra
        else:
            rep[step] = ra  # above the cut: groups stay apart
    return np.array([find(i) for i in range(n)])


def merge_close_modules(x: np.ndarray, labels: np.ndarray, cut_height: float) -> np.ndarray:
    """WGCNA mergeCloseModules (one set, iterate, relabel): merge modules
    whose eigengenes have 1 - cor below ``cut_height`` (average linkage)."""
    def dissimilarity(labels, levels):
        me = _eigengenes(x, labels, levels)
        return 1.0 - _cor(me, me)

    return _merge_modules(labels, dissimilarity, cut_height)


def _merge_modules(labels: np.ndarray, dissimilarity, cut_height: float) -> np.ndarray:
    """The mergeCloseModules loop: cluster the modules by
    ``dissimilarity(labels, levels)`` (modules x modules; average linkage),
    merge those joined below ``cut_height`` into the first of each group,
    repeat until nothing merges; modules relabelled by size."""
    labels = labels.copy()
    while True:
        levels = np.unique(labels[labels != 0])
        if len(levels) < 2:
            break
        diss = np.ascontiguousarray(dissimilarity(labels, levels))
        np.fill_diagonal(diss, 0.0)
        merge, height = require_extension().module.hclust_average(diss)
        groups = _cutree_height(merge, height, cut_height)
        new = labels.copy()
        for g in np.unique(groups):
            mods = levels[groups == g]
            for k in mods[1:]:
                new[labels == k] = mods[0]
        if len(np.unique(new[new != 0])) == len(levels):
            break
        labels = new
    return _relabel_by_size(labels)


def _check_module_options(power_method, min_size, merge_cut) -> None:
    """The options :func:`tom_modules` and consensus_tom_modules share."""
    if power_method not in (None, "pearson", "spearman", "bicor"):
        raise GCENInputError("power_method must be pearson, spearman or bicor.")
    if min_size < 2:
        raise GCENInputError("min_size must be at least 2.")
    if not 0 <= merge_cut <= 1:
        raise GCENInputError("merge_cut must be between 0 and 1.")


def wgcna_post_process(x: np.ndarray, labels: np.ndarray, signed: bool, min_size: int,
                       min_core_kme: float = 0.5, min_kme_to_stay: float = 0.3,
                       reassign_threshold: float = 1e-6) -> np.ndarray:
    """blockwiseModules after cutreeDynamic, before mergeCloseModules: the
    kME checks and the reassignment (``x``: genes x samples)."""
    labels = labels.copy()
    n_samples = x.shape[1]
    min_core_kme_size = min_size / 3
    levels = [int(k) for k in np.unique(labels[labels != 0])]
    if not levels:
        return labels
    me = _eigengenes(x, labels, levels)
    keep = []
    for i, k in enumerate(levels):
        genes = np.flatnonzero(labels == k)
        kme = _cor(x[genes], me[i:i + 1])[:, 0]
        if not signed:
            kme = np.abs(kme)
        if (kme > min_core_kme).sum() < min_core_kme_size:
            labels[genes] = 0
            continue
        if (kme < min_kme_to_stay).any():
            labels[genes[kme < min_kme_to_stay]] = 0
            if (labels[genes] > 0).sum() < min_size:
                labels[genes] = 0
                continue
        keep.append(i)
    if not (labels > 0).any():
        return labels
    me = me[keep]
    label_index = [levels[i] for i in keep]

    assigned = np.flatnonzero(labels != 0)
    prop = labels[assigned].copy()
    kme = _cor(x[assigned], me)
    if not signed:
        kme = np.abs(kme)
    delete = []
    for mod, k in enumerate(label_index):
        mod_genes = np.flatnonzero(prop == k)
        kme_mod = kme[mod_genes, mod]
        kme_best = kme[mod_genes].max(axis=1)
        cand = kme_mod < kme_best
        if not cand.any():
            continue
        p_mod = _fisher_p(kme_mod[cand], n_samples)
        which_best = kme[mod_genes[cand]].argmax(axis=1) + 1  # column number, as WGCNA
        p_best = _fisher_p(kme_best[cand], n_samples)
        with np.errstate(divide="ignore", invalid="ignore"):
            ratio = p_best / p_mod
        reassign = np.isfinite(ratio) & (ratio < reassign_threshold)
        if reassign.any():
            labels[assigned[mod_genes[cand][reassign]]] = which_best[reassign]
            if len(mod_genes) - reassign.sum() < min_size:
                delete.append(k)
    if delete:
        labels[np.isin(labels, delete)] = 0
    return labels


def tom_modules(expr, method: str = "pearson", signed: bool = False, power="auto",
                min_size: int = 20, deep_split: int = 2, merge_cut: float = 0.15,
                detect_cut: float = 0.995, pam_stage: bool = True,
                min_kme_to_stay: float = 0.3, reassign_threshold: float = 1e-6,
                power_method: str | None = None, power_r2: float | None = None,
                min_power: float = 0, powers=None,
                log: str | None = None, threads: int = 2,
                gene_names=None, layer=None, use_highly_variable: bool = False,
                sample_names=None) -> TomModules:
    """WGCNA ``blockwiseModules`` modules (see the module docstring).

    ``expr``: genes x samples (DataFrame, ExpressionMatrix, AnnData, or an
    array with ``gene_names``). ``signed``: ((1 + r) / 2)^power adjacency and
    signed kME (WGCNA ``networkType = "signed"``, ``TOMType = "signed"``),
    else ``|r|^power``. Defaults are blockwiseModules' (``minModuleSize`` 20,
    ``deepSplit`` 2, ``mergeCutHeight`` 0.15, ``detectCutHeight`` 0.995).
    Memory: one genes x genes double matrix (the correlations, the
    adjacency and 1 - TOM are built in it; the clustering works in it too)
    plus about 1/16 of that while the TOM is built.
    ``sample_names`` (array input): names of the columns, else taken from
    the DataFrame / AnnData, or sample_1, sample_2, ...

    ``power="auto"`` picks the power as pickSoftThreshold (first power with
    R^2 > 0.85, Pearson / Spearman as ``method``). ``power_method``
    ("pearson", "spearman" or "bicor"), ``power_r2``, ``min_power`` and
    ``powers`` change that: the lowest of ``powers`` above ``min_power``
    whose fit R^2 is at least ``power_r2``, from correlations by
    ``power_method``; if none, the WGCNA FAQ power. ``tom_modules(expr,
    **HDWGCNA)`` reproduces hdWGCNA's ConstructNetwork with its
    TestSoftPowers power on the same datExpr. ``pam_stage``,
    ``min_kme_to_stay`` and ``reassign_threshold`` (0: none) are the
    cutreeDynamic / blockwiseModules options of those names.
    """
    from gcen2._build_file import _tom_metadata, _tom_power
    from gcen2.network import _prepare_matrix

    samples = _sample_names(expr, sample_names)
    matrix, genes = _prepare_matrix(expr, method, log, gene_names, layer=layer,
                                    use_highly_variable=use_highly_variable)
    if samples is None:
        samples = [f"sample_{i + 1}" for i in range(matrix.shape[1])]
    elif len(samples) != matrix.shape[1]:
        raise GCENInputError("sample_names must give one name per column of expr.")
    p = _tom_power(True, power, method)
    custom_pick = p == 0 and (power_method not in (None, method) or power_r2 is not None
                              or min_power or powers is not None)
    _check_module_options(power_method, min_size, merge_cut)
    genes = list(genes)
    sd = matrix.std(axis=1)
    good = np.isfinite(matrix).all(axis=1) & (sd > 0)
    labels = np.zeros(len(genes), dtype=np.int64)
    unmerged = labels.copy()
    if good.sum() < 3:
        raise GCENInputError("fewer than 3 genes with finite, non-constant expression.")
    x = np.ascontiguousarray(matrix[good])
    core = require_extension()
    rule = None
    if custom_pick:
        p, meta, rule = _pick_power(x, power_method or method, signed, power_r2, min_power,
                                    powers, threads)
    diss, info = core.module.tom_dissimilarity(x, method, signed, p, [], threads)
    if not custom_pick:
        meta = _tom_metadata(info)
    merge, height = core.module.hclust_average_inplace(diss)
    with _blas_threads(threads):
        lab = cutree_hybrid(merge, height, diss, min_cluster_size=min_size,
                            deep_split=deep_split, cut_height=detect_cut,
                            pam_stage=pam_stage)
        del diss
        lab = wgcna_post_process(x, lab, signed, min_size, min_kme_to_stay=min_kme_to_stay,
                                 reassign_threshold=reassign_threshold)
        unmerged[good] = _relabel_by_size(lab)
        labels[good] = merge_close_modules(x, lab, merge_cut)
        me_labels, me, var = module_eigengenes(x, labels[good])
    res = TomModules(genes=genes, labels=labels, unmerged=unmerged, power=meta["power"],
                     power_estimated=meta["estimated"], fits=meta["fits"], power_rule=rule,
                     samples=samples,
                     me_labels=me_labels, eigengenes=me, var_explained=var, signed=signed,
                     _x=x, _good=good)
    # the gene tree (R's hclust convention, leaves = the good genes in order) and
    # how the TOM was made, for figures (dendrogram, hub networks, gene UMAP)
    res.metadata["tree"] = {"merge": np.asarray(merge), "height": np.asarray(height)}
    res.metadata["tom"] = {"method": method, "signed": bool(signed), "power": float(meta["power"])}
    return res


def soft_power_table(res: TomModules):
    """The scale-free fit of every soft power tried (pickSoftThreshold's
    table: power, r2, slope, mean_k, median_k, max_k), with the power used
    marked; a ``set`` column for consensus networks. Empty when the power was
    given."""
    import pandas as pd

    md = res.metadata
    if md.get("consensus"):
        frames = []
        for s in md["sets"]:
            f = pd.DataFrame(md["set_fits"][s])
            if f.empty:
                continue
            f.insert(0, "set", s)
            f["selected"] = f["power"] == md["set_powers"][s]
            frames.append(f)
        return pd.concat(frames, ignore_index=True) if frames else pd.DataFrame()
    f = pd.DataFrame(res.fits)
    if not f.empty:
        f["selected"] = f["power"] == res.power
    return f


def dendrogram_dict(res: TomModules) -> dict | None:
    """The gene tree as a JSON-ready dict: genes (the leaves, in the order of
    the merge numbers), merge (R hclust convention), height, and each leaf's
    module after and before merging (0-based as modules.tsv, -1 unassigned)."""
    tree = res.metadata.get("tree")
    if tree is None:
        return None
    good = np.flatnonzero(res._good)
    return {"genes": [res.genes[i] for i in good],
            "merge": np.asarray(tree["merge"]).astype(int).tolist(),
            "height": [float(h) for h in tree["height"]],
            "modules": (res.labels[good] - 1).astype(int).tolist(),
            "unmerged": (res.unmerged[good] - 1).astype(int).tolist()}


def hub_tom(res: TomModules, n_hubs: int = 25, threads: int = 2):
    """The TOM between every gene in a module (rows) and the ``n_hubs`` genes
    of highest kME of each module (columns), recomputed from the result's
    expression with its power: the input of the hub-gene networks and the
    module UMAP. Not for consensus networks (their TOM is a combination of
    the sets' TOMs)."""
    import pandas as pd

    md = res.metadata
    if md.get("consensus") or "tom" not in md:
        raise GCENInputError("the hub TOM is available for single networks only.")
    t = md["tom"]
    tom, _ = require_extension().module.tom_dissimilarity(
        res._x, t["method"], t["signed"], t["power"], [], threads, similarity=True)
    good = np.flatnonzero(res._good)
    lab = res.labels[good]
    kme, _ = res.kme()
    col = {k: i for i, k in enumerate(res.me_labels)}
    rows = np.flatnonzero(lab > 0)
    hubs = []
    for k in sorted(set(lab[lab > 0].tolist())):
        idx = np.flatnonzero(lab == k)
        score = kme[good[idx], col[k]]
        order = idx[np.lexsort((np.array([res.genes[good[i]] for i in idx]), -score))]
        hubs.extend(order[:n_hubs].tolist())
    names = [res.genes[i] for i in good]
    out = pd.DataFrame(tom[np.ix_(rows, hubs)], index=[names[i] for i in rows],
                       columns=[names[i] for i in hubs])
    out.index.name = "gene"
    return out


def gene_significance(res: TomModules, traits):
    """WGCNA gene significance: correlation of every gene with every trait
    (pairwise complete observations) and its Student p-value; genes x
    (GS.<trait>, p.GS.<trait>). ``traits``: samples x numeric traits."""
    import pandas as pd

    if not isinstance(traits, pd.DataFrame):
        raise GCENInputError("traits must be a DataFrame with samples as rows.")
    traits = traits.copy()
    traits.index = traits.index.map(str)
    t = traits.reindex(res.samples).to_numpy(dtype=np.float64).T
    gs = np.full((len(res.genes), t.shape[0]), np.nan)
    p = np.full_like(gs, np.nan)
    good = np.flatnonzero(res._good)
    r, pv, _ = cor_pairwise(res._x, t)
    gs[good], p[good] = r, pv
    cols = [str(c) for c in traits.columns]
    return pd.DataFrame(np.hstack([gs, p]), index=pd.Index(res.genes, name="gene"),
                        columns=[f"GS.{c}" for c in cols] + [f"p.GS.{c}" for c in cols])


def write_module_files(res: TomModules, modules=None, eigengenes=None, kme=None, hubs=None,
                       n_hubs: int = 10, trait_table=None, trait_cor=None, soft_power=None,
                       dendrogram=None, tom_hubs=None, n_tom_hubs: int = 25,
                       gs_table=None, gs=None, threads: int = 2) -> None:
    """Write :func:`tom_modules` results as TSV files (each path optional).

    Modules are numbered 0..K-1 by size (``gcen2 module`` format: module,
    gene; unassigned genes left out); eigengene ``ME<k>`` is module k and
    ``ME_unassigned`` WGCNA's ME0 (the unassigned genes). ``eigengenes``:
    samples x eigengenes; ``kme``: gene, module (NA = unassigned), kME and
    Student p-value for every eigengene; ``hubs``: :meth:`TomModules.hub_genes`;
    ``trait_cor``: ``trait_table`` (:meth:`TomModules.trait_cor`) renamed.
    """
    import pandas as pd

    fmt = "%.10g"
    names = [f"ME{k - 1}" if k else "ME_unassigned" for k in res.me_labels]
    if modules is not None:
        from gcen2._matrix_io import write_modules

        write_modules(modules, res.modules())
    if eigengenes is not None:
        pd.DataFrame(res.eigengenes.T, index=pd.Index(res.samples, name="sample"),
                     columns=names).to_csv(eigengenes, sep="\t", float_format=fmt, lineterminator="\n")
    if kme is not None:
        k, p = res.kme()
        df = pd.DataFrame(np.hstack([k, p]), index=pd.Index(res.genes, name="gene"),
                          columns=[f"k{n}" for n in names] + [f"p.k{n}" for n in names])
        df.insert(0, "module", [str(v - 1) if v else "NA" for v in res.labels])
        df.to_csv(kme, sep="\t", float_format=fmt, na_rep="NA", lineterminator="\n")
    if hubs is not None:
        res.hub_genes(n_hubs).to_csv(hubs, sep="\t", index=False, float_format=fmt, lineterminator="\n")
    if trait_cor is not None:
        rename = {f"ME{k}": n for k, n in zip(res.me_labels, names)}
        trait_table.assign(module=trait_table["module"].map(rename)).to_csv(
            trait_cor, sep="\t", index=False, float_format=fmt, na_rep="NA", lineterminator="\n")
    if soft_power is not None:
        sp = soft_power_table(res)
        if not sp.empty:
            sp.to_csv(soft_power, sep="\t", index=False, float_format=fmt, lineterminator="\n")
    if dendrogram is not None:
        import json

        d = dendrogram_dict(res)
        if d is not None:
            with open(dendrogram, "w", encoding="utf-8", newline="\n") as fh:
                json.dump(d, fh)
    if tom_hubs is not None:
        hub_tom(res, n_tom_hubs, threads).to_csv(tom_hubs, sep="\t", float_format=fmt, lineterminator="\n")
    if gs is not None:
        gs_table.to_csv(gs, sep="\t", float_format=fmt, na_rep="NA", lineterminator="\n")


def _pick_power(x, method, signed, r2_cut, min_power, powers, threads):
    """(power, metadata, rule text) for a non-default power rule."""
    from gcen2._build_file import _FIT_COLUMNS

    r2_cut = 0.85 if r2_cut is None else float(r2_cut)
    pw = [float(v) for v in (powers if powers is not None else
                             list(range(1, 11)) + list(range(12, 21, 2)))]
    if not pw or min(pw) <= 0:
        raise GCENInputError("powers must be positive.")
    fits = require_extension().correlation.soft_threshold(x, method, signed, pw, threads)
    ok = [f[0] for f in fits if f[0] > min_power and f[1] >= r2_cut]
    n = x.shape[1]
    if ok:
        p, est = min(ok), True
    else:
        u = 9 if n < 20 else 8 if n < 30 else 7 if n < 40 else 6
        p, est = float(2 * u if signed else u), False
        import warnings

        warnings.warn(f"no soft-threshold power above {min_power:g} reached R^2 >= "
                      f"{r2_cut:g}; using the WGCNA FAQ power {p:g}", stacklevel=3)
    rule = f"{method}, lowest power" + (f" above {min_power:g}" if min_power else "") \
        + f" with R^2 >= {r2_cut:g}"
    meta = {"power": p, "estimated": est, "fits": [dict(zip(_FIT_COLUMNS, f)) for f in fits]}
    return p, meta, rule


def _sample_names(expr, sample_names):
    if sample_names is not None:
        return [str(s) for s in sample_names]
    from gcen2._anndata import is_anndata

    inner = getattr(expr, "data", None) if type(expr).__name__ == "ExpressionMatrix" else expr
    if is_anndata(inner):
        return [str(s) for s in inner.obs_names]
    cols = getattr(inner, "columns", None)
    return None if cols is None else [str(s) for s in cols]
