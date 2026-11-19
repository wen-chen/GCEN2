"""Harmony (Korsunsky et al. 2019) batch correction of an embedding, as the R
package harmony 2.0.5 runs it (``RunHarmony`` defaults), in double precision.

hdWGCNA's harmonized module eigengenes (``ModuleEigengenes(group.by.vars =
...)``) run Harmony on each module's principal components; see
:func:`gcen2.single_cell.cell_eigengenes`.

The algorithm and every default follow harmony 2.0.5: cosine-normalized
cells; ``nclust = min(round(N / 30), 100)`` soft clusters, initialized by
harmony's k-means (random distance-weighted seeding, then ten Lloyd steps);
each of at most ``max_iter`` (10) rounds runs ``max_iter_cluster`` (4) passes
of the diversity-penalized soft clustering (``theta`` 2, ``sigma`` 0.1;
cells updated in random blocks of 5%) and then the mixture-of-experts ridge
correction (``lambda`` per cluster = ``alpha`` (0.2) x the expected number of
cells of each batch; batches with less than 1e-5 of a cluster left out of
its correction); it stops when the objective improves by less than
``epsilon_harmony`` (0.01). R's harmony computes in single precision and
draws its random numbers from R's generator, so results agree with it
closely but not exactly; ``seed`` makes this version reproducible.
"""

from __future__ import annotations

import math

import numpy as np

from gcen2.exceptions import GCENInputError


def _normalize_cols(m: np.ndarray) -> np.ndarray:
    n = np.linalg.norm(m, axis=0)
    n[n == 0] = 1.0
    return m / n


def _kmeans_centers(x: np.ndarray, k: int, rng: np.random.Generator) -> np.ndarray:
    """harmony's kmeans_centers: seeds drawn with probability falling with
    the cosine distance to random cells, then ten single Lloyd steps
    (Armadillo kmeans, keep_existing)."""
    n = x.shape[1]
    idx = np.floor(rng.random(k) * (n - 1)).astype(np.int64)
    y = x[:, idx].copy()
    used = set()
    for i in range(k):
        dist = np.abs(2 * (1 - y[:, i] @ x))
        with np.errstate(divide="ignore"):
            prob = -np.log(rng.random(n)) / dist
        j = int(np.argmin(prob))
        while j in used:
            prob[j] = prob.max()
            j = int(np.argmin(prob))
        used.add(j)
        y[:, i] = x[:, j]
    for _ in range(10):
        d2 = (y * y).sum(axis=0)[:, None] - 2 * (y.T @ x)
        lab = np.argmin(d2, axis=0)
        for c in range(k):
            members = lab == c
            if members.any():
                y[:, c] = x[:, members].mean(axis=1)
    return y


class _Harmony:
    def __init__(self, z, phi, b_vec, theta, sigma, lam, alpha, k, max_iter_kmeans,
                 epsilon_kmeans, epsilon_harmony, block_size, batch_prop_cutoff, rng):
        self.z_orig = z
        self.z_corr = _normalize_cols(z)
        self.phi = phi                      # B x N, 0/1
        self.n = z.shape[1]
        self.b = phi.shape[0]
        self.b_vec = b_vec
        self.bounds = np.cumsum(b_vec)
        self.pr_b = phi.sum(axis=1) / self.n
        self.theta = theta
        self.sigma = sigma
        self.lam = lam                      # None: alpha * E
        self.alpha = alpha
        self.k = k
        self.max_iter_kmeans = max_iter_kmeans
        self.eps_kmeans = epsilon_kmeans
        self.eps_harmony = epsilon_harmony
        self.block_size = 0.2 if self.n < 40 else block_size
        self.cutoff = batch_prop_cutoff
        self.rng = rng
        self.batch_of = [np.flatnonzero(phi[b]) for b in range(self.b)]
        self.obj_kmeans: list[float] = []
        self.obj_harmony: list[float] = []

    def _soft_assign(self) -> None:
        self.dist = 2 * (1 - self.y.T @ self.z_corr)
        r = np.exp(-self.dist / self.sigma[:, None])
        self.r = r / r.sum(axis=0)
        self.e = self.r.sum(axis=1)[:, None] * self.pr_b[None, :]
        self.o = self.r @ self.phi.T

    def init(self) -> None:
        self.y = _normalize_cols(_kmeans_centers(self.z_corr, self.k, self.rng))
        self._soft_assign()
        self._objective()
        self.obj_harmony.append(self.obj_kmeans[-1])

    def _objective(self) -> None:
        r = self.r
        kmeans_error = float((r * self.dist).sum())
        with np.errstate(divide="ignore", invalid="ignore"):
            ent = np.where(r > 0, r * np.log(r), 0.0)
        entropy = float((ent * self.sigma[:, None]).sum())
        cross = float(((r * self.sigma[:, None])
                       * ((self.theta[None, :] * np.log((self.o + self.e + 1)
                                                        / (2 * self.e + 1))) @ self.phi)).sum())
        self.obj_kmeans.append((kmeans_error + entropy + cross) * 2000 / self.n)

    def _update_r(self) -> None:
        n = self.n
        order = self.rng.permutation(n)
        n_blocks = math.ceil(1.0 / self.block_size)
        per = int(n * self.block_size)
        for i in range(n_blocks):
            lo = i * per
            hi = n if i == n_blocks - 1 else (i + 1) * per
            cells = order[lo:hi]
            rc = self.r[:, cells]
            phic = self.phi[:, cells]
            self.e -= rc.sum(axis=1)[:, None] * self.pr_b[None, :]
            self.o -= rc @ phic.T
            rc = np.exp(-self.dist[:, cells] / self.sigma[:, None])
            rc /= rc.sum(axis=0)
            rc *= (((2 * self.e + 1) / (self.o + self.e + 1)) ** self.theta[None, :]) @ phic
            rc /= rc.sum(axis=0)
            self.r[:, cells] = rc
            self.e += rc.sum(axis=1)[:, None] * self.pr_b[None, :]
            self.o += rc @ phic.T

    def cluster(self) -> None:
        if len(self.obj_harmony) != 1:
            self.z_corr = _normalize_cols(self.z_corr)
            self._soft_assign()
        for it in range(self.max_iter_kmeans):
            self._update_r()
            self._objective()
            if it > 3:  # window of 3, as harmony
                old = sum(self.obj_kmeans[-2 - i] for i in range(3))
                new = sum(self.obj_kmeans[-1 - i] for i in range(3))
                if abs(old - new) / abs(old) < self.eps_kmeans:
                    break
        self.obj_harmony.append(self.obj_kmeans[-1])

    def correct(self) -> None:
        z = self.z_orig.copy()
        sizes = self.phi.sum(axis=1)
        n_cov = len(self.b_vec)
        cov_of = np.searchsorted(self.bounds, np.arange(self.b), side="right")
        for k in range(self.k):
            avg_r = self.o[k] / sizes
            rep = avg_r > self.cutoff
            levels = np.bincount(cov_of[rep], minlength=n_cov)
            keep = [b for b in range(self.b) if rep[b] and levels[cov_of[b]] > 1]
            if len(keep) == self.b:
                cells = np.arange(self.n)
            else:
                if not (levels > 1).any():
                    continue
                cells = np.unique(np.concatenate([self.batch_of[b] for b in keep]))
            rk = self.r[k, cells]
            phi_moe = np.vstack([np.ones(cells.size), self.phi[keep][:, cells]])
            phi_rk = phi_moe * rk[None, :]
            lam = (np.concatenate([[0.0], self.alpha * self.e[k, keep]]) if self.lam is None
                   else np.concatenate([[0.0], self.lam[keep]]))
            cov = phi_rk @ phi_moe.T + np.diag(lam)
            w = np.linalg.inv(cov) @ (phi_rk @ self.z_orig[:, cells].T)   # (1 + B') x d
            self.y[:, k] = w[0]
            w[0] = 0.0
            z[:, cells] -= w.T @ phi_rk
        self.z_corr = z
        self.y = _normalize_cols(self.y)

    def converged(self) -> bool:
        old, new = self.obj_harmony[-2], self.obj_harmony[-1]
        return (old - new) / abs(old) < self.eps_harmony


def harmony_integrate(embedding, batches, theta=2.0, sigma: float = 0.1, lambda_=None,
                      nclust: int | None = None, max_iter: int = 10, early_stop: bool = True,
                      alpha: float = 0.2, tau: float = 0.0, block_size: float = 0.05,
                      max_iter_cluster: int = 4, epsilon_cluster: float = 1e-3,
                      epsilon_harmony: float = 1e-2, batch_prop_cutoff: float = 1e-5,
                      seed: int = 0) -> np.ndarray:
    """Harmony-corrected ``embedding`` (cells x dimensions), as R's
    ``RunHarmony(embedding, batches)`` with harmony 2.0.5's defaults.

    ``batches``: one label per cell, or a list of such vectors (several
    covariates; ``theta`` and ``lambda_`` then one value or one per
    covariate). ``lambda_`` None estimates it per cluster (``alpha`` x the
    expected cells of each batch), as harmony 2 does by default. Returns the
    corrected cells x dimensions array.
    """
    z = np.asarray(embedding, dtype=np.float64)
    if z.ndim != 2:
        raise GCENInputError("embedding must be cells x dimensions.")
    n = z.shape[0]
    if n < 6:
        raise GCENInputError("Harmony needs at least 6 cells.")
    covs = batches if (isinstance(batches, (list, tuple)) and len(batches)
                       and np.ndim(batches[0]) == 1 and len(batches[0]) == n) else [batches]
    phis, b_vec = [], []
    for c in covs:
        c = np.asarray(c).astype(str)
        if c.shape != (n,):
            raise GCENInputError("batches must give one label per cell.")
        levels = np.unique(c)
        phis.append((c[None, :] == levels[:, None]).astype(np.float64))
        b_vec.append(len(levels))
    phi = np.vstack(phis)
    if nclust is None:
        nclust = min(int(np.round(n / 30)), 100)
    nclust = max(int(nclust), 1)
    th = np.atleast_1d(np.asarray(theta, dtype=np.float64))
    if th.size == 1:
        th = np.repeat(th, len(covs))
    if th.size != len(covs):
        raise GCENInputError("give theta for each batch variable.")
    theta_vec = np.repeat(th, b_vec)
    if tau > 0:
        theta_vec = theta_vec * (1 - np.exp(-(phi.sum(axis=1) / (nclust * tau)) ** 2))
    lam = None
    if lambda_ is not None:
        lv = np.atleast_1d(np.asarray(lambda_, dtype=np.float64))
        if (lv <= 0).any():
            raise GCENInputError("lambda must be positive.")
        if lv.size == 1:
            lv = np.repeat(lv, len(covs))
        lam = np.repeat(lv, b_vec)
    h = _Harmony(z.T.copy(), phi, b_vec, theta_vec, np.full(nclust, float(sigma)), lam, alpha,
                 nclust, max_iter_cluster, epsilon_cluster,
                 epsilon_harmony if early_stop else -np.inf, block_size, batch_prop_cutoff,
                 np.random.default_rng(seed))
    h.init()
    for _ in range(max_iter):
        h.cluster()
        h.correct()
        if h.converged():
            break
    return h.z_corr.T.copy()
