"""Consensus modules (WGCNA blockwiseConsensusModules), R's RNG and Harmony."""

import numpy as np
import pandas as pd
import pytest

from gcen2._rrng import RRandom
from gcen2.consensus import (
    HDWGCNA_CONSENSUS,
    consensus_tom_modules,
    r_quantile7,
    r_quantile8,
)
from gcen2.core import require_extension as _core
from gcen2.exceptions import GCENInputError
from gcen2.harmony import harmony_integrate
from gcen2.tom_modules import HDWGCNA, tom_modules

# --- R's random number generator -------------------------------------------

def test_rrandom_matches_r():
    # set.seed(54321); runif(3) / sample(1000, 5) / sample(5e7, 20000)
    r = RRandom(54321)
    assert [f"{r.unif_rand():.17g}" for _ in range(3)] == [
        "0.42900781612843275", "0.49843042111024261", "0.17669233866035938"]
    assert RRandom(54321).sample(1000, 5).tolist() == [468, 922, 316, 575, 878]
    big = RRandom(54321).sample(int(5e7), 20000)  # R's hash algorithm
    assert big[:5].tolist() == [30637978, 20661823, 10892682, 49982738, 41659770]
    assert big[-3:].tolist() == [6602040, 31526101, 33264739]
    assert RRandom(7).sample(100, 10).tolist() == [42, 83, 31, 92, 66, 15, 90, 8, 67, 88]


def test_rrandom_sample_errors():
    with pytest.raises(ValueError):
        RRandom(1).sample(3, 4)


# --- quantiles as R / WGCNA --------------------------------------------------

def test_r_quantiles():
    x = np.array([0.12, 0.5, 0.33, 0.91, 0.07, 0.64, 0.28, 0.45, 0.8, 0.19])
    assert r_quantile8(x, 0.95) == 0.91
    assert r_quantile8(x, 0.3) == 0.22899999999999998
    m = np.array([[1.0, 4.0, np.nan], [3.0, 2.0, 5.0]])
    np.testing.assert_allclose(r_quantile7(m, 0.3, axis=1),
                               [np.quantile([1.0, 4.0], 0.3), np.quantile([3, 2, 5], 0.3)])


def test_row_quantile_is_wgcna_rowquantilec():
    # WGCNA rowQuantileC(m, 0.3), printed with 17 digits
    m = np.array([[0.3, 0.1, 0.7], [0.9, 0.2, 0.4], [0.5, 0.5, 0.6], [0.11, 0.83, 0.27]])
    out = _core().module.row_quantile([m[:, j].copy() for j in range(3)], 0.3, 1)
    assert [f"{v:.17g}" for v in out] == ["0.22", "0.32000000000000001", "0.5",
                                         "0.20600000000000002"]
    m2 = np.array([[0.3, 0.1, 0.7, 0.2, 0.9], [0.15, 0.95, 0.35, 0.55, 0.45]])
    out = _core().module.row_quantile([m2[:, j].copy() for j in range(5)], 0.3, 1)
    assert [f"{v:.17g}" for v in out] == ["0.22000000000000003", "0.36999999999999994"]


# --- consensus modules ------------------------------------------------------

def _sets(n_sets=3, n_genes=240, n_samples=40, seed=0):
    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(n_genes)]
    out = []
    for s in range(n_sets):
        f = rng.normal(size=(4, n_samples))
        load = np.zeros((n_genes, 4))
        for k in range(4):
            load[k * 50:(k + 1) * 50, k] = 1.0
        x = load @ f + 0.7 * rng.normal(size=(n_genes, n_samples))
        out.append(pd.DataFrame(x, index=genes,
                                columns=[f"s{s}_{j}" for j in range(n_samples)]))
    return out


def test_one_set_is_tom_modules():
    x = _sets(1)[0]
    opts = dict(HDWGCNA, power=6)
    for k in ("power_method", "power_r2", "min_power", "powers"):
        opts.pop(k)
    a = tom_modules(x, **opts)
    b = consensus_tom_modules([x], **opts)
    assert (a.labels == b.labels).all()
    assert a.labels.max() >= 3


def test_consensus_finds_shared_modules():
    sets = _sets()
    res = consensus_tom_modules(sets, set_names=["a", "b", "c"], power=6,
                                **{k: v for k, v in HDWGCNA_CONSENSUS.items()
                                   if k not in ("power_method", "power_r2", "min_power",
                                                "powers")})
    md = res.metadata
    assert md["sets"] == ["a", "b", "c"]
    assert md["calibration_powers"]["a"] == 1.0
    assert res.labels.max() == 4
    for k in range(4):  # each planted block is one module
        assert len(set(res.labels[k * 50:(k + 1) * 50])) == 1
    assert (res.labels[200:] == 0).mean() >= 0.8
    assert len(res.samples) == 120 and res.eigengenes.shape[1] == 120
    assert set(md["set_eigengenes"]) == {"a", "b", "c"}
    assert md["set_eigengenes"]["b"].shape[0] == 40
    # a set given with its genes in another order gives the same modules
    shuffled = [sets[0], sets[1].iloc[::-1], sets[2]]
    res2 = consensus_tom_modules(shuffled, power=6, **{
        k: v for k, v in HDWGCNA_CONSENSUS.items()
        if k not in ("power_method", "power_r2", "min_power", "powers")})
    assert (res2.labels == res.labels).all()


def test_consensus_powers_and_errors():
    sets = _sets(2, n_genes=150)
    res = consensus_tom_modules(sets, power=[5, 7], signed=True, min_size=30)
    assert list(res.metadata["set_powers"].values()) == [5.0, 7.0]
    with pytest.warns(UserWarning, match="FAQ power"):  # 40 noisy samples: no fit
        res = consensus_tom_modules(sets, signed=True, min_size=30, power_method="bicor",
                                    power_r2=0.8, min_power=3)
    assert len(res.metadata["set_powers"]) == 2
    with pytest.raises(GCENInputError):
        consensus_tom_modules(sets, power=[5, 6, 7])
    with pytest.raises(GCENInputError):
        consensus_tom_modules(sets, consensus_quantile=1.5)
    with pytest.raises(GCENInputError):
        consensus_tom_modules([sets[0], sets[1].iloc[:100]])
    with pytest.raises(GCENInputError):
        consensus_tom_modules(sets, set_names=["a", "a"])


# --- Harmony ----------------------------------------------------------------

def _batched(n=1200, seed=0):
    rng = np.random.default_rng(seed)
    grp = rng.integers(0, 3, n)
    bat = rng.integers(0, 2, n)
    z = rng.normal(size=(n, 8))
    z[:, 0] += grp * 4
    z[:, 1] += bat * 3
    return z, grp, bat


def test_harmony_removes_batch_and_keeps_groups():
    z, grp, bat = _batched()
    h = harmony_integrate(z, bat, seed=1)
    shift = abs(h[bat == 1, 1].mean() - h[bat == 0, 1].mean())
    assert shift < 0.3 * abs(z[bat == 1, 1].mean() - z[bat == 0, 1].mean())
    assert np.corrcoef(h[:, 0], grp)[0, 1] > 0.9
    np.testing.assert_array_equal(h, harmony_integrate(z, bat, seed=1))


def test_harmony_options_and_errors():
    z, _, bat = _batched(300)
    other = np.where(np.arange(300) % 3 == 0, "x", "y")
    h = harmony_integrate(z, [bat, other], theta=[2, 1], lambda_=1.0, nclust=5, seed=0)
    assert h.shape == z.shape
    with pytest.raises(GCENInputError):
        harmony_integrate(z[:5], bat[:5])
    with pytest.raises(GCENInputError):
        harmony_integrate(z, bat[:10])
    with pytest.raises(GCENInputError):
        harmony_integrate(z, bat, lambda_=-1)
