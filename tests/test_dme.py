"""hdWGCNA FindDMEs / FindAllDMEs and ModuleTraitCorrelation (reference
values from R: wilcox.test, Seurat FindMarkers, Hmisc::rcorr)."""

import numpy as np
import pandas as pd
import pytest

from gcen2.exceptions import GCENInputError
from gcen2.single_cell import find_all_dmes, find_dmes, module_trait_cor
from gcen2.single_cell.dme import _bh, _wilcox_p

X1 = [0.5, 1.2, 0, 0.3, 2.1, 0.8, 1.5, 0, 0.9, 1.1]
X2 = [0.1, 0, 0.2, 0.4, 0, 0.6, 0.05, 0.3]


def test_wilcox_as_r():
    assert _wilcox_p(np.array(X1), np.array(X2)) == pytest.approx(0.049305340625921366,
                                                                    rel=1e-12)
    # no ties, fewer than 50 per group: exact
    assert _wilcox_p(np.array([1.3, 2.2, 0.7, 1.9, 2.8]),
                     np.array([0.2, 0.9, 0.4, 1.1, 0.6, 0.3])) == pytest.approx(
        0.017316017316017316, rel=1e-12)


def _me():
    cells = [f"c{i}" for i in range(1, 19)]
    return pd.DataFrame({"m1": X1 + X2, "m2": X1[::-1] + X2[::-1]}, index=cells)


def test_find_dmes_as_findmarkers():
    me = _me()
    d = find_dmes(me, me.index[:10], me.index[10:])
    assert list(d.columns) == ["p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj", "module"]
    assert d.loc["m1", "p_val"] == pytest.approx(0.049305340625921366, rel=1e-12)
    assert d.loc["m1", "avg_log2FC"] == pytest.approx(2.8698958332585445, rel=1e-12)
    assert d.loc["m1", "pct.1"] == 0.8 and d.loc["m1", "pct.2"] == 0.75
    assert d.loc["m1", "p_val_adj"] == pytest.approx(0.098610681251842733, rel=1e-12)
    with pytest.raises(GCENInputError):
        find_dmes(me, me.index[:10], me.index[9:])
    with pytest.raises(GCENInputError):
        find_dmes(me, ["nope", "c1", "c2"], me.index[10:])
    with pytest.raises(GCENInputError):
        find_dmes(me, me.index[:2], me.index[10:])


def test_find_all_dmes_and_skip():
    me = _me()
    g = ["a"] * 10 + ["b"] * 6 + ["c"] * 2
    with pytest.warns(UserWarning, match="'c'"):
        d = find_all_dmes(me, g)
    assert list(pd.unique(d["group"])) == ["a", "b"]
    assert d.index[0].startswith("a.")
    only = find_dmes(me, me.index[:10], me.index[10:], only_pos=True)
    assert (only["avg_log2FC"] > 0).all()


def test_module_trait_cor_as_rcorr():
    me = _me().iloc[:10]
    t = pd.DataFrame({"t": [1, 2, 3, 4, 5, 6, np.nan, 8, 9, 10]}, index=me.index)
    m = module_trait_cor(me, t).set_index("module")
    assert m.loc["m1", "cor"] == pytest.approx(0.10073466091724519, rel=1e-12)
    assert m.loc["m1", "p"] == pytest.approx(0.79651246353623772, rel=1e-9)
    assert m.loc["m2", "cor"] == pytest.approx(-0.11013798774303998, rel=1e-12)
    assert m.loc["m2", "p"] == pytest.approx(0.77788446086844232, rel=1e-9)
    assert (m["n"] == 9).all() and (m["group"] == "all_cells").all()


def test_module_trait_cor_groups_and_categories():
    rng = np.random.default_rng(0)
    me = pd.DataFrame(rng.normal(size=(40, 2)), columns=["ME1", "ME2"],
                      index=[f"c{i}" for i in range(40)])
    t = pd.DataFrame({"sex": pd.Categorical(rng.choice(["f", "m"], 40)),
                      "stage": rng.choice(["early", "late"], 40)}, index=me.index)
    g = np.array(["x"] * 20 + ["y"] * 17 + ["z"] * 3)
    with pytest.warns(UserWarning, match="stage"):
        m = module_trait_cor(me, t, groups=g)
    assert list(pd.unique(m["group"])) == ["all_cells", "x", "y"]  # z: < 5 cells
    assert len(m) == 3 * 2 * 2
    with pytest.raises(GCENInputError):
        module_trait_cor(me, t.iloc[:10])
    with pytest.raises(GCENInputError):
        module_trait_cor(me, t, method="kendall")


def test_bh_as_p_adjust():
    p = np.array([0.01, 0.04, np.nan, 0.03, 0.5])
    np.testing.assert_allclose(_bh(p), [0.04, 0.16 / 3, np.nan, 0.16 / 3, 0.5])
