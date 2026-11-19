"""Module preservation (WGCNA modulePreservation; reference values from R
WGCNA 1.74) and module projection (hdWGCNA ProjectModules)."""

import numpy as np
import pandas as pd
import pytest

from gcen2.exceptions import GCENInputError
from gcen2.preservation import module_preservation


def _small():
    rng = np.random.default_rng(11)
    genes = [f"g{i}" for i in range(120)]

    def make(n, keep):
        f = rng.normal(size=(n, 3))
        x = rng.normal(size=(n, 120))
        for k in range(3):
            x[:, k * 25:(k + 1) * 25] += (1.2 if keep[k] else 0.0) * f[:, [k]]
        return pd.DataFrame(x, columns=genes, index=[f"s{i}" for i in range(n)])

    ref = make(60, [True, True, True])
    test = make(50, [True, True, False])
    mods = pd.Series(["blue"] * 25 + ["brown"] * 25 + ["turquoise"] * 25 + ["grey"] * 45,
                     index=genes)
    return ref, test, mods


def test_preservation_as_wgcna():
    ref, test, mods = _small()
    obs, z = module_preservation(ref, test, mods, n_permutations=10, threads=2)
    assert list(z.index) == ["blue", "brown", "gold", "grey", "turquoise"]
    # modulePreservation(..., nPermutations = 10, randomSeed = 12345, quickCor = 0)
    want_z = [11.2903929037, 2.48923609485, 10.1268013442, -1.96833768442, -1.14147957897]
    np.testing.assert_allclose(z["Zsummary.pres"], want_z, rtol=1e-6)
    want_kme = [0.00579649053947, 0.198728116576, 0.874979238336, 0.145705935346,
                0.115056501651]
    np.testing.assert_allclose(obs["cor.kME"], want_kme, rtol=1e-6)
    assert list(z["moduleSize"]) == [25, 25, 60, 45, 25]
    assert {"Zdensity.pres", "Zconnectivity.pres", "Zsummary.qual"} <= set(z.columns)
    assert "meanClusterCoeff.pres" not in obs.columns  # all missing: dropped


def test_preservation_threads_and_errors():
    ref, test, mods = _small()
    _, a = module_preservation(ref, test, mods, n_permutations=4, threads=1)
    _, b = module_preservation(ref, test, mods, n_permutations=4, threads=3)
    pd.testing.assert_frame_equal(a, b)
    with pytest.raises(GCENInputError):
        module_preservation(ref, test, mods, n_permutations=1)
    with pytest.raises(GCENInputError):
        module_preservation(ref, test, pd.Series("grey", index=mods.index))
    bad = test.copy()
    bad["g3"] = 1.0
    with pytest.raises(GCENInputError, match="constant"):
        module_preservation(ref, bad, mods, n_permutations=4)
    with pytest.raises(GCENInputError):
        module_preservation(ref, test, mods, network_type="weird")


ad = pytest.importorskip("anndata")


def _query(n=400, seed=0):
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    genes = [f"g{i}" for i in range(60)]
    f = rng.normal(size=(n, 2))
    logmu = np.log(rng.uniform(0.5, 2, 60))[None, :].repeat(n, 0)
    logmu[:, :20] += 0.8 * f[:, [0]]
    logmu[:, 20:40] += 0.8 * f[:, [1]]
    x = rng.poisson(np.exp(logmu)).astype(np.float32)
    obs = pd.DataFrame({"batch": np.array(["a", "b"])[np.arange(n) % 2],
                        "part": np.array(["x", "y"])[np.arange(n) // (n // 2)]},
                       index=[f"c{i}" for i in range(n)])
    return ad.AnnData(sp.csr_matrix(x), obs=obs, var=pd.DataFrame(index=genes))


def test_project_modules():
    from gcen2.single_cell import cell_eigengenes, project_modules

    q = _query()
    genes = list(q.var_names)
    mods = pd.DataFrame({"gene": genes[:20] + genes[20:40] + ["zz1", "zz2", "zz3", genes[50]],
                         "module": ["A"] * 20 + ["B"] * 20 + ["C"] * 4})
    with pytest.warns(UserWarning, match="C"):
        me, hme, table = project_modules(q, mods)
    assert list(me.columns) == ["A", "B"] and hme is None
    assert table.set_index("module").loc["C", "genes_found"] == 1
    want, _ = cell_eigengenes(q.X, np.arange(40), np.repeat([1, 2], 20), [1, 2])
    np.testing.assert_allclose(me.to_numpy(), want, atol=1e-12)
    sub, hsub, _ = project_modules(q, mods.iloc[:40], cells=q.obs["part"] == "x",
                                   harmonize_by="batch")
    assert sub.shape == (200, 2) and hsub.shape == (200, 2)
    with pytest.raises(GCENInputError):
        project_modules(q, mods.assign(gene=["nope"] * len(mods)))


def test_cli_project_and_preservation(tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    q = _query()
    path = tmp_path / "q.h5ad"
    q.write_h5ad(path)
    genes = list(q.var_names)
    mods = tmp_path / "modules.tsv"
    pd.DataFrame({"module": [0] * 20 + [1] * 20, "gene": genes[:40]}).to_csv(
        mods, sep="\t", index=False, lineterminator="\n")
    r = CliRunner().invoke(app, ["project-modules", "-i", str(path), "-m", str(mods),
                                 "-o", str(tmp_path / "p"), "--groupby", "part",
                                 "--group", "x"])
    assert r.exit_code == 0, r.output
    me = pd.read_csv(tmp_path / "p" / "projected_eigengenes.tsv", sep="\t", index_col=0)
    assert list(me.columns) == ["ME0", "ME1"] and len(me) == 200

    ref, test, smods = _small()
    rp, tp, mp = tmp_path / "ref.tsv", tmp_path / "test.tsv", tmp_path / "smods.tsv"
    ref.T.to_csv(rp, sep="\t", lineterminator="\n")
    test.T.to_csv(tp, sep="\t", lineterminator="\n")
    pd.DataFrame({"gene": smods.index, "module": smods.to_numpy()}).to_csv(mp, sep="\t",
                                                                          index=False, lineterminator="\n")
    r = CliRunner().invoke(app, ["module-preservation", "-r", str(rp), "-m", str(mp),
                                 "-q", str(tp), "--permutations", "10",
                                 "-o", str(tmp_path / "mp")])
    assert r.exit_code == 0, r.output
    z = pd.read_csv(tmp_path / "mp" / "preservation_Z.tsv", sep="\t", index_col=0)
    assert z.loc["blue", "Zsummary.pres"] == pytest.approx(11.2903929037, rel=1e-6)
    r = CliRunner().invoke(app, ["module-preservation", "-r", str(rp), "-m", str(mp),
                                 "-q", str(path), "--permutations", "4", "--groupby", "part",
                                 "--group", "x", "-k", "10", "--max-metacells", "50",
                                 "-o", str(tmp_path / "mp2")])
    assert r.exit_code == 0, r.output  # query metacells built from the h5ad (genes g0..g59)
    z2 = pd.read_csv(tmp_path / "mp2" / "preservation_Z.tsv", sep="\t", index_col=0)
    assert z2.loc["blue", "moduleSize"] == 25
