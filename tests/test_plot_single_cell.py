"""Single-cell and spatial module figures."""

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")

from gcen2 import plot  # noqa: E402
from gcen2.exceptions import GCENInputError  # noqa: E402


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("sc")
    rng = np.random.default_rng(0)
    n = 600
    cells = [f"c{i}" for i in range(n)]
    grp = np.array(["A", "B", "C"])[np.arange(n) % 3]
    me = pd.DataFrame({"ME0": rng.normal(size=n) + (grp == "A") * 2,
                       "ME1": rng.normal(size=n) - (grp == "B")}, index=cells)
    me.rename_axis("cell").to_csv(d / "cell_eigengenes.tsv", sep="\t", lineterminator="\n")
    xy = pd.DataFrame(rng.normal(size=(n, 2)), index=cells, columns=["x", "y"])
    xy.to_csv(d / "umap.tsv", sep="\t", lineterminator="\n")
    pd.DataFrame({"group": grp}, index=cells).to_csv(d / "obs.tsv", sep="\t", lineterminator="\n")
    dm = pd.DataFrame({"group": ["A", "A", "B", "B"], "module": ["ME0", "ME1", "ME0", "ME1"],
                       "p_val": [1e-20, 0.3, 0.2, 1e-5], "avg_log2FC": [1.5, -0.1, 0.2, -0.8],
                       "pct.1": [0.9, 0.5, 0.5, 0.3], "pct.2": [0.4, 0.5, 0.5, 0.6],
                       "p_val_adj": [2e-20, 0.6, 0.4, 2e-5]})
    dm.to_csv(d / "dmes.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"module": [0] * 30 + [1] * 20, "gene": [f"g{i}" for i in range(50)]}).to_csv(
        d / "modules.tsv", sep="\t", index=False, lineterminator="\n")
    z = pd.DataFrame({"moduleSize": [30, 20, 50, 50], "Zsummary.pres": [15.0, 3.0, 1.0, 0.5]},
                     index=pd.Index(["ME0", "ME1", "gold", "grey"], name="module"))
    z.to_csv(d / "z.tsv", sep="\t", lineterminator="\n")
    pd.DataFrame({"moduleSize": [30, 20, 50, 50], "medianRank.pres": [1, 2, 3, 4]},
                 index=z.index).to_csv(d / "obs_pres.tsv", sep="\t", lineterminator="\n")
    se = pd.DataFrame({"set": ["s1"] * 20 + ["s2"] * 20, "sample": range(40),
                       "ME_unassigned": rng.normal(size=40), "ME0": rng.normal(size=40)})
    se.to_csv(d / "set_me.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"cell_type": ["T", "B"], "dir": ["T", "B"], "n_cells": [500, 300],
                  "n_metacells": [100, 60], "n_genes": [3000, 2500], "power": [6, 8],
                  "modules": [5, 3], "genes_in_modules": [900, 500], "reassigned": [3, 1]}
                 ).to_csv(d / "celltypes.tsv", sep="\t", index=False, lineterminator="\n")
    return d, pd.Series(grp, index=cells)


def test_figures(data, tmp_path):
    d, g = data
    f = ["png"]
    figs = [
        plot.me_embedding(d / "cell_eigengenes.tsv", d / "umap.tsv", output=tmp_path / "a",
                          formats=f),
        plot.spatial_me(d / "cell_eigengenes.tsv", d / "umap.tsv", output=tmp_path / "b",
                        formats=f),
        plot.me_dotplot(d / "cell_eigengenes.tsv", g, output=tmp_path / "c", formats=f),
        plot.me_violin(d / "cell_eigengenes.tsv", g, output=tmp_path / "d", formats=f),
        plot.me_radar(d / "cell_eigengenes.tsv", g, output=tmp_path / "e", formats=f),
        plot.dme_volcano(d / "dmes.tsv", output=tmp_path / "f", formats=f),
        plot.dme_lollipop(d / "dmes.tsv", modules=d / "modules.tsv", output=tmp_path / "g",
                          formats=f),
        plot.preservation(d / "z.tsv", d / "obs_pres.tsv", output=tmp_path / "h", formats=f),
        plot.preservation_lollipop(d / "z.tsv", output=tmp_path / "i", formats=f),
        plot.set_eigengenes(d / "set_me.tsv", output=tmp_path / "j", formats=f),
        plot.celltype_overview(d / "celltypes.tsv", output=tmp_path / "k", formats=f),
    ]
    assert len(list(tmp_path.glob("*.png"))) == len(figs)
    assert figs[0].axes[0].get_title() == "ME0 (turquoise)"
    assert figs[1].axes[0].yaxis_inverted()  # image coordinates
    titles = [a.get_title() for a in figs[5].axes if a.get_title()]
    assert titles == ["A", "B"]


def test_embedding_options_and_errors(data):
    d, g = data
    fig = plot.me_embedding(d / "cell_eigengenes.tsv", d / "umap.tsv", modules=["1"])
    assert [a.get_title() for a in fig.axes if a.get_title()] == ["ME1 (blue)"]
    with pytest.raises(GCENInputError):
        plot.me_embedding(d / "cell_eigengenes.tsv", d / "umap.tsv", modules=["ME7"])
    with pytest.raises(GCENInputError):
        plot.me_radar(d / "cell_eigengenes.tsv", g.replace("C", "B"))
    with pytest.raises(GCENInputError):
        plot.me_dotplot(d / "cell_eigengenes.tsv", g.to_numpy()[:10])


def test_cli(data, tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    d, _ = data
    cases = [
        ["me-embedding", "-i", str(d / "cell_eigengenes.tsv"), "-c", str(d / "umap.tsv")],
        ["me-dotplot", "-i", str(d / "cell_eigengenes.tsv"), "-g", str(d / "obs.tsv"),
         "--groupby", "group"],
        ["me-radar", "-i", str(d / "cell_eigengenes.tsv"), "-g", str(d / "obs.tsv"),
         "--groupby", "group"],
        ["dme-volcano", "-i", str(d / "dmes.tsv")],
        ["dme-lollipop", "-i", str(d / "dmes.tsv"), "-m", str(d / "modules.tsv")],
        ["preservation", "-i", str(d / "z.tsv"), "--observed", str(d / "obs_pres.tsv")],
        ["preservation", "-i", str(d / "z.tsv"), "--lollipop"],
        ["set-eigengenes", "-i", str(d / "set_me.tsv")],
        ["celltype-overview", "-i", str(d / "celltypes.tsv")],
    ]
    for i, args in enumerate(cases):
        r = CliRunner().invoke(app, ["plot", *args, "-o", str(tmp_path / f"p{i}"),
                                     "--format", "pdf"])
        assert r.exit_code == 0, (args, r.output)
        assert (tmp_path / f"p{i}.pdf").exists()
    r = CliRunner().invoke(app, ["plot", "me-dotplot", "-i", str(d / "cell_eigengenes.tsv"),
                                 "-g", str(d / "obs.tsv"), "--groupby", "nope", "-o",
                                 str(tmp_path / "x")])
    assert r.exit_code == 1


def test_celltype_order_same_in_overview_and_metacells(tmp_path):
    """Cell types top to bottom in the order of celltypes.tsv, in both figures."""
    import numpy as np
    import pandas as pd

    order = ["Mono", "NK", "B"]  # not sorted by cells
    pd.DataFrame({"cell_type": order, "dir": order, "n_cells": [30, 10, 20],
                  "n_metacells": [3, 1, 2], "n_genes": [100] * 3, "modules": [2] * 3,
                  "genes_in_modules": [50] * 3}).to_csv(tmp_path / "celltypes.tsv", sep="\t",
                                                       index=False, lineterminator="\n")
    cells = [f"c{i}" for i in range(60)]
    k = 0
    for t, n in zip(order, (3, 1, 2)):
        (tmp_path / t).mkdir()
        pd.DataFrame({"metacell": np.repeat([f"{t}{i}" for i in range(n)], 3),
                      "cell": cells[k:k + 3 * n], "center": np.tile([True, False, False], n)}
                     ).to_csv(tmp_path / t / "metacell_members.tsv", sep="\t", index=False, lineterminator="\n")
        k += 3 * n
    xy = pd.DataFrame(np.random.default_rng(0).normal(size=(60, 2)), index=cells,
                      columns=["x", "y"])

    def top_to_bottom(ax):
        ticks = {lab.get_text(): y for lab, y in zip(ax.get_yticklabels(), ax.get_yticks())}
        return sorted(ticks, key=lambda c: -ax.transData.transform((0, ticks[c]))[1])

    assert top_to_bottom(plot.celltype_overview(tmp_path / "celltypes.tsv").axes[0]) == order
    assert top_to_bottom(plot.metacells(tmp_path, xy).axes[1]) == order
