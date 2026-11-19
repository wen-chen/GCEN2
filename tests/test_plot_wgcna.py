"""WGCNA / hdWGCNA figures from the files GCEN2 writes."""

import json

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")

from gcen2 import plot  # noqa: E402
from gcen2.exceptions import GCENInputError  # noqa: E402
from gcen2.tom_modules import gene_significance, tom_modules, write_module_files  # noqa: E402


@pytest.fixture(scope="module")
def run(tmp_path_factory):
    d = tmp_path_factory.mktemp("wgcna")
    rng = np.random.default_rng(0)
    x = rng.normal(size=(240, 40))
    for k in range(3):
        x[k * 60:(k + 1) * 60] += rng.normal(size=40) * 1.5
    genes = [f"g{i}" for i in range(240)]
    samples = [f"s{i}" for i in range(40)]
    with pytest.warns(UserWarning):
        res = tom_modules(x, gene_names=genes, sample_names=samples, signed=True, min_size=20,
                          power_method="pearson", power_r2=0.95, min_power=3)
    traits = pd.DataFrame({"t1": x[:60].mean(axis=0), "t2": rng.normal(size=40)},
                          index=samples)
    write_module_files(res, modules=d / "modules.tsv", eigengenes=d / "me.tsv",
                       kme=d / "kme.tsv", hubs=d / "hubs.tsv", trait_table=res.trait_cor(traits),
                       trait_cor=d / "tc.tsv", soft_power=d / "sp.tsv",
                       dendrogram=d / "tree.json", tom_hubs=d / "tomh.tsv",
                       gs_table=gene_significance(res, traits), gs=d / "gs.tsv")
    pd.DataFrame(x, index=genes, columns=samples).to_csv(d / "expr.tsv", sep="\t", lineterminator="\n")
    traits.to_csv(d / "traits.tsv", sep="\t", lineterminator="\n")
    return d, res


def test_written_files(run):
    d, res = run
    sp = pd.read_csv(d / "sp.tsv", sep="\t")
    assert sp["selected"].sum() == 1 and sp.loc[sp["selected"], "power"].iat[0] == res.power
    tree = json.loads((d / "tree.json").read_text(encoding="utf-8"))
    assert len(tree["genes"]) == len(tree["merge"]) + 1 == len(tree["modules"])
    tomh = pd.read_csv(d / "tomh.tsv", sep="\t", index_col=0)
    assert tomh.shape[1] == 25 * int(res.labels.max())
    gs = pd.read_csv(d / "gs.tsv", sep="\t", index_col=0)
    assert list(gs.columns) == ["GS.t1", "GS.t2", "p.GS.t1", "p.GS.t2"]
    assert gs["GS.t1"].iloc[:60].mean() > 0.3


def test_figures(run, tmp_path):
    d, _ = run
    f = ["png"]
    figs = {
        "soft": plot.soft_power(d / "sp.tsv", output=tmp_path / "soft", formats=f),
        "tree": plot.dendrogram(d / "tree.json", output=tmp_path / "tree", formats=f),
        "me": plot.eigengene_network(d / "me.tsv", traits=d / "traits.tsv",
                                     output=tmp_path / "me", formats=f),
        "trait": plot.module_trait(d / "tc.tsv", annotate="stars", output=tmp_path / "tr",
                                   formats=f),
        "kme": plot.kme(d / "kme.tsv", output=tmp_path / "kme", formats=f),
        "mmgs": plot.mm_gs(d / "kme.tsv", d / "gs.tsv", "ME0", "t1", output=tmp_path / "mmgs",
                           formats=f),
        "hub": plot.hub_network(d / "tomh.tsv", d / "kme.tsv", "0", output=tmp_path / "hub",
                                formats=f),
        "huball": plot.hub_network_all(d / "tomh.tsv", d / "kme.tsv",
                                       output=tmp_path / "huball", formats=f),
        "tom": plot.tom_heatmap(d / "expr.tsv", d / "modules.tsv", power=6, n_genes=100,
                                output=tmp_path / "tom", formats=f),
        "samples": plot.sample_tree(d / "expr.tsv", traits=d / "traits.tsv", cut_height=10,
                                    output=tmp_path / "samples", formats=f),
        "sizes": plot.module_sizes(d / "modules.tsv", output=tmp_path / "sizes", formats=f),
    }
    figs["umap"] = plot.module_umap(d / "tomh.tsv", d / "kme.tsv", n_hubs=5,
                                    output=tmp_path / "umap", formats=f)
    for name in figs:
        assert figs[name] is not None
    assert len(list(tmp_path.glob("*.png"))) == len(figs)
    assert "t1" in [t.get_text() for t in figs["me"].axes[1].get_xticklabels()]
    assert figs["kme"].axes[0].get_title() == "ME0 (turquoise)"
    assert figs["mmgs"].axes[0].get_title().startswith("r = ")
    hub_texts = [t.get_text() for t in figs["hub"].axes[0].texts]
    assert len(hub_texts) == 25


def test_errors(run):
    d, _ = run
    with pytest.raises(GCENInputError):
        plot.soft_power(pd.DataFrame())
    with pytest.raises(GCENInputError):
        plot.mm_gs(d / "kme.tsv", d / "gs.tsv", "ME9", "t1")
    with pytest.raises(GCENInputError):
        plot.module_trait(d / "tc.tsv", annotate="numbers")


def test_cli(run, tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    d, _ = run
    r = CliRunner().invoke(app, ["plot", "hub-network", "-t", str(d / "tomh.tsv"), "-k",
                                 str(d / "kme.tsv"), "-o", str(tmp_path / "hub"), "--format",
                                 "png"])
    assert r.exit_code == 0, r.output
    assert len(list(tmp_path.glob("hub_ME*.png"))) == 3
    for args in (["soft-power", "-i", str(d / "sp.tsv")],
                 ["dendrogram", "-i", str(d / "tree.json")],
                 ["module-trait", "-i", str(d / "tc.tsv")],
                 ["sample-tree", "-e", str(d / "expr.tsv")]):
        r = CliRunner().invoke(app, ["plot", *args, "-o", str(tmp_path / args[0]),
                                     "--format", "svg"])
        assert r.exit_code == 0, r.output
        assert (tmp_path / f"{args[0]}.svg").exists()


def test_module_cli_outputs(tmp_path, run):
    from typer.testing import CliRunner

    from gcen2.cli import app

    d, _ = run
    r = CliRunner().invoke(app, ["module", "-e", str(d / "expr.tsv"), "-o",
                                 str(tmp_path / "m.tsv"), "--signed", "--power", "6",
                                 "--min-size", "20", "--soft-power", str(tmp_path / "sp.tsv"),
                                 "--dendrogram", str(tmp_path / "t.json"), "--tom-hubs",
                                 str(tmp_path / "th.tsv"), "--traits", str(d / "traits.tsv"),
                                 "--gs", str(tmp_path / "gs.tsv")])
    assert r.exit_code == 0, r.output
    assert (tmp_path / "t.json").exists() and (tmp_path / "th.tsv").exists()
    assert (tmp_path / "gs.tsv").exists()
    r = CliRunner().invoke(app, ["module", "-e", str(d / "expr.tsv"), "-o",
                                 str(tmp_path / "m2.tsv"), "--gs", str(tmp_path / "gs2.tsv")])
    assert r.exit_code == 1
