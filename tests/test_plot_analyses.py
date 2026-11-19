"""Figures of GCEN2's other analyses: QC, networks, enrichment, RWR,
differential co-expression, trajectory and fusion."""

import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("matplotlib")

from gcen2 import plot  # noqa: E402
from gcen2.exceptions import GCENInputError  # noqa: E402


@pytest.fixture(scope="module")
def data(tmp_path_factory):
    d = tmp_path_factory.mktemp("an")
    rng = np.random.default_rng(1)
    genes = [f"g{i}" for i in range(80)]
    samples = [f"s{i}" for i in range(12)]
    x = rng.lognormal(2, 1, size=(80, 12))
    pd.DataFrame(x, index=genes, columns=samples).to_csv(d / "expr.tsv", sep="\t", lineterminator="\n")
    pd.DataFrame(x[:50] * 1.1, index=genes[:50], columns=samples).to_csv(d / "kept.tsv",
                                                                        sep="\t", lineterminator="\n")
    a, b = np.triu_indices(40, 1)
    keep = rng.random(a.size) < 0.25
    net = pd.DataFrame({"gene_a": np.array(genes)[a[keep]], "gene_b": np.array(genes)[b[keep]],
                        "correlation": rng.uniform(0.6, 1, keep.sum()),
                        "p_value": rng.uniform(0, 1e-3, keep.sum())})
    net.to_csv(d / "net.network", sep="\t", index=False, lineterminator="\n")
    shuf = net.copy()
    shuf["gene_b"] = rng.permutation(shuf["gene_b"].to_numpy())
    shuf = shuf[shuf["gene_a"] != shuf["gene_b"]]
    shuf.to_csv(d / "shuf.network", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"module": [0] * 20 + [1] * 20, "gene": genes[:40]}).to_csv(
        d / "modules.txt", sep="\t", index=False, lineterminator="\n")
    term = pd.DataFrame({"term_id": [f"GO:{i}" for i in range(6)],
                         "name": [f"process {i}" for i in range(6)],
                         "namespace": ["biological_process"] * 3 + ["cellular_component"] * 3,
                         "enrichment": [3.0] * 6, "study_count": [5, 4, 3, 6, 2, 3],
                         "study_n": [20] * 6, "pop_count": [10] * 6, "pop_n": [100] * 6,
                         "p_value": [1e-8, 1e-6, 1e-4, 1e-5, 0.01, 0.02]})
    term.to_csv(d / "enrich.tsv", sep="\t", index=False, lineterminator="\n")
    (d / "annot").mkdir()
    term.to_csv(d / "annot" / "module_0.go", sep="\t", index=False, lineterminator="\n")
    term.iloc[::-1].to_csv(d / "annot" / "module_1.go", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"gene": genes[:30], "probability": np.linspace(1, 0.01, 30)}).to_csv(
        d / "rwr.tsv", sep="\t", index=False, lineterminator="\n")
    (d / "seeds.txt").write_text("g0\ng1\n", encoding="utf-8")
    diff = pd.DataFrame({"gene_a": ["g0", "g1", "g2", "g3", "g4"],
                         "gene_b": ["g5", "g6", "g7", "g8", "g9"],
                         "cor_a": [0.9, 0.8, np.nan, 0.85, 0.95],
                         "cor_b": [0.5, 0.82, 0.9, np.nan, 0.6],
                         "z_score": [4, 0.1, np.nan, np.nan, 5],
                         "p_value": [1e-4, 0.9, np.nan, np.nan, 1e-6],
                         "category": ["shared_changed", "shared_stable", "B_only", "A_only",
                                      "shared_changed"]})
    diff.to_csv(d / "diff.tsv", sep="\t", index=False, lineterminator="\n")
    t = d / "traj"
    t.mkdir()
    pd.DataFrame({"window": range(4), "start": [0, 1, 2, 3], "end": [2, 3, 4, 5],
                  "center": [1, 2, 3, 4], "n_cells": 50, "n_samples": 10,
                  "node_count": [30, 25, 40, 35], "edge_count": [100, 80, 150, 120]}).to_csv(
        t / "windows.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"gene_a": ["a", "b"], "gene_b": ["c", "d"],
                  "pattern": ["gained", "stable"]}).to_csv(t / "edge_dynamics.tsv", sep="\t",
                                                           index=False, lineterminator="\n")
    pd.DataFrame(rng.integers(0, 20, size=(30, 4)), index=genes[:30],
                 columns=[f"degree_{i}" for i in range(4)]).rename_axis("gene").to_csv(
        t / "degree_dynamics.tsv", sep="\t", lineterminator="\n")
    pd.DataFrame({"gene_a": ["g0", "g1", "g2"], "gene_b": ["g3", "g4", "g5"],
                  "score": [1.0, 0.5, 0.4], "n_layers": [2, 1, 1],
                  "layers": ["L1,L2", "L1", "L2"], "w_L1": [1.0, 0.5, np.nan],
                  "w_L2": [1.0, np.nan, 0.4]}).to_csv(d / "fused.tsv", sep="\t", index=False, lineterminator="\n")
    return d


def test_figures(data, tmp_path):
    d = data
    f = ["png"]
    figs = [
        plot.degree_distribution(d / "net.network", output=tmp_path / "a", formats=f),
        plot.edge_weights(d / "net.network", output=tmp_path / "b", formats=f),
        plot.network_view(d / "net.network", modules=d / "modules.txt", output=tmp_path / "c",
                          formats=f),
        plot.module_size_hist(d / "modules.txt", output=tmp_path / "d", formats=f),
        plot.null_distribution(d / "net.network", [d / "shuf.network"], output=tmp_path / "e",
                               formats=f),
        plot.expression_distribution(d / "expr.tsv", d / "kept.tsv", output=tmp_path / "f",
                                     formats=f),
        plot.mean_variance(d / "expr.tsv", kept=d / "kept.tsv", output=tmp_path / "g",
                           formats=f),
        plot.sample_correlation(d / "expr.tsv", output=tmp_path / "h", formats=f),
        plot.enrich_bar(d / "enrich.tsv", output=tmp_path / "i", formats=f),
        plot.enrich_dot(d / "annot", output=tmp_path / "j", formats=f),
        plot.enrich_heatmap(d / "annot", output=tmp_path / "k", formats=f),
        plot.rwr_rank(d / "rwr.tsv", seeds=d / "seeds.txt", output=tmp_path / "l", formats=f),
        plot.rwr_subnetwork(d / "net.network", d / "rwr.tsv", d / "seeds.txt",
                            output=tmp_path / "m", formats=f),
        plot.diff_overview(d / "diff.tsv", output=tmp_path / "n", formats=f),
        plot.diff_network(d / "diff.tsv", output=tmp_path / "o", formats=f),
        plot.trajectory_overview(d / "traj", output=tmp_path / "p", formats=f),
        plot.fusion_overview(d / "fused.tsv", output=tmp_path / "q", formats=f),
    ]
    assert len(list(tmp_path.glob("*.png"))) == len(figs)
    assert figs[0].axes[1].get_title().startswith("scale-free R²")
    labels = [t.get_text() for t in figs[8].axes[0].get_yticklabels()]
    assert labels[0] == "process 0"  # smallest p first


def test_errors(data):
    d = data
    with pytest.raises(GCENInputError):
        plot.edge_weights(d / "net.network", column="tom")
    with pytest.raises(GCENInputError):
        plot.trajectory_overview(d)
    with pytest.raises(GCENInputError):
        plot.fusion_overview(d / "diff.tsv")
    with pytest.raises(GCENInputError):
        plot.rwr_rank(pd.DataFrame({"gene": [], "probability": []}))


def test_cli(data, tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    d = data
    cases = [["degree-distribution", "-i", str(d / "net.network")],
             ["network-view", "-i", str(d / "net.network"), "-m", str(d / "modules.txt")],
             ["null-distribution", "-i", str(d / "net.network"), "-s",
              str(d / "shuf.network")],
             ["expression-qc", "-e", str(d / "expr.tsv"), "--normalized", str(d / "kept.tsv")],
             ["enrich-dot", "-i", str(d / "annot"), "--top", "2"],
             ["rwr-rank", "-i", str(d / "rwr.tsv"), "-g", str(d / "seeds.txt")],
             ["diff-overview", "-i", str(d / "diff.tsv")],
             ["trajectory", "-i", str(d / "traj")],
             ["fusion", "-i", str(d / "fused.tsv")]]
    for i, args in enumerate(cases):
        r = CliRunner().invoke(app, ["plot", *args, "-o", str(tmp_path / f"p{i}"), "--format",
                                     "png"])
        assert r.exit_code == 0, (args, r.output)
        assert (tmp_path / f"p{i}.png").exists()


def test_headless_with_unreachable_display(tmp_path):
    """No backend chosen and a DISPLAY nobody answers: figures are still
    drawn at once (Agg), instead of matplotlib probing the display."""
    env = {k: v for k, v in os.environ.items() if k != "MPLBACKEND"}
    env["DISPLAY"] = "localhost:97.0"
    code = ("import pandas as pd\nfrom gcen2 import plot\n"
            "plot.module_sizes(pd.Series(['ME0', 'ME1', 'ME0'], index=list('abc')), "
            f"output=r'{tmp_path / 'x.png'}')\n")
    subprocess.run([sys.executable, "-c", code], env=env, check=True, timeout=60)
    assert (tmp_path / "x.png").exists()


def test_merge_overlap_counts(tmp_path):
    def net(pairs, name):
        f = tmp_path / f"{name}.network"
        pd.DataFrame([(a, b, 0.8, 0.0) for a, b in pairs],
                     columns=["gene_a", "gene_b", "correlation", "p_value"]).to_csv(
            f, sep="\t", index=False, lineterminator="\n")
        return f

    a = net([("x", "y"), ("x", "z"), ("y", "z")], "a")
    b = net([("y", "x"), ("z", "w")], "b")  # (y, x) is (x, y)
    c = net([("x", "y"), ("w", "v")], "c")
    m = net([("x", "y")], "m")
    fig = plot.merge_overlap([a, b, c], merged=m)
    texts = [t.get_text() for t in fig.axes[0].texts]
    assert "1\n(1 kept)" in texts            # x-y in all three, kept
    assert "2\n(0 kept)" in texts            # x-z, y-z only in a
    assert fig.axes[0].get_title().startswith("5 distinct edges in 3 networks; 1 in")
    fig4 = plot.merge_overlap([a, b, c, m], names=list("abcm"))
    assert len(fig4.axes) == 2               # UpSet: bars and matrix
    with pytest.raises(GCENInputError):
        plot.merge_overlap([a])
    with pytest.raises(GCENInputError):
        plot.merge_overlap([a, b], names=["only one"])


def test_extract_overview(data):
    d = data
    fig = plot.extract_overview(d / "net.network", ["g0", "g1", "g2", "missing"], top=3)
    venn, bars = fig.axes[:2]
    assert venn.get_title() == "3 of 4 query genes found"
    assert len(bars.get_yticklabels()) == 3
    net = pd.read_csv(d / "net.network", sep="\t")
    deg = pd.concat([net["gene_a"], net["gene_b"]]).value_counts()
    top = bars.get_yticklabels()[0].get_text()
    width = sum(p.get_width() for p in bars.patches if abs(p.get_y() + p.get_height() / 2) < 1e-9)
    assert width == deg[top]


def test_moran_and_metacells(tmp_path):
    rng = np.random.default_rng(0)
    t = pd.DataFrame({"gene": [f"g{i}" for i in range(200)], "morans_i": rng.normal(0.05, 0.1, 200),
                      "mean": rng.lognormal(0, 1, 200), "detected": rng.random(200)})
    assert len(plot.moran(t).axes) == 2
    t["morans_i_smoothed"] = t["morans_i"] + 0.3
    assert len(plot.moran(t).axes) == 3
    cells = [f"c{i}" for i in range(60)]
    xy = pd.DataFrame(rng.normal(size=(60, 2)), index=cells, columns=["x", "y"])
    mem = pd.DataFrame({"metacell": np.repeat([f"m{i}" for i in range(6)], 5),
                        "cell": cells[:30], "center": np.tile([True] + [False] * 4, 6),
                        "cell_type": ["A"] * 15 + ["B"] * 15})
    fig = plot.metacells(mem, xy, n_cells={"A": 30, "B": 30})
    assert fig.axes[0].get_title() == "6 metacells of 30 cells"
    assert "3 (50% of cells)" in [x.get_text().strip() for x in fig.axes[1].texts]
    with pytest.raises(GCENInputError):
        plot.metacells(mem.drop(columns="cell"), xy)


def test_new_cli_commands(tmp_path, data):
    from typer.testing import CliRunner

    from gcen2.cli import app

    d = data
    pd.DataFrame({"gene": ["g0", "g1"], "morans_i": [0.4, -0.1], "mean": [1.0, 2.0],
                  "detected": [1.0, 0.5]}).to_csv(tmp_path / "mi.tsv", sep="\t", index=False, lineterminator="\n")
    cases = [["merge-overlap", "-i", str(d / "net.network"), "-i", str(d / "shuf.network"),
              "--name", "A", "--name", "B"],
             ["extract-overview", "-i", str(d / "net.network"), "-g", str(d / "seeds.txt")],
             ["moran", "-i", str(tmp_path / "mi.tsv")]]
    for i, args in enumerate(cases):
        r = CliRunner().invoke(app, ["plot", *args, "-o", str(tmp_path / f"n{i}"), "--format",
                                     "png"])
        assert r.exit_code == 0, (args, r.output)
        assert (tmp_path / f"n{i}.png").exists()


def test_enrichment_figures_from_annotate_go_results(data):
    """The dict annotate_go returns (EnrichmentResult values, module_<k> keys,
    modules without terms as empty tables) plots as the annotate directory."""
    from gcen2.models import EnrichmentResult

    d = data
    term = pd.read_csv(d / "enrich.tsv", sep="\t")
    def er(t, k):
        return EnrichmentResult(results=t, ontology_type="GO", gene_or_module=k)

    res = {"module_0": er(term, "module_0"), "module_1": er(pd.DataFrame(), "module_1"),
           "module_2": er(term.iloc[::-1], "module_2")}
    fig = plot.enrich_dot(res, top=2)
    labels = [t.get_text() for t in fig.axes[0].get_xticklabels()]
    assert labels == ["ME0 (turquoise)", "ME2 (brown)"]  # module 1: no terms
    assert plot.enrich_heatmap(res, top=2).axes


def test_enrichment_colour_bar_names_the_p_column(data):
    """Tables with an fdr column are coloured by it, and the colour bar says so."""
    term = pd.read_csv(data / "enrich.tsv", sep="\t")
    res = {"0": term.assign(fdr=term["p_value"] * 2), "1": term.assign(fdr=term["p_value"])}
    for fig in (plot.enrich_dot(res, top=2), plot.enrich_heatmap(res, top=2)):
        assert fig.axes[-1].get_ylabel() == "−log10 fdr"
    assert plot.enrich_dot({"0": term}, top=2).axes[-1].get_ylabel() == "−log10 p_value"
