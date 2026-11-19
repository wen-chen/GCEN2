"""Public functions that no other test runs: the gcen2.plotting compatibility
API, PDF reports, and small helpers of io, models, config and the package."""

import numpy as np
import pandas as pd
import pytest

from gcen2.exceptions import GCENIOError
from gcen2.models import ModuleSet


@pytest.fixture
def mpl():
    matplotlib = pytest.importorskip("matplotlib")
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    yield plt
    plt.close("all")


class TestLegacyPlotting:
    def test_network_interactive(self, tmp_path, small_network, monkeypatch):
        pytest.importorskip("pyvis")
        from gcen2 import plotting

        monkeypatch.chdir(tmp_path)

        mods = ModuleSet(modules=[{"GeneA", "GeneB", "GeneC"}, {"GeneD", "GeneE"}])
        for color_by, m in (("module", mods), ("degree", None)):
            out = tmp_path / f"{color_by}.html"
            plotting.network_interactive(small_network, output=out, color_by=color_by,
                                         modules=m)
            html = out.read_text(encoding="utf-8")
            assert "GeneA" in html and "GeneE" in html
        assert not (tmp_path / "lib").exists()  # self-contained HTML, nothing else written

    def test_expression_violin_and_heatmap(self, mpl, tmp_path, small_expr):
        pytest.importorskip("seaborn")
        from gcen2 import plotting

        fig = plotting.expression_violin(small_expr, output=tmp_path / "v.png")
        assert (tmp_path / "v.png").stat().st_size > 1000 and fig.axes
        mods = ModuleSet(modules=[{"Gene0", "Gene1", "Gene2"}, {"Gene5", "Gene6"}])
        fig = plotting.heatmap(small_expr, mods, output=tmp_path / "h.png")
        assert (tmp_path / "h.png").stat().st_size > 1000

    def test_eigengene_and_trait_heatmaps(self, mpl, tmp_path):
        from gcen2 import plotting

        rng = np.random.default_rng(0)
        me = pd.DataFrame(rng.normal(size=(12, 3)), columns=["ME0", "ME1", "ME_unassigned"])
        fig = plotting.eigengene_heatmap(me, output=tmp_path / "e.png")
        texts = [t.get_text() for t in fig.axes[0].texts]
        assert len(texts) == 4  # 2 x 2 correlations; ME_unassigned left out
        assert "1.00" in texts
        tc = pd.DataFrame({"module": ["ME0", "ME0", "ME1", "ME1"],
                           "trait": ["age", "bmi"] * 2, "cor": [0.5, -0.2, 0.1, 0.9],
                           "p": [0.01, 0.3, 0.7, 1e-6]})
        fig = plotting.module_trait_heatmap(tc, output=tmp_path / "t.png")
        assert any("0.90" in t.get_text() for t in fig.axes[0].texts)

    def test_spatial_eigengenes(self, mpl, tmp_path):
        from gcen2 import plotting

        rng = np.random.default_rng(1)
        coords = rng.random((50, 2)) * 100
        me = pd.DataFrame(rng.normal(size=(50, 2)), columns=["ME0", "ME1"])
        fig = plotting.spatial_eigengenes(coords, me, output=tmp_path / "s.png")
        assert [a.get_title() for a in fig.axes if a.get_title()][:2] == ["ME0", "ME1"]

    def test_diff_network(self, mpl, tmp_path, small_network):
        from gcen2 import plotting
        from gcen2.differential import compare
        from gcen2.models import Network

        b = Network(edges=small_network.edges.assign(
            correlation=[0.2, 0.85, 0.8, 0.95, 0.7]), is_signed=False, method="pearson")
        diff = compare(small_network, b, n_a=30, n_b=30)
        fig = plotting.diff_network(diff, output=tmp_path / "d.png")
        assert (tmp_path / "d.png").stat().st_size > 1000 and fig.axes


def test_pdf_report(tmp_path, small_network):
    pytest.importorskip("jinja2")
    try:
        import weasyprint  # noqa: F401
    except (ImportError, OSError):  # not installed, or no Pango
        pytest.skip("weasyprint (with its system libraries) not available")
    import gcen2

    out = tmp_path / "r.pdf"
    gcen2.report.generate(network=small_network, output=out)
    data = out.read_bytes()
    assert data[:5] == b"%PDF-" and len(data) > 5000


def test_read_h5ad(tmp_path):
    ad = pytest.importorskip("anndata")
    from gcen2.io import read_h5ad

    x = np.arange(12, dtype=np.float32).reshape(3, 4)  # cells x genes
    a = ad.AnnData(X=x, obs=pd.DataFrame(index=["c0", "c1", "c2"]),
                   var=pd.DataFrame(index=["g0", "g1", "g2", "g3"]))
    a.write_h5ad(tmp_path / "a.h5ad")
    expr = read_h5ad(tmp_path / "a.h5ad")
    assert list(expr.gene_names) == ["g0", "g1", "g2", "g3"]
    assert list(expr.sample_names) == ["c0", "c1", "c2"]
    assert expr.is_single_cell and expr.data.shape == (3, 4)  # the AnnData itself
    np.testing.assert_array_equal(np.asarray(expr.data.X), x)
    with pytest.raises(GCENIOError):
        read_h5ad(tmp_path / "missing.h5ad")


def test_write_csv(tmp_path):
    from gcen2.io import write_csv

    df = pd.DataFrame({"a": [1, 2]}, index=["x", "y"])
    write_csv(df, tmp_path / "t.csv")
    assert (tmp_path / "t.csv").read_text(encoding="utf-8").splitlines() == [",a", "x,1", "y,2"]
    with pytest.raises(GCENIOError):
        write_csv(df, tmp_path / "no" / "dir" / "t.csv")


def test_csv_to_tsv_round_trip_with_quoted_commas(tmp_path):
    """The migration guide's replacement for v1 csv_to_tsv / tsv_to_csv:
    quoted fields with commas stay one field; header-less v1 matrices too."""
    from gcen2.io import read_csv, read_tsv, write_csv, write_tsv

    (tmp_path / "h.csv").write_text('gene,s1,s2\n"A,1",1,2\nB,3.5,4\n', encoding="utf-8")
    m = read_csv(tmp_path / "h.csv")
    assert list(m.data.index) == ["A,1", "B"] and list(m.data.columns) == ["s1", "s2"]
    write_tsv(m.data, tmp_path / "h.tsv")
    assert (tmp_path / "h.tsv").read_text(encoding="utf-8").splitlines()[1].split("\t")[0] == "A,1"
    back = read_tsv(tmp_path / "h.tsv")
    write_csv(back.data, tmp_path / "back.csv")
    assert read_csv(tmp_path / "back.csv").data.equals(m.data)
    (tmp_path / "v1.csv").write_text('"G,1",1.5,2\nG2,3,4\n', encoding="utf-8")  # v1: no header line
    v1 = read_csv(tmp_path / "v1.csv")
    assert list(v1.data.index) == ["G,1", "G2"] and v1.data.shape == (2, 2)


def test_network_to_adjacency(small_network):
    from gcen2.exceptions import GCENInputError
    from gcen2.models import Network

    e = small_network.edges
    order = list(small_network.to_networkx().nodes)  # default order: first appearance
    adj = small_network.to_adjacency()
    assert adj.shape == (5, 5) and adj.nnz == 2 * len(e)
    dense = adj.toarray()
    np.testing.assert_array_equal(dense, dense.T)
    for a, b, r in e[["gene_a", "gene_b", "correlation"]].itertuples(index=False):
        assert dense[order.index(a), order.index(b)] == r  # weighted by the correlation

    neg = Network(edges=e.assign(correlation=-e["correlation"]), is_signed=False,
                  method="pearson")
    assert neg.to_adjacency().toarray().min() == -0.9  # sign kept
    np.testing.assert_array_equal(neg.to_adjacency(absolute=True).toarray(), dense)
    np.testing.assert_array_equal(small_network.to_adjacency(weight=None).toarray(),
                                  (dense != 0).astype(float))
    np.testing.assert_array_equal(small_network.to_adjacency(weight="p_value").toarray() != 0,
                                  dense != 0)
    # chosen gene order; unknown genes give empty rows; other genes' edges left out
    sub = small_network.to_adjacency(genes=["GeneC", "GeneA", "GeneZ"]).toarray()
    np.testing.assert_array_equal(sub, [[0, 0.85, 0], [0.85, 0, 0], [0, 0, 0]])
    # repeated edges are summed, a self-loop counted once
    rep = Network(edges=pd.DataFrame({"gene_a": ["x", "y", "x"], "gene_b": ["y", "x", "x"],
                                      "correlation": [0.5, 0.25, 0.1]}),
                  is_signed=False, method="pearson")
    np.testing.assert_allclose(rep.to_adjacency().toarray(), [[0.1, 0.75], [0.75, 0]])
    with pytest.raises(GCENInputError, match="tom"):
        small_network.to_adjacency(weight="tom")


def test_config_temp_dir(tmp_path):
    import tempfile

    from gcen2.config import GCENConfig

    assert GCENConfig().get_temp_dir() == tempfile.gettempdir()
    assert GCENConfig(temp_dir=str(tmp_path)).get_temp_dir() == str(tmp_path)


def test_every_subpackage_is_an_attribute_of_gcen2():
    """After a bare ``import gcen2`` (fresh interpreter), every public
    submodule is reachable as gcen2.<name> and listed by dir()."""
    import subprocess
    import sys
    from pathlib import Path

    import gcen2

    root = Path(gcen2.__file__).parent
    public = sorted(p.stem if p.is_file() else p.name for p in root.iterdir()
                    if not p.name.startswith(("_", ".")) and p.name not in ("cli", "templates")
                    and (p.suffix == ".py" or (p / "__init__.py").is_file()))
    code = ("import gcen2, sys\n"
            f"names = {public!r}\n"
            "bad = [n for n in names if not hasattr(gcen2, n) or n not in dir(gcen2)]\n"
            "print(bad); sys.exit(1 if bad else 0)\n")
    r = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8")
    assert r.returncode == 0, r.stdout + r.stderr
    assert "plot" in public and "single_cell" in public
    assert set(public) <= set(dir(gcen2))


def test_missing_output_directories_are_created(tmp_path):
    """Outputs in a directory that does not exist yet: it is created (no
    traceback), for every command that writes files."""
    pytest.importorskip("anndata")
    import inspect

    from typer.testing import CliRunner

    from gcen2.cli import app, plot_app
    from tests.golden.golden import INPUTS

    for info in [*app.registered_commands, *plot_app.registered_commands]:
        params = inspect.signature(info.callback).parameters
        if "output" in params:
            assert hasattr(info.callback, "__wrapped__"), info.name or info.callback.__name__
    i, n = INPUTS, tmp_path / "new" / "dir"
    runner = CliRunner()
    cases = [["normalize", "-i", f"{i}/bulk.tsv", "-o", f"{n}/a/norm.tsv"],
             ["build", "-i", f"{i}/bulk_small.tsv", "-o", f"{n}/b/x.network", "-c", "0.9",
              "-l", "log2"],
             ["module", "-e", f"{i}/set1.tsv", "-o", f"{n}/c/m.tsv", "-l", "log2", "--power",
              "6", "--me", f"{n}/d/me.tsv"],
             ["enrich", "-e", f"{i}/query.txt", "-b", f"{i}/background.txt", "-g",
              f"{i}/go.obo", "-a", f"{i}/gene_go.assoc", "-o", f"{n}/e/enrich.tsv"]]
    for args in cases:
        r = runner.invoke(app, args)
        assert r.exit_code == 0, (args, r.output)
    for f in ("a/norm.tsv", "b/x.network", "c/m.tsv", "d/me.tsv", "e/enrich.tsv"):
        assert (n / f).is_file(), f


def test_cli_annotations_work_on_python_39():
    """typer evaluates the annotations of every command when the CLI is built:
    `X | None` there breaks every gcen2 command on Python 3.9 (use Optional)."""
    import ast
    from pathlib import Path

    import gcen2.cli

    bad = []
    for path in sorted(Path(gcen2.cli.__file__).parent.glob("*.py")):
        for fn in ast.walk(ast.parse(path.read_text(encoding="utf-8"))):
            if isinstance(fn, ast.FunctionDef):
                for a in [*fn.args.args, *fn.args.kwonlyargs]:
                    if a.annotation is not None and any(
                            isinstance(n, ast.BinOp) and isinstance(n.op, ast.BitOr)
                            for n in ast.walk(a.annotation)):
                        bad.append(f"{path.name}: {fn.name}({a.arg}) line {a.lineno}")
    assert not bad, bad
