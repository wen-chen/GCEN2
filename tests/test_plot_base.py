"""The plotting infrastructure: style, module colours, saving, CLI."""

import pandas as pd
import pytest

pytest.importorskip("matplotlib")

from gcen2.exceptions import GCENInputError  # noqa: E402
from gcen2.plot import (  # noqa: E402
    figsize,
    module_color,
    module_label,
    module_sizes,
    save_figure,
    set_defaults,
    style,
)


@pytest.fixture(autouse=True)
def _defaults():
    yield
    set_defaults(palette="wgcna", formats=("pdf", "png", "svg"))


def test_sizes():
    assert figsize("single")[0] == pytest.approx(89 / 25.4)
    assert figsize("double", height=50) == pytest.approx((183 / 25.4, 50 / 25.4))
    assert figsize(100, aspect=10)[1] == pytest.approx(170 / 25.4)  # capped
    with pytest.raises(GCENInputError):
        figsize("huge")


def test_module_colors():
    assert module_color("ME0") == "#40E0D0"  # turquoise, as WGCNA labels2colors(1)
    assert module_color(1) == "#0000FF"
    assert module_color("brown") == "#A52A2A"
    assert module_color("MEyellow") == "#FFFF00"
    assert module_color("grey") == module_color("ME_unassigned") == "#BEBEBE"
    assert module_color("gold") == "#FFD700"
    assert module_label("ME2") == "ME2 (brown)"
    assert module_label("ME_unassigned") == "unassigned"
    set_defaults(palette="okabe-ito")
    assert module_color("ME0") == "#E69F00"
    assert module_label("ME2") == "ME2"
    with pytest.raises(GCENInputError):
        set_defaults(palette="rainbow")
    with pytest.raises(GCENInputError):
        set_defaults(formats=["jpg"])


def test_style_is_local():
    import matplotlib as mpl

    before = mpl.rcParams["font.size"]
    with style():
        assert mpl.rcParams["font.size"] == 7
        assert mpl.rcParams["pdf.fonttype"] == 42
    assert mpl.rcParams["font.size"] == before


def test_save_formats(tmp_path):
    mods = pd.Series(["ME0"] * 5 + ["ME1"] * 3 + ["grey"] * 2, index=[f"g{i}" for i in range(10)])
    fig = module_sizes(mods, output=tmp_path / "a")
    assert sorted(p.name for p in tmp_path.iterdir()) == ["a.pdf", "a.png", "a.svg"]
    assert len(fig.axes[0].patches) == 2  # unassigned left out
    save_figure(fig, tmp_path / "b.tiff")
    assert (tmp_path / "b.tiff").exists()
    paths = save_figure(fig, tmp_path / "c", formats=["svg"])
    assert [p.name for p in paths] == ["c.svg"]
    svg = (tmp_path / "c.svg").read_text(encoding="utf-8")
    assert "<text" in svg  # text kept as text
    pdf = (tmp_path / "a.pdf").read_bytes()
    assert b"/FontFile2" in pdf  # TrueType embedded: editable text


def test_cli(tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    m = tmp_path / "modules.tsv"
    pd.DataFrame({"module": [0, 0, 1], "gene": ["a", "b", "c"]}).to_csv(m, sep="\t", index=False, lineterminator="\n")
    r = CliRunner().invoke(app, ["plot", "module-sizes", "-m", str(m), "-o",
                                 str(tmp_path / "s"), "--format", "png", "--size", "120"])
    assert r.exit_code == 0, r.output
    assert (tmp_path / "s.png").exists() and not (tmp_path / "s.pdf").exists()
    r = CliRunner().invoke(app, ["plot", "module-sizes", "-m", str(tmp_path / "no.tsv"), "-o",
                                 str(tmp_path / "x")])
    assert r.exit_code == 1
