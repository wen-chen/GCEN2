"""docs/cli.rst is generated from the commands (tools/make_cli_docs.py) and
must be up to date: a changed command or option needs the file regenerated."""

import importlib.util
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def _tool():
    spec = importlib.util.spec_from_file_location("make_cli_docs",
                                                  ROOT / "tools" / "make_cli_docs.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_cli_reference_up_to_date():
    text = _tool().render()
    assert (ROOT / "docs" / "cli.rst").read_text(encoding="utf-8") == text, \
        "docs/cli.rst is out of date: run python tools/make_cli_docs.py"


def test_cli_reference_lists_every_command():
    import typer

    from gcen2.cli import app, plot_app

    text = (ROOT / "docs" / "cli.rst").read_text(encoding="utf-8")
    for name in set(typer.main.get_command(app).commands) - {"plot"}:
        assert f"\ngcen2 {name}\n" in text, name
    for name in typer.main.get_command(plot_app).commands:
        assert f"\ngcen2 plot {name}\n" in text, name
    assert " \n" not in text  # no trailing spaces (pre-commit would strip them)
