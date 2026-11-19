"""Golden-output regression suite (tests/golden): every analysis command on
small fixed inputs must give the frozen outputs. This is the baseline for
performance work; see tests/golden/cases.py and make_golden.py."""

import json

import pytest

pytest.importorskip("anndata")

from gcen2.core import _HAVE_CORE  # noqa: E402
from tests.golden.cases import STEPS  # noqa: E402
from tests.golden.golden import EXPECTED, compare, read_expected, run_steps  # noqa: E402

pytestmark = pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")

MANIFEST = json.loads((EXPECTED / "manifest.json").read_text(encoding="utf-8"))


@pytest.fixture(scope="module")
def produced(tmp_path_factory):
    out = tmp_path_factory.mktemp("golden")
    return out, run_steps(out)


def test_every_step_has_expected_outputs():
    assert [s[0] for s in STEPS] == list(MANIFEST)
    assert all(MANIFEST[name] for name in MANIFEST), "a step writes no output"


@pytest.mark.parametrize("step", [s[0] for s in STEPS])
def test_step_outputs(produced, step):
    out, manifest = produced
    assert manifest[step] == MANIFEST[step], "different set of output files"
    diffs = []
    for rel in MANIFEST[step]:
        d = compare((out / rel).read_text(encoding="utf-8"), read_expected(rel), rel)
        if d:
            diffs.append(f"{rel}: {d}")
    assert not diffs, "\n".join(diffs)


def test_every_command_is_in_the_golden_suite():
    """A new analysis command fails this test until it has a golden step."""
    import typer

    from gcen2.cli import app

    commands = set(typer.main.get_command(app).commands) - {"plot", "report"}
    assert commands <= {args[0] for _, args, _ in STEPS}, sorted(
        commands - {args[0] for _, args, _ in STEPS})
