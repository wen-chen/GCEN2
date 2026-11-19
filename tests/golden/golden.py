"""Run the golden steps (cases.STEPS) and compare their outputs."""

from __future__ import annotations

import gzip
import json
import math
import os
from pathlib import Path

from tests.golden.cases import SKIP_SUFFIXES, STEPS

HERE = Path(__file__).parent
INPUTS = HERE / "inputs"
EXPECTED = HERE / "expected"
RTOL, ATOL = 1e-9, 1e-12
# outputs whose numbers depend on the library's summation beyond RTOL: the
# maximum adjacency ratio (MAR) statistics of module preservation for the grey
# module (tiny adjacencies; NumPy 1.24 vs 2.x differ by ~2e-7 relative)
RTOL_BY_PREFIX = {"preservation/": 1e-6}


def _files(root: Path) -> dict:
    return {p.relative_to(root).as_posix(): p.stat().st_mtime_ns for p in root.rglob("*")
            if p.is_file() and not p.name.endswith(SKIP_SUFFIXES)}


def run_steps(out: Path, inputs: Path = INPUTS, only=None) -> dict:
    """Run every step (in order) into ``out``; return {step: [files it wrote]}
    and, for steps with a stdout file, write it."""
    from typer.testing import CliRunner

    from gcen2.cli import app

    os.environ.setdefault("MPLBACKEND", "Agg")
    out.mkdir(parents=True, exist_ok=True)
    manifest = {}
    runner = CliRunner()
    for name, args, stdout in STEPS:
        before = _files(out)
        argv = [a.replace("{I}", str(inputs)).replace("{O}", str(out)) for a in args]
        r = runner.invoke(app, argv, catch_exceptions=False)
        if r.exit_code != 0:
            raise RuntimeError(f"golden step {name} failed:\n{r.output}")
        if stdout is not None:
            # absolute paths of this run are not part of the result
            (out / stdout).write_text(r.output.replace(str(out), "{O}")
                                      .replace(str(inputs), "{I}"), encoding="utf-8")
        after = _files(out)
        manifest[name] = sorted(f for f, t in after.items() if before.get(f) != t)
        if only is not None and name == only:
            break
    return manifest


def _num(tok: str):
    try:
        return float(tok)
    except ValueError:
        return None


def _same_token(a: str, b: str) -> bool:
    if a == b:
        return True
    x, y = _num(a), _num(b)
    if x is None or y is None:
        return False
    if math.isnan(x) and math.isnan(y):
        return True
    return math.isclose(x, y, rel_tol=RTOL, abs_tol=ATOL)


def _same_json(a, b, path="") -> str | None:
    if isinstance(a, dict) and isinstance(b, dict):
        if a.keys() != b.keys():
            return f"{path}: keys {sorted(a)} != {sorted(b)}"
        for k in a:
            d = _same_json(a[k], b[k], f"{path}.{k}")
            if d:
                return d
        return None
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return f"{path}: length {len(a)} != {len(b)}"
        for i, (x, y) in enumerate(zip(a, b)):
            d = _same_json(x, y, f"{path}[{i}]")
            if d:
                return d
        return None
    if isinstance(a, (int, float)) and isinstance(b, (int, float)) and \
            not isinstance(a, bool) and not isinstance(b, bool):
        return None if _same_token(repr(float(a)), repr(float(b))) else f"{path}: {a} != {b}"
    return None if a == b else f"{path}: {a!r} != {b!r}"


def read_expected(rel: str) -> str:
    with gzip.open(EXPECTED / (rel + ".gz"), "rt", encoding="utf-8") as fh:
        return fh.read()


def compare(got_text: str, want_text: str, rel: str) -> str | None:
    """None if the same, else a description of the first difference."""
    if rel.endswith(".json"):
        return _same_json(json.loads(got_text), json.loads(want_text))
    global RTOL
    rtol = next((t for p, t in RTOL_BY_PREFIX.items() if rel.startswith(p)), None)
    if rtol is not None and rtol > RTOL:
        saved, RTOL = RTOL, rtol
        try:
            return compare(got_text, want_text, "")
        finally:
            RTOL = saved
    got, want = got_text.splitlines(), want_text.splitlines()
    if len(got) != len(want):
        return f"{len(got)} lines, expected {len(want)}"
    for i, (g, w) in enumerate(zip(got, want), 1):
        if g == w:
            continue
        gt, wt = g.split(), w.split()
        if len(gt) != len(wt) or not all(_same_token(x, y) for x, y in zip(gt, wt)):
            return f"line {i}:\n  got:      {g[:200]}\n  expected: {w[:200]}"
    return None


def assert_same_file(got: Path, want: Path, rtol: float = 1e-12) -> None:
    """Two outputs of different code paths (compiled kernel and the pandas
    reference): byte for byte with pandas 3; with pandas < 3 the reference
    path's float parsing, groupby sums and float formatting differ in the
    last digit, so numbers are compared to ``rtol`` (text still exactly)."""
    import pandas as pd

    a, b = Path(got).read_bytes(), Path(want).read_bytes()
    if int(pd.__version__.split(".")[0]) >= 3:
        assert a == b
        return
    global RTOL
    saved, RTOL = RTOL, rtol
    try:
        d = compare(a.decode(), b.decode(), str(got))
    finally:
        RTOL = saved
    assert d is None, d
