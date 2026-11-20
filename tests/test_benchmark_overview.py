"""benchmarks/README.md's overview table and benchmarks/overview.tsv match the
committed benchmark results (regenerate with `python benchmarks/overview.py`)."""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.skipif(not (ROOT / "benchmarks" / "overview.py").exists(),
                    reason="benchmarks/ is not in the source distribution")
def test_overview_up_to_date():
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks/overview.py"), "--check"],
                       capture_output=True, text=True, encoding="utf-8")
    assert r.returncode == 0, r.stdout + r.stderr
