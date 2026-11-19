"""benchmarks/README.md's overview table and benchmarks/overview.tsv match the
committed benchmark results (regenerate with `python benchmarks/overview.py`)."""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_overview_up_to_date():
    r = subprocess.run([sys.executable, str(ROOT / "benchmarks/overview.py"), "--check"],
                       capture_output=True, text=True, encoding="utf-8")
    assert r.returncode == 0, r.stdout + r.stderr
