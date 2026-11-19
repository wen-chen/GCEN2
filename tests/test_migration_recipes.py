"""The Python code that docs/migration.rst gives for GCEN v1 utilities without
a gcen2 command gives the values of the v1 tools.

The code blocks are taken from the guide itself, so the guide cannot drift from
what is tested. Expected values: GCEN v1 0.6.3's data_stat on the sample
matrix, and its generate_expr_matrix_from_{rsem,stringtie} (FPKM and -t TPM)
on the small inputs in tests/data/v1_util/ (made up for this test).
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
DATA = Path(__file__).resolve().parent / "data" / "v1_util"


def _guide_code() -> dict:
    """Execute the python code blocks of the migration guide; their names."""
    text = (ROOT / "docs" / "migration.rst").read_text(encoding="utf-8")
    blocks = re.findall(r"\.\. code-block:: python\n\n((?:   .*\n|\n)+)", text)
    namespace: dict = {}
    for block in blocks:
        code = "\n".join(line[3:] for line in block.splitlines())
        if "def data_stat" in code or "def matrix_from_" in code:
            exec(compile(code, "migration.rst", "exec"), namespace)
    return namespace


@pytest.fixture(scope="module")
def guide():
    ns = _guide_code()
    for name in ("data_stat", "matrix_from_rsem", "matrix_from_stringtie"):
        assert name in ns, f"{name} not found in docs/migration.rst"
    return ns


def _numbers(text: str) -> list[float]:
    """The values printed, without the percentile labels ("5%:")."""
    text = re.sub(r"\d+%:", "", text)
    return [float(x) for x in re.findall(r"-?\d+(?:\.\d+)?(?:e[-+]?\d+)?", text)]


SAMPLE = ROOT / "data" / "sample" / "gene_expr.tsv"


@pytest.mark.skipif(not SAMPLE.exists(), reason="GCEN v1 sample matrix not present "
                    "(data/sample is not in the source distribution)")
def test_data_stat_matches_v1(guide, capsys):
    guide["data_stat"](SAMPLE)
    ours = _numbers(capsys.readouterr().out)
    v1 = _numbers((DATA / "v1_data_stat.txt").read_text(encoding="utf-8"))
    assert len(ours) == len(v1) == 16
    np.testing.assert_allclose(ours, v1, rtol=1e-12)


@pytest.mark.parametrize("tool, value, expected", [
    ("matrix_from_rsem", "FPKM", "v1_rsem_fpkm.tsv"),
    ("matrix_from_rsem", "TPM", "v1_rsem_tpm.tsv"),
    ("matrix_from_stringtie", "FPKM", "v1_stringtie_fpkm.tsv"),
    ("matrix_from_stringtie", "TPM", "v1_stringtie_tpm.tsv"),
])
def test_expression_matrix_matches_v1(guide, monkeypatch, tool, value, expected):
    monkeypatch.chdir(DATA)  # the sample lists name files relative to it
    sample_list = "rsem_sample.txt" if "rsem" in tool else "stringtie_sample.txt"
    ours = guide[tool](sample_list, value=value).sort_index()
    v1 = pd.read_csv(DATA / expected, sep="\t", index_col=0).sort_index()
    assert list(ours.columns) == list(v1.columns) == ["S1", "S2", "S3"]
    assert list(ours.index) == list(v1.index)
    # v1 prints 6 significant digits
    np.testing.assert_allclose(ours.to_numpy(), v1.to_numpy(), rtol=1e-5, atol=1e-6)
