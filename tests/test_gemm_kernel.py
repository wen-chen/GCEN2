"""Tests for the matrix-product (gemm) Pearson / Spearman kernel of
network.build / build_to_file / build_from_file."""

import hashlib
import os
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

from gcen2.core import _HAVE_CORE
from gcen2.exceptions import GCENInputError

pytestmark = pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")


def _expr(n_genes=300, n_samples=200, seed=0):
    rng = np.random.default_rng(seed)
    f = rng.normal(size=(8, n_samples))
    x = rng.normal(size=(n_genes, 8)) @ f + rng.normal(size=(n_genes, n_samples)) * 2
    x[5] = 3.0             # constant gene: r = 0
    x[7, ::3] = 0.0        # ties for Spearman
    return pd.DataFrame(x, index=[f"g{i}" for i in range(n_genes)],
                        columns=[f"s{j}" for j in range(n_samples)])


def _key(net):
    e = net.edges
    return list(zip(e["gene_a"], e["gene_b"]))


@pytest.mark.parametrize("method", ["pearson", "spearman"])
def test_same_edges_as_pairwise(method):
    from gcen2.network import build

    df = _expr()
    a = build(df, method=method, cor=0.3, pval=0.01, kernel="pairwise", threads=2)
    b = build(df, method=method, cor=0.3, pval=0.01, kernel="gemm", threads=2)
    assert _key(a) == _key(b) and len(a.edges) > 1000
    np.testing.assert_allclose(b.edges["correlation"], a.edges["correlation"], rtol=0,
                               atol=1e-13)
    np.testing.assert_allclose(b.edges["p_value"], a.edges["p_value"], rtol=1e-9)
    full_a = build(df, method=method, all_edges=True, kernel="pairwise").edges
    full_b = build(df, method=method, all_edges=True, kernel="gemm").edges
    assert len(full_b) == 300 * 299 // 2 and _key_df(full_a) == _key_df(full_b)
    const = full_b[(full_b["gene_a"] == "g5") | (full_b["gene_b"] == "g5")]
    assert (const["correlation"] == 0).all()
    np.testing.assert_allclose(full_b["correlation"], full_a["correlation"], atol=1e-13)


def _key_df(e):
    return list(zip(e["gene_a"], e["gene_b"]))


def test_independent_of_threads_and_tiles():
    from gcen2.network import build

    df = _expr(n_genes=500, n_samples=150, seed=1)
    ref = build(df, method="pearson", all_edges=True, kernel="gemm", threads=1).edges
    for t in (3, 8):
        e = build(df, method="pearson", all_edges=True, kernel="gemm", threads=t).edges
        assert e["correlation"].to_numpy().tobytes() == ref["correlation"].to_numpy().tobytes()


def test_auto_and_errors():
    from gcen2.core import _gcen_core
    from gcen2.network import GEMM_MIN_SAMPLES, build

    assert GEMM_MIN_SAMPLES == _gcen_core.correlation.GEMM_MIN_SAMPLES
    for ns, kernel in ((GEMM_MIN_SAMPLES - 1, "pairwise"), (GEMM_MIN_SAMPLES, "gemm")):
        df = _expr(n_genes=80, n_samples=ns, seed=2)
        auto = build(df, method="pearson", all_edges=True).edges["correlation"].to_numpy()
        want = build(df, method="pearson", all_edges=True,
                     kernel=kernel).edges["correlation"].to_numpy()
        assert auto.tobytes() == want.tobytes()
    with pytest.raises(GCENInputError, match="kernel"):
        build(_expr(n_genes=10), kernel="blas")


def test_streamed_and_file_paths(tmp_path):
    from gcen2._build_file import build_from_file
    from gcen2.network import build, build_to_file

    df = _expr(n_genes=200, n_samples=120, seed=3)
    want = tmp_path / "want.network"
    build(df, method="spearman", cor=0.3, kernel="gemm").to_edgelist(want)
    got = tmp_path / "got.network"
    build_to_file(df, got, method="spearman", cor=0.3, kernel="gemm",
                  max_pairs_per_block=1000)
    assert got.read_bytes() == want.read_bytes()
    src = tmp_path / "x.tsv"
    df.to_csv(src, sep="\t", float_format="%.17g", lineterminator="\n")
    got2 = tmp_path / "got2.network"
    build_from_file(src, got2, method="spearman", cor=0.3, kernel="gemm")
    assert got2.read_bytes() == want.read_bytes()
    with pytest.raises(GCENInputError, match="kernel"):
        build_from_file(src, got2, kernel="x")


def test_same_on_baseline_simd():
    """The AVX and baseline product kernels give the same bits."""
    code = (
        "import hashlib, numpy as np\n"
        "from gcen2.core import _gcen_core as c\n"
        "rng = np.random.default_rng(4); x = rng.normal(size=(150, 333))\n"
        "r = c.correlation.network_edges(x, 'pearson', 3, 0.0, 1.0, False, False, True,"
        " kernel='gemm')[1]\n"
        "print(hashlib.sha256(r.tobytes()).hexdigest())\n")
    out = {}
    for simd in ("", "baseline"):
        env = {**os.environ, "GCEN2_SIMD": simd}
        out[simd] = subprocess.run([sys.executable, "-c", code], env=env, check=True,
                                   capture_output=True, text=True, encoding="utf-8").stdout.strip()
    from gcen2.core import _gcen_core as c

    x = np.random.default_rng(4).normal(size=(150, 333))
    r = c.correlation.network_edges(x, "pearson", 1, 0.0, 1.0, False, False, True,
                                    kernel="gemm")[1]
    assert out[""] == out["baseline"] == hashlib.sha256(r.tobytes()).hexdigest()


def test_cli_kernel(tmp_path):
    from typer.testing import CliRunner

    from gcen2.cli import app

    df = _expr(n_genes=60, n_samples=100, seed=5)
    src = tmp_path / "x.tsv"
    df.to_csv(src, sep="\t", float_format="%.17g", lineterminator="\n")
    outs = {}
    for k in ("pairwise", "gemm"):
        outs[k] = tmp_path / f"{k}.network"
        r = CliRunner().invoke(app, ["build", "-i", str(src), "-o", str(outs[k]), "-c", "0.3",
                                     "--kernel", k])
        assert r.exit_code == 0, r.output
    a = pd.read_csv(outs["pairwise"], sep="\t")
    b = pd.read_csv(outs["gemm"], sep="\t")
    assert list(a.iloc[:, 0]) == list(b.iloc[:, 0]) and len(a) > 0
    bad = CliRunner().invoke(app, ["build", "-i", str(src), "-o", str(outs["gemm"]),
                                   "--kernel", "x"])
    assert bad.exit_code == 1 and "kernel" in bad.output
