"""gcen2._build_file.build_from_file: the numpy-free C++ build path must
write the same file as read_matrix + build_to_file, and fall back to that
path for inputs outside the plain subset."""

import subprocess
import sys

import numpy as np
import pytest

from gcen2._build_file import build_from_file
from gcen2._matrix_io import read_matrix
from gcen2.core import _HAVE_CORE
from gcen2.exceptions import GCENInputError, GCENIOError
from gcen2.network import build_to_file

pytestmark = pytest.mark.skipif(not _HAVE_CORE, reason="compiled extension not available")


def _core():
    from gcen2.core import _gcen_core
    return _gcen_core


def _reference(src, out, header=None, **kw):
    m, g, _ = read_matrix(src, header=header)
    return build_to_file(m, out, gene_names=g, **kw)


def _write_matrix(path, n_genes=60, n_samples=8, seed=0, header=True, eol="\n",
                  ragged=False, names=None):
    rng = np.random.default_rng(seed)
    base = rng.normal(size=(6, n_samples))
    vals = base[rng.integers(0, 6, n_genes)] + 0.3 * rng.normal(size=(n_genes, n_samples))
    names = names or [f"G{i}" for i in range(n_genes)]
    lines = []
    if header:
        lines.append("gene\t" + "\t".join(f"s{j}" for j in range(n_samples)))
    for i, row in enumerate(vals):
        fields = [repr(float(v)) if i % 3 else f"{v:.4f}" for v in row]
        if ragged and i % 7 == 3:
            fields = fields[:-2]
        lines.append(names[i] + "\t" + "\t".join(fields))
    path.write_bytes((eol.join(lines) + eol).encode())
    return path


def _same(tmp_path, src, header=None, expect_cpp=True, **kw):
    a, b = tmp_path / "a.net", tmp_path / "b.net"
    ra = _reference(src, a, header=header, **kw)
    rb = build_from_file(src, b, header=header, **kw)
    assert ra == rb
    assert a.read_bytes() == b.read_bytes()
    if expect_cpp is not None:
        has_header = header if header is not None else __import__(
            "gcen2._matrix_io", fromlist=["x"])._detect_header(src, "\t")
        cpp = _core().correlation.build_network_file(
            str(src), str(tmp_path / "c.net"), bool(has_header), kw.get("method", "spearman"),
            kw.get("cor", 0.1), kw.get("pval", 0.001), kw.get("signed", False),
            kw.get("all_edges", False), kw.get("threads", 2),
            kw.get("max_pairs_per_block", 500_000), kw.get("fdr", False))
        assert (cpp is not None) == expect_cpp
    return ra


class TestSameAsNumpyPath:
    @pytest.mark.parametrize("method", ["spearman", "pearson"])
    @pytest.mark.parametrize("opts", [
        dict(cor=0.5, pval=0.05),
        dict(cor=0.5, pval=0.05, signed=True),
        dict(cor=0.0, pval=1.0),
        dict(all_edges=True),
        dict(cor=0.3, pval=0.2, max_pairs_per_block=100, threads=3),
        dict(cor=0.5, pval=0.05, fdr=True),
        dict(cor=0.3, pval=0.2, signed=True, fdr=True, max_pairs_per_block=100, threads=3),
        dict(all_edges=True, fdr=True),
    ])
    def test_options(self, tmp_path, method, opts):
        src = _write_matrix(tmp_path / "m.tsv")
        info = _same(tmp_path, src, method=method, **opts)
        assert info["edge_count"] > 0

    @pytest.mark.parametrize("header", [True, False])
    @pytest.mark.parametrize("eol", ["\n", "\r\n"])
    def test_header_eol_ragged(self, tmp_path, header, eol):
        src = _write_matrix(tmp_path / "m.tsv", header=header, eol=eol, ragged=True)
        # cor=0 keeps NaN-free pairs only; padded rows give NaN correlations
        _same(tmp_path, src, cor=0.0, pval=1.0, method="pearson")
        _same(tmp_path, src, all_edges=True, method="pearson")

    def test_header_one_field_longer_and_blank_lines(self, tmp_path):
        src = _write_matrix(tmp_path / "m.tsv")
        text = src.read_text(encoding="utf-8").splitlines()
        text[0] = "\t".join(text[0].split("\t")[1:])  # header names samples only
        text.insert(5, "")
        text.insert(9, "  \t ")
        src.write_text("\n".join(text), encoding="utf-8")  # no final newline
        _same(tmp_path, src, header=True, cor=0.5, pval=0.05)

    def test_header_longer_than_rows(self, tmp_path):
        # extra header fields widen the matrix; rows are padded with NaN
        src = _write_matrix(tmp_path / "m.tsv")
        lines = src.read_text(encoding="utf-8").splitlines()
        lines[0] += "\textra\tmore"
        src.write_text("\n".join(lines) + "\n", encoding="utf-8")
        _same(tmp_path, src, header=True, all_edges=True)

    def test_quoted_and_unicode_gene_names(self, tmp_path):
        names = [f'"Q{i}"' if i % 2 else f"基因-{i}" for i in range(40)]
        src = _write_matrix(tmp_path / "m.tsv", n_genes=40, names=names)
        _same(tmp_path, src, cor=0.5, pval=0.05)

    def test_number_forms(self, tmp_path):
        rows = ["g1\t1\t2.5\t-3e2\t.5\t4.\t1E-3", "g2\t2\t3.5\t-2e2\t.75\t5.\t2E-3",
                "g3\t-0\t0.0\t1e+2\t0.25\t3\t5e-1", "g4\t7\t1\t2\t3\t4\t5"]
        src = tmp_path / "m.tsv"
        src.write_text("\n".join(rows) + "\n", encoding="utf-8")
        _same(tmp_path, src, header=False, all_edges=True)

    def test_file_without_edges(self, tmp_path):
        src = _write_matrix(tmp_path / "m.tsv", n_genes=10)
        info = _same(tmp_path, src, cor=0.9999, pval=1e-12)
        assert info == {"node_count": 0, "edge_count": 0}


class TestBlocks:
    """Blocks run on one reused thread pool; rows longer than a work unit
    (256 pairs) are split across threads."""

    @pytest.mark.parametrize("opts", [
        dict(all_edges=True, max_pairs_per_block=1000, threads=8),
        dict(cor=0.4, pval=0.01, max_pairs_per_block=1, threads=3),
        dict(cor=0.4, pval=0.01, max_pairs_per_block=50_000, threads=5, method="pearson"),
    ])
    def test_many_blocks_and_threads(self, tmp_path, opts):
        src = _write_matrix(tmp_path / "m.tsv", n_genes=700, n_samples=6, seed=3)
        assert _same(tmp_path, src, **opts)["edge_count"] > 0

    def test_default_block_size_per_path(self, tmp_path, monkeypatch):
        import gcen2._build_file as bf
        import gcen2.network as net

        corr = _core().correlation
        seen = {}
        real_cpp, real_np = corr.build_network_file, net.build_to_file

        def cpp(*a, **kw):
            seen["cpp"] = a[9]
            return real_cpp(*a, **kw)

        def numpy_path(*a, **kw):
            seen["numpy"] = kw["max_pairs_per_block"]
            return real_np(*a, **kw)

        monkeypatch.setattr(corr, "build_network_file", cpp)
        monkeypatch.setattr(net, "build_to_file", numpy_path)
        src = _write_matrix(tmp_path / "m.tsv")
        build_from_file(src, tmp_path / "a.net", fdr=True)  # FDR on the C++ path too
        build_from_file(src, tmp_path / "b.net", log="log2")
        assert seen == {"cpp": bf.MAX_PAIRS, "numpy": bf.MAX_PAIRS}
        build_from_file(src, tmp_path / "c.net", max_pairs_per_block=123)
        assert seen["cpp"] == 123

    def test_zero_max_pairs_rejected(self, tmp_path):
        src = _write_matrix(tmp_path / "m.tsv")
        with pytest.raises(ValueError):
            build_from_file(src, tmp_path / "a.net", max_pairs_per_block=0)


class TestFallback:
    @pytest.mark.parametrize("value", ["NA", "", "nan", "inf", "+1.5", " 1.5", '"1.5"',
                                       "1,5", "0x10", "1e400", "1_000"])
    def test_non_plain_values_use_numpy_path(self, tmp_path, value):
        src = _write_matrix(tmp_path / "m.tsv", n_genes=20)
        lines = src.read_text(encoding="utf-8").splitlines()
        f = lines[4].split("\t")
        f[2] = value
        lines[4] = "\t".join(f)
        src.write_text("\n".join(lines) + "\n", encoding="utf-8")
        try:
            _reference(src, tmp_path / "ref.net", all_edges=True)
        except Exception as exc:  # the numpy path rejects it: so must build_from_file
            with pytest.raises(type(exc)):
                build_from_file(src, tmp_path / "x.net", all_edges=True)
            return
        _same(tmp_path, src, all_edges=True, expect_cpp=False)

    def test_invalid_utf8_falls_back(self, tmp_path):
        src = _write_matrix(tmp_path / "m.tsv", n_genes=10)
        src.write_bytes(src.read_bytes().replace(b"G3\t", b"G\xff3\t"))
        with pytest.raises(UnicodeDecodeError):
            build_from_file(src, tmp_path / "x.net")
        assert _core().correlation.build_network_file(
            str(src), str(tmp_path / "c.net"), True, "spearman", 0.1, 0.001, False, False,
            1, 1000) is None
        assert not (tmp_path / "c.net").exists()

    def test_unicode_whitespace_line_falls_back(self, tmp_path):
        src = _write_matrix(tmp_path / "m.tsv", n_genes=10)
        src.write_text(src.read_text(encoding="utf-8") + " \n", encoding="utf-8")  # blank for Python's str.strip()
        _same(tmp_path, src, all_edges=True, expect_cpp=False)

    def test_errors_match(self, tmp_path):
        with pytest.raises(GCENIOError):
            build_from_file(tmp_path / "missing.tsv", tmp_path / "x.net")
        two = tmp_path / "two.tsv"
        two.write_text("g1\t1\t2\ng2\t2\t1\n", encoding="utf-8")
        with pytest.raises(GCENInputError, match="3 samples"):
            build_from_file(two, tmp_path / "x.net")
        bad = _write_matrix(tmp_path / "bad.tsv", n_genes=5)
        lines = bad.read_text(encoding="utf-8").splitlines()
        lines[0] = "gene\ts0"  # header shorter than the rows
        bad.write_text("\n".join(lines) + "\n", encoding="utf-8")
        with pytest.raises(GCENIOError, match="header"):
            build_from_file(bad, tmp_path / "x.net", header=True)
        with pytest.raises(GCENInputError, match="Unknown method"):
            build_from_file(bad, tmp_path / "x.net", method="kendall")
        good = _write_matrix(tmp_path / "good.tsv", n_genes=5)
        with pytest.raises(OSError):
            build_from_file(good, tmp_path / "no_dir" / "x.net")

    @pytest.mark.parametrize("opts", [dict(fdr=True, pval=0.5, cor=0.2), dict(log="log2")])
    def test_fdr_and_log_use_numpy_path(self, tmp_path, opts):
        src = _write_matrix(tmp_path / "m.tsv")
        src.write_text(src.read_text(encoding="utf-8").replace("-", ""), encoding="utf-8")  # log needs x + 1 > 0
        _same(tmp_path, src, expect_cpp=None, **opts)


def test_cli_build_does_not_import_numpy(tmp_path):
    src = _write_matrix(tmp_path / "m.tsv")
    out = tmp_path / "o.net"
    code = (
        "import sys\n"
        "from typer.testing import CliRunner\n"
        "from gcen2.cli import app\n"
        f"r = CliRunner().invoke(app, ['build', '-i', {str(src)!r}, '-o', {str(out)!r},"
        " '-c', '0.5', '-p', '0.05'])\n"
        "assert r.exit_code == 0, r.output\n"
        "print(sorted(m for m in ('numpy', 'pandas', 'scipy') if m in sys.modules))\n")
    res = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8",
                         check=True)
    assert res.stdout.strip().splitlines()[-1] == "[]"
    ref = tmp_path / "ref.net"
    _reference(src, ref, cor=0.5, pval=0.05)
    assert out.read_bytes() == ref.read_bytes()
