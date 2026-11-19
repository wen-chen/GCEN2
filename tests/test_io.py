"""Tests for gcen2.io — file I/O utilities."""

import numpy as np
import pandas as pd
import pytest

from gcen2.exceptions import GCENIOError
from gcen2.io import read_csv, read_tsv, write_tsv
from gcen2.models import ExpressionMatrix


@pytest.fixture
def sample_df():
    return pd.DataFrame(
        {"S1": [1.0, 2.0, 3.0], "S2": [4.0, 5.0, 6.0]},
        index=["GeneA", "GeneB", "GeneC"],
    )


class TestReadTSV:
    def test_reads_correctly(self, tmp_path, sample_df):
        p = tmp_path / "test.tsv"
        sample_df.to_csv(p, sep="\t", lineterminator="\n")
        result = read_tsv(p)
        assert isinstance(result, ExpressionMatrix)
        pd.testing.assert_frame_equal(result.data, sample_df)

    def test_gene_and_sample_names(self, tmp_path, sample_df):
        p = tmp_path / "test.tsv"
        sample_df.to_csv(p, sep="\t", lineterminator="\n")
        result = read_tsv(p)
        assert list(result.gene_names) == list(sample_df.index)
        assert list(result.sample_names) == list(sample_df.columns)

    def test_missing_file_raises(self):
        with pytest.raises(GCENIOError):
            read_tsv("/nonexistent/file.tsv")


class TestReadCSV:
    def test_handles_quoted_fields(self, tmp_path):
        """Ensure proper CSV parsing (fixes GCEN v1 naive replacement bug)."""
        content = '"gene","S1","S2"\n"Gene,A",1.0,2.0\n"GeneB",3.0,4.0\n'
        p = tmp_path / "test.csv"
        p.write_text(content, encoding="utf-8")
        result = read_csv(p)
        assert "Gene,A" in result.gene_names  # quoted comma preserved

    def test_missing_file_raises(self):
        with pytest.raises(GCENIOError):
            read_csv("/nonexistent/file.csv")


class TestWriteTSV:
    def test_round_trip(self, tmp_path, sample_df):
        p = tmp_path / "out.tsv"
        write_tsv(sample_df, p)
        result = read_tsv(p)
        pd.testing.assert_frame_equal(result.data, sample_df)


class TestHeaderDetection:
    """Regression: header-less GCEN v1 matrices must not lose their first gene."""

    V1_ROWS = "G1\t0.5\t1.25\t3\nG2\t2\t0\t1.5\nG3\t4\t5\t6\n"

    def test_headerless_v1_tsv_keeps_first_gene(self, tmp_path):
        path = tmp_path / "v1.tsv"
        path.write_text(self.V1_ROWS, encoding="utf-8")
        expr = read_tsv(path)
        assert list(expr.data.index) == ["G1", "G2", "G3"]
        assert expr.data.shape == (3, 3)
        assert expr.data.loc["G1"].tolist() == [0.5, 1.25, 3.0]

    def test_headerless_v1_csv_keeps_first_gene(self, tmp_path):
        path = tmp_path / "v1.csv"
        path.write_text(self.V1_ROWS.replace("\t", ","), encoding="utf-8")
        expr = read_csv(path)
        assert list(expr.data.index) == ["G1", "G2", "G3"]

    def test_hash_header_is_header(self, tmp_path):
        path = tmp_path / "hash.tsv"
        path.write_text("#gene\tS1\tS2\tS3\n" + self.V1_ROWS, encoding="utf-8")
        expr = read_tsv(path)
        assert list(expr.data.columns) == ["S1", "S2", "S3"]
        assert list(expr.data.index) == ["G1", "G2", "G3"]

    def test_pandas_header_still_read(self, tmp_path, sample_df):
        path = tmp_path / "pd.tsv"
        sample_df.to_csv(path, sep="\t", lineterminator="\n")
        expr = read_tsv(path)
        assert list(expr.data.columns) == ["S1", "S2"]

    def test_numeric_sample_names_with_explicit_header(self, tmp_path):
        path = tmp_path / "num.tsv"
        path.write_text("gene\t1\t2\t3\n" + self.V1_ROWS, encoding="utf-8")
        # ambiguous: integer-only first line is read as data, with a warning
        with pytest.warns(UserWarning, match="header=True"):
            assert "gene" in read_tsv(path).data.index
        assert list(read_tsv(path, header=True).data.columns) == ["1", "2", "3"]

    def test_decimal_first_row_no_warning(self, tmp_path):
        import warnings
        path = tmp_path / "v1.tsv"
        path.write_text(self.V1_ROWS, encoding="utf-8")
        with warnings.catch_warnings():
            warnings.simplefilter("error")
            read_tsv(path)


class TestReadGeneList:
    def test_v1_seed_file(self, tmp_path):
        from gcen2.io import read_gene_list
        path = tmp_path / "seeds.list"
        path.write_text("#gene\tprobability\nG1\t1.0\n\nG2\t0.5\n  G3  \n", encoding="utf-8")
        assert read_gene_list(path) == ["G1", "G2", "G3"]

    def test_missing_file(self):
        from gcen2.io import read_gene_list
        with pytest.raises(GCENIOError):
            read_gene_list("/nonexistent/genes.txt")


class TestReadNetworkGenesOnly:
    @pytest.mark.parametrize("header", [
        "", "#node1\tnode2\tcorrelation\tp-value\n", "gene_a\tgene_b\tcorrelation\tp_value\n"])
    def test_same_genes_as_full_read(self, tmp_path, header):
        from gcen2.io import read_network
        p = tmp_path / "n.network"
        p.write_text(header + "A\tB\t0.9\t1e-5\nB\tC\t-0.95\t2e-6\nA\tC\t0.91\t3e-5\n", encoding="utf-8")
        full = read_network(p)
        genes = read_network(p, genes_only=True)
        assert list(genes.edges.columns) == ["gene_a", "gene_b"]
        assert str(genes.edges["gene_a"].dtype) == "category"
        assert genes.edges["gene_a"].astype(str).tolist() == full.edges["gene_a"].tolist()
        assert genes.edges["gene_b"].astype(str).tolist() == full.edges["gene_b"].tolist()
        assert len(full.edges) == 3 and genes.node_count == 3
        assert list(full.edges.columns)[:4] == ["gene_a", "gene_b", "correlation", "p_value"]


class TestReadMatrix:
    """gcen2._matrix_io.read_matrix (pandas-free, used by ``gcen2 build``)
    must read the same matrix as read_tsv."""

    @pytest.mark.parametrize("text", [
        "\tS1\tS2\tS3\nG1\t1.5\t2\t3e-2\nG2\t0.1\t0.2\t0.3\n",       # pandas header
        "S1\tS2\tS3\nG1\t1.5\t2\t3\nG2\t4\t5\t6\n",                # header without gene column
        "#gene\tS1\tS2\tS3\nG1\t1\t2\t3\nG2\t4\t5\t6\n",          # v1 # header
        "G1\t1.5\t2.25\t3.1\nG2\t4.5\t5.5\t6.5\n",                    # header-less v1
        "gene\tS1\tS2\tS3\r\nG1\t1\tNA\t3\r\n\r\nG2\t\tnan\t6\r\n",  # CRLF, NA, blank line
        'gene\tS1\tS2\tS3\n"G1"\t"1.5"\t2\t3\nG2\t4\t5\n',          # quotes, short row
    ])
    def test_same_as_read_tsv(self, tmp_path, text):
        from gcen2._matrix_io import read_matrix
        p = tmp_path / "m.tsv"
        p.write_bytes(text.encode())
        ref = read_tsv(p)
        values, genes, samples = read_matrix(p)
        want = np.ascontiguousarray(ref.data.to_numpy(dtype=np.float64))
        np.testing.assert_array_equal(values, want)
        assert values.flags.c_contiguous and values.dtype == np.float64
        assert genes == [str(g) for g in ref.data.index]
        assert samples == [str(s) for s in ref.data.columns]

    def test_bad_value_raises(self, tmp_path):
        from gcen2._matrix_io import read_matrix
        p = tmp_path / "m.tsv"
        p.write_text("gene\tS1\tS2\nG1\t1\tabc\n", encoding="utf-8")
        with pytest.raises(GCENIOError, match="G1"):
            read_matrix(p)

    def test_missing_file(self, tmp_path):
        from gcen2._matrix_io import read_matrix
        with pytest.raises(GCENIOError):
            read_matrix(tmp_path / "none.tsv")


def test_build_path_does_not_import_pandas():
    """``import gcen2`` is lazy and the build path needs no pandas/scipy."""
    import subprocess
    import sys
    code = ("import sys, gcen2; from gcen2._matrix_io import read_matrix; "
            "from gcen2.network import build_to_file; import gcen2.cli; "
            "print(sorted(m for m in ('pandas', 'scipy') if m in sys.modules))")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8", check=True)
    assert out.stdout.strip() == "[]"
