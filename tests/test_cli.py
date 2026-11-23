"""End-to-end CLI tests for gcen2.cli.

Requires typer. Skipped when not installed.
Tests use typer's CliRunner for in-process invocation.
"""

import numpy as np
import pandas as pd
import pytest

try:
    from typer.testing import CliRunner

    from gcen2.cli import app
    HAS_TYPER = True
except ImportError:
    HAS_TYPER = False

requires_typer = pytest.mark.skipif(not HAS_TYPER, reason="typer not installed")


@pytest.fixture
def runner():
    return CliRunner()


@pytest.fixture
def sample_tsv(tmp_path):
    """Write a small expression TSV and return its path."""
    rng = np.random.default_rng(0)
    data = rng.exponential(scale=10.0, size=(15, 8))
    genes = [f"Gene{i}" for i in range(15)]
    samples = [f"S{j}" for j in range(8)]
    df = pd.DataFrame(data, index=genes, columns=samples)
    path = tmp_path / "expr.tsv"
    df.to_csv(path, sep="\t", lineterminator="\n")
    return path


@requires_typer
@pytest.mark.parametrize("flag", ["--version", "-V"])
def test_version_option(runner, flag):
    import gcen2

    result = runner.invoke(app, [flag])
    assert result.exit_code == 0
    assert result.output.strip() == f"gcen2 {gcen2.__version__}"


@requires_typer
class TestNormalizeCommand:
    def test_tmm_produces_output(self, runner, sample_tsv, tmp_path):
        out = tmp_path / "norm.tsv"
        result = runner.invoke(app, ["normalize", "-i", str(sample_tsv),
                                     "-o", str(out), "-m", "tmm"])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_unknown_method_fails(self, runner, sample_tsv, tmp_path):
        out = tmp_path / "norm.tsv"
        result = runner.invoke(app, ["normalize", "-i", str(sample_tsv),
                                     "-o", str(out), "-m", "invalid_method"])
        assert result.exit_code != 0 or "error" in result.output.lower()

    def test_missing_input_flag_fails(self, runner, tmp_path):
        out = tmp_path / "norm.tsv"
        result = runner.invoke(app, ["normalize", "-o", str(out)])
        assert result.exit_code != 0


@requires_typer
class TestFilterCommand:
    def test_filter_produces_output(self, runner, sample_tsv, tmp_path):
        out = tmp_path / "filtered.tsv"
        result = runner.invoke(app, ["filter", "-i", str(sample_tsv),
                                     "-o", str(out), "-p", "0.75"])
        assert result.exit_code == 0, result.output
        assert out.exists()
        # Filtered file should have fewer rows than original
        orig = pd.read_csv(sample_tsv, sep="\t", index_col=0)
        filt = pd.read_csv(out, sep="\t", index_col=0)
        assert len(filt) <= len(orig)


@requires_typer
class TestStatsCommand:
    def test_stats_output(self, runner, tmp_path):
        edges = pd.DataFrame({
            "gene_a": ["A", "A", "B"],
            "gene_b": ["B", "C", "C"],
            "correlation": [0.9, 0.8, 0.7],
            "p_value": [0.001, 0.002, 0.003],
        })
        net_path = tmp_path / "net.tsv"
        edges.to_csv(net_path, sep="\t", index=False, lineterminator="\n")
        result = runner.invoke(app, ["stats", "-i", str(net_path)])
        assert result.exit_code == 0, result.output
        assert "Nodes" in result.output
        assert "Edges" in result.output

    def test_stats_without_pandas(self, tmp_path):
        import subprocess
        import sys

        net_path = tmp_path / "net.tsv"
        net_path.write_text(V1_NETWORK, encoding="utf-8")
        code = (
            "import sys\n"
            "from typer.testing import CliRunner\n"
            "from gcen2.cli import app\n"
            f"r = CliRunner().invoke(app, ['stats', '-i', {str(net_path)!r}, '-t', '3'])\n"
            "assert r.exit_code == 0, r.output\n"
            "print(r.output)\n"
            "print(sorted(m for m in ('pandas', 'scipy', 'networkx') if m in sys.modules))\n")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8",
                             check=True).stdout
        assert out.strip().splitlines()[-1] == "[]"
        assert "Nodes:                  4" in out
        assert "Transitivity:           0.6000" in out

    def test_merge_without_pandas(self, tmp_path):
        import subprocess
        import sys

        net_path = tmp_path / "net.tsv"
        net_path.write_text(V1_NETWORK, encoding="utf-8")
        out_path = tmp_path / "m.tsv"
        code = (
            "import sys\n"
            "from typer.testing import CliRunner\n"
            "from gcen2.cli import app\n"
            f"r = CliRunner().invoke(app, ['merge', '-i', {str(net_path)!r}, '-i', "
            f"{str(net_path)!r}, '-o', {str(out_path)!r}, '-c', '1.65'])\n"
            "assert r.exit_code == 0, r.output\n"
            "print(r.output)\n"
            "print(sorted(m for m in ('pandas', 'scipy', 'networkx') if m in sys.modules))\n")
        out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, encoding="utf-8",
                             check=True).stdout
        assert out.strip().splitlines()[-1] == "[]"
        assert "4 nodes, 2 edges" in out
        assert out_path.read_text(encoding="utf-8") == ("gene_a\tgene_b\tcorrelation\tdirection\n"
                                        "A\tB\t1.8\t1\nC\tD\t1.7\t1\n")


V1_NETWORK = (
    "#Node1\tNode2\tCorrelation\tPvalue\tFDR\n"
    "A\tB\t0.9\t0.001\n"
    "A\tC\t0.8\t0.002\n"
    "B\tC\t0.7\t0.003\n"
    "C\tD\t0.85\t0.001\n"
)


@requires_typer
class TestV1NetworkFormat:
    """Regression: CLI commands must accept GCEN v1 ``.network`` files.

    Several commands used to read networks with a bare ``pd.read_csv`` and
    failed with ``KeyError: 'gene_a'`` on v1 headers (``#Node1`` / ``Node2``).
    """

    @pytest.fixture
    def v1_net(self, tmp_path):
        path = tmp_path / "v1.network"
        path.write_text(V1_NETWORK, encoding="utf-8")
        return path

    def test_stats(self, runner, v1_net):
        result = runner.invoke(app, ["stats", "-i", str(v1_net)])
        assert result.exit_code == 0, result.output
        assert "Nodes" in result.output

    def test_module(self, runner, v1_net, tmp_path):
        out = tmp_path / "modules.tsv"
        result = runner.invoke(app, ["module", "-i", str(v1_net), "-o", str(out)])
        assert result.exit_code == 0, result.output

    def test_rwr(self, runner, v1_net, tmp_path):
        seeds = tmp_path / "seeds.txt"
        seeds.write_text("A\n", encoding="utf-8")
        out = tmp_path / "rwr.tsv"
        result = runner.invoke(
            app, ["rwr", "-n", str(v1_net), "-g", str(seeds), "-o", str(out)])
        assert result.exit_code == 0, result.output
        assert out.exists()


@requires_typer
class TestRWRCommand:
    def test_rwr_produces_output(self, runner, tmp_path):
        edges = pd.DataFrame({
            "gene_a": ["A", "A", "B", "C"],
            "gene_b": ["B", "C", "C", "D"],
            "correlation": [0.9, 0.8, 0.7, 0.85],
            "p_value": [0.001] * 4,
        })
        net_path = tmp_path / "net.tsv"
        edges.to_csv(net_path, sep="\t", index=False, lineterminator="\n")

        seeds_path = tmp_path / "seeds.txt"
        seeds_path.write_text("A\n", encoding="utf-8")

        out_path = tmp_path / "rwr_out.tsv"
        result = runner.invoke(app, [
            "rwr", "-n", str(net_path), "-g", str(seeds_path),
            "-o", str(out_path)
        ])
        assert result.exit_code == 0, result.output
        assert out_path.exists()


@requires_typer
class TestReportCommand:
    def test_report_creates_html(self, runner, tmp_path):
        try:
            import jinja2  # noqa: F401
        except ImportError:
            pytest.skip("Jinja2 not installed")

        edges = pd.DataFrame({
            "gene_a": ["A", "B"],
            "gene_b": ["B", "C"],
            "correlation": [0.9, 0.8],
            "p_value": [0.001, 0.002],
        })
        net_path = tmp_path / "net.tsv"
        edges.to_csv(net_path, sep="\t", index=False, lineterminator="\n")
        out_path = tmp_path / "report.html"

        result = runner.invoke(app, [
            "report", "-n", str(net_path), "-o", str(out_path)
        ])
        assert result.exit_code == 0, result.output


@requires_typer
class TestMergeCommand:
    def test_merge_two_networks(self, runner, tmp_path):
        net1 = tmp_path / "net1.network"
        net2 = tmp_path / "net2.network"
        out = tmp_path / "merged.network"
        net1.write_text("gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t0.9\t0.001\nA\tC\t0.85\t0.001\n", encoding="utf-8")
        net2.write_text("gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t0.8\t0.001\nB\tC\t0.88\t0.001\n", encoding="utf-8")
        result = runner.invoke(app, [
            "merge",
            "-i", str(net1), "-i", str(net2),
            "-o", str(out),
            "-c", "0.0",
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        df = pd.read_csv(out, sep="\t")
        assert "gene_a" in df.columns
        assert "gene_b" in df.columns
        assert len(df) > 0


@requires_typer
class TestExtractCommand:
    def test_extract_subnetwork(self, runner, tmp_path):
        net_file = tmp_path / "net.network"
        genes_file = tmp_path / "genes.txt"
        out = tmp_path / "sub.network"
        net_file.write_text(
            "gene_a\tgene_b\tcorrelation\tp_value\n"
            "A\tB\t0.9\t0.001\nA\tC\t0.85\t0.001\nD\tE\t0.88\t0.001\n", encoding="utf-8"
        )
        genes_file.write_text("A\nB\n", encoding="utf-8")
        result = runner.invoke(app, [
            "extract", "-i", str(net_file), "-g", str(genes_file), "-o", str(out)
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        df = pd.read_csv(out, sep="\t")
        # D-E edge 不在 A/B 子网络里
        gene_pairs = set(zip(df["gene_a"], df["gene_b"]))
        assert ("D", "E") not in gene_pairs
        assert len(df) > 0


@requires_typer
class TestDiffCommand:
    def test_diff_two_networks(self, runner, tmp_path):
        net_a = tmp_path / "netA.network"
        net_b = tmp_path / "netB.network"
        out = tmp_path / "diff.tsv"
        net_a.write_text(
            "gene_a\tgene_b\tcorrelation\tp_value\n"
            "A\tB\t0.9\t0.001\nA\tC\t0.85\t0.001\n", encoding="utf-8"
        )
        net_b.write_text(
            "gene_a\tgene_b\tcorrelation\tp_value\n"
            "A\tB\t0.5\t0.001\nB\tC\t0.88\t0.001\n", encoding="utf-8"
        )
        result = runner.invoke(app, [
            "diff",
            "--network-a", str(net_a),
            "--network-b", str(net_b),
            "-o", str(out),
            "--n-a", "30", "--n-b", "30",
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        df = pd.read_csv(out, sep="\t")
        assert "category" in df.columns
        valid = {"A_only", "B_only", "shared_changed", "shared_stable"}
        assert set(df["category"]).issubset(valid)
        assert len(df) > 0


@requires_typer
class TestShuffleCommand:
    def test_shuffle_preserves_edge_count(self, runner, tmp_path):
        import itertools
        net_file = tmp_path / "net.network"
        out = tmp_path / "shuffled.network"
        lines = ["gene_a\tgene_b\tcorrelation\tp_value"]
        genes = [f"G{i}" for i in range(10)]
        for a, b in list(itertools.combinations(genes, 2))[:10]:
            lines.append(f"{a}\t{b}\t0.9\t0.001")
        net_file.write_text("\n".join(lines) + "\n", encoding="utf-8")
        result = runner.invoke(app, [
            "shuffle", "-i", str(net_file), "-o", str(out), "--swap-factor", "10"
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        orig = pd.read_csv(net_file, sep="\t")
        shuf = pd.read_csv(out, sep="\t")
        assert len(shuf) == len(orig)


@requires_typer
class TestSpatialBuildCommand:
    def test_filtered_mode(self, runner, tmp_path):
        """spatial-build --mode filtered should produce a .network file."""
        pytest.importorskip("anndata")
        import anndata
        import numpy as np
        import pandas as pd

        # Create minimal h5ad
        rng = np.random.default_rng(0)
        X = rng.random((20, 30)).astype(np.float32)
        coords = np.array([[i, j] for i in range(4) for j in range(5)], dtype=float)
        obs = pd.DataFrame(index=[f"s{i}" for i in range(20)])
        var = pd.DataFrame(index=[f"G{i}" for i in range(30)])
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        adata.obsm["spatial"] = coords
        h5ad = tmp_path / "data.h5ad"
        adata.write_h5ad(h5ad)
        out = tmp_path / "out.network"

        result = runner.invoke(app, [
            "spatial-build", "-i", str(h5ad), "-o", str(out),
            "--mode", "filtered", "--n-neighbors", "4",
            "--cor-cutoff", "0.0", "--p-cutoff", "1.0",
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()

    def test_regions_mode_creates_dir(self, runner, tmp_path):
        pytest.importorskip("anndata")
        import anndata
        import numpy as np
        import pandas as pd

        rng = np.random.default_rng(1)
        X = rng.random((20, 30)).astype(np.float32)
        coords = np.array([[i, j] for i in range(4) for j in range(5)], dtype=float)
        obs = pd.DataFrame({"region": ["A"] * 10 + ["B"] * 10},
                           index=[f"s{i}" for i in range(20)])
        var = pd.DataFrame(index=[f"G{i}" for i in range(30)])
        adata = anndata.AnnData(X=X, obs=obs, var=var)
        adata.obsm["spatial"] = coords
        h5ad = tmp_path / "data.h5ad"
        adata.write_h5ad(h5ad)
        out_dir = tmp_path / "regions"

        result = runner.invoke(app, [
            "spatial-build", "-i", str(h5ad), "-o", str(out_dir),
            "--mode", "regions", "--region-key", "region",
            "--n-neighbors", "4", "--cor-cutoff", "0.0", "--p-cutoff", "1.0",
        ])
        assert result.exit_code == 0, result.output
        assert out_dir.is_dir()
        network_files = list(out_dir.glob("*.network"))
        assert len(network_files) == 2


class TestAnnotateV1Layout:
    """`gcen2 annotate` accepts v1 options and writes v1-named files."""

    FIX = __import__("pathlib").Path(__file__).parent / "data" / "v1_annotate"

    @pytest.mark.parametrize("mode", ["network", "module"])
    def test_files_match_v1_names(self, runner, tmp_path, mode):
        src = ["-n", str(self.FIX / "net.network")] if mode == "network" else \
              ["-m", str(self.FIX / "modules.txt")]
        out = tmp_path / "ann"
        res = runner.invoke(app, ["annotate", "-g", str(self.FIX / "go.obo"),
                                  "-a", str(self.FIX / "gene_go.assoc"), *src,
                                  "-o", str(out), "-p", "1.0"])
        assert res.exit_code == 0, res.output
        got = {p.name for p in out.glob("*.go")}
        expected = {p.name for p in (self.FIX / f"v1_{mode}").glob("*.go")}
        assert got == expected
        if mode == "module":
            assert {p.name for p in out.glob("*.txt")} == {n[:-3] + ".txt" for n in expected}
        df = pd.read_csv(next(iter(out.glob("*.go"))), sep="\t")
        assert df["name"].ne("").all() and "namespace" in df

    def test_requires_source(self, runner, tmp_path):
        res = runner.invoke(app, ["annotate", "-g", str(self.FIX / "go.obo"),
                                  "-a", str(self.FIX / "gene_go.assoc"), "-o", str(tmp_path)])
        assert res.exit_code != 0


class TestReadModules:
    def test_v1_and_v2_formats(self, tmp_path):
        from gcen2.io import read_modules
        v1 = tmp_path / "v1.txt"
        v1.write_text("A\tB\tC\nD\tE\n", encoding="utf-8")
        v2 = tmp_path / "v2.tsv"
        v2.write_text("module\tgene\n0\tA\n0\tB\n0\tC\n1\tD\n1\tE\n", encoding="utf-8")
        assert read_modules(v1).modules == [{"A", "B", "C"}, {"D", "E"}]
        assert read_modules(v2).modules == [{"A", "B", "C"}, {"D", "E"}]


@requires_typer
@pytest.mark.parametrize("args", [
    ["build", "-i", "{m}/expr.tsv", "-o", "{o}/x.network"],
    ["module", "-i", "{m}/x.network", "-o", "{o}/m.tsv"],
    ["merge", "-i", "{m}/a.network", "-i", "{m}/b.network", "-o", "{o}/m.network"],
    ["consensus-modules", "-e", "{m}/a.tsv", "-e", "{m}/b.tsv", "-o", "{o}/c"],
    ["spatial-build", "-i", "{m}/s.h5ad", "-o", "{o}/s.network"],
    ["trajectory", "-i", "{m}/expr.tsv", "--pseudotime", "{m}/pt.tsv", "-o", "{o}/t"],
    ["fuse", "-i", "{m}/a.network", "-i", "{m}/b.network", "-o", "{o}/f.tsv"],
    ["report", "-m", "{m}/modules.tsv", "-o", "{o}/r.html"],
    ["plot", "kme", "-i", "{m}/kme.tsv", "-o", "{o}/kme"],
    ["plot", "degree-distribution", "-i", "{m}/x.network", "-o", "{o}/deg"],
])
def test_missing_input_is_an_error_message(runner, tmp_path, args):
    """Every command reports a missing input file as ``ERROR: ...`` with
    exit code 1, not as a traceback (gcen2.cli._common.wrap_command)."""
    pytest.importorskip("anndata")
    args = [a.format(m=tmp_path / "missing", o=tmp_path / "out") for a in args]
    res = runner.invoke(app, args)
    assert res.exit_code == 1, res.output
    assert "ERROR:" in res.output and "not found" in res.output.lower(), res.output
    assert res.exception is None or isinstance(res.exception, SystemExit)


@requires_typer
def test_warnings_are_printed_once_as_warning(runner, tmp_path, monkeypatch):
    """User warnings of a command are printed as ``WARNING: ...``."""
    import warnings

    import gcen2.network

    def stats_with_warning(path, threads=2):
        warnings.warn("something to know", UserWarning, stacklevel=2)
        raise gcen2.GCENInputError("and then it failed")

    monkeypatch.setattr(gcen2.network, "stats_from_file", stats_with_warning)
    res = runner.invoke(app, ["stats", "-i", str(tmp_path / "x.network")])
    assert res.exit_code == 1
    assert res.output.count("WARNING: something to know") == 1
    assert "ERROR: and then it failed" in res.output
    assert res.output.index("WARNING") < res.output.index("ERROR")


@requires_typer
def test_python_m_gcen2_cli():
    """``python -m gcen2.cli`` runs the same CLI as the gcen2 script."""
    import subprocess
    import sys

    r = subprocess.run([sys.executable, "-m", "gcen2.cli", "--help"], capture_output=True,
                       text=True, encoding="utf-8")
    assert r.returncode == 0, r.stderr
    assert "celltype-modules" in r.stdout and "plot" in r.stdout
