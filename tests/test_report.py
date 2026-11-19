"""Tests for gcen2.report — HTML report generation.

Requires Jinja2. Skipped when not installed.
"""

import pytest

try:
    import jinja2  # noqa: F401
    HAS_JINJA = True
except ImportError:
    HAS_JINJA = False

requires_jinja = pytest.mark.skipif(not HAS_JINJA, reason="Jinja2 not installed")


@requires_jinja
class TestGenerateReport:
    def test_generates_html_file(self, tmp_path, small_network):
        from gcen2.report import generate
        out = tmp_path / "report.html"
        generate(network=small_network, output=out)
        assert out.exists()
        content = out.read_text(encoding="utf-8")
        assert "<html" in content.lower()

    def test_html_contains_version(self, tmp_path):
        import gcen2
        from gcen2.report import generate
        out = tmp_path / "report.html"
        generate(output=out)
        content = out.read_text(encoding="utf-8")
        assert gcen2.__version__ in content

    def test_html_has_network_section(self, tmp_path, small_network):
        from gcen2.report import generate
        out = tmp_path / "report.html"
        generate(network=small_network, output=out)
        content = out.read_text(encoding="utf-8")
        assert "Network Overview" in content

    def test_html_has_module_section(self, tmp_path, small_network):
        from gcen2.models import ModuleSet
        from gcen2.report import generate
        modules = ModuleSet(
            modules=[{"GeneA", "GeneB", "GeneC"}, {"GeneD", "GeneE"}]
        )
        out = tmp_path / "report.html"
        generate(network=small_network, modules=modules, output=out)
        content = out.read_text(encoding="utf-8")
        assert "Module Results" in content

    def test_html_has_rwr_section(self, tmp_path, small_network):
        import pandas as pd

        from gcen2.models import RWRResult
        from gcen2.report import generate
        rwr = RWRResult(
            ranked_genes=pd.DataFrame({"gene": ["GeneD", "GeneE"],
                                       "probability": [0.05, 0.03]}),
            seed_genes=["GeneA"],
        )
        out = tmp_path / "report.html"
        generate(network=small_network, rwr=rwr, output=out)
        content = out.read_text(encoding="utf-8")
        assert "RWR" in content

    def test_pdf_requires_weasyprint(self, tmp_path, small_network, monkeypatch):
        """Without weasyprint, or without its system libraries: a clear error."""
        import builtins
        import sys

        from gcen2.exceptions import GCENInputError
        from gcen2.report import generate

        out = tmp_path / "report.pdf"
        monkeypatch.setitem(sys.modules, "weasyprint", None)  # import -> ImportError
        with pytest.raises(GCENInputError, match="pip install weasyprint"):
            generate(network=small_network, output=out)
        monkeypatch.delitem(sys.modules, "weasyprint")
        real_import = builtins.__import__

        def no_pango(name, *args, **kwargs):
            if name == "weasyprint":
                raise OSError("cannot load library 'libpango-1.0-0'")
            return real_import(name, *args, **kwargs)

        monkeypatch.setattr(builtins, "__import__", no_pango)
        with pytest.raises(GCENInputError, match="system libraries"):
            generate(network=small_network, output=out)
        assert not out.exists()

    def test_missing_template_raises(self, tmp_path, small_network):
        from gcen2.exceptions import GCENIOError
        from gcen2.report import generate

        with pytest.raises(GCENIOError, match="template not found"):
            generate(network=small_network,
                     template="/nonexistent/template.html.jinja2",
                     output=tmp_path / "out.html")


@requires_jinja
class TestReportFromFile:
    """network.summary_from_file (C++ adjacency, no edge table) gives the same
    report as read_network."""

    _TEXT = (
        "#Node1\tNode2\tCorrelation\tPvalue\n"
        "A\tB\t0.9\t1e-3\n"
        "B\tA\t0.8\t1e-3\n"          # reversed duplicate
        "A\tB\t0.7\t1e-3\n"          # duplicate
        "\"C\"\tA\t0.6\t1e-3\n"      # quoted name
        "C\tC\t0.5\t1e-3\n"          # self-loop: degree + 2
        "D\tE\t0.5\t1e-3\n"
        "101\tA\t-0.6\t1e-3\n"
        "B\tC\t0.6\t1e-3\n"
    )

    def _file(self, tmp_path, text=None):
        p = tmp_path / "net.network"
        p.write_text(self._TEXT if text is None else text, encoding="utf-8")
        return p

    def _require_core(self):
        from gcen2.core import _gcen_core
        if _gcen_core is None:
            pytest.skip("compiled extension not available")

    def test_stats_and_degrees_match_networkx(self, tmp_path):
        self._require_core()
        pytest.importorskip("networkx")
        from gcen2.io import read_network
        from gcen2.network import NetworkSummary, summary_from_file

        path = self._file(tmp_path)
        net = read_network(path)
        s = summary_from_file(path)
        assert isinstance(s, NetworkSummary)
        assert s.stats == net.stats()
        assert sorted(s.degrees.tolist()) == sorted(d for _, d in net.to_networkx().degree())
        assert (s.method, s.is_signed) == (net.method, net.is_signed)

    @pytest.mark.parametrize("text", [None, "gene_a\tgene_b\tcorrelation\tp_value\n"])
    def test_same_report_bytes(self, tmp_path, text):
        self._require_core()
        from gcen2.io import read_network
        from gcen2.network import summary_from_file
        from gcen2.report import generate

        path = self._file(tmp_path, text)
        generate(network=read_network(path), output=tmp_path / "ref.html")
        generate(network=summary_from_file(path), output=tmp_path / "fast.html")
        assert (tmp_path / "fast.html").read_bytes() == (tmp_path / "ref.html").read_bytes()

    def test_fallback_without_extension(self, tmp_path, monkeypatch):
        import gcen2.core
        from gcen2.models import Network
        from gcen2.network import summary_from_file

        monkeypatch.setattr(gcen2.core, "_gcen_core", None)
        assert isinstance(summary_from_file(self._file(tmp_path)), Network)


@requires_jinja
class TestWgcnaSections:
    """WGCNA-style module sections (gcen2 module -e / celltype-modules outputs)."""

    @pytest.fixture(scope="class")
    def outputs(self, tmp_path_factory):
        pytest.importorskip("matplotlib")
        from pathlib import Path

        from typer.testing import CliRunner

        from gcen2.cli import app

        data = Path(__file__).parent / "data" / "wgcna_tom"
        d = tmp_path_factory.mktemp("wgcna")
        files = {k: d / f"{k}.tsv" for k in ("m", "me", "kme", "hubs", "tc", "tomh")}
        files["tree"] = d / "tree.json"
        r = CliRunner().invoke(app, [
            "module", "-e", str(data / "modexpr.tsv"), "-o", str(files["m"]), "--signed",
            "--power", "6", "--min-size", "10", "--merge-cut", "0.25", "--me", str(files["me"]),
            "--kme", str(files["kme"]), "--hubs", str(files["hubs"]), "--n-hubs", "4",
            "--traits", str(data / "traits.tsv"), "--trait-cor", str(files["tc"]),
            "--dendrogram", str(files["tree"]), "--tom-hubs", str(files["tomh"])])
        assert r.exit_code == 0, r.output
        return files

    def test_section_hubs(self, outputs):
        import pandas as pd

        from gcen2.report import wgcna_section

        sec = wgcna_section(outputs["m"], hubs=outputs["hubs"], n_hubs=3)
        hubs = pd.read_csv(outputs["hubs"], sep="\t")
        assert sec["n_modules"] == hubs["module"].nunique() == len(sec["rows"])
        first = sec["rows"][0]
        assert first["hubs"] == list(hubs[hubs["module"] == 0].sort_values("rank")["gene"])[:3]
        # from kME when no hub file: the same genes
        sec2 = wgcna_section(outputs["m"], kme=outputs["kme"], n_hubs=4)
        assert [r["hubs"] for r in sec2["rows"]] == [
            list(hubs[hubs["module"] == m]["gene"]) for m in sorted(hubs["module"].unique())]
        assert wgcna_section(outputs["m"])["rows"][0]["hubs"] is None

    def test_report_with_figures(self, tmp_path, outputs):
        from gcen2.report import generate, wgcna_section

        sec = wgcna_section(outputs["m"], eigengenes=outputs["me"], hubs=outputs["hubs"],
                            trait_cor=outputs["tc"])
        out = tmp_path / "r.html"
        generate(output=out, wgcna=sec)
        html = out.read_text(encoding="utf-8")
        assert "WGCNA-style Modules" in html and "Hub genes" in html
        # module sizes, eigengene network, module-trait: inline vector SVG
        assert html.count('<figure class="gcen2-figure">') == 3
        assert html.count("<svg") == 3 and "Figure 3." in html
        assert "weight" in html

    def test_cli(self, tmp_path, outputs):
        from typer.testing import CliRunner

        from gcen2.cli import app

        out = tmp_path / "r.html"
        r = CliRunner().invoke(app, ["report", "-m", str(outputs["m"]), "--me",
                                     str(outputs["me"]), "--trait-cor", str(outputs["tc"]),
                                     "-o", str(out)])
        assert r.exit_code == 0, r.output
        assert out.read_text(encoding="utf-8").count('<figure class="gcen2-figure">') == 3
        r = CliRunner().invoke(app, ["report", "--me", str(outputs["me"]), "-o", str(out)])
        assert r.exit_code == 1

    def test_all_figures_and_figures_dir(self, tmp_path, outputs):
        from typer.testing import CliRunner

        from gcen2.cli import app

        out = tmp_path / "r.html"
        figs = tmp_path / "figs"
        r = CliRunner().invoke(app, [
            "report", "-m", str(outputs["m"]), "--me", str(outputs["me"]), "--kme",
            str(outputs["kme"]), "--hubs", str(outputs["hubs"]), "--trait-cor",
            str(outputs["tc"]), "--dendrogram", str(outputs["tree"]), "--tom-hubs",
            str(outputs["tomh"]), "-o", str(out), "--figures-dir", str(figs)])
        assert r.exit_code == 0, r.output
        html = out.read_text(encoding="utf-8")
        # sizes, dendrogram, eigengene network, module-trait, kME, hub network
        assert html.count('<figure class="gcen2-figure">') == 6
        assert "Figure 6." in html and "<img src=\"data:image/png" not in html
        assert len(list(figs.glob("*.pdf"))) == 6 and len(list(figs.glob("*.svg"))) == 6
        r = CliRunner().invoke(app, ["report", "--dendrogram", str(outputs["tree"]), "-o",
                                     str(out)])
        assert r.exit_code == 1

    def test_celltype_dir(self, tmp_path, outputs):
        import shutil

        import pandas as pd

        from gcen2.report import celltype_module_sections

        for t in ("T1", "T 2"):
            sub = tmp_path / t.replace(" ", "_")
            sub.mkdir()
            shutil.copy(outputs["m"], sub / "modules.tsv")
            shutil.copy(outputs["me"], sub / "eigengenes.tsv")
            shutil.copy(outputs["hubs"], sub / "hub_genes.tsv")
        pd.DataFrame({"cell_type": ["T1", "T 2"], "dir": ["T1", "T_2"], "n_cells": [500, 300],
                      "n_metacells": [100, 60], "n_genes": [200, 200], "power": [6.0, 6.0],
                      "modules": [7, 7], "genes_in_modules": [168, 168]}).to_csv(
            tmp_path / "celltypes.tsv", sep="\t", index=False, lineterminator="\n")
        secs = celltype_module_sections(tmp_path)
        assert [s["title"] for s in secs] == ["T1", "T 2"]
        assert "500 cells" in secs[0]["note"] and secs[1]["rows"][0]["hubs"]


@requires_jinja
def test_report_merge_extract_moran_metacells(tmp_path):
    """The figures of merge, extract, spatial-build --moran and the metacells
    of celltype-modules, in their own report sections, numbered in order."""
    pytest.importorskip("matplotlib")
    import re

    import numpy as np
    import pandas as pd
    from typer.testing import CliRunner

    from gcen2.cli import app

    rng = np.random.default_rng(0)
    genes = [f"g{i}" for i in range(30)]

    def net(name, n):
        a, b = np.triu_indices(30, 1)
        k = rng.choice(a.size, n, replace=False)
        f = tmp_path / f"{name}.network"
        pd.DataFrame({"gene_a": np.array(genes)[a[k]], "gene_b": np.array(genes)[b[k]],
                      "correlation": rng.uniform(0.6, 1, n), "p_value": 1e-5}).to_csv(
            f, sep="\t", index=False, lineterminator="\n")
        return f

    n1, n2, merged = net("a", 120), net("b", 100), net("m", 40)
    (tmp_path / "q.txt").write_text("g0\ng1\nmissing\n", encoding="utf-8")
    pd.DataFrame({"gene": genes, "morans_i": rng.normal(0.1, 0.1, 30),
                  "mean": rng.lognormal(0, 1, 30), "detected": rng.random(30)}).to_csv(
        tmp_path / "mi.tsv", sep="\t", index=False, lineterminator="\n")
    ct = tmp_path / "ct"
    (ct / "A").mkdir(parents=True)
    pd.DataFrame({"module": [1] * 15 + [2] * 15, "gene": genes}).to_csv(
        ct / "A" / "modules.tsv", sep="\t", index=False, lineterminator="\n")
    cells = [f"c{i}" for i in range(40)]
    pd.DataFrame({"metacell": np.repeat(["m0", "m1", "m2", "m3"], 5), "cell": cells[:20],
                  "center": np.tile([True] + [False] * 4, 4)}).to_csv(
        ct / "A" / "metacell_members.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame({"cell_type": ["A"], "dir": ["A"], "n_cells": [40], "n_metacells": [4],
                  "n_genes": [30], "power": [6.0], "modules": [2],
                  "genes_in_modules": [30]}).to_csv(ct / "celltypes.tsv", sep="\t", index=False, lineterminator="\n")
    pd.DataFrame(rng.normal(size=(40, 2)), index=cells, columns=["x", "y"]).to_csv(
        tmp_path / "umap.tsv", sep="\t", lineterminator="\n")
    out = tmp_path / "r.html"
    r = CliRunner().invoke(app, [
        "report", "-n", str(n1), "--merge-input", str(n1), "--merge-input", str(n2),
        "--merged", str(merged), "--extract-genes", str(tmp_path / "q.txt"),
        "--moran", str(tmp_path / "mi.tsv"), "--celltype-modules", str(ct),
        "--coords", str(tmp_path / "umap.tsv"), "-o", str(out),
        "--figures-dir", str(tmp_path / "fig")])
    assert r.exit_code == 0, r.output
    html = out.read_text(encoding="utf-8")
    for h in ("2b. Merged and Extracted Networks", "2c. Spatial Autocorrelation"):
        assert h in html
    caps = re.findall(r"Figure (\d+)\.</b> ([^<]{0,60})", html)
    assert [int(n) for n, _ in caps] == list(range(1, len(caps) + 1))
    first = [c for _, c in caps]
    assert first[1].startswith("Edges shared by the merged")
    assert first[2].startswith("Query genes of the extracted")
    assert first[3].startswith("Spatial autocorrelation")
    assert first[4].startswith("Cell types at a glance")
    assert first[5].startswith("Metacells on the cell embedding")
    for name in ("merge_overlap", "extract_overview", "moran", "celltypes_metacells"):
        assert (tmp_path / "fig" / f"{name}.svg").is_file()
    # options that need another one
    for args in (["--merged", str(merged)], ["--extract-genes", str(tmp_path / "q.txt")],
                 ["--coords", str(tmp_path / "umap.tsv")]):
        r = CliRunner().invoke(app, ["report", *args, "-o", str(tmp_path / "x.html")])
        assert r.exit_code == 1 and "ERROR" in r.output, args
