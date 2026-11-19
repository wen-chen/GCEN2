"""Integration tests for the complete GCEN2 pipeline.

Runs the full recommended pipeline (normalize → filter → build → module →
annotate → rwr) on the example data in the repository (data/example/bulk:
the GCEN v1 sample data) and verifies:
  1. Each step produces output of the expected shape / type.
  2. Module identification returns non-empty modules within the network's genes.
  3. RWR converges and returns a valid probability distribution.

The edge-by-edge comparison with the network GCEN v1 built
(``gene_co_expr.network``, 186 MB, not in the repository) runs when the
full GCEN v1 sample package is in data/sample.
"""

from __future__ import annotations

import os

import pytest

ROOT = os.path.join(os.path.dirname(__file__), "..")
EXAMPLE_DIR = os.path.join(ROOT, "data", "example", "bulk")
V1_DIR = os.path.join(ROOT, "data", "sample")  # local copy of the GCEN v1 package
HAS_V1_NETWORK = os.path.isfile(os.path.join(V1_DIR, "gene_co_expr.network"))

from gcen2 import annotation, module, network, preprocessing, rwr
from gcen2.io import read_network, read_tsv

# ---------------------------------------------------------------------------
# Helper
# ---------------------------------------------------------------------------


def _sample(filename: str) -> str:
    return os.path.join(EXAMPLE_DIR, filename)


@pytest.fixture(scope="module")
def expr():
    return read_tsv(_sample("gene_expr.tsv")).data


@pytest.fixture(scope="module")
def normed(expr):
    return preprocessing.normalize(expr, method="tmm")


@pytest.fixture(scope="module")
def filtered(normed):
    return preprocessing.filter_genes(normed, percent_mean=0.75)


@pytest.fixture(scope="module")
def net(filtered):
    """Spearman network of the filtered sample data (|r| >= 0.9, p <= 0.001)."""
    return network.build(filtered, method="spearman", cor=0.9, pval=0.001, threads=4)


@pytest.fixture(scope="module")
def files(tmp_path_factory, normed, net):
    d = tmp_path_factory.mktemp("integration")
    normed.to_csv(d / "gene_expr_norm.tsv", sep="\t", lineterminator="\n")
    net.to_edgelist(str(d / "net.network"))
    return d


# ---------------------------------------------------------------------------
# Core pipeline
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestBulkPipeline:
    """Full normalize → filter → build → module → annotate → rwr pipeline."""

    def test_read_expression_matrix(self, expr):
        assert expr.shape == (32940, 15)
        assert not expr.empty

    def test_normalize_tmm_preserves_shape(self, expr, normed):
        # Shape is preserved
        assert normed.shape == expr.shape
        # No negative values after TMM
        assert (normed.values >= 0).all()

    def test_normalize_all_methods(self, expr):
        for method in ("upqt", "median", "deseq", "tmm"):
            normed = preprocessing.normalize(expr, method=method)
            assert normed.shape == expr.shape, f"{method} changed shape"
            assert not normed.isnull().all().all(), f"{method} produced all-NaN output"

    def test_filter_reduces_genes(self, normed, filtered):
        # Filtering to top 75% should reduce gene count but keep some
        assert 0 < filtered.shape[0] < normed.shape[0]

    def test_normalized_file_loadable(self, files, normed):
        back = read_tsv(str(files / "gene_expr_norm.tsv")).data
        assert back.shape == normed.shape

    def test_network_build_from_filtered_data(self, net, filtered):
        assert net.edge_count > 1000
        genes = set(net.edges["gene_a"]) | set(net.edges["gene_b"])
        assert genes <= set(filtered.index)
        assert (net.edges["correlation"].abs() >= 0.9).all()
        assert (net.edges["p_value"] <= 0.001).all()

    @pytest.mark.skipif(not HAS_V1_NETWORK, reason="GCEN v1 reference network not present "
                        "(data/sample/gene_co_expr.network of the GCEN v1 0.6.3 package)")
    def test_network_matches_v1_reference(self):
        """Rebuild the full GCEN v1 reference network and compare edge by edge.

        ``gene_co_expr.network`` was produced by GCEN v1 from the raw,
        header-less ``gene_expr.tsv`` with Spearman, |r| >= 0.9, p <= 0.001,
        unsigned (33k genes, 3.6M edges). v1 prints values with 6 decimals.
        This file lacks every |r| == 0.9 edge, which a fresh v1 0.6.3 build
        does output (and matches v2 exactly), so those are allowed as extras.
        """
        ref = read_network(os.path.join(V1_DIR, "gene_co_expr.network")).edges
        expr = read_tsv(_sample("gene_expr.tsv")).data
        cor = 0.9
        got = network.build(
            expr, method="spearman", cor=cor, pval=0.001, signed=False, threads=4
        ).edges

        merged = ref.merge(got, on=["gene_a", "gene_b"], how="left",
                           suffixes=("_v1", "_v2"), indicator=True)
        missing = merged["_merge"] != "both"
        assert not missing.any(), (
            f"{missing.sum()} v1 edges not reproduced, e.g.\n{merged[missing].head()}"
        )
        assert (merged["correlation_v1"] - merged["correlation_v2"]).abs().max() < 1e-6
        rel_p = (merged["p_value_v1"] - merged["p_value_v2"]).abs() / merged["p_value_v1"]
        assert rel_p.max() < 1e-5

        # Extra v2 edges may only sit exactly on the |r| == cutoff boundary,
        # which the v1 0.6.3 binary rounded just below the threshold.
        extra = got.merge(ref[["gene_a", "gene_b"]], how="left", indicator=True)
        extra = extra[extra["_merge"] == "left_only"]
        assert (extra["correlation"].abs() == cor).all()

    def test_network_stats(self, net):
        stats = network.stats(net)
        assert stats.node_count > 0
        assert stats.edge_count > 0
        assert 0.0 <= stats.density <= 1.0
        assert stats.average_degree > 0

    def test_module_identification(self, net):
        modules = module.identify(net, similarity=0.5, threads=2)
        assert modules.module_count > 0
        # All module genes must come from the network
        net_genes = set(net.edges["gene_a"]) | set(net.edges["gene_b"])
        for m in modules.modules:
            assert m.issubset(net_genes), (
                f"Module contains genes not in network: {m - net_genes}"
            )

    def test_go_annotation(self, net):
        obo = _sample("go-basic.obo.gz")
        assoc = _sample("gene_go.assoc.gz")
        results = annotation.annotate_go(net, obo, assoc, pval=0.05, threads=1)
        assert len(results) > 0
        for er in results.values():
            assert er.p_cutoff == 0.05

    def test_rwr_prediction(self, net):
        from gcen2.io import read_gene_list

        seed_genes = read_gene_list(_sample("rwr_seed_genes.list"))
        assert seed_genes
        # Use first 5 seeds
        result = rwr.predict(net, seed_genes[:5], gamma=0.5, prob_cutoff=0.0)
        assert result.steps > 0, "RWR did not converge"
        assert len(result.ranked_genes) > 0
        # Probabilities should be non-negative
        assert (result.ranked_genes["probability"] >= 0).all()

    def test_enrichment_standalone(self):
        query_file = _sample("enrichment_gene.list")
        bg_file = _sample("background_gene.list")
        obo = _sample("go-basic.obo.gz")
        assoc = _sample("gene_go.assoc.gz")
        with open(query_file, encoding="utf-8") as f:
            query_genes = [line.strip() for line in f if line.strip()]
        with open(bg_file, encoding="utf-8") as f:
            bg_genes = [line.strip() for line in f if line.strip()]
        ontology = annotation.load_go(obo)
        assoc_dict = annotation.load_association(assoc, ontology)
        result = annotation.enrich(query_genes, bg_genes, ontology, assoc_dict, pval=0.05)
        assert len(result.results) > 0


# ---------------------------------------------------------------------------
# Network utilities
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestNetworkUtilities:

    def test_network_merge(self, filtered):
        # networks of two halves of the samples
        net1 = network.build(filtered.iloc[:3000, :8], method="pearson", cor=0.9, pval=0.01)
        net2 = network.build(filtered.iloc[:3000, 7:], method="pearson", cor=0.9, pval=0.01)
        merged = network.merge([net1, net2], cor_cutoff=0.0)
        assert merged.edge_count > 0
        # Direction column should be present (GCEN2 v1 bug fix)
        assert "direction" in merged.edges.columns, (
            "merge() must output a 'direction' column (GCEN2 v1 bug fix)"
        )
        assert set(merged.edges["direction"].unique()).issubset({1, -1, "1", "-1"})

    def test_network_extract(self, net):
        with open(_sample("gene_list.txt"), encoding="utf-8") as f:
            genes = [line.strip() for line in f if line.strip()]
        sub = network.extract(net, genes[:10])
        assert sub.edge_count <= net.edge_count

    def test_network_shuffle_preserves_degrees(self, net):
        shuffled = network.shuffle(net, swap_factor=1, seed=42)
        # Same edge count after shuffle
        assert shuffled.edge_count == net.edge_count
        # Same node set
        orig_nodes = set(net.edges["gene_a"]) | set(net.edges["gene_b"])
        shuf_nodes = set(shuffled.edges["gene_a"]) | set(shuffled.edges["gene_b"])
        assert orig_nodes == shuf_nodes


# ---------------------------------------------------------------------------
# CLI integration
# ---------------------------------------------------------------------------

@pytest.mark.integration
class TestCLIPipeline:

    def _get_runner(self):
        try:
            from typer.testing import CliRunner

            from gcen2.cli import app
            return CliRunner(), app
        except ImportError:
            pytest.skip("typer not installed")

    def test_cli_normalize(self, tmp_path):
        runner, app = self._get_runner()
        output = str(tmp_path / "normed.tsv")
        result = runner.invoke(
            app, ["normalize", "-i", _sample("gene_expr.tsv"), "-o", output, "-m", "tmm"]
        )
        assert result.exit_code == 0, result.output
        assert os.path.isfile(output)

    def test_cli_filter(self, tmp_path, files):
        runner, app = self._get_runner()
        output = str(tmp_path / "filtered.tsv")
        result = runner.invoke(
            app,
            ["filter", "-i", str(files / "gene_expr_norm.tsv"), "-o", output, "-p", "0.75"],
        )
        assert result.exit_code == 0, result.output
        assert os.path.isfile(output)

    def test_cli_stats(self, tmp_path, files):
        runner, app = self._get_runner()
        result = runner.invoke(app, ["stats", "-i", str(files / "net.network")])
        assert result.exit_code == 0, result.output
        # Output should contain numeric stats
        assert any(c.isdigit() for c in result.output)

    def test_cli_full_pipeline(self, tmp_path):
        """normalize → filter → build end-to-end via CLI."""
        runner, app = self._get_runner()
        normed = str(tmp_path / "normed.tsv")
        filtered = str(tmp_path / "filtered.tsv")
        net_out = str(tmp_path / "net.tsv")

        r1 = runner.invoke(
            app, ["normalize", "-i", _sample("gene_expr.tsv"), "-o", normed, "-m", "tmm"]
        )
        assert r1.exit_code == 0, r1.output

        r2 = runner.invoke(app, ["filter", "-i", normed, "-o", filtered, "-p", "0.75"])
        assert r2.exit_code == 0, r2.output

        r3 = runner.invoke(
            app,
            ["build", "-i", filtered, "-o", net_out, "-m", "spearman", "-c", "0.8", "-p", "0.001"],
        )
        assert r3.exit_code == 0, r3.output
        assert os.path.isfile(net_out)
