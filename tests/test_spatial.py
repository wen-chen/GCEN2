"""Tests for gcen2.spatial module."""
from __future__ import annotations

import numpy as np
import pytest

try:
    import anndata  # noqa: F401
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False

requires_single_cell = pytest.mark.skipif(
    not HAS_ANNDATA, reason="anndata not installed"
)


@pytest.fixture
def grid_coords():
    """4×5 grid of spots, shape (20, 2)."""
    return np.array([[i, j] for i in range(4) for j in range(5)], dtype=float)


@pytest.fixture
def mock_adata(grid_coords):
    """20 spots, 50 genes, two regions A/B."""
    anndata = pytest.importorskip("anndata")
    import pandas as pd

    rng = np.random.default_rng(42)
    X = rng.random((20, 50)).astype(np.float32)
    obs = pd.DataFrame(
        {"region": ["A"] * 10 + ["B"] * 10},
        index=[f"s{i}" for i in range(20)],
    )
    var = pd.DataFrame(index=[f"G{i}" for i in range(50)])
    adata = anndata.AnnData(X=X, obs=obs, var=var)
    adata.obsm["spatial"] = grid_coords
    return adata


# ── Task 1 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestBuildNeighborGraph:
    def test_knn_shape(self, grid_coords):
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, n_neighbors=4)
        assert W.shape == (20, 20)

    def test_knn_row_nnz(self, grid_coords):
        """Each row should have at least n_neighbors non-zeros (symmetry may add more)."""
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, n_neighbors=4)
        for i in range(20):
            assert W.getrow(i).nnz >= 4, f"row {i} has {W.getrow(i).nnz} neighbors, expected >= 4"

    def test_radius_distant_zero(self, grid_coords):
        """Spots farther than radius should have weight 0."""
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, radius=1.5)
        # spot 0=(0,0) and spot 19=(3,4) are distance ~5, should be 0
        assert W[0, 19] == 0.0

    def test_radius_near_nonzero(self, grid_coords):
        """Adjacent spots (dist=1) should have nonzero weight with radius=2."""
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, radius=2.0)
        # spot 0=(0,0), spot 1=(0,1), dist=1 → should be neighbors
        assert W[0, 1] > 0.0

    def test_both_raises(self, grid_coords):
        from gcen2.exceptions import GCENInputError
        from gcen2.spatial.neighbors import build_neighbor_graph
        with pytest.raises(GCENInputError, match="not both"):
            build_neighbor_graph(grid_coords, n_neighbors=4, radius=1.0)

    def test_neither_raises(self, grid_coords):
        from gcen2.exceptions import GCENInputError
        from gcen2.spatial.neighbors import build_neighbor_graph
        with pytest.raises(GCENInputError, match="Must specify"):
            build_neighbor_graph(grid_coords, n_neighbors=None, radius=None)

    def test_weights_positive(self, grid_coords):
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, n_neighbors=4)
        assert (W.data > 0).all()

    def test_n_neighbors_zero_raises(self, grid_coords):
        from gcen2.exceptions import GCENInputError
        from gcen2.spatial.neighbors import build_neighbor_graph
        with pytest.raises(GCENInputError, match="n_neighbors must be >= 1"):
            build_neighbor_graph(grid_coords, n_neighbors=0)

    def test_radius_negative_raises(self, grid_coords):
        from gcen2.exceptions import GCENInputError
        from gcen2.spatial.neighbors import build_neighbor_graph
        with pytest.raises(GCENInputError, match="radius must be non-negative"):
            build_neighbor_graph(grid_coords, radius=-1.0)


# ── Task 2 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestBuildSpatialFiltered:
    def test_returns_network(self, mock_adata):
        from gcen2.models import Network
        from gcen2.spatial.core import build_spatial_filtered
        net = build_spatial_filtered(mock_adata, n_neighbors=4, cor_cutoff=0.0, p_cutoff=1.0)
        assert isinstance(net, Network)

    def test_fewer_edges_than_unfiltered(self, mock_adata):
        """Spatial filtering should remove some edges vs. unfiltered build."""
        import pandas as pd

        import gcen2
        from gcen2.spatial.core import _morans_i, build_spatial_filtered
        from gcen2.spatial.neighbors import build_neighbor_graph

        X = mock_adata.X if not hasattr(mock_adata.X, "toarray") else mock_adata.X.toarray()
        expr_df = pd.DataFrame(X.T, index=mock_adata.var_names,
                               columns=mock_adata.obs_names)
        unfiltered = gcen2.network.build(expr_df, cor=0.0, pval=1.0)
        filtered = build_spatial_filtered(mock_adata, n_neighbors=4, cor_cutoff=0.3, p_cutoff=1.0,
                                          moran_fdr=1.0)
        assert filtered.edge_count <= unfiltered.edge_count

        # Verify every edge in filtered network connects two spatially autocorrelated genes
        W = build_neighbor_graph(mock_adata.obsm["spatial"], n_neighbors=4)
        mi = _morans_i(X, W)
        gene_names = list(mock_adata.var_names)
        gene_to_mi = dict(zip(gene_names, mi))
        if filtered.edge_count > 0:
            for _, row in filtered.edges.iterrows():
                assert gene_to_mi.get(row["gene_a"], 0) > 0
                assert gene_to_mi.get(row["gene_b"], 0) > 0

    def test_moran_fdr_drops_genes_without_a_pattern(self):
        """Two co-expressed genes with a spatial gradient and two co-expressed
        noise genes: I > 0 alone may keep the noise pair, the test does not."""
        import anndata as ad

        from gcen2.spatial.core import build_spatial_filtered

        rng = np.random.default_rng(3)
        n = 300
        coords = rng.uniform(0, 10, (n, 2))
        grad = coords[:, 0]
        noise = rng.normal(size=n)
        X = np.column_stack([grad + rng.normal(0, 0.5, n), grad + rng.normal(0, 0.5, n),
                             noise + rng.normal(0, 0.3, n), noise + rng.normal(0, 0.3, n)])
        a = ad.AnnData(X, obsm={"spatial": coords})
        a.var_names = ["s1", "s2", "n1", "n2"]
        kw = dict(n_neighbors=6, cor_cutoff=0.5, p_cutoff=1.0, fraction=0, normalize=False)
        pairs = {tuple(sorted(e)) for e in
                 build_spatial_filtered(a, **kw).edges[["gene_a", "gene_b"]].to_numpy()}
        assert pairs == {("s1", "s2")}
        with pytest.raises(Exception, match="moran_fdr"):
            build_spatial_filtered(a, moran_fdr=0, **kw)

    def test_moran_test_matches_permutations(self):
        from gcen2.spatial.core import _morans_i, _morans_test
        from gcen2.spatial.neighbors import build_neighbor_graph

        rng = np.random.default_rng(1)
        n = 200
        W = build_neighbor_graph(rng.uniform(0, 10, (n, 2)), n_neighbors=6)
        X = np.column_stack([rng.poisson(2, n), rng.exponential(1, n) ** 2]).astype(float)
        I, z, p = _morans_test(X, W)
        perm = np.array([_morans_i(X[rng.permutation(n)], W) for _ in range(2000)])
        np.testing.assert_allclose((I + 1 / (n - 1)) / z, perm.std(axis=0), rtol=0.08)
        assert ((p > 0) & (p < 1)).all()

    def test_missing_spatial_raises(self, mock_adata):
        from gcen2.exceptions import GCENInputError
        from gcen2.spatial.core import build_spatial_filtered
        del mock_adata.obsm["spatial"]
        with pytest.raises(GCENInputError, match="spatial"):
            build_spatial_filtered(mock_adata, n_neighbors=4)

    def test_empty_network_early_return(self, mock_adata):
        """Very high cor_cutoff -> empty initial network -> empty filtered network."""
        from gcen2.models import Network
        from gcen2.spatial.core import build_spatial_filtered
        net = build_spatial_filtered(mock_adata, n_neighbors=4, cor_cutoff=2.0, p_cutoff=1.0)
        assert isinstance(net, Network)
        assert net.edge_count == 0


# ── Task 3 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestBuildSpatialWeighted:
    def test_returns_network(self, mock_adata):
        from gcen2.models import Network
        from gcen2.spatial.core import build_spatial_weighted
        net = build_spatial_weighted(mock_adata, n_neighbors=4, cor_cutoff=0.0, p_cutoff=1.0)
        assert isinstance(net, Network)

    def test_positive_edge_count(self, mock_adata):
        from gcen2.spatial.core import build_spatial_weighted
        net = build_spatial_weighted(mock_adata, n_neighbors=4, cor_cutoff=0.0, p_cutoff=1.0)
        assert net.edge_count > 0

    def test_bandwidth_param(self, mock_adata):
        """Explicit bandwidth should not raise."""
        from gcen2.spatial.core import build_spatial_weighted
        net = build_spatial_weighted(mock_adata, n_neighbors=4, bandwidth=2.0,
                                     cor_cutoff=0.0, p_cutoff=1.0)
        assert net.edge_count > 0

    def test_missing_spatial_raises(self, mock_adata):
        from gcen2.exceptions import GCENInputError
        from gcen2.spatial.core import build_spatial_weighted
        del mock_adata.obsm["spatial"]
        with pytest.raises(GCENInputError, match="spatial"):
            build_spatial_weighted(mock_adata, n_neighbors=4)


# ── Task 4 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestBuildSpatialRegions:
    def test_manual_region_key(self, mock_adata):
        from gcen2.models import Network
        from gcen2.spatial.core import build_spatial_regions
        result = build_spatial_regions(
            mock_adata, region_key="region", n_neighbors=4,
            cor=0.0, pval=1.0,
        )
        assert isinstance(result, dict)
        assert set(result.keys()) == {"A", "B"}
        for net in result.values():
            assert isinstance(net, Network)

    def test_region_key_missing_raises(self, mock_adata):
        from gcen2.exceptions import GCENInputError
        from gcen2.spatial.core import build_spatial_regions
        with pytest.raises(GCENInputError, match="not found in adata.obs"):
            build_spatial_regions(mock_adata, region_key="nonexistent", n_neighbors=4)

    def test_small_region_skipped(self, mock_adata):
        """Regions with fewer than 3 spots should be skipped with a warning."""
        mock_adata.obs["tiny"] = ["X"] * 19 + ["Y"]   # Y has only 1 spot
        from gcen2.spatial.core import build_spatial_regions
        with pytest.warns(UserWarning, match="only"):
            result = build_spatial_regions(
                mock_adata, region_key="tiny", n_neighbors=4,
                cor=0.0, pval=1.0,
            )
        assert "Y" not in result

    def test_auto_regions(self, mock_adata):
        """Auto KMeans: result should have n_regions keys."""
        pytest.importorskip("sklearn")
        from gcen2.spatial.core import build_spatial_regions
        result = build_spatial_regions(
            mock_adata, region_key=None, n_regions=3, n_neighbors=4,
            cor=0.0, pval=1.0,
        )
        assert isinstance(result, dict)
        assert len(result) >= 1


# ── Task 5 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestSpatialPublicAPI:
    def test_import_from_gcen2_spatial(self):
        from gcen2 import spatial
        assert hasattr(spatial, "build_spatial_filtered")
        assert hasattr(spatial, "build_spatial_weighted")
        assert hasattr(spatial, "build_spatial_regions")

    def test_import_via_gcen2(self):
        import gcen2
        assert hasattr(gcen2, "spatial")


@requires_single_cell
class TestSpatialAutocorrelation:
    def test_values_and_smoothing(self, mock_adata):
        from gcen2.spatial import spatial_autocorrelation
        from gcen2.spatial.core import _morans_i
        from gcen2.spatial.neighbors import build_neighbor_graph

        a = mock_adata.copy()
        a.X[:, 0] = a.obsm["spatial"][:, 0] * 3 + 1  # a gradient: strongly autocorrelated
        t = spatial_autocorrelation(a, fraction=0, normalize=False, n_neighbors=4,
                                    smoothed=True)
        assert list(t.columns) == ["gene", "morans_i", "pval", "fdr", "mean", "detected",
                                   "morans_i_smoothed"]
        assert t["gene"].iloc[0] == "G0"
        assert t["fdr"].iloc[0] < 1e-3 and (t["fdr"] >= t["pval"]).all()
        assert t["morans_i"].is_monotonic_decreasing
        W = build_neighbor_graph(a.obsm["spatial"], n_neighbors=4)
        X = np.asarray(a.X, dtype=float)
        want = dict(zip(a.var_names, _morans_i(X, W)))
        np.testing.assert_allclose(t["morans_i"], t["gene"].map(want), rtol=1e-12)
        # smoothing averages neighbours: autocorrelation goes up
        assert (t["morans_i_smoothed"] > t["morans_i"]).mean() > 0.9

    def test_cli_moran(self, tmp_path, mock_adata):
        import pandas as pd
        from typer.testing import CliRunner

        from gcen2.cli import app

        mock_adata.write_h5ad(tmp_path / "a.h5ad")
        for mode, cols in (("filtered", 6), ("weighted", 7)):
            out = tmp_path / f"{mode}.tsv"
            r = CliRunner().invoke(app, ["spatial-build", "-i", str(tmp_path / "a.h5ad"),
                                         "-o", str(tmp_path / f"{mode}.network"), "-m", mode,
                                         "-c", "0.3", "-p", "1", "--n-neighbors", "4",
                                         "--moran", str(out)])
            assert r.exit_code == 0, r.output
            t = pd.read_csv(out, sep="\t")
            assert t.shape == (50, cols)
            assert "Moran's I of 50 genes" in r.output
