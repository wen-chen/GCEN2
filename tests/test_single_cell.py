"""Tests for single-cell functionality (pseudo-bulk, cell-type networks).

Requires anndata. Skipped when anndata is not installed.
"""

import numpy as np
import pandas as pd
import pytest

from gcen2.single_cell.core import build_by_celltype

try:
    import anndata as ad
    HAS_ANNDATA = True
except ImportError:
    HAS_ANNDATA = False

requires_sc = pytest.mark.skipif(not HAS_ANNDATA, reason="anndata not installed")


def _make_adata(n_cells=50, n_genes=20, n_celltypes=3, seed=42):
    """Create a small synthetic AnnData for testing."""
    rng = np.random.default_rng(seed)
    X = rng.poisson(lam=5.0, size=(n_cells, n_genes)).astype(float)
    cell_types = [f"type_{i % n_celltypes}" for i in range(n_cells)]
    obs = pd.DataFrame({"cell_type": cell_types},
                       index=[f"cell_{i}" for i in range(n_cells)])
    var = pd.DataFrame(index=[f"gene_{j}" for j in range(n_genes)])
    return ad.AnnData(X=X, obs=obs, var=var)


@requires_sc
class TestPseudobulk:
    def test_returns_dataframe(self):
        from gcen2.single_cell import pseudobulk
        adata = _make_adata()
        result = pseudobulk(adata, groupby="cell_type", min_cells=1)
        assert isinstance(result, pd.DataFrame)

    def test_columns_are_cell_types(self):
        from gcen2.single_cell import pseudobulk
        adata = _make_adata(n_celltypes=3)
        result = pseudobulk(adata, groupby="cell_type", min_cells=1)
        assert set(result.columns) == {"type_0", "type_1", "type_2"}

    def test_rows_are_genes(self):
        from gcen2.single_cell import pseudobulk
        adata = _make_adata(n_genes=20)
        result = pseudobulk(adata, groupby="cell_type", min_cells=1)
        assert len(result) == 20

    def test_sum_method(self):
        from gcen2.single_cell import pseudobulk
        adata = _make_adata()
        result_sum = pseudobulk(adata, groupby="cell_type", method="sum", min_cells=1)
        result_mean = pseudobulk(adata, groupby="cell_type", method="mean", min_cells=1)
        # sum should be >= mean for positive expression values
        assert (result_sum.values >= result_mean.values - 1e-9).all()

    def test_min_cells_filter(self):
        import warnings

        from gcen2.exceptions import GCENInputError
        from gcen2.single_cell import pseudobulk
        adata = _make_adata(n_cells=30, n_celltypes=3)
        # Every group has < 100 cells: each is skipped with a warning and,
        # since nothing is left to aggregate, a GCENInputError is raised.
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            with pytest.raises(GCENInputError, match="min_cells"):
                pseudobulk(adata, groupby="cell_type", min_cells=100)
        assert any("min_cells" in str(w.message) for w in caught)

    def test_min_cells_partial_filter(self):
        import warnings

        from gcen2.single_cell import pseudobulk
        # 31 cells over 3 types -> sizes 11/10/10; only type_0 survives
        adata = _make_adata(n_cells=31, n_celltypes=3)
        sizes = adata.obs["cell_type"].value_counts()
        threshold = int(sizes.max())
        with warnings.catch_warnings(record=True):
            warnings.simplefilter("always")
            result = pseudobulk(adata, groupby="cell_type", min_cells=threshold)
        assert list(result.columns) == ["type_0"]

    def test_invalid_groupby_raises(self):
        from gcen2.exceptions import GCENInputError
        from gcen2.single_cell import pseudobulk
        adata = _make_adata()
        with pytest.raises(GCENInputError):
            pseudobulk(adata, groupby="nonexistent_key")

    def test_sparse_matrix(self):
        import scipy.sparse as sp

        from gcen2.single_cell import pseudobulk
        adata = _make_adata()
        adata.X = sp.csr_matrix(adata.X)
        result = pseudobulk(adata, groupby="cell_type", min_cells=1)
        assert isinstance(result, pd.DataFrame)


@requires_sc
class TestBuildByCelltype:
    def test_returns_dict_of_networks(self):
        from gcen2.models import Network
        from gcen2.single_cell import build_by_celltype
        adata = _make_adata(n_cells=60, n_genes=15, n_celltypes=2)
        try:
            networks = build_by_celltype(
                adata, celltype_key="cell_type",
                min_cells=5, pseudo_replicates=5, strategy="pseudo_replicates",
                method="pearson", cor=0.0, pval=1.0, threads=1
            )
        except Exception as e:
            # Skip if C++ extension not available
            if "_gcen_core" in str(e) or "not built" in str(e):
                pytest.skip("C++ extension not compiled")
            raise
        assert isinstance(networks, dict)
        for net in networks.values():
            assert isinstance(net, Network)

    def test_use_highly_variable_requires_annotation(self):
        from gcen2.exceptions import GCENInputError
        from gcen2.single_cell import build_by_celltype
        adata = _make_adata()
        with pytest.raises(GCENInputError):
            build_by_celltype(adata, celltype_key="cell_type",
                              use_highly_variable=True)


@pytest.mark.requires_single_cell
class TestChunkedCorrelation:
    """Chunked computation gives same results as full computation."""

    def test_chunked_matches_full(self):
        """max_genes_per_chunk < n_genes should give same edge count as no chunking."""
        pytest.importorskip("anndata")
        import anndata as ad

        rng = np.random.default_rng(42)
        n_genes, n_cells, n_groups = 60, 300, 3
        X = np.abs(rng.standard_normal((n_cells, n_genes)))
        obs = pd.DataFrame(
            {'celltype': [f'CT{i % n_groups}' for i in range(n_cells)]},
            index=[f'cell{i}' for i in range(n_cells)]
        )
        var = pd.DataFrame(index=[f'G{i}' for i in range(n_genes)])
        adata = ad.AnnData(X=X, obs=obs, var=var)

        kwargs = dict(
            celltype_key='celltype',
            pseudo_replicates=5,
            strategy='pseudo_replicates',
            seed=42,
            cor=0.3,
            pval=0.1,
            threads=1
        )

        # Full computation
        net_full = build_by_celltype(adata, **kwargs)

        # Chunked computation (force chunking with small chunk_size)
        net_chunk = build_by_celltype(adata, max_genes_per_chunk=20, **kwargs)

        assert set(net_full.keys()) == set(net_chunk.keys())
        for ct in net_full:
            full_count = net_full[ct].edge_count
            chunk_count = net_chunk[ct].edge_count
            assert full_count == chunk_count, (
                f"CT={ct}: full={full_count} edges, chunked={chunk_count} edges"
            )

    def test_chunked_no_op_when_small(self):
        """When n_genes <= max_genes_per_chunk, behaviour is identical."""
        pytest.importorskip("anndata")
        import anndata as ad

        rng = np.random.default_rng(7)
        n_genes, n_cells = 10, 100
        X = np.abs(rng.standard_normal((n_cells, n_genes)))
        obs = pd.DataFrame(
            {'celltype': ['A'] * 50 + ['B'] * 50},
            index=[f'c{i}' for i in range(n_cells)]
        )
        var = pd.DataFrame(index=[f'G{i}' for i in range(n_genes)])
        adata = ad.AnnData(X=X, obs=obs, var=var)

        kwargs = dict(celltype_key='celltype', pseudo_replicates=4, strategy='pseudo_replicates', seed=0, cor=0.3, pval=0.1, threads=1)
        net_normal = build_by_celltype(adata, **kwargs)
        net_big_chunk = build_by_celltype(adata, max_genes_per_chunk=1000, **kwargs)

        for ct in net_normal:
            assert net_normal[ct].edge_count == net_big_chunk[ct].edge_count


def _sparse_adata(n_cells=120, n_genes=25, seed=0):
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    X = rng.poisson(0.8, (n_cells, n_genes)).astype(np.float32)
    obs = pd.DataFrame({"cell_type": [f"t{i % 3}" for i in range(n_cells)]},
                       index=[f"c{i}" for i in range(n_cells)])
    var = pd.DataFrame({"highly_variable": [i % 2 == 0 for i in range(n_genes)]},
                       index=[f"g{j}" for j in range(n_genes)])
    adata = ad.AnnData(X=sp.csr_matrix(X), obs=obs, var=var)
    adata.layers["counts"] = sp.csr_matrix(X * 3)
    return adata, X


@requires_sc
class TestAnnDataInput:
    def test_adata_matrix_subsets(self):
        from gcen2.single_cell import adata_matrix

        adata, X = _sparse_adata()
        m, genes = adata_matrix(adata)
        assert m.flags.c_contiguous and m.dtype == np.float64
        np.testing.assert_array_equal(m, X.T)
        m, genes = adata_matrix(adata, layer="counts", use_highly_variable=True,
                                groupby="cell_type", groups=["t0", "t2"])
        cells = np.array([i % 3 in (0, 2) for i in range(len(X))])
        np.testing.assert_array_equal(m, (X * 3)[cells][:, ::2].T)
        assert genes == [f"g{j}" for j in range(0, 25, 2)]
        m1, _ = adata_matrix(adata, groupby="cell_type", groups="t1")
        assert m1.shape == (25, 40)

    @pytest.mark.parametrize("kwargs, match", [
        ({"layer": "nope"}, "Layer"),
        ({"groupby": "nope", "groups": "a"}, "groupby"),
        ({"groupby": "cell_type"}, "groups is required"),
        ({"groupby": "cell_type", "groups": "zz"}, "No cells"),
    ])
    def test_adata_matrix_errors(self, kwargs, match):
        from gcen2.exceptions import GCENInputError
        from gcen2.single_cell import adata_matrix

        with pytest.raises(GCENInputError, match=match):
            adata_matrix(_sparse_adata()[0], **kwargs)

    def test_build_sparse_equals_dataframe(self):
        import gcen2

        adata, X = _sparse_adata()
        df = pd.DataFrame(X.T.astype(float), index=adata.var_names, columns=adata.obs_names)
        for method in ("pearson", "spearman"):
            got = gcen2.network.build(adata, method=method, cor=0.0, pval=1.0)
            ref = gcen2.network.build(df, method=method, cor=0.0, pval=1.0)
            pd.testing.assert_frame_equal(got.edges, ref.edges)
        em = gcen2.models.ExpressionMatrix(adata, np.array(adata.var_names),
                                           np.array(adata.obs_names), True)
        pd.testing.assert_frame_equal(gcen2.network.build(em, cor=0.0, pval=1.0).edges,
                                      gcen2.network.build(df, cor=0.0, pval=1.0).edges)

    def test_build_layer_hvg_and_to_file(self, tmp_path):
        import gcen2

        adata, X = _sparse_adata()
        sub = pd.DataFrame((X * 3)[:, ::2].T.astype(float),
                           index=adata.var_names[::2], columns=adata.obs_names)
        got = gcen2.network.build(adata, layer="counts", use_highly_variable=True,
                                  cor=0.0, pval=1.0)
        pd.testing.assert_frame_equal(got.edges,
                                      gcen2.network.build(sub, cor=0.0, pval=1.0).edges)
        gcen2.network.build_to_file(adata, tmp_path / "a.tsv", layer="counts",
                                    use_highly_variable=True, cor=0.0, pval=1.0)
        gcen2.network.build_to_file(sub, tmp_path / "b.tsv", cor=0.0, pval=1.0)
        assert (tmp_path / "a.tsv").read_bytes() == (tmp_path / "b.tsv").read_bytes()

    def test_layer_needs_anndata(self):
        import gcen2
        from gcen2.exceptions import GCENInputError

        df = pd.DataFrame(np.random.default_rng(0).random((5, 6)))
        with pytest.raises(GCENInputError, match="AnnData"):
            gcen2.network.build(df, layer="x")

    def test_cli_build_h5ad(self, tmp_path):
        typer_testing = pytest.importorskip("typer.testing")
        from gcen2.cli import app

        adata, X = _sparse_adata()
        adata.write_h5ad(tmp_path / "x.h5ad")
        cells = np.array([i % 3 == 1 for i in range(len(X))])
        sub = pd.DataFrame(X[cells].T.astype(float), index=adata.var_names,
                           columns=adata.obs_names[cells])
        sub.to_csv(tmp_path / "x.tsv", sep="\t", lineterminator="\n")
        runner = typer_testing.CliRunner()
        r = runner.invoke(app, ["build", "-i", str(tmp_path / "x.h5ad"), "-o",
                                str(tmp_path / "a.net"), "--groupby", "cell_type",
                                "--group", "t1", "-c", "0.2", "-p", "0.5"])
        assert r.exit_code == 0, r.output
        assert "40 cells x 25 genes" in r.output
        r = runner.invoke(app, ["build", "-i", str(tmp_path / "x.tsv"), "-o",
                                str(tmp_path / "b.net"), "-c", "0.2", "-p", "0.5"])
        assert r.exit_code == 0, r.output
        assert (tmp_path / "a.net").read_bytes() == (tmp_path / "b.net").read_bytes()
        r = runner.invoke(app, ["build", "-i", str(tmp_path / "x.tsv"), "-o",
                                str(tmp_path / "c.net"), "--hvg"])
        assert r.exit_code == 1 and "h5ad" in r.output
        r = runner.invoke(app, ["build", "-i", str(tmp_path / "x.h5ad"), "-o",
                                str(tmp_path / "c.net"), "--groupby", "cell_type"])
        assert r.exit_code == 1
