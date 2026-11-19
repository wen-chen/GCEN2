# Spatial Co-Expression Network Analysis — Implementation Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 实现 `gcen2/spatial/` 模块，支持三种空间转录组共表达网络分析方法，并提供 `gcen2 spatial-build` CLI 命令。

**Architecture:** 独立模块目录 `gcen2/spatial/`，`neighbors.py` 提供共用邻居图构建，`core.py` 实现三种建网方法，均通过 `gcen2/spatial/__init__.py` 对外暴露。三种方法共用同一 `build_neighbor_graph()` 入口，返回现有 `Network` 或 `dict[str, Network]`，不引入新数据类。

**Tech Stack:** Python 3.12, anndata, numpy, scipy (KDTree + sparse), gcen2.network.build, typer (CLI)

**设计文档：** `docs/plans/2026-03-16-spatial-design.md`

---

## 背景知识（执行前必读）

### 现有相关 API

| 函数 | 文件 | 说明 |
|------|------|------|
| `gcen2.network.build(expr_df, ...)` | `gcen2/network.py` | 标准共表达网络构建，输入 genes×samples DataFrame |
| `gcen2.single_cell.pseudobulk(adata, groupby)` | `gcen2/single_cell/core.py` | 按分组聚合 AnnData → DataFrame |
| `Network.to_edgelist(path)` | `gcen2/models.py` | 输出网络为 TSV |
| `_require_anndata()` | `gcen2/single_cell/core.py` | 检查 anndata 是否安装，未安装时抛 GCENInputError |
| `GCENInputError` | `gcen2/exceptions.py` | 输入类型/值错误 |

### AnnData 空间数据格式
```python
adata.X              # (n_spots, n_genes) 表达矩阵
adata.obs            # spot 元数据 DataFrame，含 region 等 obs 列
adata.var_names      # 基因名列表
adata.obsm['spatial']  # (n_spots, 2) 空间坐标，float64
adata.layers['xxx']  # 可选额外层
```

### 三种方法简介

1. **空间过滤**：标准建网 → 计算每个基因的 Moran's I（空间自相关）→ 保留两端基因均有空间信号（Moran's I > 0）的边
2. **空间加权**：用空间邻居均值平滑表达矩阵（dropout 降噪）→ 用平滑后的矩阵建标准网络
3. **分区域建网**：按区域切分 spots → 每个区域聚合 → 各自建网 → 返回字典

### 测试守卫装饰器
```python
# 与 test_single_cell.py 相同模式
requires_single_cell = pytest.mark.requires_single_cell
# 在 conftest.py 中已注册 requires_single_cell marker
```

---

## Task 1：空间邻居图 `gcen2/spatial/neighbors.py`

**Files:**
- Create: `gcen2/spatial/__init__.py`（空文件，Task 5 再填充）
- Create: `gcen2/spatial/neighbors.py`
- Test: `tests/test_spatial.py`

### Step 1: 创建 `tests/test_spatial.py`，写失败测试

```python
"""Tests for gcen2.spatial module."""
from __future__ import annotations

import numpy as np
import pytest

requires_single_cell = pytest.mark.requires_single_cell


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
        """Each row should have exactly n_neighbors non-zeros."""
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, n_neighbors=4)
        for i in range(20):
            assert W.getrow(i).nnz == 4, f"row {i} has {W.getrow(i).nnz} neighbors, expected 4"

    def test_radius_distant_zero(self, grid_coords):
        """Spots farther than radius should have weight 0."""
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, radius=1.5)
        # spot (0,0) and spot (3,4) are distance ~5, should be 0
        assert W[0, 19] == 0.0

    def test_radius_near_nonzero(self, grid_coords):
        """Adjacent spots (dist=1) should have nonzero weight with radius=2."""
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, radius=2.0)
        # spot 0 = (0,0), spot 1 = (0,1), dist=1 → should be neighbors
        assert W[0, 1] > 0.0

    def test_both_raises(self, grid_coords):
        from gcen2.spatial.neighbors import build_neighbor_graph
        from gcen2.exceptions import GCENInputError
        with pytest.raises(GCENInputError, match="not both"):
            build_neighbor_graph(grid_coords, n_neighbors=4, radius=1.0)

    def test_neither_raises(self, grid_coords):
        from gcen2.spatial.neighbors import build_neighbor_graph
        from gcen2.exceptions import GCENInputError
        with pytest.raises(GCENInputError, match="Must specify"):
            build_neighbor_graph(grid_coords, n_neighbors=None, radius=None)

    def test_weights_positive(self, grid_coords):
        from gcen2.spatial.neighbors import build_neighbor_graph
        W = build_neighbor_graph(grid_coords, n_neighbors=4)
        assert (W.data > 0).all()
```

### Step 2: 确认测试失败

```bash
cd /mnt/d/code_base/gcen2
python3 -m pytest tests/test_spatial.py::TestBuildNeighborGraph -v --tb=short 2>&1 | tail -10
```
期望：`ModuleNotFoundError: No module named 'gcen2.spatial'`

### Step 3: 创建 `gcen2/spatial/__init__.py`（占位）

```python
"""Spatial transcriptomics co-expression network analysis."""
```

### Step 4: 实现 `gcen2/spatial/neighbors.py`

```python
"""Spatial neighbor graph construction (KNN or radius-based)."""
from __future__ import annotations

import numpy as np
import scipy.sparse as sp

from gcen2.exceptions import GCENInputError


def build_neighbor_graph(
    coords: np.ndarray,
    n_neighbors: int | None = 6,
    radius: float | None = None,
) -> sp.csr_matrix:
    """Build a Gaussian-weighted spatial neighbor graph.

    Parameters
    ----------
    coords:
        Array of shape ``(n_spots, 2)`` with 2-D spatial coordinates.
    n_neighbors:
        Number of nearest neighbors (KNN mode). Mutually exclusive with *radius*.
    radius:
        Distance cutoff (radius mode). Mutually exclusive with *n_neighbors*.

    Returns
    -------
    scipy.sparse.csr_matrix of shape ``(n_spots, n_spots)`` with Gaussian
    edge weights ``exp(-d² / (2σ²))``, where σ = median neighbor distance.
    """
    if n_neighbors is not None and radius is not None:
        raise GCENInputError("Specify either n_neighbors or radius, not both")
    if n_neighbors is None and radius is None:
        raise GCENInputError("Must specify either n_neighbors or radius")

    from scipy.spatial import KDTree

    coords = np.asarray(coords, dtype=float)
    n = len(coords)
    tree = KDTree(coords)

    if n_neighbors is not None:
        # Query k+1 because the first result is the point itself (distance 0)
        dists, indices = tree.query(coords, k=n_neighbors + 1)
        dists = dists[:, 1:]      # drop self (col 0)
        indices = indices[:, 1:]  # drop self (col 0)

        sigma = float(np.median(dists[dists > 0])) or 1.0
        weights = np.exp(-(dists.ravel() ** 2) / (2 * sigma ** 2))
        rows = np.repeat(np.arange(n), n_neighbors)
        cols = indices.ravel()

    else:  # radius mode
        pairs = tree.query_pairs(radius, output_type="ndarray")
        if len(pairs) == 0:
            return sp.csr_matrix((n, n), dtype=float)

        i_idx, j_idx = pairs[:, 0], pairs[:, 1]
        dists_pair = np.linalg.norm(coords[i_idx] - coords[j_idx], axis=1)
        sigma = float(np.median(dists_pair)) or 1.0
        w = np.exp(-(dists_pair ** 2) / (2 * sigma ** 2))

        # Symmetric: add both directions
        rows = np.concatenate([i_idx, j_idx])
        cols = np.concatenate([j_idx, i_idx])
        weights = np.concatenate([w, w])

    W = sp.csr_matrix((weights, (rows, cols)), shape=(n, n), dtype=float)
    return W
```

### Step 5: 运行测试确认通过

```bash
python3 -m pytest tests/test_spatial.py::TestBuildNeighborGraph -v --tb=short 2>&1 | tail -15
```
期望：7 passed（或因无 anndata 而 skipped）

### Step 6: Commit

```bash
cd /mnt/d/code_base/gcen2
git add gcen2/spatial/__init__.py gcen2/spatial/neighbors.py tests/test_spatial.py
git commit -m "feat(spatial): add build_neighbor_graph (KNN + radius, Gaussian weights)"
```

---

## Task 2：空间过滤建网 `build_spatial_filtered`

**Files:**
- Create: `gcen2/spatial/core.py`
- Modify: `tests/test_spatial.py`（追加测试类）

### Step 1: 追加测试到 `tests/test_spatial.py`

```python
# ── Task 2 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestBuildSpatialFiltered:
    def test_returns_network(self, mock_adata):
        from gcen2.spatial.core import build_spatial_filtered
        from gcen2.models import Network
        net = build_spatial_filtered(mock_adata, n_neighbors=4, cor_cutoff=0.0, p_cutoff=1.0)
        assert isinstance(net, Network)

    def test_fewer_edges_than_unfiltered(self, mock_adata):
        """Spatial filtering should remove some edges vs. unfiltered build."""
        import gcen2
        from gcen2.spatial.core import build_spatial_filtered
        import pandas as pd

        # Build unfiltered reference network
        X = mock_adata.X if not hasattr(mock_adata.X, "toarray") else mock_adata.X.toarray()
        expr_df = pd.DataFrame(X.T, index=mock_adata.var_names,
                               columns=mock_adata.obs_names)
        unfiltered = gcen2.network.build(expr_df, cor_cutoff=0.0, p_cutoff=1.0)
        filtered = build_spatial_filtered(mock_adata, n_neighbors=4, cor_cutoff=0.0, p_cutoff=1.0)
        # Spatial filter can only remove edges, not add
        assert filtered.edge_count <= unfiltered.edge_count

    def test_missing_spatial_raises(self, mock_adata):
        from gcen2.spatial.core import build_spatial_filtered
        from gcen2.exceptions import GCENInputError
        del mock_adata.obsm["spatial"]
        with pytest.raises(GCENInputError, match="spatial"):
            build_spatial_filtered(mock_adata, n_neighbors=4)
```

### Step 2: 确认测试失败

```bash
python3 -m pytest tests/test_spatial.py::TestBuildSpatialFiltered -v --tb=short 2>&1 | tail -10
```
期望：`ModuleNotFoundError: No module named 'gcen2.spatial.core'`

### Step 3: 创建 `gcen2/spatial/core.py` 实现 `build_spatial_filtered`

```python
"""Spatial co-expression network analysis: three build methods."""
from __future__ import annotations

import warnings
from typing import Any

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.models import Network
from gcen2.spatial.neighbors import build_neighbor_graph


def _require_anndata() -> None:
    try:
        import anndata  # noqa: F401
    except ImportError as exc:
        raise GCENInputError(
            "anndata is required for spatial analysis: pip install gcen2[single-cell]"
        ) from exc


def _get_expression(adata: Any, layer: str | None) -> np.ndarray:
    """Return dense (n_spots, n_genes) expression array."""
    import scipy.sparse as sp

    X = adata.layers[layer] if layer else adata.X
    if sp.issparse(X):
        X = X.toarray()
    return np.asarray(X, dtype=float)


def _get_coords(adata: Any) -> np.ndarray:
    if "spatial" not in adata.obsm:
        raise GCENInputError(
            "No spatial coordinates found in adata.obsm['spatial']. "
            "Make sure you have loaded a spatial AnnData object."
        )
    return np.asarray(adata.obsm["spatial"], dtype=float)


def _morans_i(X: np.ndarray, W: "scipy.sparse.csr_matrix") -> np.ndarray:
    """Compute Moran's I for each gene (column of X).

    Parameters
    ----------
    X : (n_spots, n_genes) expression matrix (mean-centered per gene)
    W : (n_spots, n_spots) spatial weight matrix

    Returns
    -------
    Array of shape (n_genes,) with Moran's I values.
    """
    import scipy.sparse as sp

    n = X.shape[0]
    w_sum = W.sum()
    if w_sum == 0:
        return np.zeros(X.shape[1])

    z = X - X.mean(axis=0)          # mean-center each gene
    z2_sum = (z ** 2).sum(axis=0)   # sum of squared deviations per gene
    z2_sum = np.where(z2_sum == 0, 1e-12, z2_sum)

    # Wz = W @ z  → spatial lag
    Wz = W @ z                       # (n_spots, n_genes)
    numerator = (z * Wz).sum(axis=0)
    I = (n / w_sum) * (numerator / z2_sum)
    return np.asarray(I).ravel()


def build_spatial_filtered(
    adata: Any,
    *,
    layer: str | None = None,
    n_neighbors: int | None = 6,
    radius: float | None = None,
    cor_cutoff: float = 0.6,
    p_cutoff: float = 0.05,
    **build_kwargs,
) -> Network:
    """Build a co-expression network filtered by spatial autocorrelation.

    Constructs a standard co-expression network, then removes edges where
    either gene lacks significant spatial autocorrelation (Moran's I ≤ 0).

    Parameters
    ----------
    adata:
        AnnData with ``obsm['spatial']`` coordinates.
    layer:
        Expression layer to use (``adata.X`` if None).
    n_neighbors:
        KNN neighbours. Mutually exclusive with *radius*.
    radius:
        Radius cutoff. Mutually exclusive with *n_neighbors*.
    cor_cutoff:
        Minimum absolute correlation for an edge.
    p_cutoff:
        Maximum p-value for an edge.

    Returns
    -------
    Network with spatially unsupported edges removed.
    """
    _require_anndata()
    import gcen2

    coords = _get_coords(adata)
    X = _get_expression(adata, layer)          # (n_spots, n_genes)
    gene_names = list(adata.var_names)

    # Build standard network
    expr_df = pd.DataFrame(X.T, index=gene_names, columns=list(adata.obs_names))
    net = gcen2.network.build(expr_df, cor_cutoff=cor_cutoff, p_cutoff=p_cutoff,
                               **build_kwargs)

    if net.edge_count == 0:
        return net

    # Build spatial neighbor graph and compute Moran's I per gene
    W = build_neighbor_graph(coords, n_neighbors=n_neighbors, radius=radius)
    mi = _morans_i(X, W)                       # (n_genes,)
    gene_to_mi = dict(zip(gene_names, mi))

    # Keep only edges where both genes have positive Moran's I
    edges = net.edges.copy()
    mask = (
        edges["gene_a"].map(gene_to_mi).fillna(0) > 0
    ) & (
        edges["gene_b"].map(gene_to_mi).fillna(0) > 0
    )
    return Network(edges=edges[mask].reset_index(drop=True))
```

### Step 4: 运行测试确认通过

```bash
python3 -m pytest tests/test_spatial.py::TestBuildSpatialFiltered -v --tb=short 2>&1 | tail -10
```
期望：3 passed

### Step 5: Commit

```bash
git add gcen2/spatial/core.py tests/test_spatial.py
git commit -m "feat(spatial): add build_spatial_filtered (Moran's I edge filter)"
```

---

## Task 3：空间加权建网 `build_spatial_weighted`

**Files:**
- Modify: `gcen2/spatial/core.py`（追加函数）
- Modify: `tests/test_spatial.py`（追加测试类）

### Step 1: 追加测试到 `tests/test_spatial.py`

```python
# ── Task 3 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestBuildSpatialWeighted:
    def test_returns_network(self, mock_adata):
        from gcen2.spatial.core import build_spatial_weighted
        from gcen2.models import Network
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
        assert net.edge_count >= 0
```

### Step 2: 确认测试失败

```bash
python3 -m pytest tests/test_spatial.py::TestBuildSpatialWeighted -v --tb=short 2>&1 | tail -10
```
期望：`ImportError: cannot import name 'build_spatial_weighted'`

### Step 3: 追加 `build_spatial_weighted` 到 `gcen2/spatial/core.py`

在文件末尾追加：

```python
def build_spatial_weighted(
    adata: Any,
    *,
    layer: str | None = None,
    n_neighbors: int | None = 6,
    radius: float | None = None,
    bandwidth: float | None = None,
    cor_cutoff: float = 0.6,
    p_cutoff: float = 0.05,
) -> Network:
    """Build a co-expression network using spatially smoothed expression.

    Smooths the expression matrix by replacing each spot's expression with a
    weighted average of its spatial neighbors (Gaussian kernel). The smoothed
    matrix is then used for standard co-expression network construction.

    Parameters
    ----------
    adata:
        AnnData with ``obsm['spatial']`` coordinates.
    layer:
        Expression layer to use (``adata.X`` if None).
    n_neighbors:
        KNN neighbours. Mutually exclusive with *radius*.
    radius:
        Radius cutoff. Mutually exclusive with *n_neighbors*.
    bandwidth:
        Gaussian kernel bandwidth. Auto-estimated from median neighbor
        distance if None.
    cor_cutoff:
        Minimum absolute correlation for an edge.
    p_cutoff:
        Maximum p-value for an edge.

    Returns
    -------
    Network built on spatially smoothed expression values.
    """
    _require_anndata()
    import gcen2
    import scipy.sparse as sp

    coords = _get_coords(adata)
    X = _get_expression(adata, layer)          # (n_spots, n_genes)
    gene_names = list(adata.var_names)

    # Build spatial weight matrix
    W = build_neighbor_graph(coords, n_neighbors=n_neighbors, radius=radius)

    # Override bandwidth if provided by reconstructing weights
    if bandwidth is not None:
        W = W.copy().astype(float)
        # Re-weight using provided bandwidth
        W_coo = W.tocoo()
        if len(W_coo.data) > 0:
            dists = np.sqrt(
                ((coords[W_coo.row] - coords[W_coo.col]) ** 2).sum(axis=1)
            )
            W_coo.data[:] = np.exp(-(dists ** 2) / (2 * bandwidth ** 2))
            W = W_coo.tocsr()

    # Row-normalise: each spot's weight vector sums to 1 (include self)
    n = X.shape[0]
    # Add identity so each spot includes itself in the weighted average
    I_sp = sp.eye(n, format="csr")
    W_self = W + I_sp
    row_sums = np.asarray(W_self.sum(axis=1)).ravel()
    row_sums = np.where(row_sums == 0, 1.0, row_sums)
    D_inv = sp.diags(1.0 / row_sums)
    W_norm = D_inv @ W_self                    # row-normalised

    # Smooth expression
    X_smooth = W_norm @ X                      # (n_spots, n_genes)

    # Build standard co-expression network on smoothed values
    expr_df = pd.DataFrame(X_smooth.T, index=gene_names, columns=list(adata.obs_names))
    return gcen2.network.build(expr_df, cor_cutoff=cor_cutoff, p_cutoff=p_cutoff)
```

### Step 4: 运行测试确认通过

```bash
python3 -m pytest tests/test_spatial.py::TestBuildSpatialWeighted -v --tb=short 2>&1 | tail -10
```
期望：3 passed

### Step 5: Commit

```bash
git add gcen2/spatial/core.py tests/test_spatial.py
git commit -m "feat(spatial): add build_spatial_weighted (Gaussian smoothing)"
```

---

## Task 4：分区域建网 `build_spatial_regions`

**Files:**
- Modify: `gcen2/spatial/core.py`（追加函数）
- Modify: `tests/test_spatial.py`（追加测试类）

### Step 1: 追加测试到 `tests/test_spatial.py`

```python
# ── Task 4 tests ──────────────────────────────────────────────────────────────

@requires_single_cell
class TestBuildSpatialRegions:
    def test_manual_region_key(self, mock_adata):
        from gcen2.spatial.core import build_spatial_regions
        from gcen2.models import Network
        result = build_spatial_regions(
            mock_adata, region_key="region", n_neighbors=4,
            cor_cutoff=0.0, p_cutoff=1.0,
        )
        assert isinstance(result, dict)
        assert set(result.keys()) == {"A", "B"}
        for net in result.values():
            assert isinstance(net, Network)

    def test_region_key_missing_raises(self, mock_adata):
        from gcen2.spatial.core import build_spatial_regions
        from gcen2.exceptions import GCENInputError
        with pytest.raises(GCENInputError, match="not found in adata.obs"):
            build_spatial_regions(mock_adata, region_key="nonexistent", n_neighbors=4)

    def test_small_region_skipped(self, mock_adata):
        """Regions with fewer than 3 spots should be skipped with a warning."""
        import pandas as pd
        mock_adata.obs["tiny"] = ["X"] * 19 + ["Y"]   # Y has only 1 spot
        from gcen2.spatial.core import build_spatial_regions
        with pytest.warns(UserWarning, match="only"):
            result = build_spatial_regions(
                mock_adata, region_key="tiny", n_neighbors=4,
                cor_cutoff=0.0, p_cutoff=1.0,
            )
        assert "Y" not in result

    def test_auto_regions(self, mock_adata):
        """Auto Leiden: result should have n_regions keys."""
        pytest.importorskip("sklearn")
        from gcen2.spatial.core import build_spatial_regions
        result = build_spatial_regions(
            mock_adata, region_key=None, n_regions=3, n_neighbors=4,
            cor_cutoff=0.0, p_cutoff=1.0,
        )
        assert isinstance(result, dict)
        assert len(result) >= 1   # at least some regions
```

### Step 2: 确认测试失败

```bash
python3 -m pytest tests/test_spatial.py::TestBuildSpatialRegions -v --tb=short 2>&1 | tail -10
```
期望：`ImportError: cannot import name 'build_spatial_regions'`

### Step 3: 追加 `build_spatial_regions` 到 `gcen2/spatial/core.py`

```python
def build_spatial_regions(
    adata: Any,
    *,
    region_key: str | None = None,
    n_regions: int = 5,
    layer: str | None = None,
    n_neighbors: int | None = 6,
    radius: float | None = None,
    **build_kwargs,
) -> dict[str, Network]:
    """Build one co-expression network per spatial region.

    Partitions spots into regions (manual annotation or automatic Leiden
    clustering on the spatial neighbor graph), then builds a separate
    co-expression network for each region.

    Parameters
    ----------
    adata:
        AnnData with ``obsm['spatial']`` coordinates.
    region_key:
        Column in ``adata.obs`` to use as region labels. If None, performs
        automatic clustering using KMeans into *n_regions* clusters.
    n_regions:
        Number of regions for automatic clustering (used only when
        *region_key* is None).
    layer:
        Expression layer to use (``adata.X`` if None).
    n_neighbors:
        KNN neighbours for spatial graph. Mutually exclusive with *radius*.
    radius:
        Radius cutoff. Mutually exclusive with *n_neighbors*.

    Returns
    -------
    dict mapping region name → Network. Regions with fewer than 3 spots
    are skipped with a warning.
    """
    _require_anndata()
    import gcen2
    import scipy.sparse as sp

    _get_coords(adata)  # validates obsm['spatial'] exists

    # Determine region labels
    if region_key is not None:
        if region_key not in adata.obs.columns:
            raise GCENInputError(
                f"region_key '{region_key}' not found in adata.obs. "
                f"Available columns: {list(adata.obs.columns)}"
            )
        labels = adata.obs[region_key].astype(str).tolist()
    else:
        # Automatic clustering via KMeans (sklearn)
        try:
            from sklearn.cluster import KMeans
        except ImportError as exc:
            raise GCENInputError(
                "scikit-learn is required for automatic region clustering: "
                "pip install scikit-learn"
            ) from exc
        coords = _get_coords(adata)
        km = KMeans(n_clusters=n_regions, random_state=42, n_init=10)
        cluster_ids = km.fit_predict(coords)
        labels = [f"region_{c}" for c in cluster_ids]

    X = _get_expression(adata, layer)   # (n_spots, n_genes)
    gene_names = list(adata.var_names)
    spot_names = list(adata.obs_names)
    unique_regions = sorted(set(labels))

    result: dict[str, Network] = {}
    for region in unique_regions:
        mask = [l == region for l in labels]
        n_spots = sum(mask)
        if n_spots < 3:
            warnings.warn(
                f"Region '{region}' has only {n_spots} spots (< 3), skipping.",
                UserWarning,
                stacklevel=2,
            )
            continue

        X_region = X[mask, :]          # (n_region_spots, n_genes)
        spot_subset = [s for s, m in zip(spot_names, mask) if m]
        expr_df = pd.DataFrame(
            X_region.T, index=gene_names, columns=spot_subset
        )
        net = gcen2.network.build(expr_df, **build_kwargs)
        result[region] = net

    return result
```

### Step 4: 运行测试确认通过

```bash
python3 -m pytest tests/test_spatial.py::TestBuildSpatialRegions -v --tb=short 2>&1 | tail -15
```
期望：4 passed（`test_auto_regions` 若无 sklearn 则 skip）

### Step 5: Commit

```bash
git add gcen2/spatial/core.py tests/test_spatial.py
git commit -m "feat(spatial): add build_spatial_regions (manual + auto KMeans)"
```

---

## Task 5：模块集成 `__init__.py` + `gcen2` 包注册

**Files:**
- Modify: `gcen2/spatial/__init__.py`（填充公开 API）
- Modify: `gcen2/__init__.py`（注册 spatial 模块）
- Modify: `tests/test_spatial.py`（追加集成测试）

### Step 1: 追加测试到 `tests/test_spatial.py`

```python
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
```

### Step 2: 确认测试失败

```bash
python3 -m pytest tests/test_spatial.py::TestSpatialPublicAPI -v --tb=short 2>&1 | tail -10
```

### Step 3: 填充 `gcen2/spatial/__init__.py`

```python
"""Spatial transcriptomics co-expression network analysis."""

from gcen2.spatial.core import (
    build_spatial_filtered,
    build_spatial_weighted,
    build_spatial_regions,
)

__all__ = [
    "build_spatial_filtered",
    "build_spatial_weighted",
    "build_spatial_regions",
]
```

### Step 4: 修改 `gcen2/__init__.py`，在 `from gcen2 import ...` 行加入 `spatial`

找到第 17 行：
```python
from gcen2 import preprocessing, network, module, annotation, rwr, differential, plotting, report
```
改为：
```python
from gcen2 import preprocessing, network, module, annotation, rwr, differential, plotting, report, spatial
```

并在 `__all__` 列表中加入 `"spatial"`（在 `"differential"` 之后）。

### Step 5: 运行测试

```bash
python3 -m pytest tests/test_spatial.py -v --tb=short 2>&1 | tail -15
```
期望：全部 passed 或 skipped（无 failed）

### Step 6: Commit

```bash
git add gcen2/spatial/__init__.py gcen2/__init__.py tests/test_spatial.py
git commit -m "feat(spatial): expose public API, register gcen2.spatial module"
```

---

## Task 6：CLI `gcen2 spatial-build` + 一致性检查器更新

**Files:**
- Modify: `gcen2/cli.py`（追加 `cmd_spatial_build`）
- Modify: `tests/test_cli.py`（追加 `TestSpatialBuildCommand`）
- Modify: `tools/check_doc_consistency.py`（更新 C2 期望命令集和 C3 文件列表）

### Step 1: 追加测试到 `tests/test_cli.py`

```python
@requires_typer
class TestSpatialBuildCommand:
    def test_filtered_mode(self, runner, tmp_path):
        """spatial-build --mode filtered should produce a .network file."""
        pytest.importorskip("anndata")
        import anndata, numpy as np, pandas as pd

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
        import anndata, numpy as np, pandas as pd

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
```

### Step 2: 确认测试失败

```bash
python3 -m pytest tests/test_cli.py::TestSpatialBuildCommand -v --tb=short 2>&1 | tail -10
```

### Step 3: 追加 `cmd_spatial_build` 到 `gcen2/cli.py`

在 `@app.command("report")` 定义**之前**插入：

```python
@app.command("spatial-build")
def cmd_spatial_build(
    input: Path = typer.Option(..., "-i", "--input",
                               help="Input AnnData file (.h5ad) with obsm['spatial']."),
    output: Path = typer.Option(..., "-o", "--output",
                               help="Output .network file (filtered/weighted) or directory (regions)."),
    mode: str = typer.Option("filtered", "-m", "--mode",
                             help="Analysis mode: filtered | weighted | regions."),
    n_neighbors: int | None = typer.Option(6, "--n-neighbors",
                                           help="KNN neighbours. Use --radius to switch to radius mode."),
    radius: float | None = typer.Option(None, "--radius",
                                        help="Radius cutoff (overrides --n-neighbors)."),
    bandwidth: float | None = typer.Option(None, "--bandwidth",
                                           help="Gaussian bandwidth for weighted mode (auto if omitted)."),
    region_key: str | None = typer.Option(None, "--region-key",
                                          help="obs column for region labels (regions mode)."),
    n_regions: int = typer.Option(5, "--n-regions",
                                  help="Number of auto clusters when --region-key is omitted."),
    cor_cutoff: float = typer.Option(0.6, "-c", "--cor-cutoff", help="Minimum absolute correlation."),
    p_cutoff: float = typer.Option(0.05, "-p", "--p-cutoff", help="Maximum p-value."),
) -> None:
    """Spatial transcriptomics co-expression network analysis."""
    try:
        import anndata as ad
    except ImportError:
        typer.echo("ERROR: anndata is required: pip install gcen2[single-cell]", err=True)
        raise typer.Exit(1)

    import gcen2

    # Resolve neighbor args: radius overrides n_neighbors
    nb_kwargs: dict = {}
    if radius is not None:
        nb_kwargs["n_neighbors"] = None
        nb_kwargs["radius"] = radius
    else:
        nb_kwargs["n_neighbors"] = n_neighbors

    adata = ad.read_h5ad(input)
    typer.echo(f"Loaded {adata.n_obs} spots × {adata.n_vars} genes from {input}")

    if mode == "filtered":
        net = gcen2.spatial.build_spatial_filtered(
            adata, cor_cutoff=cor_cutoff, p_cutoff=p_cutoff, **nb_kwargs
        )
        net.to_edgelist(output)
        typer.echo(f"Filtered network: {net.node_count} nodes, {net.edge_count} edges → {output}")

    elif mode == "weighted":
        kw = dict(cor_cutoff=cor_cutoff, p_cutoff=p_cutoff, **nb_kwargs)
        if bandwidth is not None:
            kw["bandwidth"] = bandwidth
        net = gcen2.spatial.build_spatial_weighted(adata, **kw)
        net.to_edgelist(output)
        typer.echo(f"Weighted network: {net.node_count} nodes, {net.edge_count} edges → {output}")

    elif mode == "regions":
        output.mkdir(parents=True, exist_ok=True)
        kw: dict = dict(cor_cutoff=cor_cutoff, p_cutoff=p_cutoff, **nb_kwargs)
        if region_key is not None:
            kw["region_key"] = region_key
        else:
            kw["region_key"] = None
            kw["n_regions"] = n_regions
        nets = gcen2.spatial.build_spatial_regions(adata, **kw)
        for region_name, net in nets.items():
            out_file = output / f"{region_name}.network"
            net.to_edgelist(out_file)
        typer.echo(f"Regions network: {len(nets)} regions → {output}/")
        for rname, net in sorted(nets.items()):
            typer.echo(f"  {rname}: {net.node_count} nodes, {net.edge_count} edges")

    else:
        typer.echo(f"ERROR: Unknown mode '{mode}'. Choose: filtered | weighted | regions", err=True)
        raise typer.Exit(1)
```

### Step 4: 运行 CLI 测试

```bash
python3 -m pytest tests/test_cli.py::TestSpatialBuildCommand -v --tb=short 2>&1 | tail -15
```
期望：2 passed（或因无 anndata/typer 而 skip）

### Step 5: 更新一致性检查器

在 `tools/check_doc_consistency.py` 中：

**C2** — `doc_commands` 集合中加入 `"spatial-build"`：
```python
doc_commands = {
    "normalize", "filter", "build", "stats",
    "module", "annotate", "rwr", "enrich", "report",
    "merge", "extract", "shuffle", "diff", "spatial-build",
}
```

**C1** — `check_c1_api()` 中加入 spatial 模块检查：
```python
("gcen2/spatial/core.py", [
    "build_spatial_filtered", "build_spatial_weighted", "build_spatial_regions",
]),
("gcen2/spatial/neighbors.py", ["build_neighbor_graph"]),
```

**C3** — `check_c3_files()` 中加入：
```python
"gcen2/spatial/__init__.py",
"gcen2/spatial/neighbors.py",
"gcen2/spatial/core.py",
"tests/test_spatial.py",
```

### Step 6: 运行全量验证

```bash
python3 tools/check_doc_consistency.py 2>&1 | tail -5
python3 -m pytest tests/ --ignore=tests/test_integration.py -q --tb=short 2>&1 | tail -5
```
期望：`差异: 0 项` + `0 failed`

### Step 7: Commit

```bash
git add gcen2/cli.py tests/test_cli.py tools/check_doc_consistency.py
git commit -m "feat(spatial): add CLI spatial-build command, update consistency checker"
```

---

## 完成标准

1. `python3 -m pytest tests/test_spatial.py -v` — 全部 PASSED 或 SKIPPED（0 FAILED）
2. `python3 -m pytest tests/test_cli.py::TestSpatialBuildCommand -v` — PASSED 或 SKIPPED
3. `python3 tools/check_doc_consistency.py` — `差异: 0 项`
4. `python3 -m pytest tests/ --ignore=tests/test_integration.py -q` — 0 failed, 0 errors
5. `gcen2 spatial-build --help` 输出包含 `filtered | weighted | regions`
