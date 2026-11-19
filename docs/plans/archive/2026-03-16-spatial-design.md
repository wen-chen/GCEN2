# Spatial Co-Expression Network Analysis — Design Document

**Date:** 2026-03-16
**Status:** Approved, ready for implementation

---

## Goal

Add spatial transcriptomics co-expression network analysis to GCEN2 as an independent module `gcen2/spatial/`, supporting three complementary methods and full CLI integration.

## Scope

**In scope (v1):**
- AnnData + `obsm['spatial']` 2D coordinate input (Visium, MERSCOPE, Xenium)
- Three analysis modes: spatial filtering, spatial-weighted correlation, region-based networks
- KNN and radius-based spatial neighbor graphs (both supported, user-selected)
- CLI command `gcen2 spatial-build` with `--mode` switch
- Output: `Network` object (single network) or `dict[str, Network]` (regions)

**Out of scope (future):**
- H&E image registration
- 3D spatial coordinates
- Spatial visualization (separate plotting module concern)
- Trajectory-based spatial dynamics

---

## Architecture

### Module Structure

```
gcen2/spatial/
├── __init__.py      # Public API: build_spatial_filtered, build_spatial_weighted, build_spatial_regions
├── neighbors.py     # build_neighbor_graph() — KNN + radius, returns scipy sparse weight matrix
└── core.py          # Three build functions, call neighbors.py + gcen2.network.build
```

**Key decisions:**
- Independent module (not under `single_cell/`) — spatial is a distinct analysis paradigm
- No new dataclass — returns existing `Network` or `dict[str, Network]`
- `anndata` / `scipy` / `numpy` only for core; `squidpy` optional (not required)
- Guard with `_require_anndata()` pattern (same as `single_cell/core.py`)

---

## API Design

### `neighbors.py`

```python
def build_neighbor_graph(
    coords: np.ndarray,           # shape (n_spots, 2)
    n_neighbors: int | None = 6,  # KNN — Visium default 6 (honeycomb)
    radius: float | None = None,  # radius cutoff; mutually exclusive with n_neighbors
) -> scipy.sparse.csr_matrix:    # (n_spots, n_spots) Gaussian-weighted adjacency
```

- `n_neighbors` and `radius` are mutually exclusive; both None → `GCENInputError`
- Edge weight = Gaussian kernel: `exp(-d² / (2 * bandwidth²))` where `bandwidth = median(neighbor distances)`
- Returns symmetric CSR matrix

### `core.py` — Three Methods

#### Method 1: `build_spatial_filtered`
```python
def build_spatial_filtered(
    adata,
    *,
    layer: str | None = None,
    n_neighbors: int | None = 6,
    radius: float | None = None,
    cor_cutoff: float = 0.6,
    p_cutoff: float = 0.05,
    **build_kwargs,
) -> Network:
```
**Pipeline:** standard `network.build(expr)` → filter edges to retain only pairs whose expressing spots are spatially adjacent.

#### Method 2: `build_spatial_weighted`
```python
def build_spatial_weighted(
    adata,
    *,
    layer: str | None = None,
    n_neighbors: int | None = 6,
    radius: float | None = None,
    bandwidth: float | None = None,   # None = auto (median neighbor distance)
    cor_cutoff: float = 0.6,
    p_cutoff: float = 0.05,
) -> Network:
```
**Pipeline:** build Gaussian spatial weight vector per spot → compute spatially weighted Pearson correlation matrix (NumPy) → threshold → return `Network`.

#### Method 3: `build_spatial_regions`
```python
def build_spatial_regions(
    adata,
    *,
    region_key: str | None = None,   # obs column; None = auto Leiden clustering
    n_regions: int = 5,              # used only when region_key is None
    layer: str | None = None,
    n_neighbors: int | None = 6,
    radius: float | None = None,
    **build_kwargs,
) -> dict[str, Network]:
```
**Pipeline:** partition spots by `region_key` (or auto Leiden) → pseudobulk per region → call `network.build()` per region → return `{region_name: Network}`.

---

## CLI Design

Single command `gcen2 spatial-build`, `--mode` selects method:

```bash
# Spatial filtering
gcen2 spatial-build -i data.h5ad -o net.network --mode filtered --n-neighbors 6

# Spatial weighted correlation
gcen2 spatial-build -i data.h5ad -o net.network --mode weighted --bandwidth 50.0

# Region-based (manual annotation)
gcen2 spatial-build -i data.h5ad -o net_regions/ --mode regions --region-key cell_type

# Region-based (auto Leiden)
gcen2 spatial-build -i data.h5ad -o net_regions/ --mode regions --n-regions 5
```

**Output:**
- `filtered` / `weighted`: single `.network` TSV file
- `regions`: directory, one `{region_name}.network` file per region

---

## Error Handling

| Scenario | Exception | Message |
|---------|-----------|---------|
| `adata.obsm['spatial']` missing | `GCENInputError` | "No spatial coordinates found in adata.obsm['spatial']" |
| Both `n_neighbors` and `radius` specified | `GCENInputError` | "Specify either n_neighbors or radius, not both" |
| Both `n_neighbors` and `radius` are None | `GCENInputError` | "Must specify either n_neighbors or radius" |
| `region_key` not in `adata.obs` | `GCENInputError` | "region_key '{key}' not found in adata.obs" |
| Region has < 3 spots | `warnings.warn` | "Region '{r}' has only {n} spots, skipping" |
| `anndata` not installed | `GCENInputError` | "anndata is required: pip install gcen2[single-cell]" |

---

## Testing Strategy

**File:** `tests/test_spatial.py`

**Fixture:**
```python
@pytest.fixture
def mock_adata():
    """20 spots, 50 genes, 2D grid (4×5), two regions A/B"""
    # X: random (20, 50), coords: grid, obs.region: A×10 + B×10
```

**Test cases:**

| Test | Validates |
|------|-----------|
| `test_neighbor_graph_knn` | KNN=6 → sparse matrix, 6 non-zeros per row |
| `test_neighbor_graph_radius` | spots > radius → weight = 0 |
| `test_neighbor_graph_both_raises` | both args → `GCENInputError` |
| `test_filtered_returns_network` | returns `Network`, edge_count ≤ unfiltered |
| `test_weighted_returns_network` | returns `Network`, edge_count > 0 |
| `test_regions_returns_dict` | returns `{"A": Network, "B": Network}` |
| `test_regions_auto_leiden` | `region_key=None` → dict with `n_regions` keys |
| `test_missing_spatial_raises` | no `obsm['spatial']` → `GCENInputError` |

All tests guarded with `@pytest.mark.requires_single_cell`.

---

## Dependencies

| Package | Required | Usage |
|---------|----------|-------|
| `anndata` | Yes (optional extra) | AnnData input |
| `numpy` | Yes | Weighted correlation computation |
| `scipy` | Yes | Sparse neighbor graph, KD-tree for KNN |
| `scikit-learn` | Optional | Leiden clustering for auto-regions |
| `squidpy` | No | Not required; spatial viz is out of scope |

No new entries needed in `pyproject.toml` — all core deps already in `single-cell` extra.
