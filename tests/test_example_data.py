"""data/example: present, as described in its README, and readable (GO and
association files gzip-compressed)."""

import gzip
from pathlib import Path

import pytest

EX = Path(__file__).resolve().parents[1] / "data" / "example"


def test_gzip_text_files_read_as_plain(tmp_path):
    from gcen2.annotation import load_association, load_go
    from gcen2.io import read_gene_list

    data = Path(__file__).parent / "data" / "v1_annotate"
    for name in ("go.obo", "gene_go.assoc"):
        with open(data / name, "rb") as fi, gzip.open(tmp_path / f"{name}.gz", "wb") as fo:
            fo.write(fi.read())
    go, go_gz = load_go(data / "go.obo"), load_go(tmp_path / "go.obo.gz")
    assert {k: (t.name, t.namespace, t.parents) for k, t in go.terms.items()} == \
        {k: (t.name, t.namespace, t.parents) for k, t in go_gz.terms.items()}
    assert load_association(data / "gene_go.assoc", go) == \
        load_association(tmp_path / "gene_go.assoc.gz", go_gz)
    # any name: recognised by the gzip header
    (tmp_path / "genes.txt").write_bytes(gzip.compress(b"g1\ng2\n"))
    assert read_gene_list(tmp_path / "genes.txt") == ["g1", "g2"]


def test_bulk_example():
    from gcen2.annotation import load_association, load_go
    from gcen2.io import read_tsv

    d = EX / "bulk"
    expr = read_tsv(d / "gene_expr.tsv")
    assert (len(expr.gene_names), len(expr.sample_names)) == (32940, 15)
    go = load_go(d / "go-basic.obo.gz")
    assert len(go.terms) > 40000
    assoc = load_association(d / "gene_go.assoc.gz", go)
    assert len(set(assoc) & set(expr.gene_names)) > 10000
    for f in ("K2ko.tsv", "gene_kegg.assoc", "rwr_seed_genes.list", "gene_list.txt",
              "enrichment_gene.list", "background_gene.list"):
        assert (d / f).stat().st_size > 0, f


def test_single_cell_and_spatial_examples():
    ad = pytest.importorskip("anndata")
    import numpy as np

    sc = ad.read_h5ad(EX / "single_cell" / "pbmc_subset.h5ad")
    assert sc.shape == (1800, 6873)
    assert sc.obs["cell_type"].value_counts().to_dict() == {"CD14 Mono": 800, "B": 600,
                                                            "NK": 400}
    assert {"X_pca", "X_umap"} <= set(sc.obsm)
    assert np.issubdtype(sc.X.dtype, np.integer)  # raw counts
    sp = ad.read_h5ad(EX / "spatial" / "mouse_brain_anterior_subset.h5ad")
    assert sp.shape == (2695, 2000)
    assert {"array_row", "array_col", "region"} <= set(sp.obs.columns)
    assert sp.obsm["spatial"].shape == (2695, 2)
