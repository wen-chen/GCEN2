# GCEN2 example data

Small data sets for trying every part of GCEN2 (and for the README and the
tutorials). About 27 MB in all.

## `bulk/`: GCEN v1 sample data

From the GCEN v1 0.6.3 release package
(https://github.com/wen-chen/GCEN/releases/tag/v0.6.3): 32,940 genes x 15
samples, with GO and KEGG annotations.

| file | content |
|---|---|
| `gene_expr.tsv` | expression matrix, genes x samples (GCEN v1 format: no header line) |
| `go-basic.obo.gz` | Gene Ontology (data-version releases/2018-10-07), gzip-compressed; GCEN2 reads it as it is |
| `gene_go.assoc.gz` | gene -> GO terms (gzip-compressed) |
| `K2ko.tsv`, `gene_kegg.assoc` | KEGG: K number -> pathway, gene -> K numbers |
| `rwr_seed_genes.list` | seed genes for `gcen2 rwr` |
| `gene_list.txt` | genes for `gcen2 extract` |
| `enrichment_gene.list`, `background_gene.list` | query and background for `gcen2 enrich` |

```bash
cd data/example/bulk
gcen2 normalize -i gene_expr.tsv -o norm.tsv -m tmm
gcen2 filter    -i norm.tsv -o filt.tsv -p 0.75
gcen2 build     -i filt.tsv -o net.network -m spearman -c 0.8 -p 0.001 -f -t 8
gcen2 module    -i net.network -o module.txt -s 0.5 -t 8
gcen2 annotate  -g go-basic.obo.gz -a gene_go.assoc.gz -m module.txt -o annot/ -t 8
gcen2 annotate  -k K2ko.tsv -a gene_kegg.assoc -m module.txt -o annot_kegg/ -t 8
gcen2 rwr       -n net.network -g rwr_seed_genes.list -o rwr.tsv
gcen2 enrich    -e enrichment_gene.list -b background_gene.list -g go-basic.obo.gz \
                -a gene_go.assoc.gz -o enrich.tsv
gcen2 module    -e filt.tsv -o wgcna_modules.tsv --signed -l log2 --me me.tsv -t 8
```

The last command (WGCNA modules of all ~24,000 genes) takes about 2 minutes
and 5 GB of memory; the others take seconds.

## `single_cell/pbmc_subset.h5ad`: 10x PBMC 10k

1,800 cells of three cell types (CD14 Mono 800, B 600, NK 400, randomly
chosen, seed 0) x the 6,873 genes detected in at least 10 % of them; raw
counts (integers) in `X`; `obs`: `cell_type`, `cluster`, `total_counts`;
`obsm`: `X_pca` (of all 10,685 cells) and `X_umap`.

Source: 10x Genomics PBMC 10k v3 (`pbmc_10k_v3`, filtered feature-barcode
matrix), licensed under CC BY 4.0. Cell types, PCA and the clustering are
those of the GCEN2 hdWGCNA benchmark (`benchmarks/hdwgcna_pbmc/01_preprocess.py`).

```bash
gcen2 celltype-modules -i pbmc_subset.h5ad --groupby cell_type -o modules/ -t 8 --write-tom
gcen2 plot metacells -i modules/ -c pbmc_subset.h5ad -o metacells
gcen2 plot me-embedding -i modules/B/cell_eigengenes.tsv -c pbmc_subset.h5ad -o b_eigengenes
gcen2 report --celltype-modules modules/ --coords pbmc_subset.h5ad -o report.html
```

## `spatial/mouse_brain_anterior_subset.h5ad`: 10x Visium mouse brain

All 2,695 spots of the anterior sagittal section x the 2,000 genes detected in
the most spots; raw counts in `X`; `obs`: `array_row`, `array_col` (Visium
grid) and `region` (k-means, 12 clusters, seed 0, of 30 PCs of the
log-normalised highly variable genes); `obsm["spatial"]`: spot coordinates.

Source: 10x Genomics `V1_Mouse_Brain_Sagittal_Anterior` (Space Ranger 1.1.0),
licensed under CC BY 4.0. Prepared as in the GCEN2 hdWGCNA spatial benchmark
(`benchmarks/hdwgcna_spatial/01_prepare.py`).

```bash
gcen2 spatial-modules -i mouse_brain_anterior_subset.h5ad --groupby region -o spatial_modules/ -t 8 --plot
gcen2 spatial-build   -i mouse_brain_anterior_subset.h5ad -m weighted -o smoothed.network \
                      -c 0.6 -t 8 --moran moran.tsv
gcen2 plot moran -i moran.tsv -o moran
```

## Licences

GCEN v1 sample data: GPL-3.0 (GCEN v1). Gene Ontology: CC BY 4.0
(http://geneontology.org/docs/go-citation-policy/). KEGG pathway data: KEGG is
free for academic use only (https://www.kegg.jp/kegg/legal.html); the KEGG
files are those distributed with GCEN v1. 10x Genomics data sets: CC BY 4.0.
