# Case study: what GCEN2 finds in GTEx whole blood

The other benchmarks show that GCEN2 gives the reference tools' results. This
one asks what the results say about the biology, on the GTEx v10 whole-blood
data of `../wgcna_gtex` (803 donors, 11,944 genes, top 5 expression PCs
regressed out):

1. **Modules.** `gcen2 module -e` (the same 23 modules as WGCNA's
   blockwiseModules, `../wgcna_gtex`), GO enrichment, hub genes and the
   correlation of each module eigengene with the donor's sex, age and death
   circumstances (Hardy scale), checked in two independent halves of the
   donors.
2. **Gene function prediction.** Random walk with restart on the
   co-expression network (`gcen2 build`, `gcen2 rwr`), cross-validated on
   GO biological processes against controls for node degree.

## Data

Besides `../wgcna_gtex` (work/expr_11944_pc5.tsv, work/half_A.txt,
half_B.txt) and the GTEx subject phenotypes in `data/`
(`GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt`, see `../dgca_gtex`),
the Gene Ontology release 2026-08-05 in `data/go/`:

| file | URL | SHA-256 |
|---|---|---|
| `go-basic.obo` | `https://release.geneontology.org/2026-08-05/ontology/go-basic.obo` | `b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc` |
| `HUMAN-uniprot.gaf.gz` | `https://release.geneontology.org/2026-08-05/annotations/gaf/HUMAN-uniprot.gaf.gz` | `a0afba19dfb1f8fa996bc1bdcd61fd0c9bd4cf0d2bf2509d09ac86993c0e70a2` |

Annotations with a NOT qualifier, inferred from expression (IEP, HEP) or
without data (ND) are left out, so that co-expression is not scored against
annotations derived from co-expression. 10,140 of the 11,944 genes are
annotated.

## Run

```bash
cd <gtex>
python -I .../01_prepare.py     # work/case/: GO association, donor traits
bash .../02_run.sh              # runs/case/: modules, network, shuffled network
python -I .../03_modules.py     # results/module_*.tsv, trait_halves.tsv
python -I .../04_rwr_cv.py      # results/rwr_cv*.tsv (about 4 minutes)
python -I .../05_figures.py     # results/fig_*
```

## Results

### 1. Modules are the cell types and programmes of blood

The 23 modules (4,431 genes; `results/module_summary.tsv`, all terms with
FDR <= 0.05 in `results/module_go.tsv`) are, by their hub genes and top GO
biological processes (fold enrichment, FDR: the `fdr` column of
`gcen2.annotation.enrich`, Benjamini-Hochberg over every GO term annotated
among the 11,944 genes):

| module | genes | hub genes | what it is (top GO process) |
|---|---|---|---|
| ME2 | 438 | ANK1, EIF2AK1, MARCHF8, XPO7 | erythroid cells (heme B biosynthesis, 22x, 5e-7) |
| ME7 | 166 | MPIG6B, ITGB3, NRGN, SPARC | platelets (cell adhesion, 5.2x, 3e-6) |
| ME11 | 140 | CEBPE, CEACAM8, PGLYRP1, CD24 | immature neutrophils, granule proteins (killing of cells of another organism, 24x, 2e-7) |
| ME3 | 414 | CD68, CSF1R, NAGA, ADAP2 | monocytes (endocytosis, 3.6x, 2e-7) |
| ME8 | 162 | TBX21, KLRD1, PYHIN1, ADGRG1 | NK and cytotoxic T cells (immune system process, 3.4x, 2e-14) |
| ME18 | 68 | LEF1, TCF7, CCR7, TRABD2A | naive T cells (regulation of T cell activation, 7.6x, 2e-4) |
| ME19 | 66 | TRAC, CD5, TRBC2, CD6 | T cells (T cell receptor signalling, 20x, 1e-8) |
| ME16 | 85 | CD19, MS4A1, CD22, BANK1 | B cells (B cell activation, 18x, 4e-9) |
| ME6 | 214 | IGKC, MZB1, IGLC1 | plasma cells (immunoglobulin-mediated immune response, 37x, 2e-60) |
| ME13 | 116 | IFIT3, PARP9, EPSTI1, DDX58 | type I interferon response (response to virus, 17x, 3e-35) |
| ME10 | 143 | TOP2A, FOXM1, KIF11 | proliferating cells (mitotic cell cycle, 11x, 6e-40) |
| ME14 | 107 | HSPH1, HSP90AA1, FKBP4 | heat-shock response (response to unfolded protein, 47x, 4e-20) |
| ME9, ME12, ME0, ME21 | | RPL*, MRPL*, COX7C, MT-* | cytoplasmic and mitochondrial translation, respiration |
| ME20 | 51 | CRYAB, DES, TNNC1 | skeletal-muscle genes (muscle contraction, 26x, 8e-14): not a blood cell type, probably traces of muscle tissue in some samples |

The naive T-cell module (`results/fig_hub_ME18`) has the canonical naive and
central-memory genes among its 25 top-kME genes: LEF1, TCF7, CCR7, IL7R, CD27,
CD28, BCL11B. ME15 (pseudogenes of ribosomal proteins) and ME22 (histone
genes) are technical or annotation artefacts.

**Age.** Five modules correlate with the donor's age, each in both
independent halves of the donors (`results/trait_halves.tsv`, p < 0.05 in
each; all 803 donors in `results/module_summary.tsv`):

| module | r (all) | p (all) | r, half A / B |
|---|---|---|---|
| ME18 naive T cells | −0.17 | 7e-7 | −0.16 / −0.19 |
| ME16 B cells | −0.16 | 4e-6 | −0.15 / −0.17 |
| ME3 monocytes | +0.14 | 5e-5 | +0.14 / +0.15 |
| ME12 mitochondrial translation | +0.13 | 2e-4 | +0.13 / +0.13 |
| ME8 NK / cytotoxic T cells | +0.13 | 3e-4 | +0.16 / +0.10 |

These are the known shifts of the ageing immune system: fewer naive T cells
(thymic involution) and B cells, relatively more myeloid cells and more
cytotoxic effector cells. Women have higher B-cell (ME16, r = +0.12, p = 9e-4)
and T-cell module eigengenes (ME18, ME19; r = +0.10, +0.11), in line with the
known sex differences in adaptive immunity; these are positive in both halves
but reach p < 0.05 in only one of them, so they are weaker evidence than the
age associations.

**Death circumstances.** No module associates consistently with the Hardy
scale or death on a ventilator: the 5 regressed PCs carried that variation
(PC2 of the uncorrected data, 26% of the variance, r = −0.72 with death on a
ventilator; PC1, 36%, r = +0.35; `results/pc_traits.tsv`). The correction that
removes the confounding (`../wgcna_gtex`) also removes this biology; age,
less correlated with the PCs (|r| <= 0.28), survives it.

### 2. RWR predicts gene function beyond node degree

2,144 GO biological processes with 20–200 genes in the network (11,236
genes, 533,232 edges; 8,992 genes with a BP annotation), 3-fold
cross-validation (`results/rwr_cv.tsv`, `results/rwr_cv_summary.tsv`,
`results/fig_rwr_cv`):

| score | median AUROC | terms with AUROC >= 0.7 | terms above the shuffled network |
|---|---|---|---|
| RWR (`gcen2 rwr`) | 0.607 | 322 | 2,016 |
| RWR, \|r\| weights | 0.608 | 323 | 2,019 |
| neighbour voting | 0.578 | 229 | 1,789 |
| RWR on a degree-preserving shuffle | 0.515 | 31 | — |
| degree only | 0.518 | 38 | 1,215 |

RWR beats the shuffled network on 94% of the terms and the degree alone on
95% (median gain in AUROC over the shuffle 0.087; Wilcoxon signed-rank p
below double precision for both), and neighbour voting
(p = 6e-137): the walk uses who is connected to whom, not only which genes
are hubs (Gillis & Pavlidis, PLoS One 2011). The best-predicted terms are
blood functions — cellular response to type I interferon (AUROC 0.95),
immunoglobulin-mediated immune response (0.94), B cell-mediated immunity
(0.90), T cell selection (0.89) — and core cellular machinery (translation,
oxidative phosphorylation, the TCA cycle); terms of processes not active in
blood are near chance (median of terms named neuro- or synaptic: 0.55–0.56).
For a few terms the degree alone already does well (e.g. "negative regulation
of syncytium formation", 0.91 by degree), because their genes are the
highly connected ribosomal proteins; the controls flag such cases.
Weighting the walk by |r| changes little (median AUROC +0.001, p = 7e-6),
because every edge has r >= 0.4.

The cross-validation (19,296 walks, one at a time) takes about 4 minutes.

## Files

| file | contents |
|---|---|
| `01_prepare.py` | GO association (symbols), donor traits, background genes |
| `02_run.sh` | `gcen2 module -e` with eigengenes, kME, hubs, traits; `gcen2 build`; `gcen2 shuffle` |
| `03_modules.py` | GO enrichment per module, trait correlations in each half, PCs vs traits |
| `04_rwr_cv.py` | cross-validated RWR against neighbour voting, a shuffled network and degree |
| `05_figures.py` | figure panels with `gcen2.plot` |
| `results/` | the tables and figures above |

