# GCEN2 multi-layer networks vs RandomWalkRestartMH and R

Two ways to combine gene networks of different kinds:

| GCEN2 | compared with | what it gives |
|---|---|---|
| `gcen2 rwr-multiplex` | RandomWalkRestartMH 1.13.2 (Valdeolivas et al., Bioinformatics 2019: RWR-M, the authors' implementation) | genes ranked by a random walk with restart over all layers at once |
| `gcen2 fuse` | an R implementation of the same definition (`03_fuse.R`, data.table); there is no published package for it | one network, every edge scored by the mean of its normalised weights over the layers |

**Data.** Three layers over 18,082 human genes:

| layer | source | edges |
|---|---|---|
| co-expression | GTEx whole blood, 11,944 genes, 5 PCs removed (`../wgcna_gtex`); `gcen2 build -m pearson -c 0.5 -p 0.001 -f -s`; weight r | 142,757 |
| protein interaction | STRING v12 human, combined score >= 700; weight score / 1000 | 236,930 |
| regulation | TRRUST v2 human, transcription factor – target pairs | 8,300 |

Seeds: nine interferon-response genes (`IFIT1 IFIT3 ISG15 MX1 OAS1 RSAD2
IFI44L STAT1 IRF7`).

**Settings.** RandomWalkRestartMH's defaults: restart 0.7, delta 0.5,
unweighted layers, geometric mean over the layers; GCEN2 `--gamma 0.7
--delta 0.5` (`--transition supra` is the default). Fusion: weights divided
by each layer's largest, equal layer weights.

## Run

Data in `data/`: TRRUST v2,
`https://www.grnpedia.org/trrust/data/trrust_rawdata.human.tsv`
(sha256 `9b909319ccc8e36588b5a1bd3640e0dfba03936c6fc154cc4304d742bec11bed`);
STRING v12 as in [`../hdwgcna_pbmc`](../hdwgcna_pbmc/README.md).

```bash
cd <data directory>
python .../01_prepare.py <GTEx>/work/expr_11944_pc5.tsv <STRING directory> \
    data/trrust_rawdata.human.tsv          # work/{coexpression,ppi,regulation}.tsv
bash .../run_all.sh                        # both tools, /usr/bin/time
python .../04_compare.py                   # results/compare.tsv
```

## Results

`results/compare.tsv` (times: mean of two runs):

| | GCEN2 | other | agreement |
|---|---|---|---|
| RWR-M | **0.49 s, 0.18 GB** | RandomWalkRestartMH: 28.3 s, 0.37 GB | all 18,073 non-seed genes, scores within 1.5e-12, the same top 100 |
| fusion | **0.42 s, 0.19 GB** | R data.table: 0.57 s, 0.15 GB | the same 376,750 edges, scores within 6e-16, the same layer counts |

RandomWalkRestartMH builds the supra-adjacency matrix as a sparse `Matrix`
and iterates in R; GCEN2 iterates in C++. The fusion is a join and a sum, which
data.table does about as fast as GCEN2.

## Notes

- **Transition matrix.** RandomWalkRestartMH puts the layers, scaled by
  1 − delta, and delta / (L − 1) between the copies of a gene into one
  supra-adjacency matrix and normalises its columns. A gene of degree d then
  jumps to another layer with probability delta / ((1 − delta) d + delta),
  so hubs mostly stay in their layer. Before this benchmark GCEN2 normalised
  each layer first (every gene jumps with probability delta): the same top 10
  shared 1 gene. GCEN2 now follows RandomWalkRestartMH by default and keeps
  the per-layer version as `--transition layer`.
- **Weights.** RandomWalkRestartMH rescales the weights of a weighted layer
  to `(w − min)/(max − min) · (1 − a) + a` with `a = min/max`, which is
  w / max; GCEN2 `--weighted` does the same (checked against
  RandomWalkRestartMH in `tests/test_multilayer.py`, weighted and not).
- **Gene names.** RandomWalkRestartMH removes the first `_1` from the
  reported names (`gsub("_1", "", ...)`): `hCG_1984214` becomes `hCG984214`.
  `04_compare.py` matches those two genes after the same substitution.
