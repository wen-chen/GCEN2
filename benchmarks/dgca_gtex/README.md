# GCEN2 vs DGCA: differential co-expression

Which gene pairs are correlated differently in two conditions? DGCA
(Differential Gene Correlation Analysis, McKenzie et al., BMC Syst Biol 2016)
tests every pair with Fisher's z: the two correlations are z-transformed and
their difference divided by its standard error. `gcen2 diff` applies the same
test to the edges of two networks; built with `-c 0 -p 1`, the networks hold
every gene pair, as DGCA reports.

| DGCA 1.0.3 | GCEN2 |
|---|---|
| `ddcorAll(corrType = "pearson" / "spearman", nPerms = 0, adjust = "BH")` | `gcen2 build -m pearson / spearman -c 0 -p 1` for each condition, then `gcen2 diff --n-a --n-b -m pearson / spearman -t 8` |

Both use the standard error `sqrt(v/(n_a - 3) + v/(n_b - 3))` of the
difference of the z-transformed correlations, `v = 1` for Pearson and
`v = 1.06` for Spearman (Fieller et al. 1957), and adjust the p-values of all
pairs by Benjamini-Hochberg (DGCA `pValDiff_adj`, GCEN2 `fdr`).

**Data.** GTEx v10 whole blood (log2(TPM + 1), the 2,000 and 5,000 most
variable genes of `../wgcna_gtex`), split by sex: 541 male (A) and 262 female
(B) donors.

**Settings.** Pearson and Spearman correlation; 8 threads for both (DGCA through
`WGCNA::allowWGCNAThreads`). DGCA's `corr_cutoff` is set to 0.9999, the
clipping GCEN2 applies before the z transform (DGCA's default is 0.99).
DGCA's z is B − A, GCEN2's A − B; they are compared with the sign flipped.

## Run

Data: `../wgcna_gtex` (its `01_prepare.py` writes `work/expr_5000.tsv`) and
the GTEx subject phenotypes,
`https://storage.googleapis.com/adult-gtex/annotations/v10/metadata-files/GTEx_Analysis_v10_Annotations_SubjectPhenotypesDS.txt`
(sha256 `2313caeaaf3718ac8aead8c6b9eab92b3f613509cfe1a6f7b6286310a7200e73`),
in `data/`.

```bash
cd <GTEx data directory>
python .../01_prepare.py      # work/dc_<N>_{male,female}.tsv
bash .../run_all.sh           # both tools, /usr/bin/time
python .../03_compare.py      # results/compare.tsv
```

## Results

`results/compare.tsv` (times: mean of two runs):

| | Pearson, 2,000 genes | Pearson, 5,000 | Spearman, 2,000 | Spearman, 5,000 |
|---|---|---|---|---|
| gene pairs (both tools) | 1,999,000 | 12,497,500 | 1,999,000 | 12,497,500 |
| max \|Δr\| | 5e-07 | 5e-07 | 5e-07 | 5e-07 |
| max \|Δz\|, DGCA's formula on GCEN2's r | 7e-15 | 9e-15 | 6e-15 | 7e-15 |
| max relative Δp, ΔFDR | 1.2e-03 | 1.7e-03 | 1.1e-03 | 1.1e-03 |
| pairs with p <= 0.05 that differ | 4 | 21 | 2 | 19 |
| pairs with FDR <= 0.05: GCEN2 / DGCA (differing) | 11,777 / 11,777 (0) | 132,270 / 132,270 (0) | 18,556 / 18,556 (0) | 266,175 / 266,173 (4) |
| GCEN2: time, peak memory | 2.1 s, 0.40 GB | 9.0 s, 1.27 GB | 2.2 s, 0.40 GB | 8.6 s, 1.27 GB |
| DGCA: time, peak memory | 4.6 s, 0.96 GB | 23.2 s, 4.03 GB | 258 s, 0.96 GB | 1,582 s, 4.33 GB |

The two tools compute the same test. GCEN2 writes correlations with six
decimals (the GCEN v1 network format), which is the whole difference: on the
same rounded correlations DGCA's formula gives GCEN2's z to 1e-14, and the
few pairs classified differently lie on the 0.05 boundary.

DGCA's Spearman runs used one CPU core throughout (its Pearson runs use the
8 WGCNA threads); GCEN2 ranks each gene once and correlates the ranks with
its Pearson kernel. Of GCEN2's 9 s on 5,000 genes, `gcen2 diff` (reading
the two 12.5M-edge files, the test, BH and writing, 8 threads) takes 5.3-5.9
s and the two `build` runs 1.5-1.7 s each (`results/time.txt`). In practice
the networks are built with a cutoff and are much smaller.

## Differences

- **Clipping.** DGCA caps r at `corr_cutoff` from above only; GCEN2 clips to
  ±0.9999.
- **Beyond the test.** DGCA also classifies each pair by the sign and
  significance of its correlation in each condition (`+/0`, `-/+`, ...) and
  offers permutation-based FDR and module-level tests; GCEN2 classifies the
  edges by network membership (`A_only`, `B_only`, `shared_changed`,
  `shared_stable`).
- **What is tested.** DGCA tests every pair of the input genes; `gcen2 diff`
  the pairs that are edges of both networks, so `-c 0 -p 1` is needed for
  the same set (the FDR is over the tested pairs).
