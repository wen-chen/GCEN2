# Thread scaling and parameter sensitivity

Two questions the comparisons with the reference tools leave open: how the
main commands scale with the number of threads (`-t`), and how much the
results depend on the parameters. Data of the other benchmarks:
GTEx v10 whole blood (`../wgcna_gtex`, `../dgca_gtex`) and the Human Lung
Cell Atlas (`../hdwgcna_hlca`); STRING v12 gold standard of `../hdwgcna_pbmc`.

## Run

```bash
# on an idle machine; about 20 minutes
GTEX_DIR=<gtex> HLCA_DIR=<hlca> bash .../01_threads.sh <out>
python -I .../02_summary.py <out>                 # results/threads.tsv, fig_threads
cd <gtex>
GCEN2_STRING_GOLD=<pbmc>/work/string_gold.tsv python -I .../03_sensitivity.py
```

## Thread scaling

Intel Core i7-13700 (8 performance cores with 2 threads each, 8 efficiency
cores: 16 cores, 24 threads), 123 GB RAM, idle (load 0.18 at the start,
2026-11-17); mean of two runs (`results/threads.tsv`, `results/fig_threads`):

| command | 1 thread | 2 | 4 | 8 | 16 | speed-up at 8 | peak memory |
|---|---|---|---|---|---|---|---|
| `build` with FDR, GTEx 11,944 genes × 803 samples | 8.3 s | 5.0 s | 3.2 s | 2.6 s | 2.8 s | 3.2× | 0.55–0.66 GB |
| `module -e`, same matrix | 45.4 s | 25.7 s | 14.6 s | 9.0 s | 9.9 s | 5.0× | 1.41–1.50 GB |
| `diff`, two networks of 12.5M pairs | 12.9 s | 9.1 s | 6.8 s | 5.6 s | 4.8 s | 2.3× | 1.09–1.51 GB |
| `celltype-modules --backed`, 78,816 alveolar macrophages of the whole HLCA | 143 s | 101 s | 70 s | 53 s | 50 s | 2.7× | 4.2 GB |

- The compute kernels run in parallel; what remains is serial: reading and
  parsing the input, writing, and the steps not parallelised. Fitted to
  Amdahl's law at 8 threads, the serial part is 8% of the one-thread time of
  `module -e`, 21% of `build`, 28% of `celltype-modules` and 35% of `diff`
  (which mostly reads and writes 12.5M-line files).
- Beyond 8 threads the work moves onto the second threads of the performance
  cores and onto the efficiency cores: `build` and `module -e` are slightly
  slower with 16 than with 8, `diff` and `celltype-modules` slightly faster.
  8 threads, the setting of the benchmarks, is close to the best on this CPU.
- Memory hardly depends on the thread count (each thread's buffers are small
  next to the data).
- `celltype-modules` at 1 and 2 threads varied between the two runs
  (136 / 151 s, 113 / 89 s): with `--backed` the expression is read from disk
  as needed, so the page cache matters.

## Parameter sensitivity

GTEx whole blood, 5,000 genes, top 5 PCs regressed out (the main data of
`../wgcna_gtex`); one parameter at a time from the defaults used elsewhere.
STRING precision: the fraction of gene pairs (edges, or pairs within a
module) that are STRING interactions; the baseline over all pairs is 0.07%.

**Correlation cutoff** (`gcen2 build -p 0.001 -f -s -c r`, then `gcen2 module
-s 0.5`; `results/sensitivity_network.tsv`):

| r >= | edges | genes | edge STRING precision | GCEN modules (genes) | within-module enrichment |
|---|---|---|---|---|---|
| 0.3 | 420,249 | 4,950 | 0.6% | 99 (1,868) | 8× |
| 0.4 | 140,064 | 4,719 | 1.4% | 123 (1,495) | 30× |
| 0.5 | 44,411 | 3,954 | 3.5% | 192 (1,289) | 78× |
| 0.6 | 14,692 | 2,391 | 8.4% | 291 (1,212) | 137× |
| 0.7 | 5,191 | 1,135 | 17.3% | 200 (783) | 283× |
| 0.8 | 1,595 | 475 | 30.0% | 96 (371) | 510× |

A higher cutoff trades coverage for precision smoothly: each step of 0.1
divides the edges by about 3 and multiplies the share of STRING
interactions by 1.7–2.5; no setting is a cliff.

**GCEN module similarity** (`gcen2 module -s`, r >= 0.4 network;
`results/sensitivity_gcen_modules.tsv`): 0.3 merges most genes into one
module of 2,575 genes (enrichment 1.8×), 0.4 still one of 1,067 (5.7×); from
0.5 (the GCEN v1 default) modules are compact: 123 modules, 1,495 genes, 30×;
0.6 and 0.7 keep fewer genes in tighter modules (68×, 137×).

**WGCNA-style modules** (`gcen2 module -e --signed`; defaults: soft power
auto (10 here), deepSplit 2, minimum size 30, merge cut 0.25: 19 modules,
3,062 genes, 7.6×; ARI with the default over all 5,000 genes, unassigned as
one group, and over the genes in modules in both runs;
`results/sensitivity_wgcna_modules.tsv`):

| parameter | values | modules | genes in modules | within-module enrichment | ARI with default: all genes | ARI: genes in modules in both |
|---|---|---|---|---|---|---|
| soft power | 6 / 8 / 12 / 14 / 16 | 32 / 22 / 16 / 13 / 13 | 4,053 / 3,873 / 1,896 / 1,414 / 1,140 | 8× / 6× / 16× / 19× / 31× | 0.31 / 0.43 / 0.45 / 0.30 / 0.23 | 0.67 / 0.80 / 0.91 / 0.80 / 0.89 |
| deepSplit | 0 / 1 / 3 / 4 | 13 / 14 / 31 / 35 | 2,099 / 2,128 / 3,814 / 3,888 | 9× / 9× / 9× / 10× | 0.53 / 0.54 / 0.52 / 0.48 | 0.71 / 0.72 / 0.75 / 0.70 |
| minimum size | 20 / 50 / 100 | 23 / 14 / 8 | 3,492 / 2,902 / 2,019 | 7.5× / 6.4× / 8.6× | 0.75 / 0.89 / 0.51 | 0.91 / 0.89 / 0.79 |
| merge cut | 0.15 / 0.35 | 19 / 19 | 3,062 / 3,062 | 7.6× / 7.6× | 1.0 / 1.0 | 1.0 / 1.0 |

The soft power and deepSplit change the modules most: a higher power keeps
fewer genes in modules, and those are more often STRING partners (as a
higher correlation cutoff does); deepSplit sets how finely the tree is cut.
Most of the change is in which genes are assigned at all: among the genes
in modules in both runs the ARI is 0.67–0.91. The merge cut makes no
difference between 0.15 and 0.35 on these data. These are properties of WGCNA's method, which
GCEN2 reproduces (`../wgcna_gtex`); `--power auto` (pickSoftThreshold)
removes the most influential choice from the user.

## Files

| file | contents |
|---|---|
| `01_threads.sh` | the four commands with `-t` 1–16, twice, under `/usr/bin/time` |
| `02_summary.py` | means, speed-up, efficiency; `results/threads.tsv`, `fig_threads` |
| `03_sensitivity.py` | the parameter grid; `results/sensitivity_*.tsv` |
