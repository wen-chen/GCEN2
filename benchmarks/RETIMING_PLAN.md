# 统一重测计划（功能冻结后执行）

**已于 2026-11-17 执行**（记录见第 5 节）。之前提交的计时是开发过程中分几次测的，测量时本机常有其他任务在跑，会影响运行时间和内存；现在各 README 和总表中的数字都来自这次重测。

执行脚本：`/ssd1/gcen-retime/retime.sh`（第 0–10 项，第 7 项除外）、`retime_fix.sh`（修复后重测 `gcen2 module -e`）和 `retime_item7.sh`（PBMC 逐对内核计时与第 7 项），日志在 `/ssd1/gcen-retime/logs/`。

这次重测只重新测时间和峰值内存。结果的一致性（相同的边、模块、ARI）已经验证并提交，除非代码有改动，否则不用重新验证。

## 0. 本机的目录和环境

后面的命令都引用这几个变量。换一台机器时只改这里（数据的下载方法见各基准的 README，环境见 `environment-r.yml`、`requirements-python.txt`）：

```bash
export GCEN2_REPO=/ssd1/gcen2                  # 本仓库
export GCEN_V1=/ssd1/gcen-v1                   # GCEN v1 0.6.3（含 bin/）
export GTEX_DIR=/ssd1/gcen-bulk                # GTEx 数据目录（data/、work/）
export PBMC_DIR=/ssd1/gcen-pbmc                # PBMC 数据目录
export SPATIAL_DIR=/ssd1/gcen-spatial          # 空间数据目录
export HLCA_DIR=/ssd1/gcen-hlca                # HLCA 数据目录（consensus/ 在其下）
export DC_DIR=/ssd1/gcen-dc                    # 多层网络和轨迹基准的数据目录
# 本机：hdwgcna 环境的 R；DGCA 和 RandomWalkRestartMH 装在单独的库 /ssd1/gcen-rlib
export RSCRIPT="env R_LIBS=/ssd1/gcen-rlib /home/chenwen/program/miniforge3/envs/hdwgcna/bin/Rscript"
export GCEN2_STRING_GOLD=$PBMC_DIR/work/string_gold.tsv   # GTEx 基准用
B=$GCEN2_REPO/benchmarks
```

## 1. 什么时候开始

满足以下全部条件后才开始：

1. **功能冻结**：不再改 `gcen2/` 和 `src/` 的代码。记下冻结时的提交号，所有 GCEN2 计时都在这个提交上测：
   ```bash
   git -C $GCEN2_REPO rev-parse --short HEAD    # 写进下面的记录表
   export PATH=$GCEN2_REPO/.venv/bin:$PATH
   pip install --no-build-isolation -e ".[dev]"  # 确保已安装的扩展模块对应这个提交
   python -m pytest -q                           # 全部通过再测
   python -m pytest -q tests/test_golden.py      # 黄金输出：每个命令的结果与冻结的一致
   ```
2. **机器空闲**：除系统进程外，没有占用 CPU、内存或磁盘 I/O 的任务。
   ```bash
   uptime                                   # load average 应 < 1
   ps -eo pid,pcpu,pmem,etime,args --sort=-pcpu | head -15
   free -g
   ```
   不要结束用户自己的其他任务。如果有任务在跑，就等它结束，或者先问用户。
3. 重测期间不在本机启动其他计算任务，包括 Claude 会话里的测试和构建。

## 2. 测量约定

- 计时一律用 `/usr/bin/time -f "%e s %M KB"`，测量范围是整个命令或脚本，包括读入数据。
- **线程数**：与 GCEN v1 比较时用 4 线程（`--threads 4`），其余都用 8 线程（`-t 8`），两个工具一致。
- **重复次数**：
  - GCEN2 每项至少跑 2 次，取平均值或中位数（与各 README 的写法保持一致）。
  - R 工具（WGCNA、hdWGCNA）每项至少跑 1 次，时间允许就跑 2 次。
  - 第 8–10 项的 `run_all.sh` 默认每个工具跑 2 次（`REPS=2`），比较脚本取平均值。
- **记录负载**：每组开始前和结束后各记一次 `uptime`，写进下面的记录表。
- **只跑一个**：任何时候只运行一个基准，不并行。

## 3. 重测项目和顺序

数据都已在本机，不需要重新下载。时间是根据以前的运行估计的。命令中的 `.../` 即 `$B/`。

| # | 基准 | 数据目录 | 命令（在数据目录下执行） | 预计 | 更新的结果文件 |
|---|---|---|---|---|---|
| 1 | GCEN v1 vs GCEN2 | `$GCEN_V1`（v1 程序）、v1 样例数据 | `python benchmarks/compare_v1_v2.py --v1 $GCEN_V1 --sizes 5000 10000 20000 0 --threads 4 --mem-limit-gb 32 --out <临时目录>`（在仓库根目录执行） | 20 分钟 | `benchmarks/results.json`、`summary.tsv` |
| 2 | WGCNA，GTEx 全血 | `$GTEX_DIR` | `bash .../wgcna_gtex/run_all.sh` 跑 2 次（每次覆盖 `runs/time.txt`，先另存第一次的）；再跑 `python -I .../06_tom_modules.py`、`python -I .../04_compare.py` | 15 分钟 | `wgcna_gtex/results/time*.txt`、`time.tsv`、`tom_modules.tsv` |
| 3 | hdWGCNA 单细胞，PBMC 10k | `$PBMC_DIR` | `bash .../hdwgcna_pbmc/run_all.sh`（24 组，两个工具）；再跑 `python -I .../05_tom_modules.py`、`python -I .../04_compare.py`、`python -I .../07_cell_level.py runs/hdwgcna` | 45 分钟 | `hdwgcna_pbmc/results/gcen2_time*.txt`、`hdwgcna_time.txt`、`tom_modules.tsv` |
| 4 | hdWGCNA 空间，Visium | `$SPATIAL_DIR` | `bash .../hdwgcna_spatial/run_all.sh`（两个样本，各跑 2 次，含 `04_compare.py`） | 20 分钟 | `hdwgcna_spatial/results/compare.tsv`、`time.txt` |
| 5 | hdWGCNA 空间，Slide-seqV2 | `$SPATIAL_DIR` | `bash .../hdwgcna_spatial/run_slideseq.sh`（两个工具各 2 次，含种子对照和 `05_compare_slideseq.py`） | 15 分钟 | `slideseqv2_compare.tsv`、`time_slideseqv2.txt` |
| 5b | hdWGCNA 空间，Visium HD 16 µm | `$SPATIAL_DIR` | `bash .../hdwgcna_spatial/run_visium_hd.sh`（两个工具各 2 次，含种子对照和 `05_compare_slideseq.py runs visium_hd time_visium_hd.txt`） | 60 分钟 | `visium_hd_compare.tsv`、`time_visium_hd.txt`、`visium_hd_hdwgcna_timing.tsv` |
| 6 | hdWGCNA 单细胞，HLCA 肺图谱 | `$HLCA_DIR` | `bash .../hdwgcna_hlca/run_all.sh`（5 万细胞子集与全图谱；GCEN2 各 2 次，hdWGCNA 各 1 次），再跑 GCEN2 种子 1、2（见 `hdwgcna_hlca/README.md`）和 `04_compare.py` | 65 分钟 | `hdwgcna_hlca/results/*` |
| 6b | hdWGCNA 共识网络 + harmonized MEs + DMEs + 模块-性状相关，HLCA NK 细胞 | `$HLCA_DIR/consensus` | `bash .../hdwgcna_consensus/run_all.sh`（两个工具各 2 次，含种子 1、2 和 `03_compare.py`） | 20 分钟 | `hdwgcna_consensus/results/*` |
| 6c | hdWGCNA 模块投射与保守性，HLCA NK 细胞 | `$HLCA_DIR/consensus` | `bash .../hdwgcna_preservation/run_all.sh`（两个工具各 2 次，500 次置换；hdWGCNA 每次约 37 分钟） | 90 分钟 | `hdwgcna_preservation/results/*` |
| 8 | DGCA，差异共表达，GTEx 全血 | `$GTEX_DIR` | `bash .../dgca_gtex/run_all.sh`（Pearson 和 Spearman，各 2 次），再跑 `python -I .../dgca_gtex/03_compare.py` | 70 分钟（DGCA 的 Spearman 是单线程，5,000 个基因每次约 30 分钟） | `dgca_gtex/results/compare.tsv`、`time.txt` |
| 9 | RandomWalkRestartMH（多层 RWR）和 R 网络融合 | `$DC_DIR` | `bash .../multilayer_string/run_all.sh`，再跑 `python -I .../multilayer_string/04_compare.py` | 2 分钟 | `multilayer_string/results/compare.tsv`、`time.txt` |
| 10 | R 实现的轨迹网络，paul15 | `$DC_DIR` | `bash .../trajectory_paul15/run_all.sh`，再跑 `python -I .../trajectory_paul15/03_compare.py` | 1 分钟 | `trajectory_paul15/results/compare.tsv`、`time.txt` |
| 7 | HLCA 的 hdWGCNA 种子对照 | `$HLCA_DIR` | `$RSCRIPT .../02_hdwgcna.R work/hlca_<size>.h5 work/types.txt runs/hdwgcna/<size>_seed<s> <s>`，s = 1、2，两种规模 | 80 分钟 | 在 `04_compare.py` 里加上 hdWGCNA 种子间 ARI |

不含第 7 项共约 7 小时，含第 7 项约 8.5 小时（第 7 项放在最后）。

注意：
- **第 1 项**：`--out` 用临时目录，跑完后把 `results.json` 和 `summary.tsv` 复制回 `benchmarks/`。
- **第 2 项**：GTEx 的 GCEN2 时间现在是两次运行的平均值（`time_gcen2_gemm.txt`）。`run_all.sh` 每次运行都会改写 `runs/time.txt`，所以第二次运行前先把它另存。
- **第 3 项**：PBMC 表里 GCEN2 的时间（`gcen2_time_gemm.txt`，两次运行的平均值）是 6 个细胞类型单独循环跑两次 `03_gcen2.py` 得到的，不来自 `run_all.sh`。照此再跑一遍：
  ```bash
  for rep in 1 2; do for ct in "CD14 Mono" "CD4 T" "B" "CD8 T" "NK" "FCGR3A Mono"; do
    /usr/bin/time -f "${ct// /_} $rep %e %M" -a -o runs/gcen2_time_retime.txt \
      python -I .../03_gcen2.py "$ct" runs/gcen2_retime/${ct// /_} - 42 > /dev/null 2>&1
  done; done
  ```
- **第 3、6 项**：hdWGCNA 脚本（`hdwgcna_pbmc/02_hdwgcna.R`、`hdwgcna_hlca/02_hdwgcna.R`）现在包含 `ModuleEigengenes` + `ModuleConnectivity`（计时字段 `t_me`），GCEN2 的 `celltype-modules` 默认也计算细胞层面的特征基因和 kME。这样两边的"完整流程"包含同样的步骤，新数字与以前提交的数字不能直接比较，README 里要注明这一点。HLCA 的 `04_compare.py` 会用 hdWGCNA 改派后的 `modules_final.tsv` 比较完整流程。
- **第 4、5 项**：hdWGCNA 的空间计时是在 6bf2fc4 时测的，这次一起重测。
- **第 5b 项**：hdWGCNA 约 20 GB 内存，每次约 10 分钟（4 次，其中 2 次是种子对照）。数据已解压并准备好（`work/visium_hd.h5ad`、`work/visium_hd_export/`）。GCEN2 的 spot 层面特征基因和 kME 现在分块计算（Visium HD 从 11.3 GB 降到 2.3 GB），其他空间样本的结果不变（spot kME 只在第 10 位有效数字上不同）。
- **第 8–10 项**：输入已准备好（`$GTEX_DIR/work/dc_*.tsv`、`$DC_DIR/work/`），结果的一致性已在 2026-11-17 验证。`RSCRIPT` 必须能加载 DGCA 和 RandomWalkRestartMH（本机见第 0 节）。
- **第 6 项**：`hdwgcna_hlca/run_all.sh` 现在还会用 `--backed` 再跑一遍 GCEN2（`gcen2_top10_backed`、`gcen2_all_backed`），README 里两种模式都要报告。HLCA 需要约 20 GB 内存（GCEN2 13 GB，hdWGCNA 18 GB），两个工具不要同时跑。

## 4. 重测之后

1. 把新结果复制进各基准的 `results/`，按新数字更新各 README（表格、倍数、正文里引用的时间）。
2. 运行 `python benchmarks/overview.py`，重新生成总表和 `overview.tsv`。
3. 删掉 `benchmarks/README.md` 和 `hdwgcna_hlca/README.md` 里"计时是临时数字"的说明，改成本次测量的提交号、日期和负载。
4. 检查根目录 `README.md` 的 Benchmarks 小节和 `docs/design/` 里引用的时间（`grep -rn " s, \| GB" README.md docs/design`）。
5. 运行 `python benchmarks/overview.py --check`、`python tools/check_doc_consistency.py` 和 `python -m pytest -q`。
6. 在本机提交（不推送）。

## 5. 记录表

机器在开测前空闲（load 0.20），重测期间没有运行其他任务。各项依次运行，"开始前 load" 是上一项刚结束时的值（1 分钟平均），不是其他任务造成的。

| # | 基准 | 开始时间 | 开始前 load | 结束后 load | 备注 |
|---|---|---|---|---|---|
| 冻结提交 | 91bf932 | 2026-11-17 00:22 | 0.20 | 1.11 | 913 个测试通过 |
| 1 | GCEN v1 | 00:24 | 1.11 | 1.02 | |
| 2 | WGCNA GTEx | 00:40 | 1.02 | 13.02 | run_all 2 次；`module -e` 见 2f |
| 3 | PBMC | 00:46 | 13.02 | 10.69 | hdWGCNA 1 次，GCEN2 逐细胞类型 2 次；`module -e` 见 3f |
| 4 | Visium | 01:33 | 10.69 | 8.35 | 两个工具各 2 次 |
| 5 | Slide-seqV2 | 01:43 | 8.35 | 11.07 | 两个工具各 2 次；网络步骤见 5f |
| 5b | Visium HD | 01:58 | 11.07 | 16.83 | 两个工具各 2 次；网络步骤见 5bf |
| 6 | HLCA | 02:37 | 16.83 | 9.29 | GCEN2 2 次，hdWGCNA 1 次；网络步骤见 6f |
| 6b | HLCA 共识网络 | 05:00 | 9.29 | 20.68 | 两个工具各 2 次 |
| 6c | HLCA 投射与保守性 | 05:13 | 20.68 | 3.92 | 两个工具各 2 次 |
| 8 | DGCA | 06:41 | 3.92 | 1.59 | REPS=2 |
| 9 | 多层网络 | 07:45 | 1.59 | 1.21 | REPS=2 |
| 10 | 轨迹 | 07:46 | 1.21 | 1.17 | REPS=2 |
| 修复 | faed7ac | | | | 发现 `module -e` 的特征基因没有限制 scipy 自带 OpenBLAS 的线程数（在 24 个逻辑 CPU 上运行），修复后重测下面各项；其他命令不受影响（已核对 HLCA、空间、共识的 GCEN2 时间不变） |
| 2f | GTEx `06_tom_modules` × 2 | 08:02 | 0.88 | 4.32 | 另在空闲时手动测量（`results/time_gcen2_gemm.txt`）：无 FDR 与逐对内核的 `build` 各 2 次（07:50，load 0.72），`05_tom` 与 `--tom` 建网各 2 次（08:14，load 1.25） |
| 3f | PBMC `05_tom_modules` × 2 | 08:03 | 4.32 | 7.31 | |
| 5f, 5bf | Slide-seqV2、Visium HD 网络步骤 | 08:10 | 7.31 | 8.03 | |
| 6f | HLCA 网络步骤（`04_compare`） | 08:10 | 8.03 | 3.13 | |
| 3k | PBMC 逐对内核（`03_gcen2.py ... pairwise`）× 2 | 08:50 | 0.00 | 3.39 | `results/gcen2_time.txt`；输出与自动选择的内核相同（同样的边和顺序、相同模块） |
| 7 | HLCA 的 hdWGCNA 种子 1、2 | 08:51 | 3.39 | 4.32 | 两种规模，结束于 11:00；`04_compare.py` 加入 hdWGCNA 种子间与全部 9 对种子的工具间 ARI（`results/compare.tsv`、`time_hdwgcna_seeds.txt`） |
