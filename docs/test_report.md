# GCEN2 测试报告

> 测试日期：2026-11-17；版本 2.0.0.dev0（功能冻结后的代码审查、结构优化和最后检查之后）。
> 本报告是功能测试（正确性）。性能测试（与 GCEN v1、WGCNA、hdWGCNA、DGCA、RandomWalkRestartMH
> 等比较用时、内存和结果）
> 见 `benchmarks/README.md`，其中的计时于 2026-11-17 在空闲机器上统一重测（`benchmarks/RETIMING_PLAN.md`）。

## 1. 结论

| 项目 | 结果 |
|---|---|
| 测试数 | 946 |
| Python 3.14（开发环境） | 946 通过，0 失败，0 跳过（约 80 秒，24 核） |
| 从源码包安装（干净环境，现场编译 C++ 扩展） | 912 通过，1 跳过（需要本地的 GCEN v1 完整样例包，见第 5 节） |
| Python 3.9、3.10、3.11、3.12、3.13 | 结构优化之前全部通过；优化之后尚未重测（见第 4 节） |
| Python 3.9 + 最老依赖 | 同上 |
| 行覆盖率（`gcen2` 包） | 94.1%（8,330 行中 493 行未执行）；没有从未被执行的函数 |
| 代码检查 | `ruff check` 通过；`tools/check_doc_consistency.py` 153 项全部通过 |

## 2. 测试环境

- Linux x86_64（内核 7.0），GCC 15.3，24 核。
- 开发环境：Python 3.14.10，numpy 2.5.4，pandas 3.1.1，scipy 1.18.1，anndata 0.13.4，
  matplotlib 3.11.2，typer 0.29.1。
- 最老依赖：Python 3.9，numpy 1.24，pandas 2.0，scipy 1.10，anndata 0.10，matplotlib 3.7，
  typer 0.9（`pyproject.toml` 中的下限）。
- 外部参考结果：GCEN v1 0.6.3（本地编译）、R 4.4.3 + WGCNA 1.74、hdWGCNA 0.4.12、
  RandomWalkRestartMH 1.13.2。

## 3. 测试内容

### 3.1 按类别

| 类别 | 内容 | 主要文件 |
|---|---|---|
| 黄金输出（回归基准） | 每个 `gcen2` 命令（38 个步骤）在固定的小输入上运行，输出与冻结的参考一致：文本完全一致，数值相对误差 ≤ 1e-9；1 线程与 4 线程结果相同；能检出 1e-7 的改变 | `test_golden.py`、`tests/golden/` |
| 与 GCEN v1 一致 | 网络的边、模块、注释与 v1 0.6.3 的输出逐条比较；v1 的语义（符号、过滤、数字格式） | `test_integration.py`、`test_annotation.py`、`test_build_file.py`、`tests/data/v1_annotate/` |
| 与 R 工具一致 | TOM、bicor、软阈值、模块、特征基因、kME、空间 metaspot 等与 WGCNA / hdWGCNA 的参考输出比较；多层 RWR 与 RandomWalkRestartMH 的得分比较（加权与不加权，两个 delta）；差异共表达的 BH 校正与 R 的 `p.adjust` 比较。参考由 `tests/data/*/make_*.R` 生成 | `test_tom.py`、`test_tom_modules.py`、`test_spatial_modules.py`、`test_multilayer.py`、`tests/data/wgcna_tom/`、`tests/data/hdwgcna_spatial/`、`tests/data/rwrmh/` |
| 两种实现一致 | C++ 内核与 pandas / numpy 参考实现的输出逐字节相同（pandas 2 下到 1e-12）；多线程结果与单线程逐字节相同（例如 `diff` 用 1、3、4 线程和不同读块大小） | `test_network.py`、`test_build_file.py`、`test_gemm_kernel.py`、`test_multilayer.py` |
| 单元测试 | 预处理、相关、网络、模块、富集、RWR、差异、轨迹、多组学、单细胞、空间、consensus、保守性、DME、Harmony 等各函数的结果、边界情况和错误处理 | 其余 `test_*.py` |
| 命令行 | 每个命令的参数、输出文件、错误信息（缺失输入一律为 `ERROR:` 而非 traceback，警告为 `WARNING:`）；输出目录不存在时自动创建；Python 3.9 可用的类型注解；`python -m gcen2.cli` | `test_cli.py`、`test_golden.py`、`test_api_gaps.py` |
| 绘图与报告 | 全部 43 个 `gcen2 plot` 命令都生成 PDF / PNG / SVG；报告（HTML）包含全部章节和图 | `test_plot_commands.py`、`test_plot_*.py`、`test_report.py` |
| 文档与数据 | CLI 参考与代码一致；示例数据存在且可读；`gcen2` 的所有子模块可访问 | `test_cli_docs.py`、`test_example_data.py`、`test_api_gaps.py` |
| 性能比较代码 | 汇总表的生成；小规模的性能基准（pytest-benchmark） | `test_benchmark.py`、`test_benchmark_overview.py` |

### 3.2 规则（由测试强制）

- 每个 `gcen2` 命令必须有黄金输出步骤，每个 `gcen2 plot` 命令必须有测试；新增命令而没有测试时，测试失败。
- 结果与线程数、平台、依赖版本无关：排序都有明确的次序（例如得分相同时按基因名）。
- 改变结果的修改须重新生成黄金输出，并在提交说明中写明原因（`CONTRIBUTING.md`）。

### 3.3 按文件的测试数

| 文件 | 数量 | 文件 | 数量 |
|---|---|---|---|
| `test_network.py` | 133 | `test_rwr.py` | 18 |
| `test_multilayer.py` | 100 | `test_report.py` | 17 |
| `test_annotation.py` | 54 | `test_preprocessing.py` | 16 |
| `test_build_file.py` | 47 | `test_celltype_modules.py` | 16 |
| `test_plot_commands.py` | 44 | `test_spatial_modules.py` | 14 |
| `test_golden.py` | 40 | `test_api_gaps.py` | 13 |
| `test_tom_modules.py` | 39 | `test_dropout.py`、`test_correlation.py` | 各 12 |
| `test_cli.py` | 34 | `test_enrichment.py` | 13 |
| `test_module.py` | 33 | `test_plot_analyses.py` | 10 |
| `test_io.py` | 26 | `test_consensus.py`、`test_plotting.py` | 9、8 |
| `test_spatial.py`、`test_metacell.py` | 各 25 | `test_gemm_kernel.py`、`test_benchmark.py` | 各 7 |
| `test_trajectory.py` | 23 | `test_dme.py` | 6 |
| `test_tom.py` | 22 | `test_plot_wgcna.py`、`test_plot_base.py` | 各 5 |
| `test_differential.py` | 30 |
| `test_single_cell.py` | 21 | `test_preservation.py`、`test_plot_single_cell.py` | 各 4 |
| `test_integration.py` | 19 | `test_example_data.py`、`test_cli_docs.py`、`test_benchmark_overview.py` | 3、2、1 |

### 3.4 覆盖率最低的模块

| 模块 | 覆盖率 |
|---|---|
| `gcen2/single_cell/projection.py` | 80% |
| `gcen2/io.py` | 86% |
| `gcen2/rwr.py` | 86% |
| `gcen2/cli/trajectory.py`、`gcen2/cli/__init__.py` | 87% |
| `gcen2/__init__.py` | 88% |

未覆盖的行主要是错误处理分支和可选依赖缺失时的提示。

## 4. 兼容性

在 Python 3.9 至 3.13 和最老依赖上的测试于 2026-11-16 全部通过，并据此修复了若干问题
（`CHANGELOG.md` 的 “Fixed”）。之后的代码审查与结构优化（命令行拆分为 `gcen2/cli/` 包等）
只在 Python 3.14 上测试过；新代码沿用已验证的写法（命令注解用 `Optional`，有测试检查），
之后已由 CI 重跑确认：GitHub Actions 运行 44385168109（2026-11-18，提交 dc653c4 的公开版本）：从源码编译 C++ 扩展并跑全部测试——Windows Server 2022 / Python 3.12：935 通过、4 跳过；macOS 15（Apple Silicon）/ Python 3.12：934 通过、5 跳过；Linux：Python 3.9–3.14 与最老依赖全部通过；ruff 与 Sphinx 文档构建通过。跳过的测试需要不在公开仓库中的 v1 完整样例包，或在 macOS 上无法限制 Accelerate 的线程数。

## 5. 源码包

`gcen2-2.0.0.dev0.tar.gz`（26 MB，含示例数据）在一个只装了 Python 3.14 的
干净虚拟环境中用 `pip install "gcen2-2.0.0.dev0.tar.gz[dev]"` 安装，现场编译 C++ 扩展，
然后用源码包中的测试（删去包内的 `gcen2/` 目录，确保测的是安装后的版本）运行：
912 通过，1 跳过。跳过的测试需要 GCEN v1 0.6.3 的完整样例包（`data/sample/`，277 MB，
不随仓库分发），在开发环境中运行并通过。

## 6. 尚未测试

- wheel 打包（cibuildwheel）：CI 已配置，尚未运行（只在推送版本标签时运行）。macOS、Windows 上的源码编译和测试已在 CI 上通过。
- 编译器最低版本（GCC 7、Clang 5、MSVC 2017）：只用 GCC 15 编译过。
- PDF / SVG 在 Illustrator、Inkscape 中能否编辑：需要人工检查（PNG 已逐张检查）。

## 7. 如何复现

```bash
python -m pytest -q                          # 全部测试
python -m pytest -q --cov=gcen2              # 覆盖率
python -m pytest -q tests/test_golden.py     # 黄金输出
ruff check gcen2 tests tools benchmarks
python tools/check_doc_consistency.py
```
