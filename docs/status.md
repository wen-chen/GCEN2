# GCEN2 开发状态

> 更新日期：2026-11-18。
> 各项变化见 `CHANGELOG.md`，详细过程见 `git log`。

## 1. 功能

| 方面 | 状态 |
|---|---|
| GCEN v1 的全部命令（normalize、filter、build、module、annotate、enrich、rwr、merge、extract、shuffle、stats） | 完成；与 v1 0.6.3 比对：边、模块完全一致（`docs/migration.rst`）；v1 `util/` 中没有对应命令的 4 个工具，迁移指南给出与 v1 结果相同的 Python 代码 |
| WGCNA 模块（`module -e`）、consensus（`consensus-modules`）、模块保守性（`module-preservation`） | 完成；与 WGCNA 1.74 比对一致（ARI 1.0；保守性 Z 相对误差 < 7e-8） |
| 单细胞 hdWGCNA 流程（metacell、`celltype-modules`、harmony、DME、性状相关、`project-modules`） | 完成；与 hdWGCNA 0.4.12 比对一致 |
| 空间（`spatial-modules`、`spatial-build`、Moran's I） | 完成；Visium、Slide-seqV2、Visium HD 与 hdWGCNA 比对一致 |
| 差异共表达（`diff`）、多组学（`fuse`、`rwr-multiplex`）、轨迹（`trajectory`） | 完成；`diff` 与 DGCA 1.0.3、`rwr-multiplex` 与 RandomWalkRestartMH 比对一致；`fuse` 和 `trajectory` 没有已发表的对应工具，与按同一定义写的 R 实现比对一致 |
| 绘图（`gcen2 plot`，43 个命令）与报告（`gcen2 report`） | 完成；图库经过逐张视觉检查 |

## 2. 测试

- 946 个测试，在 Python 3.14 上全部通过。2026-11-18 在 Python 3.9 + 最老依赖（numpy 1.24、pandas 2.0、scipy 1.10、anndata 0.10、matplotlib 3.7、typer 0.9；与 CI 的“最低版本”任务相同）上从导出的源码安装后重测：927 通过、3 跳过（源码中不含 v1 样例数据），修正了 3 处只在旧版本上出现的测试问题（SciPy 1.10 没有 `false_discovery_control`；样例数据缺失时跳过；模块保守性中 grey 模块的 MAR 统计量随 NumPy 版本相差约 2e-7，黄金输出对该目录放宽到 1e-6）。三个目标平台都已在 CI 上验证：GitHub Actions 运行 44385168109（2026-11-18，提交 dc653c4 的公开版本）：从源码编译 C++ 扩展并跑全部测试——Windows Server 2022 / Python 3.12：935 通过、4 跳过；macOS 15（Apple Silicon）/ Python 3.12：934 通过、5 跳过；Linux：Python 3.9–3.14 与最老依赖全部通过；ruff 与 Sphinx 文档构建通过。跳过的测试需要不在公开仓库中的 v1 完整样例包，或在 macOS 上无法限制 Accelerate 的线程数。
- 黄金输出测试（`tests/golden/`）：每个命令在固定的小输入上运行，输出与冻结的参考一致；与线程数无关（1、4 线程一致）。这是性能优化的基准。
- 每个 `gcen2` 命令、每个 `gcen2 plot` 命令、每个公开函数都有测试；有测试检查新命令必须带测试。
- `gcen2` 包里没有未被测试执行的函数；行覆盖率 94.1%。
- 源码包（含示例数据）在干净环境中安装、现场编译后，包内测试 912 通过、1 跳过（需要本地的 v1 完整样例包）。
- 详见测试报告 `docs/test_report.md`。
- ruff 检查通过；`tools/check_doc_consistency.py` 通过。
- 功能冻结后做了两轮全面的代码审查与结构优化和一次最后检查（见 `docs/development_log.md` 阶段 8、9、12），结果不变（黄金输出测试通过）；有意的结果变化：`rwr-multiplex` 的默认转移矩阵改为与 RandomWalkRestartMH 相同；`diff` 增加 BH 校正（`fdr` 列）和 Spearman 方差 1.06/(n − 3)（阶段 12）。

## 3. 性能比较

`benchmarks/` 中有与 GCEN v1（8 个命令）、WGCNA（GTEx 全血）、hdWGCNA（PBMC、HLCA 58 万细胞、consensus、保守性、4 个空间数据集）、DGCA（差异共表达）、RandomWalkRestartMH（多层 RWR）以及 R 实现（网络融合、轨迹网络）的比较，每组都有脚本、原始结果和说明。

性能比较可在其他机器上复现：运行环境见 `benchmarks/environment-r.yml`（已验证可求解）和 `requirements-python.txt`，数据下载地址和校验值见各 README，脚本不含固定路径。

**计时**：2026-11-17 在空闲机器上按 `benchmarks/RETIMING_PLAN.md` 统一重测（冻结提交 91bf932；重测中发现 `gcen2 module -e` 的特征基因计算没有限制 scipy 自带 OpenBLAS 的线程数，修复后（faed7ac）重测了 `module -e` 的各项）。各 README 和总表中的数字都来自这次测量，记录表见 `RETIMING_PLAN.md` 第 5 节。

## 4. 文档

各文档的分工见 `docs/README.md`（文档地图）。2026-11-17 做了一次全面检查与结构整理：README 精简为概览，用法移入教程（新增“轨迹与多组学”）；四份技术文档各管一件事，删去相互重复的部分，并按当前代码改正了过时的描述（实现文档重写为 3.0）。

## 5. 发布目标

- **地址**：网站 https://www.biochen.org/gcen2/ ，源代码 https://github.com/wen-chen/GCEN2 （GCEN v1 保留在 https://www.biochen.org/gcen/ 和 https://github.com/wen-chen/GCEN ，网站的 `/gcen/` 作为 GCEN 系列总入口）。
- **分发方式**：源代码、PyPI（三个平台的预编译 wheel）、Bioconda（Linux、macOS；Bioconda 不为 Windows 构建）。
- **支持平台**：Linux x86-64（AMD64）、Windows x86-64（AMD64）、macOS Apple Silicon（arm64）。其他平台（Linux ARM、Intel Mac）不是目标；CI 的 wheel 构建按这三个平台配置。
- **命名**：GCEN2 之后的版本在同一仓库发布（3.x、4.x 仍叫 GCEN2）；只有不兼容的重写才换名并另开仓库。

## 6. 尚未完成或尚未验证

1. **发布**：2.0.0 已于 2026-11-20 发布到 PyPI（sdist 和 18 个 wheel，由 `release.yml` 推送 `v2.0.0` 标签发布；本地 `pip install "gcen2[all]==2.0.0"` 全部测试通过）。Bioconda 配方已提交（bioconda-recipes PR #71954，CI 全部通过，等待维护者合并），网站和 README 已按可用处理。文档发布在网站 <https://www.biochen.org/gcen2/docs/>，不使用 Read the Docs。
2. **Bioconda 合并后**：验证 `conda install -c conda-forge -c bioconda gcen2`（Linux、macOS arm64）。
3. **编译器**：文档只写已测试的版本（GCC 11、14、15，Clang 19，Apple Clang 17，MSVC 19.44 / Visual Studio 2022），不声明最低版本。
4. **PDF / SVG 的可编辑性**未在 Illustrator / Inkscape 中人工检查（PNG 已逐张检查）。
