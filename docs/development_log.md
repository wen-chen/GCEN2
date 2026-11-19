# GCEN2 开发记录

> 更新日期：2026-11-17。按阶段概述 GCEN2 的开发过程、主要决定和发现的问题。
> 每项变化见 `CHANGELOG.md`，每次提交见 `git log`，当前状态见 `docs/status.md`，
> 当时的计划和快照见 `docs/plans/archive/`。

## 阶段 0：初始实现（2026-03-15 以前，未纳入版本控制）

GCEN2 起初在版本控制之外开发：Python 包（预处理、网络构建、模块识别、功能注释、RWR、
单细胞伪重复网络、差异共表达、绘图、报告）、pybind11 绑定的 C++ 核心（相关、富集、模块、
RWR）、测试、Sphinx 文档、CI 和 bioconda 配方的骨架。2026-03-16 时约有 3,100 行
Python，127 个测试通过、2 个失败，README 为空。

## 阶段 1：空间模块与命令补全（2026-03-15 至 03-16）

- 计划：`2026-03-15-development-continuation.md`、`2026-03-15-p1-completion.md`
  （补齐 GCEN v1 的 `merge`、`extract`、`shuffle`、`diff` 命令）、
  `2026-03-15-doc-code-consistency.md`（文档与代码一致性检查器
  `tools/check_doc_consistency.py`）。
- 空间转录组模块 `gcen2.spatial`（提交 7250e95 至 d6721a0）：邻居图（KNN、半径、
  高斯权重），按 Moran's I 过滤边、空间平滑后的相关、按区域建网，`gcen2 spatial-build`。

## 阶段 2：纳入版本控制，与 GCEN v1 对齐（2026-11-06 至 11-07）

- 75fb764：把此前未提交的包、C++ 核心、测试、文档和构建配置纳入版本控制，同时修复了
  可扩展性问题（例如 C++ 相关内核对 F-order 输入给出错误的 r；RWR 改为稀疏转移矩阵；
  模块识别、统计改为稀疏算法）。
- 5bca872、8ab6bd2、2e4829e：在边的符号、注释/富集语义和数字格式上与 GCEN v1 0.6.3 一致；
  构建按块流式写出网络。
- fe01fec、b5bdd1f：与 GCEN v1 的逐命令比较（`benchmarks/compare_v1_v2.py`）；
  有意保留的差异写入 `docs/migration.rst`（例如 v1 的 `merge` 把方向拼进相关系数字段，
  是 v1 的缺陷）。

## 阶段 3：每个命令的性能（2026-11-07 至 11-13）

逐个命令去掉 pandas / scipy 的开销，把热点移到 C++（21cc1c4 至 d78164f）：
富集的超几何尾概率、模块识别的邻接构建与并查集、RWR、构建（可复用线程池、
不依赖 numpy 的读取路径）、统计与边打乱（与 numpy PCG64 相同的随机数流）、
`merge`、`extract`、`diff` 的流式 C++ 实现、报告的统计。每一步都用 v1 比较确认结果不变，
并重新测量（ff395fc、d78164f）。

## 阶段 4：新分析（2026-11-13 至 11-14）

- 轨迹网络 `gcen2 trajectory`（2ca6523）；多组学融合 `gcen2 fuse` 和多层 RWR
  `gcen2 rwr-multiplex`（3abd9d2 至 2d253aa）。
- 单细胞（9c13206 至 466922f）：AnnData 输入、kNN metacell、按细胞类型建网
  `gcen2 celltype`、考虑 dropout 的 Pearson。
- 外部比较：与 hdWGCNA 在 PBMC 10k 上（627fe4e、b3a4703，修正了比较脚本中 hdWGCNA 的
  随机种子设置），与 WGCNA 在 GTEx 全血上（ca5f8f5）。

## 阶段 5：WGCNA / hdWGCNA 工作流（2026-11-14 至 11-15）

- WGCNA 模块：TOM 边分数（0f33970）、`gcen2 module -e`（blockwiseModules：软阈值、TOM、
  层次聚类、动态剪树、合并）、特征基因、kME、性状相关（0574fef、acc633a）、
  `--preset hdwgcna`（336ab81）。
- 结果与平台无关：确定性的 AVX TOM 乘积内核（47b941f），整个 C++ 不做浮点缩并
  （FP contraction），AVX 只用于 x86（0e08592）；结果与线程数无关。
- hdWGCNA 工作流：`celltype-modules`（9d5fc65）、`spatial-modules`（6bf2fc4，Visium 网格
  metaspot；1c97268 加入无网格数据的 kNN metaspot）、细胞水平的特征基因与 kME（d37c3ae）、
  consensus 网络与 Harmony 整合（f604b26）、差异模块特征基因与性状相关（1c4f1db）、
  模块投射与保守性（f9e5962，使用与 R 相同的随机数生成器，置换与 R 一致）。
- 大数据：矩阵乘法相关内核（a2c3822）；按细胞类型逐个处理、`--backed` 读取图谱，
  不复制整个矩阵（99ad42f、73e2674、9c8a14b）。
- 外部比较：Human Lung Cell Atlas 58 万细胞（349dae1）、Visium HD（ed0ff58）、
  consensus 与保守性；所有比较汇总为一张表（064bed0）。
- 比较中发现 hdWGCNA 的一个问题：consensus 网络中，数据集不按排序给出时软阈值会配错数据集
  （`benchmarks/hdwgcna_consensus/README.md`）。
- 决定：开发期间的计时都是临时的（机器上常有其他任务），功能冻结后在空闲机器上统一重测
  （4ef48f4，`benchmarks/RETIMING_PLAN.md`）。

## 阶段 6：绘图与报告（2026-11-15 至 11-16）

- 规划见 `docs/plans/archive/plotting_plan.md`：期刊尺寸、PDF / SVG 中文字可编辑、
  每个模块在所有图中颜色一致。
- `gcen2.plot` 与 `gcen2 plot`（43 个命令；7126ba6 至 8d25b81），报告使用矢量图和
  编号图注（0e08a12）。
- 逐张视觉检查后的修正（bd99318）；metacell、Moran's I、merge、extract 的图（16a2489、
  aa64fc8）；所有图中细胞类型顺序一致（25bbdf9）。

## 阶段 7：质量与可发布性（2026-11-16）

- 黄金输出测试（9e2d62c）：每个命令在固定输入上的输出被冻结，作为此后性能优化的基准；
  每个命令、绘图命令和公开函数都有测试。
- 示例数据纳入仓库、`LICENSE`（fc72253）；ruff 检查通过（419abb8）。
- 兼容性（870c1cc）：在 Python 3.9 至 3.14 和最老依赖上验证，修复了在 Python 3.9 上
  所有命令都无法运行、旧版 typer、numpy 1、matplotlib < 3.9、anndata < 0.11 上的问题；
  RWR 中得分相同的基因改为按基因名排序，使结果与 numpy 版本无关。
- 文档（6ae4d33、fe6f236）：CHANGELOG、状态、引用、贡献指南、CLI 参考、教程，
  需求/设计/实现文档和代码说明更新到当前功能。
- 性能比较可在其他机器上复现（2024257）：运行环境文件、数据下载地址与校验值、
  脚本不含固定路径；Read the Docs 配置。
- 开发记录、测试报告、社区文件，开发计划归档（74999d8）。

## 阶段 8：功能冻结后的代码审查与结构优化（2026-11-16）

用静态分析（重复代码、死代码、复杂度）和逐模块阅读审查了全部 Python 和 C++ 代码，只改结构、
不改结果（黄金输出测试、与 GCEN v1 / WGCNA / hdWGCNA 参考结果的测试全部通过）：

- 命令行：2,264 行的 `cli.py` 拆为 `gcen2/cli/` 包（按领域分 8 个模块，共 1,956 行），
  命令按组排列；错误和警告由包装器统一输出（`ERROR:` / `WARNING:`），取代了 16 处重复的
  try/except 和 46 处“输出错误再退出”；43 个绘图命令从三种注册方式统一为一种声明式写法；删去从未生效的全局选项
  `--threads / --seed / --verbose`（文档本来就写“没有全局选项”）。
- 公共函数：C++ 扩展的加载（原 7 处各写一份）统一为 `gcen2.core.extension()` /
  `require_extension()`，并删去 `gcen2/core` 中把 `build/` 插到 `sys.path` 最前面的无效代码
  （它会遮蔽同名模块）；文件存在检查（30 余处）、TSV 字段引号（4 处）、模块文件写出、
  边文件列解析、绘图建图（9 处）、matplotlib 检查各合并为一个函数；WGCNA 与 consensus
  的模块合并循环和参数检查合并；只被测试使用的参考函数移到 `tests/reference.py`。
- C++：Fisher 检验的循环、相关内核的线程调度与排序、模块识别的逐行读取（改用已有的
  `edge_text::for_each_line`，同时修复了读到格式错误的行时文件句柄未关闭）、绑定层的 CSR
  检查和 numpy 数组转换各合并为一处（`src/bindings/binding_util.hpp`）。
- 结果：Python 16,603 → 16,071 行，C++ 7,461 → 7,417 行；`gcen2` 包语句 8,773 → 8,377，
  测试覆盖率 92.9% → 94.0%；新增测试检查所有命令对缺失输入给出 `ERROR:` 而非 traceback。

## 阶段 9：第二轮代码审查（2026-11-16）

第一轮之后 Python 已无 4 行以上的重复代码和死代码，C++ 在 `-Wall -Wextra -Wshadow` 下没有
项目代码的警告。第二轮看的是重复的实现、过长的函数和遗留的旧设计：

- 去掉 `build_by_celltype(max_genes_per_chunk=...)` 的分块建网（参数保留、不再起作用）：
  内核本身按块处理基因对，而分块会在每块内单独做 FDR 校正，与不分块的结果不同。
- WGCNA 特征基因的两份实现（模块、保守性）合为一份；hdWGCNA 式特征基因的两份 PCA 合为一份；
  anndata 检查的两份合为 `_anndata.require_anndata`。
- 长函数按步骤拆分：`consensus_tom_modules`（各数据集的 TOM、校准与 consensus）、
  `iter_modules_by_celltype`（参数检查、选细胞类型、元细胞成员表）。
- 未用到的参数、变量和再导出删去；`build_spatial_regions` 文档写明 `n_neighbors` / `radius`
  不起作用。

## 阶段 10：文档检查与结构整理（2026-11-17）

- 新增文档地图 `docs/README.md`：每份文档回答一个问题，同一件事只在一处写，其余处给出链接。
- README 从 303 行精简为概览（安装、快速入门、文档索引、性能比较汇总），各分析的用法移入教程；
  新增教程“轨迹与多组学”；图的总览移入绘图教程。README 与新教程中的命令都在示例数据上运行过。
- 四份技术文档分工：需求规格说明（2.2）只写要求与是否满足，删去项目结构、构建、示例代码和开发
  历史；设计文档（2.2）删去构建、测试、迁移三节（改为链接）；实现文档重写为 3.0（构建、实现
  规则、测试体系、发布、常见问题、新增代码的步骤）；代码说明删去与文档地图重复的文档清单。
- 按当前代码改正的过时描述：模块识别（`jaccard_modules`，并查集）、相关内核的线程模型与入口
  （`network_edges`）、富集内核、`merge`（相关系数相加，而非 Fisher's Z）、`build` 需要编译扩展、
  需求 3.2.5（分块计算已去掉）、bulk 教程中 WGCNA 的用时（约 3 秒，原写 2 分钟）。
- 发现：C++ 代码没有用到 Eigen 和 OpenMP，`strim.hpp` 未被引用；随后按决定删除（见阶段 11）。

## 阶段 11：删除不用的构建依赖（2026-11-17）

- 删除 Eigen 子模块（`third_party/eigen`）：C++ 代码早已不用它（RWR 改为 CSR、TOM 改为自己的矩阵乘法内核之后）。
  克隆不再需要 `--recursive`，CMake、conda 配方和 CI 去掉相应步骤，源码包从 29 MB 降到 26 MB（剩下的主要是示例数据）。
- 去掉 OpenMP：代码中没有 OpenMP 指令，线程都是 `std::thread`（`WorkerPool`、`parallel.hpp`）。
- 删除未被引用的 `src/util/strim.hpp` 和 `src/util/thread_pool.hpp`（GCEN v1 遗留，已被 `WorkerPool` 取代）。
- 从零重新编译后全部测试通过；扩展不再链接 OpenMP 运行库。

## 阶段 12：补充比较与最后检查（2026-11-17）

- 新增三组性能比较：差异共表达与 DGCA 1.0.3（GTEx 全血，男性对女性，1,250 万个基因对；
  同一个 Fisher z 检验，z 在相同 r 上差 1e-14）；多层 RWR 与 RandomWalkRestartMH，网络融合与
  按同一定义写的 R（data.table）实现（GTEx 共表达 + STRING + TRRUST 三层，18,082 个基因）；
  轨迹网络与 R（WGCNA 函数）实现（paul15 髓系分化，2,730 个细胞，9 个窗口，边与边的动态完全一致）。
  DGCA 已从 CRAN 撤下（2024-05），从 CRAN 存档安装；RandomWalkRestartMH 不在当前 Bioconductor 中，
  从 GitHub 固定提交安装。
- 发现：`rwr-multiplex` 的结果与 RandomWalkRestartMH 不同（前 10 名只有 1 个相同）。原因是转移矩阵：
  RandomWalkRestartMH 对整个超邻接矩阵按列归一化（度为 d 的基因跳层概率为
  delta / ((1 - delta) d + delta)），GCEN2 先在每层内归一化（每个基因都以 delta 跳层）。按决定，
  默认改为与 RandomWalkRestartMH 相同（`--transition supra`，分数差 1.5e-12），原来的做法保留为
  `--transition layer`（结果逐字节不变）；`tests/data/rwrmh/` 保存 RandomWalkRestartMH 的参考结果。
- 发现：需求规格说明写 `diff` 做 FDR 校正，实际没有（只按 p 值阈值分类）；DGCA 对 Spearman
  用方差 1.06/(n − 3)，GCEN2 用 1/(n − 3)。按决定补上：`diff` 输出 `fdr` 列（对检验过的基因对做
  Benjamini-Hochberg，`-f` 按 FDR 判定变化的边），`-m spearman` 用 1.06/(n − 3)；`celltype --diff`
  对两个 Spearman 网络自动使用。与 DGCA 的比较加上 Spearman 和 FDR，两者都一致。
- `diff` 的读写改为多线程（`-t`）：按 64 MB 整行分块，每个线程用自己的基因名表解析，按文件顺序合并
  （结果与线程数无关，出错时仍报准确的行号）；稳定并行排序；BH 校正移到 C++ 并行排序；z 检验分块
  多线程；写出用与 `fuse` 共用的 `write_blocks`。1,250 万对：13.6 s → 5.5 s（8 线程），输出逐字节不变。
- 发现：RandomWalkRestartMH 输出基因名时删去 `_1`（`hCG_1984214` 变成 `hCG984214`），比较脚本已处理。
- 新基准的 `run_all.sh` 默认每个工具跑 2 次（`REPS`），比较脚本取平均；`RETIMING_PLAN.md` 加入这三项。

## 阶段 13：统一重测（2026-11-17）

- 在空闲机器上按 `benchmarks/RETIMING_PLAN.md` 依次重测第 0–6c、8–10 项（冻结提交 91bf932，00:22–07:46，
  全部成功）；记录表见 `RETIMING_PLAN.md` 第 5 节。
- 补做两项（08:50–11:00）：PBMC 逐对内核计时（`03_gcen2.py` 加内核参数），输出与自动内核相同；
  HLCA 的 hdWGCNA 种子 1、2（第 7 项，`04_compare.py` 加入 hdWGCNA 种子间和全部 9 对种子的工具间 ARI）。
  hdWGCNA 种子间 ARI 中位数 0.64 / 0.60（5 万 / 全部细胞），与 GCEN2 种子间（0.69 / 0.61）和工具间
  （0.65 / 0.63）相当：各自 metacell 时的差异来自 metacell 抽样，而非方法。
- 发现：GTEx 上 `gcen2 module -e` 比开发时慢一倍且波动大。原因是 771998d 把特征基因改为
  `scipy.linalg.eigh`，而 scipy 自带的 OpenBLAS 在 `_blas_threads` 设定线程上限之后才第一次加载，
  上限对它无效，在 24 个逻辑 CPU 上运行。修复（faed7ac）：先加载 `scipy.linalg` 再设定上限；
  测试在新的解释器中检查。5,000 基因 3.3 → 1.7 s，模块不变。只影响没有先加载 scipy（经 scanpy、
  sklearn）的命令行 `module -e`；修复后重测了用到它的各项（GTEx、PBMC、Slide-seqV2、Visium HD、HLCA
  的网络步骤），HLCA、空间、共识的整体时间核对后不变。
- 与开发时的数字相比：WGCNA、hdWGCNA 在空闲机器上快了很多（如 GTEx 11,944 基因 WGCNA 69 → 36 s），
  所以部分倍数变小（WGCNA 模块 7.5× → 4.0×）；HLCA 的 hdWGCNA 现在包括细胞层面的特征基因
  （1,618 → 2,952 s），与 GCEN2 做的事相同；`diff` 多线程后与 DGCA 的比较 16.4 → 8.6 s。
- 各 README、总表、根 README、CHANGELOG 和设计文档中的数字改为这次测量的结果，删去“临时数字”的说明。

## 阶段 14：生物学案例（2026-11-17）

- 补生物学案例 `benchmarks/case_gtex_blood`（GTEx 全血，GO 2026-08-05 版；
  去掉 IEP、HEP、ND 证据和 NOT 限定，避免用共表达得来的注释评价共表达）。
- 模块：23 个 WGCNA 模块对应血液的细胞类型和程序（红系、血小板、未成熟中性粒、单核、NK / 细胞毒、
  初始 T、T、B、浆细胞、干扰素、增殖、热休克）；5 个模块与年龄相关且在两组对半样本中都重复
  （初始 T、B 细胞下降，单核、NK / 细胞毒、线粒体翻译上升），符合免疫衰老。死亡情况（Hardy、呼吸机）
  的信号在去掉的 5 个 PC 中（未校正数据 PC2 与呼吸机 r = −0.72）。
- RWR 功能预测的 3 折交叉验证（2,144 个 BP 条目）：AUROC 中位数 0.607，邻居投票 0.578，
  打乱网络 0.515，只用度 0.518；血液相关过程最好（I 型干扰素 0.95），神经相关接近随机。
- 修正：`plot enrich-dot`、`enrich-heatmap` 的颜色条在按 FDR 着色时仍写“−log10 p”，改为写出所用的列；
  加测试（915 个）。
- `annotate`、`enrich` 加 `fdr` 列：每个查询在背景中有注释的全部条目上做 BH（未命中的条目 p = 1，
  与 R `p.adjust` 相同）。BH 需要全部 p 值，而融合内核只保留显著的条目，所以在 C++ 内核里对每个查询
  排序本行全部 p 值后计算，再只输出显著条目的值；无扩展时的 numpy 路径相同（`_bh_by_row`）。
  其余列逐字节不变，用时不变（GTEx 网络注释 0.86 → 0.87 s）；黄金输出只为新列重新生成；
  案例改用工具的 `fdr`。917 个测试。
- 扩展性与参数敏感性 `benchmarks/scaling`（空闲机器，16:04–16:24，load 0.18）：8 线程相对 1 线程，
  `module -e` 5.0×、`build` 3.2×、`celltype-modules` 2.7×、`diff` 2.3×（Amdahl 串行部分 8–35%，
  主要是读写）；16 线程不再变快（8 个性能核）。参数：相关阈值平滑地以覆盖换精度；GCEN 模块相似度
  0.5 以下出现巨型模块；WGCNA 模块对软阈值和 deepSplit 最敏感，但两次都入模块的基因上 ARI 0.67–0.91。
- PyWGCNA 对照 `benchmarks/pywgcna_gtex`（2026-11-18，空闲机器）：PyWGCNA 2.2.1 在 pandas 3 和
  scipy 1.18 下都会出错，用其发布时（2025 年 2 月）的依赖版本单独建环境。参照是它所复现的 WGCNA
  分步流程（cutreeDynamic，pamRespectsDendro = FALSE，mergeCloseModules），不是 blockwiseModules。
  软阈值（给 WGCNA 的 RsquaredCut 0.85；PyWGCNA 默认 0.9，会选 12 而非 10）、TOM（相对差 5e-9）、
  基因树都与 WGCNA 相同；模块不同：校正数据上 ARI 0.30–0.44。逐步对照定位在 cutreeHybrid：
  树切割本身少分配一半基因（3,628 对 1,737 未分配；两者都分配的基因 ARI 0.97）；按大小重新编号时
  未分配的基因得到一个模块的编号，最大模块得到 0，于是 PyWGCNA 从不报告未分配基因，并输出小于
  最小模块大小的模块（2–20 个基因）。11,944 基因 197 s、6.9 GB，GCEN2 9.1 s、1.45 GB。
- `hdwgcna_pbmc/08_seed_ari.py`（2026-11-18）：hdWGCNA 流程在全部基因上的种子间 ARI（原来误用
  `between_tools.tsv` 与 `seed_stability.tsv`，那是 GCEN2 自己的相关网络模块，且 ARI 只算两边都在模块中的基因）；
  结果与 README 中的数字一致（工具间、GCEN2 种子间、hdWGCNA 种子间中位数 0.77 / 0.77 / 0.78）。
- 细胞作为样本的对照 `hdwgcna_pbmc/09_strategies.py`（2026-11-18）：六个细胞类型上比较 metacell、
  随机伪重复（10、50 组）、单细胞和 MAGIC 3.0.0（装在数据目录的单独路径，不进 `.venv`）。空基因起初
  在全部细胞间打乱，归一化后跟随 1 / 文库大小，使 metacell 出现假边；改为在文库大小相近的细胞（50 箱）
  内打乱。伪重复的和也要先按细胞归一化，否则核糖体基因对排在最前。结果：10 组伪重复排序与随机相当、
  两半不重复，边几乎全是预期的假阳性；MAGIC 使 19–37% 的空基因对成为边（1,000–1,700 万条边）；
  单细胞的最强基因对排序最好但边少、结构重复性低于 metacell；CD14 Mono 的 metacell 空基因假边率
  0.19%，约为名义值的 4 倍（有效样本数是近似）。
- `spatial-build` 示例 `hdwgcna_spatial/07_spatial_build.py`（DLPFC 151673）：默认阈值下网络 39 条边，
  STRING 精度 0.74；Moran's I > 0 过滤一条也不去掉（82% 的基因、69% 的空基因 I > 0，I > 0 不是检验）；
  平滑后边多 15–100 倍，精度低 3–8 倍，并让空基因得到边。
- 按对照结果改代码：`--strategy pseudo_replicates` 弃用（FutureWarning）；`zero_handling="impute"`
  警告插补后 p 值无效；`spatial-build -m filtered` 改为 Moran's I 的单侧检验（Cliff–Ord 随机化假设下的
  矩，允许非对称 kNN 权重；与置换检验的标准差相差 1–2%）加 BH 校正，`--moran-fdr` 默认 0.05，1 即旧规则
  I > 0；`--moran` 表加 `pval`、`fdr` 两列。黄金输出只为这三个文件重新生成（`sp_filtered.network`、
  `sp_moran.tsv`、`sp_moran_w.tsv`）。DLPFC 上检验保留 37% 的基因（I > 0 为 82%）、19% 的空基因，
  默认阈值下网络不变（0.3 时少 1 条边）。921 个测试。
- metacell 有效样本数改为逐基因对（2026-11-17）。诊断：CD14 Mono 上 0.19% 的假边率来自空基因的设计，
  不是 n / k：在文库大小分箱内打乱保留了 E[计数 | 文库大小]，而 RPS2 等基因的占比随细胞状态变化，
  归一化后空基因跟随状态（null_RPS2 一个就有 999 条边）。改用"占比在细胞间打乱、按文库大小泊松抽样"
  的空基因后假边率只有 2e-5。理论：与基因 y 无关的基因在 metacell 上的协方差是重叠矩阵
  C = 共享细胞数 / k，其与 y 的相关方差为 s(y) / M，s(y) = y'Cy / y'y；六个细胞类型上 M / s(y)
  与空基因相关的实测方差一致（`10_effective_samples.py`；CD14 Mono 按 s 分组 372 / 274 / 158 / 102，
  实测 351 / 260 / 160 / 110）。固定的 n / k（133）对 90% 的基因过严，对最平滑的 1% 过松。实现：
  `metacell_effective_samples`（M / s，裁剪到 [3, M]，Spearman 用秩），一对基因取较小者；
  `--effective-samples cells` 保留旧规则；细胞类型间的 Fisher z 检验仍用 n / k。结果：空基因对成边
  1.2e-5–3.3e-4（名义上限 5e-4）；PBMC 网络边数为原来的 1.0–5.6 倍；小细胞类型的模块变多、两半 ARI
  下降（B 0.87 → 0.41，FCGR3A Mono 0.79 → 0.45），因为 p ≤ 0.001 且不做 FDR 时偶然边按名义比例出现。
  PBMC 的 GCEN2 运行（24 次）、比较和计时（空闲机器，两次平均）重跑。黄金输出只为
  `celltype` 重新生成。空间示例的空基因同样改为泊松抽样：空基因 I > 0 为 53%（接近随机），检验只保留 1.3%；平滑后空基因不再得边（原来的 710 条来自文库大小的空间分布）。924 个测试。
- `gcen2 celltype` 默认开启 FDR（2026-11-18）：p 值校准后，单用 p ≤ 0.001 会让每个细胞类型
  三四千万个基因对中约 0.1% 成为偶然边，小细胞类型的网络和模块大多由它们组成。`build_by_celltype`
  新增 `fdr=True`（`network.build` 仍默认关闭），CLI 为 `--fdr/--no-fdr`。`09_strategies.py` 固定
  `fdr=False`，空基因对仍检验 p 值本身。PBMC 重跑 24 次与计时（空闲机器，两次平均）：边数
  355–122,020（原 23,245–302,654）；模块 9–87 个，覆盖 63–1,130 个基因，模块内 STRING 富集 11–208 倍；
  两半 ARI 0.62–0.93（FCGR3A Mono 两半只有 9 和 18 条边，两半都在模块中的基因只有 2 个）；
  用时 1.5–7.3 s（比 hdWGCNA 快 14–55 倍），内存最多 1.09 GB（FDR 要保留全部候选边，多约 0.4 GB）。
  黄金输出只为 `celltype` 重新生成。
- 细胞类型间差异检验改用逐对有效样本数（2026-11-18）：`write_celltype_networks` 在两个网络都有
  `gene_effective_samples` 时，每条共享边在每个细胞类型取两个基因中较小的有效样本数（不足 4 的边不检验，
  记为 `shared_stable`），否则仍用 n // k。`differential.compare` / `compare_files` 的 `n_a` / `n_b`
  可为"基因 → 样本数"映射；C++ `DiffTable` 新绑定 `names()`、`lo()`、`hi()`。校准
  （`11_diff_calibration.py`）：同一细胞类型的两半当作两个细胞类型，随机 2,000 个基因的全部 200 万对，
  p ≤ 0.05 的比例逐对为 4.1–4.9%（按有效样本数五等分各组 3.7–5.5%），n // k 只有 0.26–0.76%
  （最平滑的五分之一 0.95–2.9%，最嘈杂的 0.02–0.17%），过严，漏掉多数变化。`celltypes.tsv` 的
  `effective_samples` 和 diff 汇总的 `n_a` / `n_b` 改为中位基因的值。黄金输出只为 `celltype` 重新生成
  （中位基因 55 / 53 个有效样本，旧 n // k 为 60；变化的共享边 89 → 60）。930 个测试。
- 单精度 TOM 的影响重新检查（2026-11-18，新写 `wgcna_gtex/07_tom_precision.py`）：TOM 按单精度计算（`build --tom` 的边打分，误差最大 2.8e-7）或舍入到单精度、1 − TOM 仍用双精度时，模块完全
  相同；只有把 1 − TOM 本身舍入到 float32（1 附近间距 6e-8，几乎每对都偏 1.5e-8）才改变平均连接树，59 个基因换模块
  （树切割 ARI 0.97，最终 0.98）。原记录的"单精度 TOM 使 571 个基因换模块"来自最初实现（0574fef），现已不能复现，
  `wgcna_gtex/README.md` 和设计文档改为新结果。另核对线程数据：16 线程对 `diff` 仍有 1.2 倍提升、使 `build` 和
  `module -e` 略慢，内存随线程数变化最多 0.4 GB（原写"无进一步提升""内存不变"）。
- 迁移指南的 CSV 条目改正（2026-11-18）：GCEN2 没有 `csv_to_tsv` / `tsv_to_csv` 命令，原文写成"现在能正确处理引号
  里的逗号"。改为命令对照表加一行、"CSV files"一条说明用 `gcen2.io.read_csv` + `write_tsv` 转换；需求规格里"用
  Python 的 csv 模块"也改为实际的 pandas。补测试时发现表头检测按逗号直接切第一行：无表头的 v1 CSV 若第一个基因名是
  带逗号的引号字段，会被当成表头而丢掉第一个基因；`_detect_header` 改用 `csv.reader` 拆行。931 个测试。
- 迁移指南补上 v1 `util/` 里剩下的四个工具（2026-11-18）：`data_stat`、`generate_expr_matrix_from_rsem`、
  `generate_expr_matrix_from_stringtie` 各给一个 Python 函数（总体标准差、按 int(n·q) 取百分位、RSEM 的 FPKM/TPM 列、
  StringTie GTF 中 transcript 记录按基因求和），用本机编译的 v1 0.6.3 在样例矩阵和自造的小输入上对照，值相同（StringTie
  只差 v1 输出的 6 位有效数字）；`tests/test_migration_recipes.py` 直接执行指南里的代码块并与存下的 v1 输出比较。
  `calculate_accuracy`（邻居 GO 富集预测、命中自身注释的基因比例）不复现，只说明其算法和 GCEN v1 用的交叉验证。
  936 个测试。
- 全面检查（2026-11-18）：代码、测试、文档逐项核对。`differential.compare` 不给样本数时改用网络自己
  p 值的样本数（逐基因有效样本数优先），快速入门的例子原来按元细胞数检验。GCEN v1 没有差异共表达，CHANGELOG、
  Code Guide、status 和功能对比表改正。空间教程改为 Moran 检验（在示例数据上重跑，数字不变）。之前标为 2026-11-19 的改动日期改为实际的 11-17 / 11-18。Python 3.9 + 最老依赖
  重测，修正 3 处测试问题（见 `docs/status.md`）。937 个测试。
- 首次 GitHub CI（2026-11-18）：macOS 编译失败
  （Apple libc++ 没有浮点 `from_chars`，新增 `src/util/parse_double.hpp`，两条路径在 20 万个字符串上逐位相同）；
  Windows 编译失败（缺 `<iterator>`）；Linux 上 Spearman 网络的黄金输出在部分机器上第 6 位不同。加 CI 的
  Machine 步骤后确认：失败的机器有 AVX-512（numpy X86_V4），通过的只有 AVX2；本机用同样的库版本也通过。
  机理：NumPy 在 AVX-512 上的 `log1p` 最后一位不同，归一化值里数学上相等、浮点上略有差别的值（两次舍入）
  因而换序，Spearman 秩随之改变。按输入值确定地把 `log1p` 扰动 ±1 ulp 可在本机复现同样四个用例的差异。
  修正：对数归一化改为 x·1e4/library 一次舍入（metacell 用成员计数之和），相等的比值得到相同的值；
  平滑模式的加权和仍可能有这种情况，黄金用例 `spatial_weighted` 改用 Pearson。修正后扰动 6 个种子都不再改变
  任何黄金输出。黄金输出为 celltype、空间等步骤重新生成（最后几位变化）。
- 第三轮 CI：Linux 全部通过；Windows 又停在 MSVC 拒绝的常量 `1.0 / 0.0`（改为 `numeric_limits::infinity()`）；
  macOS 能编译，922 通过，空间三步（Moran's I、平滑、kNN metaspot）不同：测试数据在规则网格上，等距点很多，
  `cKDTree` 在不同平台上取的等距邻居不同。新增 `gcen2/_knn.py`（自身在前、按平方距离、等距按编号；第 k 与
  第 k+1 个等距时取该距离内的全部点再排），邻接图、kNN metaspot、metacell 都用它；`tests/test_knn.py` 与暴力
  排序逐行比较。macOS 上 SciPy 用 Accelerate，`threadpoolctl` 管不到，`test_blas_limit_covers_scipy` 在这种情况下跳过。
  空间三步的黄金输出重新生成。946 个测试。
- 第五轮 CI：Linux 与 macOS 15 全部通过；Windows 的"卡住"实为约 40 个测试失败，其中一个比较两个大文件的字节，
  pytest 用 difflib 生成差异报告极慢。原因是换行符：pandas `to_csv` 默认用 `os.linesep`、Python 文本模式写文件
  都会在 Windows 上写 `\r\n`，而 C++ 写 `\n`。本机把 `os.linesep` 设为 `\r\n` 即复现 41 个失败。修正：用
  语法树脚本给 `gcen2/` 与 `tests/` 中全部 109 处 `to_csv` 加 `lineterminator="\n"`，另三处文本写入用
  `newline="\n"`；模拟下全部通过。CI 的 pytest 加 `-rf`。
- 第六轮 CI：Windows 剩 23 个失败。换行符之外还有：测试与工具读写文本未指定编码（Windows 默认 cp1252：
  `cli.rst` 里的 `×` 读错、带非 ASCII 字符的输入写错）；pyvis 以默认编码写网页（`UnicodeEncodeError`）；
  测试子进程打印 `→` 到 cp1252 管道；黄金清单路径用了反斜杠；一个测试用 `write_text` 写参照文件（`\r\n`）。
  修正：包、工具、总览脚本和测试中所有文本模式读写都加 `encoding="utf-8"`（语法树脚本，按 UTF-8 字节偏移
  插入；包 3 处、工具与总览 9 处、测试 101 处）；pyvis 用 `generate_html` 后自行以 UTF-8 写出；
  `tests/conftest.py` 设 `PYTHONIOENCODING=utf-8`；清单路径 `as_posix()`。本机用包装过的 `open` 模拟
  cp1252 默认编码与 `\r\n` 换行，全部测试通过。
- 第七轮 CI：Windows 只剩 1 个失败——子进程按 UTF-8 输出后，测试用 `text=True` 读取时在 Windows 上按 cp1252
  解码，`stdout` 为 None。测试中 13 处 `subprocess.run(..., text=True)` 都加 `encoding="utf-8"`。
- 第八轮 CI（44385168109）全部通过：Windows 935 通过、macOS 15 934 通过，Linux 3.9–3.14 与最老依赖、ruff、
  Sphinx 全部通过。从首次推送起共八轮，修正的跨平台问题：Apple libc++ 的 `from_chars`、MSVC 的头文件与常量除零、
  AVX-512 上 `log1p` 末位导致的 Spearman 秩变化（归一化一次舍入）、等距 kNN 的平台相关取舍（按编号）、Windows 的
  `\r\n` 与 cp1252 默认编码（统一 `\n` 与 UTF-8）。

- 受影响的空间 benchmark 重跑（2026-11-19，空闲机器，只重跑 GCEN2 部分，hdWGCNA 输出与计时保留；旧输出备份为
  `runs/*.pre_ties`）：DLPFC、小鼠脑、Slide-seqV2、Visium HD 与 hdWGCNA 的一致性逐项不变（ARI 1.0、spot 特征基因
  5e-9、种子间 ARI 相同），模块与 hub 基因不变；用时 9.0 / 18.9 / 7.2 / 19.2 s
  （原 8.8 / 17.7 / 7.0 / 18.6 s，在重复间波动之内）。kNN 与网格 metaspot 的比较在小鼠脑上变化（knn 对网格 ARI
  0.54–0.59，原 0.44–0.57；种子间 0.62–0.74，原 0.48–0.69；最小 |r| 0.94）；`spatial-build` 只有平滑模式
  0.3 的边数 285,876 → 285,890。各 README 与总览更新。

## 主要原则

- **与 GCEN v1 兼容**：同样的命令和文件格式；结果应一致处完全一致，有意的差异写入迁移指南。
- **与原工具对照**：每个对应的工作流都与原工具（WGCNA、hdWGCNA、DGCA、RandomWalkRestartMH）在真实数据上
  比较结果和用时；没有已发表工具的分析与按同一定义写的 R 实现比较。
- **结果可复现**：与线程数、平台、依赖版本无关；排序都有明确的次序；随机过程使用固定的
  随机数流（与 numpy 或 R 相同）。
- **改动有测试把关**：黄金输出测试保证优化不改变结果；改变结果的提交须说明原因。
