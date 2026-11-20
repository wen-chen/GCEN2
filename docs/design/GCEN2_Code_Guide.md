# GCEN2 代码说明书

> 版本：2026-11-16（23 个命令和 `gcen2 plot` 的 43 个绘图命令；884 个测试；一致性检查器全部通过）。上一版：2026-03-16。
> 项目：GCEN2 — 基因共表达网络分析 Python 工具包

---

## 项目概览

GCEN2 是 GCEN v1（C++ 命令行工具）的 Python 重实现，并加入了 WGCNA / hdWGCNA 的模块分析（bulk、单细胞、空间）、轨迹与多组学网络、出版级绘图和报告。计算密集的部分（相关性、TOM、层次聚类、富集分析、模块识别、随机游走、网络文件读写与统计）用 C++ 实现，通过 pybind11 暴露为扩展模块 `_gcen_core`。Python 层提供 API、CLI、绘图和报告。

```
整体架构：

  用户代码 / CLI
      │
  gcen2 Python 包 (纯 Python API)
      │
  gcen2.core._gcen_core (C++ pybind11 扩展)
      │
  src/ C++ 算法实现（只用标准库和 pybind11）
```

---

## 目录结构

```
GCEN2/                   项目根目录
├── gcen2/               Python 包（主体）
│   ├── cli/             命令行（按领域分模块：network、modules、single_cell、spatial、
│   │                    trajectory、multilayer、report、plot；_common 为共用部分）
│   ├── core/            C++ 扩展加载层
│   ├── single_cell/     单细胞：元细胞、细胞类型网络与模块（hdWGCNA 流程）、DME、投射
│   ├── spatial/         空间转录组：spatial-build、Moran's I、空间模块（metaspot）
│   ├── plot/            出版级绘图（43 个 gcen2 plot 命令）
│   ├── tom_modules.py   WGCNA 模块
│   ├── consensus.py     WGCNA consensus 模块
│   ├── preservation.py  WGCNA 模块保守性
│   ├── harmony.py       Harmony 批次校正
│   ├── trajectory.py    轨迹（拟时序滑动窗口）共表达网络
│   ├── multilayer.py    多组学网络整合：加权融合与多层 RWR
│   └── templates/       HTML 报告模板
├── src/                 C++ 算法源码
│   ├── bindings/        pybind11 绑定层
│   └── util/            C++ 工具头文件
├── tests/               测试套件（tests/golden/：黄金输出；tests/data/：参考结果）
├── data/example/        示例数据（bulk、single_cell、spatial），纳入仓库
├── benchmarks/          与 GCEN v1、WGCNA、hdWGCNA、DGCA、RandomWalkRestartMH 等的比较
├── tools/               开发脚本（示例数据、CLI 参考、文档一致性检查）
├── docs/                所有文档（Sphinx 网站源码 + 设计文档）
│   ├── tutorials/       教程（bulk、单细胞、空间、绘图与报告）
│   ├── design/          技术设计文档（Markdown，不参与 Sphinx 构建）
│   └── plans/archive/   早期开发计划（归档，Markdown，不参与 Sphinx 构建）
├── conda-recipe/        Bioconda 打包配方
├── CMakeLists.txt       C++ 扩展构建配置
├── pyproject.toml       Python 包配置（scikit-build-core）
└── README.md、CHANGELOG.md、CONTRIBUTING.md、CODE_OF_CONDUCT.md、CITATION.cff、LICENSE
```

`data/sample/`（GCEN v1 完整样例包，277 MB）只在本地，不纳入仓库。

---

## 一、`gcen2/` — Python 主包

### `gcen2/__init__.py`
包入口文件，按需加载（PEP 562 `__getattr__`）：`import gcen2` 本身很轻，子模块和顶层函数在第一次访问时才导入。
- 设置版本号 `__version__`、作者信息
- `_SUBMODULES` 列出所有公开子模块（`preprocessing`、`network`、`module`、`tom_modules`、`consensus`、`preservation`、`annotation`、`rwr`、`differential`、`single_cell`、`harmony`、`spatial`、`trajectory`、`multilayer`、`plot`、`plotting`、`report`、`io` 等），`gcen2.<名字>` 均可直接访问（测试会在新进程中逐个检查）
- `GCENConfig`、异常类、I/O 函数按需从各模块取出
- `_HAS_CPP_EXTENSION` 标志（供用户检查 C++ 扩展是否可用）

---

### `gcen2/models.py`
**核心数据模型**，定义项目中所有数据结构（dataclass）：

| 类 | 说明 |
|----|------|
| `ExpressionMatrix` | 基因表达矩阵容器，支持 bulk（DataFrame）和单细胞（AnnData） |
| `Network` | 共表达网络，以边列表（DataFrame）存储，自动计算节点数/边数 |
| `NetworkStats` | 网络统计信息（节点数、边数、密度、平均度） |
| `ModuleSet` | 模块集合，包含所有模块（set 列表）及源网络引用 |
| `EnrichmentResult` | 富集分析结果（GO/KEGG），包含 p 值截断和 term 统计 |
| `RWRResult` | RWR 预测结果，包含排名基因表和收敛步数 |
| `DiffCoExResult` | 差异共表达分析结果 |

`Network.__post_init__` 自动从边列表推导节点集合和计数。

---

### `gcen2/preprocessing.py`
**数据预处理模块**，对应 GCEN v1 的 `data_norm` / `data_filter` 命令。

| 函数 | 说明 |
|------|------|
| `normalize(expr, method)` | 归一化基因表达，支持 5 种方法：`upqt`（上四分位数）、`median`（中位数）、`deseq`（DESeq2 几何均值）、`tmm`（TMM）、`hkg`（看家基因） |
| `filter_genes(expr, percent_mean, ...)` | 按表达均值百分位数过滤低表达基因 |

**Bug fix vs v1**：`filter_genes` 中的 `deviation()` 改用样本标准差（`ddof=1`）而非总体标准差（`ddof=0`）。

---

### `gcen2/network.py`
**网络构建与操作模块**，对应 GCEN v1 的 `network_build` / `network_merge` / `network_extract` / `network_shuffle` 命令。

| 函数 | 说明 |
|------|------|
| `build(expr, method, cor, pval, ..., layer=None, use_highly_variable=False, kernel="auto")` | 构建共表达网络，计算所有基因对的相关性（`kernel`：`pairwise` 与 GCEN v1 逐对求和相同，`gemm` 用矩阵乘法、样本多时快约 10 倍，`auto` 在 ≥ 64 个样本时用 `gemm`；`METHODS`：`pearson`、`spearman` 以及 dropout 感知的 `pearson_nonzero`、`pearson_pairwise_complete`，见 `src/correlation.cpp`），并按阈值过滤。`expr` 可以是基因×样本 DataFrame、`ExpressionMatrix` 或 AnnData（细胞×基因，稠密或稀疏，经 `_anndata.adata_matrix` 转为基因×细胞矩阵；`layer` / `use_highly_variable` 只用于 AnnData）；`build_to_file` 同样接受 AnnData |
| `stats(net)` | 计算网络统计信息（密度、平均度等） |
| `merge(networks, cor_cutoff)` | 合并多个网络（同一条边的相关系数相加，与 GCEN v1 相同），符号另列为 `direction` 列；`merge_files` 直接读写边文件（C++） |
| `extract(net, genes)` | 从大网络中提取指定基因的子网络 |
| `shuffle(net, swap_factor, seed)` | 随机打乱网络边（保度分布）用于零模型构建 |

`build()` 调用 C++ `_gcen_core.correlation.network_edges()`，`build_to_file()` 调用 `build_network_array()`（按块写文件）；两者都需要编译好的扩展（未编译时抛出 `GCENInputError`）。`extract_file`、`shuffle_file`、`stats_from_file`、`summary_from_file` 直接读写文件，无扩展时回退到 pandas。

---

### `gcen2/module.py`
**模块识别模块**，对应 GCEN v1 的 `module_identify` 命令。

| 函数 | 说明 |
|------|------|
| `identify(network, similarity, threads)` | 基于 Jaccard 相似度识别基因模块（结果与 GCEN v1 相同） |
| `modules_from_file(path, similarity, threads)` | 同上，直接读边文件（C++ 读取，不用 pandas；`gcen2 module -i` 使用） |
| `read_file_adjacency(path)` | 边文件 → 对称 CSR 邻接（C++ `read_edge_adjacency`） |

网络先转成对称 CSR 邻接，C++ `jaccard_modules` 并行计算有边相连的节点对的 Jaccard 相似度（各线程用自己的并查集），超过阈值的节点对的连通分量即模块（与 v1 的 DFS 结果相同）；单基因模块自动排除。

### `gcen2/tom_modules.py`
**WGCNA 式模块**（`gcen2 module -e`）：`tom_modules(expr, ...)` 复现 WGCNA 1.74 `blockwiseModules` 单块流程——C++ 双精度 1 − TOM（原地构建）、`hclust_average_inplace`（fastcluster 平均连接）、`cutree_hybrid`（dynamicTreeCut `cutreeHybrid` 逐行移植）、`wgcna_post_process`（kME 检查与重新分配）、`merge_close_modules`。`HDWGCNA` 预设复现 hdWGCNA `ConstructNetwork`（含 `TestSoftPowers` 的 bicor 选幂规则）。结果 `TomModules`：`modules()`、`eigengene_frame()`、`kme()`、`hub_genes(n)`、`trait_cor(traits)`；`write_module_files()` 写出模块、特征基因、kME、hub 基因、表型相关 TSV。

---

### `gcen2/annotation.py`
**GO/KEGG 功能注释与富集分析模块**，对应 GCEN v1 的 `network_go_enrich` / `network_kegg_enrich` 命令。

| 函数/类 | 说明 |
|---------|------|
| `GOTerm` | GO term 数据类（id, 名称, 命名空间, 父节点） |
| `GOOntology` | GO 本体图，支持祖先传播（`propagate()`） |
| `load_go(obo_path)` | 解析 `.obo` 文件，构建 GOOntology 对象 |
| `load_association(assoc_path, ontology)` | 解析基因-GO 关联文件（GAF 格式） |
| `enrich(query_genes, bg_genes, ontology, assoc, pval)` | 对基因集做超几何检验富集分析（调用 C++ `_gcen_core.enrichment`） |
| `annotate_go(net, obo, assoc, pval, threads)` | 对整个网络的每个基因做 GO 富集 |
| `annotate_kegg(net, k2ko, kegg_assoc, pval, threads)` | KEGG 通路富集分析 |

---

### `gcen2/rwr.py`
**随机游走预测模块（Random Walk with Restart）**，对应 GCEN v1 的 `rwr` 命令。

| 函数 | 说明 |
|------|------|
| `predict(network, seed_genes, gamma, prob_cutoff, ...)` | 以种子基因为起点做 RWR，预测功能相关基因 |

迭代公式：`P(t+1) = (1-γ)·W·P(t) + γ·P0`，收敛条件为 L1 范数 < tol。调用 C++ `_gcen_core.rwr.random_walk_restart()`。

---

### `gcen2/differential.py`
**差异共表达分析模块**（GCEN2 新增，v1 无对应命令）。

| 函数 | 说明 |
|------|------|
| `compare(network_a, network_b, method, pval, n_a, n_b, fdr=False, correlation="pearson")` | 用 Fisher's Z 变换检验两个网络共有的边的相关是否显著不同（Spearman 网络的方差为 1.06/(n−3)），p 值做 Benjamini-Hochberg 校正（`fdr` 列）；`fdr=True` 时按 FDR 判定 `shared_changed`。`n_a`/`n_b` 也可为"基因 → 有效样本数"的映射：每条边取两个基因中较小者，不足 4 的边不检验；不给时取网络元数据（`gene_effective_samples`，否则 `p_value_samples`，否则 `n_samples`，否则 30），元细胞网络因此不会按元细胞数检验 |

输出 `DiffCoExResult`，含 `z_score`、`p_value`、`fdr`、`category`（`A_only`、`B_only`、`shared_changed`、`shared_stable`）。

---

### `gcen2/plot/` — 出版级绘图
每一类分析的图（52 个函数、43 个 `gcen2 plot` 命令），统一风格：期刊版面宽度（89 / 120 / 183 mm）、5–7 pt 字号、PDF/SVG 中文字可编辑、每个模块在所有图中颜色一致（默认 WGCNA 颜色名，可选 Okabe-Ito）；默认输出 PDF、PNG、SVG。

| 文件 | 说明 |
|------|------|
| `_style.py` | 样式（`rc_params`、`style()`）、模块颜色与标签、`set_defaults`、`headless_backend`（无图形界面时用 Agg，避免 DISPLAY 失效时卡住） |
| `_save.py` | `save_figure`：同一数据输出相同文件（SVG id 固定、不写创建时间） |
| `_io.py` | 读取各类结果文件（文件或内存中的表；所有绘图函数经 `table` 读取 TSV） |
| `_common.py` | 画布与新建图（`canvas`、`figure`）、R 规约的树状图、标签避让 `repel_labels`、连通分量排布 `pack_components`、对数坐标、维恩图与 UpSet 图 |
| `wgcna.py` | 软阈值、树状图、特征基因网络、模块—表型、kME、MM-GS、hub 网络、模块 UMAP、TOM 热图、样本聚类、模块大小 |
| `single_cell.py` | 特征基因在嵌入/组织上的分布、点图、小提琴图、雷达图、DME、细胞类型总览、metacell、consensus 各数据集、保守性、Moran's I |
| `network.py`、`enrichment.py`、`function.py`、`diff.py`、`trajectory.py`、`multiomics.py`、`qc.py` | 网络（含 merge / extract）、富集、RWR、差异、轨迹、融合、表达 QC 的图 |

CLI 的绘图命令在 `gcen2/cli/plot.py` 中逐个声明（`_plot`：输入文件与选项映射到绘图函数的关键字参数），运行时才导入绘图代码。

---

### `gcen2/plotting.py`
**早期可视化 API**（保留以兼容旧代码；新代码用 `gcen2.plot`）。

| 函数 | 说明 |
|------|------|
| `network_interactive(net, ...)` | 生成交互式 HTML 网络可视化（基于 pyvis） |
| `network_static(net, ...)` | 生成静态出版质量网络图（基于 networkx + matplotlib） |
| `dotplot(result, ...)` | 富集结果气泡图（dot plot） |
| `barplot(result, ...)` | 富集结果水平柱状图 |
| `expression_boxplot(expr, ...)` | 基因表达量箱线图（跨样本对比） |
| `expression_violin(expr, ...)` | 基因表达量小提琴图 |
| `heatmap(expr, modules, ...)` | 模块内基因表达热图（clustermap） |
| `degree_distribution(net, ...)` | 节点度分布直方图 |
| `diff_network(result, ...)` | 差异共表达网络可视化（按边类型着色） |
| `trajectory_dynamics(traj, ...)` | 轨迹网络：每个窗口的边数，以及变化最大的边沿拟时序的相关系数曲线 |
| `eigengene_heatmap(me, ...)` / `module_trait_heatmap(tc, ...)` | 模块特征基因相关热图、模块—表型热图（WGCNA 式） |
| `spatial_eigengenes(coords, me, ...)` | 每个模块一幅 spot 特征基因空间图（以 0 为中心的 RdBu 色标） |

所有函数均懒加载 matplotlib / pyvis，未安装时抛出友好错误提示。`network_interactive` 输出单个自包含的 HTML 文件（`cdn_resources="in_line"`，不再在当前目录写 `lib/`）。

---

### `gcen2/report.py`
**HTML/PDF 报告生成模块**（GCEN2 新增）。

| 函数 | 说明 |
|------|------|
| `generate(network, modules, enrichment, rwr, expr, output, template, wgcna, figures_dir, celltypes, moran, merge, extract)` | 生成 HTML/PDF 报告：网络统计、merge/extract、Moran's I、模块、WGCNA 式模块章节、细胞类型总览与 metacell、富集、RWR |
| `wgcna_section(...)` / `celltype_module_sections(dir)` | 由 `module -e` / `celltype-modules` 的输出文件生成模块章节 |
| `celltype_overview(dir, coords, basis)` | 细胞类型总览与 metacell 图（需嵌入坐标） |

使用 Jinja2 模板（`templates/report_default.html.jinja2`）渲染；图由 `gcen2.plot` 绘制，以矢量 SVG 嵌入并带编号图注，`figures_dir` 时另存 PDF/PNG/SVG。PDF 需 weasyprint；缺少模板、或 weasyprint 缺系统库时给出明确错误。

---

### `gcen2/io.py`
**文件 I/O 工具**，提供统一的读写接口。

| 函数 | 说明 |
|------|------|
| `read_tsv(path)` | 读取 TSV 格式表达矩阵，返回 `ExpressionMatrix` |
| `read_csv(path)` | 读取 CSV 格式表达矩阵，返回 `ExpressionMatrix` |
| `read_h5ad(path)` | 读取单细胞 `.h5ad` 文件（需 anndata），返回 `ExpressionMatrix` |
| `read_network(path)` | 读取 GCEN 格式网络文件（`.network`），自动处理 v1/v2 列名差异，返回 `Network` |
| `write_tsv(df, path)` | 写出 TSV 文件 |
| `write_csv(df, path)` | 写出 CSV 文件 |

**v1 兼容处理**：`read_network()` 自动将 v1 的 `#Node1`/`Node2`/`Pvalue` 列名映射为 v2 的 `gene_a`/`gene_b`/`p_value`。

---

### `gcen2/cli/`
**命令行接口**，使用 Typer 框架，对应 GCEN v1 的所有命令行工具，并加入新功能的命令。完整参数见 `docs/cli.rst`（由 `tools/make_cli_docs.py` 生成）。没有全局选项。

| 文件 | 内容 |
|------|------|
| `__init__.py` | 按顺序导入各模块（即命令在 `--help` 中的顺序），用 `wrap_command` 包装每个命令 |
| `_common.py` | `app`、`plot_app`；`wrap_command`；共用函数：`fail`（输出 `ERROR:` 并以 1 退出）、`read_h5ad`、`check_groups`、`h5ad_only`、`module_options`（预设参数合并）、`echo_tom`、`read_module_table` |
| `network.py` | GCEN v1 的命令：normalize、filter、build、module、annotate、enrich、rwr、merge、extract、shuffle、stats；以及新增的 diff（v1 没有差异共表达） |
| `modules.py` | `module -e`（`wgcna_modules`）、consensus-modules、module-preservation、project-modules |
| `single_cell.py`、`spatial.py` | celltype、celltype-modules；spatial-modules、spatial-build |
| `trajectory.py`、`multilayer.py`、`report.py` | trajectory；fuse、rwr-multiplex；report |
| `plot.py` | 43 个 `gcen2 plot` 命令（`_plot` 声明式注册） |

`wrap_command` 统一处理每个命令：GCEN2 的错误（`GCENError`）输出为 `ERROR: <信息>`、退出码 1，用户警告输出为 `WARNING: <信息>`，因此命令本身只需抛出异常；先创建输出文件所在的目录；未给出的列表参数统一为 `None`（typer < 0.12 会传入空元组）。导入 CLI 不导入 numpy 和 pandas（各命令在函数内导入）。命令参数的注解用 `Optional[X]` 而不是 `X | None`（typer 会求值注解，后者在 Python 3.9 上出错；测试检查 `gcen2/cli/` 下的所有文件）。

已实现子命令：

| 命令 | 对应 v1 工具 | 说明 |
|------|-------------|------|
| `gcen2 normalize` | `data_norm` | 归一化表达矩阵 |
| `gcen2 filter` | `data_filter` | 过滤低表达基因 |
| `gcen2 build` | `network_build` | 构建共表达网络；输入为 `.h5ad` 时可用 `--layer`、`--hvg`、`--groupby`/`--group`（按 obs 列取细胞子集） |
| `gcen2 stats` | `network_stat` | 计算网络统计 |
| `gcen2 module` | `module_identify` | 模块识别 |
| `gcen2 annotate` | `network_go_enrich` / `network_kegg_enrich` | GO/KEGG 富集注释（合并为单一命令） |
| `gcen2 rwr` | `rwr` | 随机游走预测 |
| `gcen2 enrich` | （新增） | 独立富集分析（任意基因集） |
| `gcen2 report` | （新增） | 生成 HTML 分析报告 |
| `gcen2 merge` | `network_merge` | 合并多个网络 |
| `gcen2 extract` | `network_extract` | 提取子网络 |
| `gcen2 shuffle` | `network_shuffle` | 随机化网络 |
| `gcen2 diff` | （新增） | 差异共表达分析 |
| `gcen2 spatial-build` | （新增） | 空间转录组共表达网络（filtered/weighted/regions 三种模式；`--fraction`、`--normalize`、`--method`、`-t`） |
| `gcen2 spatial-modules` | （新增） | 空间样本的 WGCNA 式模块（hdWGCNA 空间流程：metaspot、TOM 模块、spot 特征基因、`--plot` 空间图）；Visium 用网格 metaspot，其他平台（`--metaspots knn`，`-k`）用空间近邻 |
| `gcen2 fuse` | （新增） | 多组学网络加权融合（共表达、PPI、调控等多层网络） |
| `gcen2 rwr-multiplex` | （新增） | 多层（multiplex）网络随机游走（RWR-M） |
| `gcen2 trajectory` | （新增） | 沿拟时序滑动窗口建网，输出每个窗口的网络与边/度的动态表 |
| `gcen2 module -e` | （新增） | WGCNA 模块（`--preset hdwgcna` 为 hdWGCNA），可写出特征基因、kME、hub、表型相关、软阈值表、树状图、hub TOM、基因显著性 |
| `gcen2 consensus-modules` | （新增） | WGCNA consensus 模块（多个数据集） |
| `gcen2 module-preservation` | （新增） | WGCNA 模块保守性（Zsummary、median rank） |
| `gcen2 celltype-modules` | （新增） | hdWGCNA 单细胞流程：每个细胞类型的元细胞、模块、细胞特征基因与 kME；`--consensus-by`、`--harmonize-by`、`--dme-by`、`--trait`、`--backed` |
| `gcen2 project-modules` | （新增） | hdWGCNA `ProjectModules`：参考模块在查询数据上的特征基因 |
| `gcen2 plot <图>` | （新增） | 43 个绘图命令，见 `gcen2/plot/` |
| `gcen2 celltype` | （新增） | 单细胞 `.h5ad` 按 obs 列为每个细胞类型建网（默认 kNN 元细胞，`--strategy`、`-k`、`--max-shared`、`--max-metacells`、`--use-rep`、`--normalize/--no-normalize`、`--layer`、`--hvg`、`--group`），写出 `<细胞类型>.network` 与 `celltypes.tsv`；`--diff` 时两两做差异共表达，写出 `diff/<A>__vs__<B>.tsv` 和 `diff/summary.tsv` |

---

### `gcen2/config.py`
**配置对象**（公开 API，保留以兼容）。

`GCENConfig` dataclass：`n_threads`（`GCEN2_THREADS` 环境变量，默认 2）、`random_seed`、`verbosity`、`temp_dir`（`get_temp_dir()`，未设置时为系统临时目录）、`figure_dpi`、`figure_format`。库函数和命令不读取它：线程数、随机种子由各函数的参数和各命令的 `-t`、`--seed` 给出。

---

### `gcen2/exceptions.py`
**异常体系**。

```
GCENError（基类）
├── GCENInputError   — 参数非法或输入文件有误
├── GCENComputeError — C++ 内核数值错误
└── GCENIOError      — 文件读写失败
```

---

### `gcen2/core/__init__.py`
**C++ 扩展加载层**。

导入 `_gcen_core`；未编译时为 `None`（`_HAVE_CORE = False`），Python 实现继续可用。其他模块统一通过 `extension()`（可无扩展的代码：返回扩展或 `None`）和 `require_extension()`（需要扩展的代码：未编译时抛出 `GCENInputError`）取用扩展；两者在调用时读取模块属性，测试用 `monkeypatch.setattr(gcen2.core, "_gcen_core", None)` 即可关闭扩展。

---

### `gcen2/_matrix_io.py`、`gcen2/_build_file.py`
`_matrix_io`：不依赖 pandas 的文本读写（`iter_lines`：二进制按行读、gzip 压缩文件按文件头自动解压）、表头判断、表达矩阵读取；各模块共用的 `require_file`（文件不存在时抛出 `GCENIOError`）、`tsv_field` / `tsv_header`（与 pandas `to_csv` 相同的引号规则）、`write_modules`（模块文件）、`correlation_specs`（merge、diff 交给 C++ 读取的边文件列）。`_build_file`：`gcen2 build` 从矩阵文件到边文件的分块流程（C++ 或 numpy 路径）。

### `gcen2/tom_modules.py`、`gcen2/consensus.py`、`gcen2/preservation.py`、`gcen2/_rrng.py`、`gcen2/harmony.py`
WGCNA / hdWGCNA 的方法：`tom_modules`（blockwiseModules：软阈值、TOM、平均连接聚类、动态剪树、特征基因合并；特征基因、kME、表型相关）；`consensus_tom_modules`（blockwiseConsensusModules）；`module_preservation`（modulePreservation）；`_rrng.RRandom`（R 的 Mersenne Twister 与 `sample()`，使置换与 R 完全相同）；`harmony_integrate`（R harmony 2.0.5 的 `RunHarmony`）。WGCNA 的特征基因只有一个实现（`tom_modules.module_eigengene`：较小 Gram 矩阵的最大特征向量），模块、consensus 和保守性共用；WGCNA 与 consensus 的模块合并共用 `_merge_modules`。`consensus_tom_modules` 按步骤分为 `_set_toms`（各数据集的 TOM 与软阈值）和 `_consensus_tom`（校准与 consensus 分位数）。hdWGCNA（Seurat）式的细胞 / spot 特征基因见 `spatial.modules.module_pc_scores`（`module_eigengene_scores` 为其第一列）。算法细节见设计文档 4.4–4.15。

### `gcen2/_anndata.py`
不导入 pandas 的 AnnData 工具（`network._prepare_matrix` 的 numpy 路径也使用）：`is_anndata(obj)` 按属性判断 AnnData；`adata_matrix(adata, layer, use_highly_variable, groupby, groups)` 先取细胞和基因子集，再把（可能稀疏的）表达矩阵直接转成基因×细胞的 C 顺序 float64 矩阵和基因名。由 `gcen2.single_cell` 再导出。`require_anndata()`：anndata 缺失时给出安装提示（单细胞和空间共用）。

### `gcen2/single_cell/` — 单细胞子包

#### `gcen2/single_cell/__init__.py`
子包入口，导入 `adata_matrix`、`is_anndata`、`metacells`、`pseudobulk`、`build_by_celltype`、`build_dropout_aware`。

#### `gcen2/single_cell/metacell.py`
`metacells(adata, groupby, k=25, max_shared=None, max_metacells=1000, min_cells=100, use_rep=None, n_pcs=30, layer, aggregate="mean", seed, threads, groups=None)`：仿 hdWGCNA `MetacellsByGroups`。`groups` 只处理这些组（结果与先取子集相同，但不复制表达矩阵；自算 PCA 时只用这些组的细胞）；聚合时只取成员细胞所在的行，不把整个矩阵转成 float64。每组内在 `obsm[use_rep]`（默认 `X_pca`，没有则自算 PCA：整数计数先按每细胞 1 万归一化并 log1p，有 `highly_variable` 时只用高变基因，稀疏矩阵隐式中心化后用 `svds`，稠密输入也转成稀疏走同一路径）上用 `cKDTree` 找 k 近邻；按随机顺序尝试每个细胞的邻域，与已选元细胞共享不超过 `max_shared`（默认 `min(10, k-1)`）个细胞才接受；成员表达取均值或总和；整数计数（`normalize=None` 自动判断）再把每个元细胞缩放到总数 1 万并 log1p（同 hdWGCNA `NormalizeMetacells`，否则测序深度会让所有基因正相关）。返回 AnnData（`obs` 含分组列、`n_cells`、`center`，`uns["metacell_members"]` 为成员细胞下标）。

#### `gcen2/single_cell/modules.py`
`modules_by_celltype(adata, celltype_key, ...)`：hdWGCNA 流程——每个细胞类型取在至少 `round(fraction × 细胞数)` 个细胞中表达的基因（`fraction_genes`，与 hdWGCNA `SetupForWGCNA(gene_select = "fraction")` 的取整一致）、在全部基因上做元细胞、再用 `tom_modules`（默认 `HDWGCNA` 预设）找模块。`iter_modules_by_celltype()`（参数相同）逐个细胞类型产出 `(类型, TomModules)`：元细胞成员一次选好（结果与一次性计算相同），表达谱到该类型时才聚合，CLI 写完一个类型就释放，内存只需容纳表达矩阵和一个类型的元细胞与 TOM；`modules_by_celltype()` 即 `dict(iter_modules_by_celltype(...))`。默认（`cell_level=True`，CLI `--cell-eigengenes`）还按 hdWGCNA `ModuleEigengenes` / `ModuleConnectivity` 在该类型的单个细胞上计算特征基因和 kME（`cell_eigengenes()`：按行分块计算，每个模块的基因取 LogNormalize 后的稠密块，用与空间流程相同的 `module_eigengene_scores`（固定列优先布局，结果与调用方的布局无关），kME 按行块累加，不需要细胞 × 基因的稠密或 float64 副本）。CLI `--backed` 时表达矩阵留在磁盘上（`read_h5ad(backed="r")`，X 须为 CSR），每个细胞类型的行只读一次，供基因筛选、元细胞和细胞层面共用，输出与读入内存时相同，`reassign=True` 时按 `ReassignModules` 改派 kME 为负的基因。`write_celltype_modules()`（接受字典或上述生成器）写出 `<类型>/modules.tsv`、`eigengenes.tsv`、`kme.tsv`、`hub_genes.tsv`（这几个基于元细胞）、`cell_eigengenes.tsv`、`cell_kme.tsv`（基于单个细胞；可选 `metacells.tsv`）和 `celltypes.tsv`（含 `reassigned` 列）。CLI：`gcen2 celltype-modules`（`--cell-eigengenes/--no-cell-eigengenes`、`--reassign/--no-reassign`）。

#### `gcen2/single_cell/dme.py`、`gcen2/single_cell/projection.py`
`find_dmes` / `find_all_dmes`：hdWGCNA 的差异模块特征基因（Seurat `FindMarkers` 的 Wilcoxon 检验）；`module_trait_cor`：模块—细胞表型相关。`project_modules`：hdWGCNA `ProjectModules`，参考模块（与查询数据重叠至少一半的基因）在查询细胞上的特征基因，可选 Harmony。`celltype-modules` 另写出 `metacell_members.tsv`（每个元细胞的细胞）、`dmes.tsv`、`trait_cor.tsv`、`cell_harmonized_eigengenes.tsv`、consensus 时的 `sets.tsv` / `set_eigengenes.tsv`。

#### `gcen2/single_cell/core.py`
**单细胞专项分析**（GCEN2 新增）：

| 函数 | 说明 |
|------|------|
| `pseudobulk(adata, groupby, layer, method, min_cells)` | 按细胞类型/cluster 聚合单细胞数据为 pseudo-bulk 样本（sum 或 mean） |
| `build_by_celltype(adata, celltype_key, ..., strategy="metacell", layer, max_genes_per_chunk, effective="overlap", fdr=True)` | 对每个细胞类型分别构建共表达网络。默认 `strategy="metacell"`：先在细胞类型内做 kNN 元细胞（`k`、`max_shared`、`max_metacells`、`use_rep`、`n_pcs` 传给 `metacells`），以元细胞为样本算相关；`"pseudo_replicates"` 为旧做法（细胞随机分成 `pseudo_replicates` 份求和，结果与以前相同），已弃用并发出 FutureWarning：PBMC 上 10 份的排序与随机相当（`benchmarks/hdwgcna_pbmc/09_strategies.py`）。元细胞重叠且每个是 k 个细胞的平均，独立样本远少于元细胞数：元细胞网络的 p 值和 FDR 由 `_effective_pvalues` 按每对基因的有效样本数（`metacell_effective_samples`，两基因中较小者；`effective="cells"` 时为 `n_cells // k`）重新计算后再按 `pval`/`fdr` 过滤；`fdr` 默认开启（与 `network.build` 不同：p 值校准后，单用 p ≤ 0.001 会让每个细胞类型三四千万个基因对中约 0.1% 成为偶然边；CLI `--no-fdr` 关闭）（内核按元细胞数预过滤的结果是其超集，BH 排名精确）；按全部元细胞计算时，模拟数据上会出现上千条模块基因与无关基因之间的假边。有效样本不足 4 的细胞类型跳过。模拟数据（4 个基因模块，3000 个细胞）上前 1200 条边命中真实模块边的比例：元细胞 99.9%，随机分份 10.8%。`max_genes_per_chunk` 不再起作用（为兼容保留）：相关内核本身按块处理基因对，而按基因分块会在每块内单独做 FDR 校正 |
| `effective_samples(net)` | 元细胞网络的一个样本数 `min(元细胞数, n_cells // k)`（随机分份为份数），用于样本过少的检查；网络没有逐基因有效样本数时（`effective="cells"`、随机分份）也用于细胞类型间的 Fisher z 检验 |
| `metacell_effective_samples(expr, members, k, method, log)` | 每个基因的有效独立元细胞数 `M / s(y)`，`s(y) = y'Cy / y'y`，`C` 为元细胞共享细胞数 / k；一对基因用两者中较小的值算 p 值（`effective="overlap"`，默认；`"cells"` 为旧的 `n_cells // k`） |
| `write_celltype_networks(networks, output_dir, diff=False, pval=0.05, threads=2)` | 写出各细胞类型网络（文件名替换为安全字符并去重）、`celltypes.tsv`；`diff=True` 时用 `differential.compare_files` 两两比较：两个网络都有逐基因有效样本数（`gene_effective_samples`）时，每条边取两个基因中较小者（各细胞类型分别取；不足 4 的边不检验，记为 `shared_stable`），否则取 `effective_samples`（不足 4 则跳过并警告）；`celltypes.tsv` 的 `effective_samples` 与 `summary.tsv` 的 `n_a`/`n_b` 为中位基因的值 |
| `build_dropout_aware(expr, method="pearson", zero_handling="both_nonzero", **build_kwargs)` | dropout 感知建网。`both_nonzero`：每对基因只用两者都非零的细胞（`network.build(method="pearson_nonzero")`，C++）；`pairwise_complete`：只去掉两者都为零的细胞（`pearson_pairwise_complete`，会发出警告）；`impute`：MAGIC 插补，发出警告：插补后的细胞不是独立样本，p 值无效（PBMC 上 12–27% 的空基因对成为边）。`expr` 可为 DataFrame 或 AnnData。模拟数据（基因模块 + 泊松计数 + 随机技术性 dropout）前 K 条边命中真实模块边的比例：随机 dropout 60–80% 且表达较高时 `both_nonzero` 0.87–0.99、普通 Pearson 0.21–0.63；零值主要来自低表达时普通 Pearson 0.76–0.995、`both_nonzero` 0.01–0.17；`pairwise_complete` 在所有场景都接近 0（剩下"一个为零、一个非零"的细胞使无关基因对的平均 r 约为 -0.43） |

---

### `gcen2/spatial/` — 空间转录组子包

设计文档：`docs/plans/archive/2026-03-16-spatial-design.md`。实现完成于 2026-03-16，22 个测试全部通过（1 个因 sklearn 未安装而跳过）。

#### `gcen2/spatial/__init__.py`
子包入口，对外暴露：`build_spatial_filtered`、`build_spatial_weighted`、`build_spatial_regions`（建网），`metaspots`、`spatial_modules`、`spot_eigengenes`、`group_trait_cor`、`write_spatial_modules`（WGCNA 式模块）。

#### `gcen2/spatial/neighbors.py`
**空间邻居图构建**：

| 函数 | 说明 |
|------|------|
| `build_neighbor_graph(coords, n_neighbors=None, radius=None)` | 基于 2D 坐标构建 KNN 或半径邻居图，返回高斯加权 scipy sparse CSR 矩阵；两参数互斥，均不传或同时传会抛 `GCENInputError` |

- **KNN 模式**：`KDTree.query(k+1)` 排除自身，高斯权重 `exp(-d²/2σ²)`，`σ` = 非零距离中位数；`W.maximum(W.T)` 对称化
- **半径模式**：`KDTree.query_pairs(radius)`，权重全为 1.0，天然对称
- 输入验证：`n_neighbors >= 1`；`radius >= 0`

#### `gcen2/spatial/core.py`
**三种空间共表达建网方法**：

三种方法都先取在至少 `round(fraction × spot 数)` 个 spot 中表达的基因（`fraction` 默认 0.05，0 为全部基因），整数计数再做 LogNormalize（`normalize=None` 自动判断），只把这些基因转成稠密矩阵。

| 函数 | 说明 |
|------|------|
| `build_spatial_filtered(adata, ..., moran_fdr=0.05)` | 标准建网后检验每个基因的 Moran's I（`_morans_test`，随机化假设下的单侧检验，BH 校正），保留两端基因均满足 I > 0 且 FDR ≤ `moran_fdr` 的边；`moran_fdr=1` 即旧规则 I > 0 |
| `build_spatial_weighted(adata, ...)` | 先用高斯邻居均值平滑表达矩阵（降低 dropout 噪声），再做标准建网；支持显式指定 `bandwidth` |
| `build_spatial_regions(adata, ...)` | 按空间区域分别建网：手动（`region_key` 指定 obs 列）或自动（KMeans 聚类坐标），< 3 个 spot 的区域发出 `UserWarning` 并跳过，返回 `dict[str, Network]` |

内部辅助函数：
- `_morans_i(X, W)` — 计算每个基因的 Moran's I 空间自相关指数，返回 shape `(n_genes,)` 数组
- `_morans_test(X, W)` — Moran's I、z 与单侧 p 值（Cliff–Ord 随机化假设下的矩，允许非对称的 kNN 权重）；`_bh(p)` — BH 校正
- `_get_coords(adata)` — 从 `adata.obsm['spatial']` 读取坐标，未找到时抛 `GCENInputError`
- `_get_expression(adata, layer, fraction, normalize)` — 基因筛选和归一化后返回密集 `(n_spots, n_genes)` float64 矩阵及基因名

#### `gcen2/spatial/modules.py`
**WGCNA 式空间模块（hdWGCNA 空间流程）**：

| 函数 | 说明 |
|------|------|
| `metaspots(adata, groupby, row_key="array_row", col_key="array_col", min_spots=50, layer, mode="sum", method="auto", spatial_key="spatial", k=7, max_shared=None, seed=42, threads=1)` | `method="grid"` 仿 hdWGCNA `MetaspotsByGroups`/`ConstructMetaspots`：每组（按名称排序；< `min_spots` 的组警告后跳过）内从最小列起每隔 4 列取中心列（第 4、8、… 列），行取与中心列同奇偶的行；中心 spot 与组内网格距离最小的 spot（六边形网格上为 4 个斜向邻居）以及同行、在组内去重列序列中相隔 2 位的 spot 相加（或平均）。与 hdWGCNA 一样，没有斜向邻居时同一个邻居可能被计入两次。`method="knn"` 用于没有 Visium 网格的平台（Slide-seq、Visium HD、Stereo-seq、成像平台的单细胞）：在 `obsm[spatial_key]` 坐标上取每个中心及其组内最近的 k − 1 个点，中心按随机顺序（`seed`）尝试，与已选 metaspot 共享不超过 `max_shared`（默认 round(2k/7)，k = 7 时为 2）个点才保留；点越小 k 应越大（10 µm 的 Slide-seq bead 约 100）。`method="auto"`：`obs` 有六边形网格（行 + 列奇偶一致）时用 grid，否则 knn（Visium HD 的正方形网格也走 knn）。返回 AnnData（索引为中心 barcode，`obs["n_spots"]`，`uns["metaspot_members"]`，`uns["metaspot_params"]`） |
| `log_normalize(X, columns)` | Seurat LogNormalize（每行按全部列的总数缩放到 1 万再 log1p），只输出 `columns` 列的稠密矩阵 |
| `spot_eigengenes(data, labels, levels)` | 仿 hdWGCNA `ModuleEigengenes`：对每个模块的基因做 Seurat ScaleData（z 分数，ddof=1，上限截断为 10），取不再中心化矩阵的第一主成分得分 u·s（由较小 Gram 矩阵的最大特征向量求得），符号与模块基因平均表达正相关 |
| `spatial_modules(adata, groupby, fraction=0.05, use_highly_variable, genes, min_spots=50, mode="sum", row_key, col_key, layer, threads, preset="hdwgcna", spot_level=True, reassign=True, method="auto", spatial_key, k=7, max_shared=None, seed=42, **tom_kwargs)` | 基因选择（`fraction_genes` 或 `genes` 列表）→ `metaspots` → LogNormalize → `tom_modules`（默认 `HDWGCNA` 预设）；`spot_level` 时在全部 spot 上计算特征基因（`metadata["spot_eigengenes"]`，spot × ME<k>）和 spot 层面 kME（`metadata["spot_kme"]`，基因与 spot 特征基因的 Pearson 相关，同 hdWGCNA `ModuleConnectivity`）；`reassign=True`（默认）时按 hdWGCNA `ReassignModules` 把对本模块 spot kME 为负的基因改派到 kME 最大的模块（全为负则归为未分配），特征基因不重算（与 hdWGCNA 相同）。需要原始计数 |
| `reassign_negative_kme(labels, kme, levels)` | 上述改派规则，返回新标签和被改派基因的下标 |
| `group_trait_cor(spot_me, groups)` | 每个 spot 特征基因与各组指示变量（0/1）的相关系数和 Student p 值（`module, trait, cor, p, n`，即 `gcen2 report --trait-cor` 格式） |
| `write_spatial_modules(res, output_dir, groups, n_hubs, write_metaspots)` | 写出 `modules.tsv`、`eigengenes.tsv`（metaspot）、`kme.tsv`、`hub_genes.tsv`、`spot_eigengenes.tsv`（ME0…）、`spot_kme.tsv`、`group_cor.tsv`、可选 `metaspots.tsv`，以及 `summary.tsv` |

DLPFC 151673 与 10x 小鼠脑两个 Visium 样本上，与 hdWGCNA 0.4.12 的基因、metaspot、datExpr（差 < 1e-14）、改派前后的模块（ARI 均为 1.0）、spot 特征基因和 spot kME 一致，见 `benchmarks/hdwgcna_spatial`。CLI：`gcen2 spatial-modules`。

---

### `gcen2/trajectory.py` — 轨迹共表达网络

需求 §3.6.1。细胞（或时间序列样本）按拟时序排序后切成等大小、可重叠的滑动窗口，每个窗口用 `network.build` 建网，再跟踪所有窗口中出现过的边。

| 符号 | 说明 |
|------|------|
| `build_trajectory(data, pseudotime, window_size=0.2, step=None, metacell_size=1, ...)` | `data` 为 AnnData（`pseudotime` 为 obs 列名或数值）或基因×样本 DataFrame；`window_size`/`step` 为细胞数（int）或比例（float），默认步长为半个窗口，末尾补一个窗口以覆盖全部细胞；`metacell_size>1` 时把窗口内相邻细胞平均成 metacell；拟时序缺失的细胞被剔除并警告 |
| `TrajectoryNetworks` | `windows`（`TrajectoryWindow` 列表）、`dynamics`（每条边在每个窗口的相关系数 `cor_k`、`n_windows`、`first`、`last`、`pattern`、`cor_change`、`trend`）、`window_table()`、`degree_dynamics()`、`write(dir)` |
| `PATTERNS` | 边模式：`stable`（所有窗口）、`gained`（首窗无、末窗有）、`lost`（首窗有、末窗无）、`transient`（仅中间窗口）、`intermittent`（首末都有但非全部） |

边动态表中的相关系数对所有窗口都计算（不只是通过阈值的窗口），与 C++ 内核结果一致（误差 < 1e-9）。

### `gcen2/multilayer.py` — 多组学网络整合

需求 §3.7.1。图层（layer）可以是任意基因网络（共表达、PPI、调控……），以 `Network`、DataFrame 或边列表文件给出；边权为 `|correlation|`、`|weight|` 或 1。

| 符号 | 说明 |
|------|------|
| `read_layer(path, weight_column=None, sep="\t")` | 读边列表文件为 `gene_a, gene_b, weight`；自动识别表头（GCEN 网络文件、STRING `protein1 protein2 combined_score` 等） |
| `read_layer_codes(path, ...)` | 同上，但返回 `(names, a, b, w)`：基因名列表与两列基因编码。首行由 pandas 读取以确定表头和权重列，其余行由 C++ `read_layer_file` 读取并编码（不为每行建字符串）；C++ 不能确定 pandas 读法的输入（引号、BOM、缺失值、`+` 开头的数字等）回退到 pandas。`fuse`/`fuse_to_file`/`rwr_multiplex` 对文件图层都走这条路径，多个文件并行读取 |
| `as_layer(layer)` | 把 Network / DataFrame / 路径统一为图层表 |
| `fuse(layers, weights=None, normalize="max", min_layers=1, min_score=0)` | 加权融合：每层边权先归一化（`max` / `rank` / `none`），融合分数为各层权重的加权平均（缺失层记 0）；输出 `gene_a, gene_b, score, n_layers, layers, w_<层名>` 的 `Network`。基因与基因对用整数编码处理；去重、跨层合并、打分、过滤和排序由 C++ `fuse_layer_codes` 完成，含 NaN 权重或超过 62 层时用 numpy 版 `_fuse_arrays_numpy`（两者结果逐值相同；分数按层顺序逐项相加，不用 BLAS 矩阵乘） |
| `fuse_to_file(layers, output, ..., threads=4)` | `gcen2 fuse` 使用：融合结果由 C++ `write_coded_table` 按基因编码直接写出（多线程，与 `fuse(...).to_edgelist()` 逐字节相同），不在 Python 中构造基因名列；无扩展时回退到 pandas。样本网络两层 650 万条边：写出 8.3 s → 1.0 s；加上 C++ 读取编码和按排序去重（numpy 2.x 的哈希 `np.unique` 较慢），整个 `gcen2 fuse` 16.8 s → 3.7 s，`gcen2 rwr-multiplex` 5.5 s → 2.1 s；去重与合并改用 C++ 后 `gcen2 fuse` 3.7 s → 2.6 s，内存 1.6 GB → 1.0 GB |
| `rwr_multiplex(layers, seed_genes, gamma=0.5, delta=0.5, weighted=False, aggregate="geometric", prob_cutoff=0.0, transition="supra", ...)` | 多层网络随机游走（RWR-M，Valdeolivas 等 2019）：每个基因每层一个节点；超邻接矩阵的对角块为 `(1 - delta)` 乘各层邻接矩阵（加权时除以该层最大权重），基因在层间的副本之间为 `delta / (L - 1)`。`transition="supra"`（默认）对整个矩阵按列归一化，与作者的 RandomWalkRestartMH 相同（`tests/data/rwrmh/` 参考结果）；`"layer"` 先在每层内归一化，每个基因都以 `delta` 概率跳层。无边的节点总是跳层；转移矩阵交给 C++ `random_walk_restart`。每层得分按几何平均 / 算术平均 / 求和聚合。单层时与 `rwr.predict` 完全一致 |

### `gcen2/templates/report_default.html.jinja2`
HTML 报告的 Jinja2 模板。包含网络统计表、模块摘要、富集结果表格、RWR 排名列表的 HTML 布局。

---

## 二、`src/` — C++ 算法实现

### `src/correlation.cpp` / `src/correlation.hpp`
**相关性计算内核**。

- 逐对内核 `correlate_rows`（与 GCEN v1 相同的逐对求和；`compute_pairwise_correlation` 为其包装），由 `network_stream.cpp` 按块调用：`network_edges`（返回边表，`network.build`）、`build_network_array` / `build_network_file`（边直接写文件，`build_to_file` / `gcen2 build`）
- 大样本快速路径：`prepare_gemm_rows`（每行减均值、除以离差平方和的平方根）+ `correlate_rows_gemm`（用 TOM 的 `gemm_nt` 按 96×256 分块求 r，再按同样的阈值、符号和 p 值过滤并按 (i, j) 排序）。求和顺序固定，结果与线程数、分块和指令集无关；r 与逐对内核相差 1e-15–1e-14。`use_gemm_kernel(kernel, n_samples)`：`"auto"` 在样本数 ≥ `kGemmMinSamples`（64）时使用；`make_correlate` 统一选择，`build`/`build_to_file`/`build_from_file` 和 CLI `build`、`spatial-build` 都有 `kernel`/`--kernel` 选项。GTEx 11,944 基因 × 803 样本 FDR 建网 13.3 → 2.7 s，DLPFC `spatial-build` 39 → 3.3 s
- 线程：基因对按行切成小块动态分给 `WorkerPool` 的线程，各线程结果合并后按 (i, j) 排序，结果与线程数无关
- Spearman 模式：先预计算秩，再做 Pearson 相关
- p 值：基于 t 分布（`df = n_samples - 2`）
- 绑定层先取出 numpy 数据，再释放 GIL 计算
- dropout 感知 Pearson（`prepare_nonzero_rows` + `correlate_rows_nonzero`，`pairwise_correlation` 的 `method` 为 `pearson_nonzero` / `pearson_pairwise_complete`）：每对基因只用部分样本，p 值按该样本数计算；不足 3 个样本或在子集上为常数时 r = 0、p = 1。`BothNonzero`（两基因都非零的样本）只遍历较稀疏基因的非零位置，另一基因的零值以 NaN 标记；`AnyNonzero`（至少一个非零，即去掉双零）用全样本中心化的点积加上双零样本数（位图 AND + popcount）精确修正，无双零时与普通 Pearson 逐位相同。逐对与 scipy `pearsonr` 核对：r 误差 < 1e-13

### `src/enrichment.cpp` / `src/enrichment.hpp`
**超几何检验（Fisher 精确检验）内核**。

- `hypergeometric_pvalue` / `batch_hypergeometric`：GCEN v1 的双侧 Fisher 检验（共用对数阶乘表；p 截断到 1）
- `hypergeom_sf`、`study_counts`、`significant_counts`：GCEN2 注释与富集所用的单侧超几何尾概率、研究计数以及“计数 + p 值 + 阈值”一次完成；`format_annotation_rows` 按 `to_csv` 格式写结果行

### `src/module_identify.cpp` / `src/module_identify.hpp`
**模块识别内核**。

- `jaccard_modules(indptr, indices, cutoff, n_threads)`：对称 CSR 邻接上，有边相连的节点对的 Jaccard 相似度
  `J(i,j) = (|N(i)∩N(j)| + a(i,j)) / (|N(i)| + |N(j)| - |N(i)∩N(j)| - a(i,j))`
  并行计算（每线程一个并查集），超过阈值的对合并为模块；每个节点一个标签（-1 为不属于任何模块）
- `read_edge_adjacency(path, ...)`：读边文件（`edge_text::for_each_line`），节点按名字排序，返回 CSR 邻接

### `src/rwr.cpp` / `src/rwr.hpp`
**随机游走内核**。

- `random_walk_restart_csr(indptr, indices, data, P0, ...)` — 在 CSR 稀疏转移矩阵上迭代至收敛（与 GCEN v1 相同的更新与停止规则；内存随边数增长，而不是 n²）
- `random_walk_restart_counts` — 在未归一化的邻接计数上迭代，列归一化在迭代中完成

### `src/table_writer.cpp` / `src/table_writer.hpp`
`write_coded_table`：把"基因编码 + 数值/整数/标签列"的边表写成制表符分隔文件，与 pandas `to_csv(sep="\t", index=False)` 逐字节一致（字符串按 QUOTE_MINIMAL 加引号，浮点数用 Python `repr`，NaN 为空字段）。按 32768 行分块多线程格式化、按顺序写出（`parallel.hpp` 的 `write_blocks`）。绑定为 `_gcen_core.network.write_coded_table`，供 `multilayer.fuse_to_file` 使用。

### `src/layer_file.cpp` / `src/layer_file.hpp`
`read_layer_file`：把一层的边列表文件读成"基因名（首次出现顺序）+ int32 编码 + |权重|"，与 `multilayer.read_layer` 的 pandas 读法一致（`pd.to_numeric` 对 `+` 开头的数字不是正确舍入，因此这类输入交回 pandas）。绑定为 `_gcen_core.network.read_layer_file`。

### `src/layer_fuse.cpp` / `src/layer_fuse.hpp`
`fuse_layer_codes`：`multilayer.fuse` 的计算部分。每层基因对编码为 `lo * n_genes + hi`，按 `lo` 计数排序后在桶内排序去重（保留最大权重，多线程按基因区间分块），再按 `max` / `rank`（与 pandas `rank(method="average", pct=True)` 相同）归一化；各层并集逐层归并，分数按层顺序逐项相加后除以层权重和（以 `-ffp-contract=off` 编译，不合并乘加），按分数降序、基因对升序多线程排序。NaN 权重或超过 62 层时返回 `ok=false`，由 Python 用 numpy 版计算。绑定为 `_gcen_core.network.fuse_layer_codes`。

### `src/tom.cpp` / `src/tom.hpp`
WGCNA 拓扑重叠：`TomScorer`（`build --tom` 的单精度边分数）、`tom_dissimilarity`（模块用的双精度 1 − TOM，在输出矩阵内原地构建）、`soft_threshold_fits` / `scale_free_fit`（pickSoftThreshold）。bicor 行由 `network_stream.cpp` 的 `tom_rows("bicor")` 准备。

### `src/tom_gemm.cpp` / `src/tom_gemm.hpp` / `src/tom_gemm_avx.cpp` / `src/tom_gemm_impl.hpp` / `src/tom_gemm_includes.hpp`
TOM 的矩阵乘法 `gemm_nt`：固定求和顺序（double 4 路、float 8 路交错累加，乘加分开舍入，再两两归约），AVX 内核（`-mavx` 编译，运行时检测，CMake 仅对 x86 目标启用）与基线内核（SSE2 / AArch64 NEON / 纯 C++，`tom_gemm_includes.hpp` 选择）结果逐位相同；`GCEN2_SIMD=baseline` 运行时强制基线，`-DGCEN2_GEMM_SCALAR` 编译时强制纯 C++。

### `src/hclust.cpp` / `src/hclust.hpp`
fastcluster 的 NN-chain 平均连接层次聚类（合并顺序与高度与 R 相同），`hclust_average` 与原地版本 `hclust_average_inplace`；以 `-ffp-contract=off` 编译。

### 网络文件与统计：`build_file`、`network_stream`、`edge_files`、`edge_table`、`edge_text.hpp`、`edge_format.hpp`、`network_stats`、`edge_shuffle`、`pcg64.hpp`、`huge_buffer.hpp`、`py_repr.hpp`
`gcen2 build` 从矩阵文件分块写出边文件（GCEN v1 数字格式）；merge / extract / diff 的边文件读写（diff 按 64 MB 整行分块多线程解析、稳定并行排序、多线程写出，并有并行排序的 BH 校正 `bh_adjust`；结果与线程数无关）；图统计；保度打乱（与 numpy PCG64 随机流相同）；大缓冲区；与 Python `repr` 一致的浮点格式。

### `src/util/parallel.hpp`
按块并行执行的辅助函数：`run_parallel`、`parallel_sort`（可选稳定排序）、`write_blocks`（多线程格式化、按顺序写出）；`layer_fuse`、`network_stream`、`edge_files`、`table_writer` 使用；相关内核的线程池是 `correlation.hpp` 中的 `WorkerPool`。

### `src/bindings/` — pybind11 绑定层

| 文件 | 说明 |
|------|------|
| `module.cpp` | 扩展主入口，定义 `_gcen_core` 模块并注册子模块（`correlation`、`enrichment`、`module`、`network`、`rwr`） |
| `correlation.cpp` | 相关性内核（pairwise / GEMM）、矩阵到边文件的建网，处理 GIL 释放和 numpy→raw pointer 转换 |
| `enrichment.cpp` | 超几何检验、批量显著性计数、注释行格式化 |
| `module_identify.cpp` | Jaccard 模块识别（稀疏 CSR 输入） |
| `tom_modules.cpp` | TOM、1 − TOM、软阈值拟合、层次聚类、行分位数 |
| `network.cpp` | 边文件读写、统计、打乱、merge / extract / diff、多组学图层读取与融合、表格写出 |
| `rwr.cpp` | CSR 稀疏随机游走 |

---

## 三、`tests/` — 测试套件

| 文件 | 测试内容 |
|------|---------|
| `conftest.py` | pytest fixtures（小型模拟网络、表达矩阵等） |
| `test_preprocessing.py`、`test_network.py`、`test_module.py`、`test_annotation.py`、`test_rwr.py`、`test_differential.py`、`test_io.py`、`test_build_file.py` | GCEN v1 命令对应的功能；C++ 路径与 Python 参考实现逐字节一致 |
| `test_correlation.py`、`test_gemm_kernel.py`、`test_enrichment.py`、`test_tom.py` | C++ 内核的精度 |
| `test_tom_modules.py`、`test_consensus.py`、`test_preservation.py` | WGCNA 模块、consensus、保守性，与 R 生成的参考结果比较（`tests/data/wgcna_tom/`）；R 随机数、Harmony |
| `test_single_cell.py`、`test_metacell.py`、`test_dropout.py`、`test_celltype_modules.py`、`test_dme.py` | 单细胞：伪 bulk、元细胞、dropout、hdWGCNA 流程、DME 与表型 |
| `test_spatial.py`、`test_spatial_modules.py` | 空间：三种建网、Moran's I；metaspot 与 spot 特征基因和 hdWGCNA 参考一致（`tests/data/hdwgcna_spatial/`） |
| `test_trajectory.py`、`test_multilayer.py` | 轨迹网络；多组学融合与多层 RWR（与 RandomWalkRestartMH 的参考结果比较，`tests/data/rwrmh/`） |
| `test_plot_base.py`、`test_plot_wgcna.py`、`test_plot_single_cell.py`、`test_plot_analyses.py`、`test_plotting.py` | 绘图样式与每个绘图函数 |
| `test_plot_commands.py` | 每个 `gcen2 plot` 命令在黄金输出上运行；新增命令未加用例时测试失败 |
| `test_golden.py` | 黄金输出：每个命令（38 步）在固定输入上的输出与 `tests/golden/expected/` 一致（文字完全相同，数值相对误差 1e-9）；新增命令未加步骤时测试失败。`python -m tests.golden.make_golden` 重新生成（只在有意改变结果时） |
| `test_report.py`、`test_cli.py`、`test_cli_docs.py`、`test_api_gaps.py` | 报告；CLI；CLI 参考手册与代码同步；兼容 API、PDF 报告、子模块访问、CLI 注解在 3.9 上可用等 |
| `test_integration.py` | 端到端集成测试（`data/example/bulk`；与 v1 网络逐边比较需本地 `data/sample`） |
| `test_example_data.py` | 示例数据与 gzip 读取 |
| `test_benchmark.py`、`test_benchmark_overview.py` | 性能微基准（需 `pytest-benchmark`）；`benchmarks/README.md` 总表与结果一致 |

运行方式：
```bash
python -m pytest -q                        # 全部（约 3–4 分钟）
python -m pytest -q tests/test_golden.py   # 黄金输出（几秒）
python -m pytest -q -m "not benchmark"     # 与 CI 相同
```

---

## 四、`data/example/` — 示例数据（纳入仓库）

共约 27 MB；来源、许可和每部分的命令见 `data/example/README.md`，`tools/make_example_data.py` 可从原始数据重新生成（结果逐字节相同）。

| 目录 | 内容 |
|------|------|
| `bulk/` | GCEN v1 0.6.3 样例数据：表达矩阵（32,940 基因 × 15 样本）、GO（`go-basic.obo.gz`、`gene_go.assoc.gz`，gzip 压缩，GCEN2 可直接读）、KEGG（`K2ko.tsv`、`gene_kegg.assoc`）、RWR 种子、提取与富集用的基因列表 |
| `single_cell/pbmc_subset.h5ad` | 10x PBMC 10k 中 3 个细胞类型的 1,800 个细胞（计数、PCA、UMAP） |
| `spatial/mouse_brain_anterior_subset.h5ad` | 10x Visium 小鼠脑前部切片（2,695 spot、2,000 基因、网格坐标、区域） |

h5ad 文件按旧版 anndata（0.10）可读的编码写出。GCEN v1 完整样例包（含 186 MB 的参考网络）只放在本地 `data/sample/`，不纳入仓库。

---

## 五、`docs/` — 文档

各文档的分工见 `docs/README.md`（文档地图）：Sphinx 用户文档（`docs/*.rst`、`docs/tutorials/`，`cli.rst` 由 `tools/make_cli_docs.py` 生成）、`docs/design/` 下的四份技术文档、项目记录（`status.md`、`development_log.md`、`test_report.md`）和 `docs/plans/archive/`（早期计划，原样保存）。

---

## 六、构建与打包配置

### `pyproject.toml`
Python 包配置（PEP 517，使用 scikit-build-core 作为构建后端）：
- 项目元数据（名称、版本、依赖）
- 核心依赖：numpy ≥ 1.24、pandas ≥ 2.0、scipy ≥ 1.10、threadpoolctl、typer ≥ 0.9（下限为实际测试过的最老版本）
- 可选依赖组：`single-cell`（anndata ≥ 0.10 / scanpy）、`visualization`（matplotlib ≥ 3.7 / seaborn / pyvis / networkx）、`report`（Jinja2/weasyprint）、`cli`（空，为兼容保留）、`all`（全部）、`dev`（pytest、sphinx、ruff 等）
- ruff 配置（与代码惯例一致：矩阵用大写名、typer 默认值、CLI 注解用 Optional）
- pytest 配置（标记注册、测试路径）
- 定义 CLI 入口点 `gcen2 = gcen2.cli:app`

### `CMakeLists.txt`
C++ 扩展构建配置：
- C++17 标准，`-O3 -Wall -Wextra -ffp-contract=off`（MSVC `/O2 /W3`）；线程用 `std::thread`
- 非 MSVC 编译器对所有文件加 `-ffp-contract=off`（不把乘加合并成 FMA），使结果不随目标平台变化（AArch64 GCC 和 `-march=native` 默认会合并），C++ 路径与 numpy 版逐位一致
- TOM 乘法的 AVX 内核 `tom_gemm_avx.cpp` 仅对 x86 目标以 `-mavx` 编译（`CMAKE_OSX_ARCHITECTURES` 须仅为 `x86_64`，编译器须支持 `-mavx`）
- 自动查找 Python、pybind11（无其他 C++ 依赖）
- 构建单一扩展模块 `_gcen_core`，链接 `src/` 下的全部内核和 `src/bindings/` 的绑定
- 可选 `-DGCEN2_NATIVE=ON` 启用 `-march=native`（本地开发时提速 10-15%；结果与默认构建逐位相同）

### `conda-recipe/`

| 文件 | 说明 |
|------|------|
| `meta.yaml` | Bioconda 打包配方（从 PyPI 的 sdist 构建；构建脚本写在 `build.script` 中，`pip install .`；构建依赖、运行依赖、测试导入）。只构建 Linux x86-64 和 macOS arm64，Windows 用 PyPI。提交在 bioconda-recipes 的 `recipes/gcen2/` |

---

## 七、CI / 开发配置

### `.github/workflows/ci.yml`
GitHub Actions CI：Linux 上 Python 3.9–3.14、macOS 15 与 Windows 上 3.12 的全部测试；最老依赖版本；Sphinx 文档（`-W`）；ruff；性能微基准。已在 GitHub 上全部通过。

### `.github/workflows/release.yml`
发布：构建 sdist 和 cp39–cp314 的 wheel（cibuildwheel；Linux x86-64、Windows AMD64、macOS arm64），每个 wheel 在不含源码的目录里跑 golden 与 CLI 测试。手动运行发布到 TestPyPI 并在三个平台上安装验证；推送 `v*` 标签（须与 `pyproject.toml` 的版本一致）发布到 PyPI。用 trusted publishing，不存 token。

### `.pre-commit-config.yaml`
pre-commit hooks：行尾空格、文件结尾、YAML/TOML 检查、大文件检查（`data/example/` 和 `tests/golden/` 除外）、ruff。

### 根目录文档与 `tools/`
`README.md`、`CHANGELOG.md`（版本变化）、`CONTRIBUTING.md`（开发、测试与提交规则）、`CITATION.cff`（引用信息）、`LICENSE`（GPL-3.0）。`tools/`：`make_example_data.py`（重新生成示例数据）、`make_cli_docs.py`（生成 `docs/cli.rst`）、`check_doc_consistency.py`（文档与代码一致性检查）。

### `benchmarks/`
与 GCEN v1（`compare_v1_v2.py`）、WGCNA（`wgcna_gtex/`）、hdWGCNA（`hdwgcna_pbmc/`、`hdwgcna_hlca/`、`hdwgcna_consensus/`、`hdwgcna_preservation/`、`hdwgcna_spatial/`）、DGCA（`dgca_gtex/`）、RandomWalkRestartMH 与 R 实现的网络融合（`multilayer_string/`）、R 实现的轨迹网络（`trajectory_paul15/`）的比较，每组有准备脚本、两边的运行脚本、比较脚本、原始结果和 README（含数据下载地址和校验值）；`overview.py` 生成总表。运行环境：`environment-r.yml`（R 4.4.3、WGCNA 1.74、Seurat 5.5.1、harmony 2.0.5 等；hdWGCNA 0.4.12、DGCA 1.0.3、RandomWalkRestartMH 不在 conda 上，按文件中的命令安装）、`requirements-python.txt`。脚本在各自的数据目录中运行，不含任何固定路径（Rscript 用 `RSCRIPT` 指定）。计时于 2026-11-17 功能冻结后在空闲机器上按 `RETIMING_PLAN.md` 统一测量。

### `website/`
项目网站 <https://www.biochen.org/gcen2/>（静态，Bootstrap 5.3.9，英文；设计与实现见 `GCEN2_Website_Design.md`、`GCEN2_Website_Implementation.md`）。`build.py` 把 `templates/`（Jinja2）与 `content.py` 渲染成 `_site/`，版本、benchmark 表、安装矩阵、命令数和图数从 `pyproject.toml`、`README.md` 和 CLI 读取；文档用 `sphinx_theme/gcen2/`（继承 basic，套上网站的导航栏和页脚）构建到 `_site/docs/`。`static/`：样式、脚本、logo、workflow 图和图库图片，以及自带的 Bootstrap（不用 CDN）。`gallery/make_gallery.py` 在示例数据上运行教程命令生成图库，`workflow/make_workflow.py` 画 workflow 图。`tests/check_site.py`（链接、图片、页面结构）和 `tests/visual_test.py`（Playwright：三种宽度 × 浅色/深色截图、溢出、控制台错误、对比度、交互）。

### `docs/requirements.txt`
文档构建的依赖：只有纯 Python 包，不编译 C++ 扩展（文档构建不需要）。文档不在 Read the Docs 上发布，由 `website/build.py` 构建到网站的 `docs/`，发布在 <https://www.biochen.org/gcen2/docs/>。

### `.gitignore`
忽略 `build/`、`*.so`、`*.egg-info`、`__pycache__` 等构建产物。

---

## 八、数据流图

### 典型 bulk 分析流程

```
gene_expr.tsv
     │
     ▼ read_tsv()
ExpressionMatrix
     │
     ▼ preprocessing.normalize(method="tmm")
ExpressionMatrix (normalized)
     │
     ▼ preprocessing.filter_genes(percent_mean=0.75)
ExpressionMatrix (filtered)
     │
     ▼ network.build(method="spearman", cor=0.8, pval=0.001)
     │   └── C++: _gcen_core.correlation.network_edges()
Network (edge list DataFrame)
     │
     ├──▼ module.identify(similarity=0.5)
     │    └── C++: _gcen_core.module.jaccard_modules()
     │   ModuleSet
     │
     ├──▼ annotation.annotate_go(obo, assoc, pval=0.05)
     │    └── C++: _gcen_core.enrichment.hypergeom_sf()
     │   dict[gene → EnrichmentResult]
     │
     └──▼ rwr.predict(seed_genes, gamma=0.5)
          └── C++: _gcen_core.rwr.random_walk_restart()
         RWRResult
              │
              ▼ report.generate()
           HTML 报告
```

### 空间转录组分析流程

```
visium_data.h5ad（含 obsm['spatial']）
     │
     ▼ gcen2.spatial.build_spatial_filtered(adata, n_neighbors=6, cor_cutoff=0.6, p_cutoff=0.05)
     │   ├── build_neighbor_graph(coords, n_neighbors=6)  → W (稀疏 CSR)
     │   ├── _network_build(expr_df, cor=0.6, pval=0.05)  → 初始网络
     │   └── _morans_test(X, W) → 去掉 I ≤ 0 或 FDR > moran_fdr 的基因的边
     ▼
Network（空间过滤后的共表达网络）

     或

     ▼ gcen2.spatial.build_spatial_weighted(adata, n_neighbors=6, bandwidth=2.0)
     │   ├── build_neighbor_graph(coords, n_neighbors=6)  → W
     │   ├── W_norm = diag(1/rowsums) @ (W + I)           → 行归一化（含自环）
     │   └── X_smooth = W_norm @ X                        → 空间平滑表达矩阵
     ▼
Network（基于空间平滑表达的共表达网络）

     或

     ▼ gcen2.spatial.build_spatial_regions(adata, region_key="tissue_layer")
     │   ├── 按 adata.obs["tissue_layer"] 切分 spots
     │   └── 每个区域独立调用 _network_build()
     ▼
dict[region_name → Network]（每个空间区域一个网络）
```
