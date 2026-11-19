# GCEN2 绘图模块规划

## 1. 目标

1. **能直接用于发表。** 默认输出符合期刊规范的矢量图：尺寸、字号、字体和配色都统一，文字在 Illustrator 或 Inkscape 中可以直接编辑。
2. **覆盖每一步分析。** 不只是 hdWGCNA 和 WGCNA 有的图要有，GCEN2 的每个分析命令只要结果能画，都要有对应的图。
3. **报告用同一套图。** `gcen2 report` 嵌入的图来自同一个绘图模块，与单独导出的图完全一样，且是矢量图。
4. **三种调用方式。**
   - Python API：每个函数返回 `Figure`，可以传入 `ax` 组合成多面板图；
   - 命令行：新增 `gcen2 plot <图名>`，直接读取各命令已写出的结果文件绘图；
   - 分析命令的 `--plot` 选项：分析完成时顺便出图。

## 2. 调研依据

- **期刊规范**（以 Nature 的作者指南为代表，其他期刊大同小异）：
  - 宽度：单栏 89 mm，1.5 栏 120–136 mm，双栏 183 mm；最大高度 170 mm。
  - 字号 5–7 pt，用无衬线字体。
  - 首选可编辑、未合并图层的矢量格式（PDF、EPS、SVG），文字不要转成轮廓。
  - 不要只靠颜色区分信息：用图例或标注，并照顾色盲读者。
- **hdWGCNA 0.4.12** 导出约 40 个绘图函数，与网络分析有关的有：
  - 软阈值和聚类树：`PlotSoftPowers`、`PlotDendrogram`；
  - 网络图：`ModuleNetworkPlot`（每个模块的 hub 基因网络，内圈 10 个、外圈 15 个）、`HubGeneNetworkPlot`（跨模块的 hub 网络）、`RunModuleUMAP` + `ModuleUMAPPlot`（在 TOM 上做基因 UMAP）；
  - 细胞层面：`ModuleFeaturePlot`、`ModuleRadarPlot`、`DotPlot`（特征基因在各组的平均值和表达比例）；
  - 模块关系：`ModuleCorrelogram`、`PlotKMEs`；
  - 下游分析：`PlotModuleTraitCorrelation`、`PlotDMEsVolcano`、`PlotDMEsLollipop`、`PlotModulePreservation`、`PlotModulePreservationLollipop`、`EnrichrBarPlot`、`EnrichrDotPlot`；
  - metacell：`DimPlotMetacells`；
  - 拓扑：`ModuleTopologyBarplot`、`ModuleTopologyHeatmap`。
  - 转录因子调控网络相关的图（`TFNetworkPlot` 等）依赖 GCEN2 没有的分析，不在本规划内。
- **WGCNA 教程**里的标准图：
  - 软阈值的无标度拟合与平均连接度（`pickSoftThreshold`）；
  - 聚类树加模块颜色条（`plotDendroAndColors`）；
  - 模块-性状热图（`labeledHeatmap`）；
  - 特征基因网络：特征基因的聚类树加相关热图（`plotEigengeneNetworks`）；
  - TOM 热图（`TOMplot`）；
  - 模块成员度与基因显著性的散点图（`verboseScatterplot`）；
  - 样本聚类树（`plotClusterTreeSamples`）。

参考资料：
- [Nature research figure guide](https://research-figure-guide.nature.com/figures/preparing-figures-our-specifications)
- [hdWGCNA network visualization vignette](https://rdrr.io/github/smorabit/hdWGCNA/f/vignettes/network_visualizations.Rmd)
- [hdWGCNA reference index](https://smorabit.github.io/hdWGCNA/reference/index.html)
- [hdWGCNA differential MEs](https://smorabit.github.io/hdWGCNA/articles/differential_MEs.html)
- [hdWGCNA module preservation vignette](https://rdrr.io/github/smorabit/hdWGCNA/f/vignettes/module_preservation.Rmd)
- [WGCNA FemaleLiver tutorial](https://rdrr.io/github/milescsmith/WGCNA/f/vignettes/FemaleLiver-02-networkConstr-auto.Rmd)

## 3. 设计原则

1. **统一样式**（`gcen2/plot/_style.py`）
   - 预设三种尺寸：`single`（89 mm）、`onehalf`（120 mm）、`double`（183 mm）；高度按图的类型取默认值。
   - 正文 7 pt，刻度 6 pt，最小 5 pt；无衬线字体，按 Arial、Helvetica、DejaVu Sans 的顺序取第一个可用的。
   - 线宽 0.5–0.75 pt；去掉上框线和右框线。
   - 配色：
     - 分类变量用 Okabe-Ito 色盲友好色板；
     - 连续变量用 viridis；
     - 有正负的量（相关系数、特征基因、log2 倍数变化）用发散色板，例如 RdBu_r，以 0 为中心，色标范围对称。
   - 用样式上下文（context manager）实现，不修改用户的全局 matplotlib 设置。
2. **模块颜色全局一致。**
   - 默认按 WGCNA 的 `standardColors` 取模块颜色（turquoise、blue、brown……），未分配基因为 grey。
   - 同一模块在所有图中颜色相同；hdWGCNA 用户看到的颜色也一致。
   - 图中一律同时标出模块名，例如"ME0 (turquoise)"，不让读者只靠颜色辨认。
3. **保存**（`save_figure`）
   - 默认同时写 PDF、PNG 和 SVG：PDF 为矢量、文字可编辑（`pdf.fonttype = 42`）；PNG 为 300 dpi；SVG 中的文字保留为文本（`svg.fonttype = none`）。
   - 可选 TIFF（600 dpi）。
   - 散点很多的图层（例如十万个细胞）只把这些点栅格化，坐标轴和文字仍保持矢量，文件不会过大。
4. **函数约定**
   - `plot_xxx(data, ..., ax=None, size="single", output=None) -> Figure`；
   - 输入直接用各命令写出的结果，可以是 DataFrame，也可以是文件路径；
   - 画多面板图时可以传入 `ax`；
   - 可选同时写出绘图所用的数据表（`--write-data`），方便作者核对或在其他软件里重画。
5. **依赖**
   - 必需：matplotlib。
   - 可选：networkx（网络布局）、umap-learn（基因 UMAP；没装时退回到 TOM 的多维尺度分析 MDS）。
   - 不引入 R 或 plotly。交互式网络沿用现有的 pyvis。
6. **兼容**：现有 `gcen2.plotting` 中的函数保留原名，改为调用新模块，旧代码照常可用。

## 4. 模块结构

```
gcen2/plot/
  __init__.py     # 公开函数；gcen2.plotting 保持兼容
  _style.py       # 尺寸、字号、色板、模块颜色、样式上下文
  _save.py        # PDF/PNG/SVG/TIFF 保存；大量散点栅格化
  _io.py          # 读取各命令的输出文件（modules.tsv、kme.tsv、eigengenes.tsv 等）
  qc.py           # 表达数据的质量控制
  network.py      # 网络统计、度分布、无标度拟合、网络图
  wgcna.py        # 软阈值、聚类树、特征基因、模块-性状、kME、TOM、hub 网络、基因 UMAP
  single_cell.py  # 细胞层面特征基因、点图、雷达图、DMEs、metacell、共识网络、投射、保守性
  spatial.py      # 空间特征基因图、Moran's I、区域相关
  enrichment.py   # GO/KEGG 富集条形图、点图、模块×通路热图
  function.py     # RWR、多层 RWR
  diff.py         # 差异共表达
  trajectory.py   # 拟时序网络动态
  multiomics.py   # 网络融合
gcen2/cli.py: gcen2 plot <图名> ...（子命令，读取结果文件）
```

## 5. 图的清单

优先级：P1 为核心图，最先做；P2 为各分析的结果图；P3 为补充图。

### 5.1 WGCNA 与 hdWGCNA 核心（`wgcna.py`），P1

| 图 | 对应 | 输入（GCEN2 输出） | 说明 |
|---|---|---|---|
| `soft_power` | PlotSoftPowers；pickSoftThreshold 图 | 软阈值拟合表（需要把各 power 的拟合结果写成文件） | 无标度拟合 R² 与平均、中位、最大连接度对 power 作图，标出所选 power 和阈值线；共识网络时每个数据集一条线 |
| `dendrogram` | PlotDendrogram；plotDendroAndColors | 层次聚类结果与模块（需要写出聚类树） | 基因聚类树加模块颜色条，可同时显示合并前和合并后两行 |
| `eigengene_network` | plotEigengeneNetworks；ModuleCorrelogram | eigengenes.tsv | 特征基因聚类树加相关热图（下三角，带数值），可以加上性状 |
| `module_trait` | labeledHeatmap；PlotModuleTraitCorrelation | trait_cor 表 | 相关系数热图，格内标 r 和显著性星号（或 p 值）；单细胞时按组分面 |
| `kme` | PlotKMEs | kme.tsv | 每个模块的 kME 排序曲线，标出前几个 hub 基因 |
| `mm_gs` | verboseScatterplot | kme.tsv 与基因显著性 | 模块成员度对基因显著性的散点图，标出 r 和 p |
| `hub_network` | ModuleNetworkPlot | TOM 或网络、kME | 单个模块的 hub 网络：内圈和外圈布局，节点大小表示 kME，边按 TOM 取最强的若干条 |
| `hub_network_all` | HubGeneNetworkPlot | 同上 | 各模块的 hub 基因合成一个网络，按模块着色 |
| `module_umap` | RunModuleUMAP + ModuleUMAPPlot | TOM、kME | 在 TOM 上做基因 UMAP：点按模块着色，大小表示 kME，hub 基因加标签 |
| `tom_heatmap` | TOMplot | TOM | 按聚类顺序排列的 TOM 热图，附模块颜色条；基因多时先抽样 |
| `module_sizes` | ModuleTopologyBarplot 一类 | modules.tsv | 各模块基因数的条形图 |
| `sample_tree` | plotClusterTreeSamples | 表达矩阵、性状 | 样本聚类树，下方附性状颜色条，用于找离群样本 |

### 5.2 单细胞与空间（`single_cell.py`、`spatial.py`），P1–P2

| 图 | 对应 | 输入 | 说明 |
|---|---|---|---|
| `me_embedding` | ModuleFeaturePlot | cell_eigengenes.tsv 与 obsm 中的 UMAP | 在细胞 UMAP 上画每个模块的特征基因，用发散色标；可选 harmonized 特征基因 |
| `me_dotplot` | DotPlot | 细胞特征基因与分组 | 行为模块、列为组：颜色表示平均特征基因，点大小表示大于 0 的细胞比例 |
| `me_violin` | VlnPlot | 同上 | 各组的特征基因小提琴图 |
| `me_radar` | ModuleRadarPlot | 同上 | 每个模块在各组的雷达图 |
| `dme_volcano` | PlotDMEsVolcano | dmes.tsv | 横轴为 log2 倍数变化，纵轴为 −log10 校正 p 值，点按模块着色并标注 |
| `dme_lollipop` | PlotDMEsLollipop | dmes.tsv | 棒棒糖图：点大小表示模块基因数，不显著的点画成 ×；按组分面 |
| `celltype_overview` | （新增） | celltypes.tsv | 每个细胞类型的模块数、基因数和 power 汇总 |
| `metacells` | DimPlotMetacells | metacell 成员与 UMAP | 在细胞 UMAP 上标出 metacell 中心 |
| `consensus_sets` | （新增） | sets.tsv、set_eigengenes.tsv | 各数据集的软阈值与校准；同一模块的特征基因在各数据集中的分布 |
| `preservation` | PlotModulePreservation | preservation_Z.tsv | Zsummary 与 medianRank 对模块大小作图，标出 Z = 2 和 Z = 10 两条阈值线 |
| `preservation_lollipop` | PlotModulePreservationLollipop | 同上 | 按 Zsummary 排序的棒棒糖图 |
| `projection` | ProjectModules 后的 FeaturePlot | projected_eigengenes.tsv | 投射后的特征基因画在查询数据的 UMAP 上 |
| `spatial_me` | （已有，重做） | spot_eigengenes.tsv 与坐标 | 空间特征基因图：多面板，色标对称，可按区域裁剪 |
| `spatial_group_cor` | （新增） | group_cor.tsv | 模块与区域的相关热图 |
| `moran` | （新增） | spatial-build 的 Moran's I 结果 | Moran's I 分布，并标出过滤阈值 |

### 5.3 GCEN2 的其他分析（P2）

| 命令 | 图 | 说明 |
|---|---|---|
| `normalize`、`filter` | `qc_boxplot`、`qc_density`、`mean_variance`、`sample_correlation` | 归一化前后的分布、均值-方差关系（标出过滤阈值）、样本相关热图 |
| `build`、`stats` | `degree_distribution`、`scale_free_fit`、`edge_weight_hist` | 度分布的双对数图和 R² 拟合，边权重分布 |
| `module` | `module_sizes`、`module_network` | Jaccard 法模块的大小分布；模块子网络图 |
| `annotate`、`enrich` | `enrich_bar`、`enrich_dot`、`enrich_heatmap` | 每个模块前 N 个通路的条形图或点图（对应 EnrichrBarPlot、EnrichrDotPlot），模块×通路的热图 |
| `rwr`、`rwr-multiplex` | `rwr_rank`、`rwr_subnetwork` | 候选基因得分排序图；种子基因与前若干个候选的子网络 |
| `extract`、`merge` | `subnetwork` | 子网络图；合并网络时各来源边数的维恩图或条形图 |
| `shuffle` | `null_distribution` | 真实统计量相对于随机网络零分布的位置 |
| `diff` | `diff_network`（重做）、`diff_volcano`、`diff_scatter` | 差异边的网络图；Δr 与显著性的火山图；两种条件下 r 的散点图 |
| `celltype` | `celltype_network_summary` | 各细胞类型网络的边数、密度，以及网络之间的重叠 |
| `trajectory` | `trajectory_dynamics`（重做）、`trajectory_heatmap` | 各时间窗内连接度或模块活性随拟时序的变化 |
| `fuse` | `fusion_weights` | 各数据层的权重；融合前后网络的重叠 |
| `spatial-build` | `moran`、`smoothing` | Moran's I 分布；平滑前后的对比 |

### 5.4 报告（P3，与 P1、P2 同步推进）

- 报告中的图全部改用新模块生成，以内嵌 SVG 的形式放入 HTML，矢量且文字可检索；PDF 报告随之也是矢量图。
- 每个章节自动选取有数据的图：
  - WGCNA 章节：软阈值、聚类树、特征基因网络、模块-性状、hub 网络；
  - 单细胞章节：特征基因的 UMAP 图和点图、DMEs、保守性；
  - 空间章节：空间特征基因图。
- 每张图下方附简短的图注，说明图中各元素的含义。
- 可选 `--figures-dir`：把报告里的每张图另存为 PDF 和 PNG，供发表直接使用。

## 6. 分析命令需要补充的输出

有些图需要的中间结果，目前还没有写出来：

1. **软阈值拟合表**：`module -e`、`celltype-modules`、`consensus-modules` 都已经计算了，只需写成 `soft_power.tsv`。
2. **基因聚类树**：hclust 的 merge 和 height，以及合并前的模块标签，写成 `dendrogram.tsv`，需加选项 `--write-tree`。
3. **TOM**：hub 网络、基因 UMAP 和 TOM 热图需要它。写出整个 TOM 太大，改为只写模块内部的 TOM（hdWGCNA 也是这样用的），格式为 npz 或 TSV 边表，需加选项 `--write-tom`。
4. **基因显著性**：模块成员度与基因显著性的散点图需要。在 `module -e --traits` 时一并计算 GS 和 p 值。
5. **细胞 UMAP 坐标**：从 h5ad 的 obsm 中读取，`gcen2 plot` 支持直接给 h5ad。

## 7. 实施步骤

每一步都包括：代码、单元测试（绘图冒烟测试，并检查图中元素的数量和标签）、文档，以及用真实基准数据生成示例图供人工检查。完成一步提交一次，只在本机提交。

| 步骤 | 内容 | 预计工作量 |
|---|---|---|
| 0 | 基础设施：`gcen2/plot` 包、样式、保存、模块颜色、IO、`gcen2 plot` 子命令框架、`gcen2.plotting` 兼容层；补充第 6 节的输出 | 1 个单元 |
| 1 | WGCNA 核心图（5.1 节的 12 张） | 2 个单元 |
| 2 | 单细胞与空间（5.2 节的 15 张），含基因 UMAP 与 hub 网络的布局 | 2 个单元 |
| 3 | 其他分析（5.3 节） | 2 个单元 |
| 4 | 报告：改用新模块的 SVG 图、自动选图、图注、`--figures-dir` | 1 个单元 |
| 5 | 示例图库：用 PBMC、HLCA、Visium、GTEx 的基准结果生成全部图，做成一个页面供审阅；README 加一节"绘图" | 1 个单元 |

验收标准：
- 每张图在单栏和双栏宽度下都没有文字重叠，字号不小于 5 pt；
- PDF 中的文字可以编辑；
- 转成灰度或模拟色盲后仍能分辨；
- 同一模块在所有图中颜色一致；
- 文档一致性检查和全部测试都通过。

## 8. 已确认的事项

1. **模块颜色**：默认用 WGCNA 的颜色名（`standardColors`：turquoise、blue……），并同时标出模块编号；另外提供 `--palette okabe-ito` 选项。
2. **默认输出格式**：PDF、PNG、SVG 三种都写；TIFF 用选项另加。
3. **图中文字**：一律用英文。

## 9. 进度

| 步骤 | 状态 | 提交 |
|---|---|---|
| 0 基础设施 | 完成 | 7126ba6 |
| 1 WGCNA 核心图（12 张）及所需输出（soft_power.tsv、dendrogram.json、tom_hubs.tsv、GS） | 完成 | e9a0f6b |
| 2 单细胞与空间（13 张） | 完成 | 451d13d、f531e86 |
| 3 其他分析（18 张） | 完成 | 8d25b81 |
| 4 报告：内嵌 SVG、图注、`--figures-dir` | 完成 | 0e08a12 |
| 5 示例图库与文档 | 完成 | 1f94e5f |
| 6 视觉检查：逐张修改排版 | 完成 | bd99318 |
| 7 补齐的三类图及所需输出 | 完成 | 16a2489 |
| 8 新图接入报告：细胞类型总览与 metacells（`--coords`）、merge / extract（`--merge-input`、`--merged`、`--extract-genes`）、Moran's I（`--moran`） | 完成 | 本次提交 |

共 43 个 `gcen2 plot` 命令，`gcen2.plot` 导出 52 个函数。

与计划相比的调整：
- **后补的图（阶段 7）**：
  - `metacells`（DimPlotMetacells）：`celltype-modules` 新写出 `metacell_members.tsv`（metacell、cell、center），图中每个 metacell 画在其细胞的平均位置，并给出每个细胞类型的 metacell 数和细胞覆盖率；
  - `moran`：`spatial-build --moran` 写出每个基因的 Moran's I（weighted 模式另有平滑后的 I），图中有 I 的分布、I 与平均表达、平滑前后对比；
  - `merge-overlap`、`extract-overview`：维恩图（2-3 个网络；更多网络时为 UpSet 图）直接从 merge / extract 的输入计算，不需要新的输出文件；merge 给出合并结果时，还标出每个区域通过阈值的边数。
- **TOM 热图**：从表达矩阵和 power 重新计算抽样基因的 TOM，不读 TOM 文件。
- **新增的两项措施**：
  - 无图形界面时自动用 Agg 后端：DISPLAY 指向无响应的 X 服务器时，matplotlib 自动选后端会卡住几分钟；
  - 输出可重复：SVG 元素 id 固定，SVG 和 PDF 不写创建时间。
