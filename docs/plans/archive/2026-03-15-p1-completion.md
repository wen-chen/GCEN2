# GCEN2 P1 Completion Plan

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 补全 GCEN2 距离可发布状态的全部 P1 缺口：4 条缺失 CLI 命令（merge/extract/shuffle/diff）、文档状态更新（3.2.5 分块计算已完成标记）、一致性检查器同步。

**Architecture:** CLI 命令严格遵循 `cli.py` 现有模式（typer 装饰器 + `gcen2.*` API 调用）；文档更新用 Edit 工具精准替换；检查器同步更新 C2 期望集。测试用 `typer.testing.CliRunner` 在 `tests/test_cli.py` 追加测试类。

**Tech Stack:** Python 3.12, typer, pytest, typer.testing.CliRunner, gcen2.network, gcen2.differential

---

## 背景知识（执行前必读）

### cli.py 命令模式
每个命令都是：`@app.command("命令名")` + `def cmd_xxx(...)` 函数体。
命令内统一用 `import gcen2` 延迟导入，用 `gcen2.read_network()` 或 `pd.read_csv` 读网络文件。

### 需要添加的 4 条命令对应 API

| CLI 命令 | Python API | 位置 |
|---------|-----------|------|
| `gcen2 merge` | `gcen2.network.merge(networks, cor_cutoff)` | `gcen2/network.py:108` |
| `gcen2 extract` | `gcen2.network.extract(network, genes)` | `gcen2/network.py:132` |
| `gcen2 shuffle` | `gcen2.network.shuffle(network, swap_factor, seed)` | `gcen2/network.py:137` |
| `gcen2 diff` | `gcen2.differential.compare(net_a, net_b, pval, n_a, n_b)` | `gcen2/differential.py:12` |

### 网络文件格式
GCEN v1/v2 网络文件是 TSV，列名经 `read_network()` 标准化为 `gene_a, gene_b, correlation, p_value`。
CLI 中读网络用 `gcen2.io.read_network(path)` 而不是裸 `pd.read_csv`（会自动处理 v1 列名）。

### DiffCoExResult 输出格式
`gcen2.differential.compare()` 返回 `DiffCoExResult`，其 `.results` 属性是 DataFrame，列：
`gene_a, gene_b, cor_a, cor_b, z_score, p_value, category`（category: A_only/B_only/shared_changed/shared_stable）。

---

## Task 1：CLI `gcen2 merge` 命令

**Files:**
- Modify: `gcen2/cli.py`（在 `cmd_stats` 前插入）
- Test: `tests/test_cli.py`（追加 `TestMergeCommand` 类）

### Step 1: 写失败测试

在 `tests/test_cli.py` 末尾追加：

```python
@requires_typer
class TestMergeCommand:
    def test_merge_two_networks(self, runner, tmp_path):
        # 创建两个小网络文件
        net1 = tmp_path / "net1.network"
        net2 = tmp_path / "net2.network"
        out = tmp_path / "merged.network"
        net1.write_text("gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t0.9\t0.001\nA\tC\t0.85\t0.001\n")
        net2.write_text("gene_a\tgene_b\tcorrelation\tp_value\nA\tB\t0.8\t0.001\nB\tC\t0.88\t0.001\n")
        result = runner.invoke(app, [
            "merge",
            "-i", str(net1), "-i", str(net2),
            "-o", str(out),
            "-c", "0.0",
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        df = pd.read_csv(out, sep="\t")
        assert "gene_a" in df.columns
        assert "gene_b" in df.columns
        assert len(df) > 0
```

### Step 2: 运行测试确认失败

```bash
cd /mnt/d/code_base/gcen2
python3 -m pytest tests/test_cli.py::TestMergeCommand -v --tb=short 2>&1 | tail -10
```

期望：`FAILED` 或 `ERROR`（命令未注册）

### Step 3: 实现命令

在 `gcen2/cli.py` 的 `cmd_stats` 函数定义（`@app.command("stats")`）**前**插入：

```python
@app.command("merge")
def cmd_merge(
    inputs: list[Path] = typer.Option(..., "-i", "--input",
                                      help="Input network files (repeat for each network)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output merged network TSV file."),
    cor_cutoff: float = typer.Option(0.5, "-c", "--cor-cutoff",
                                     help="Minimum absolute aggregated correlation."),
) -> None:
    """Merge two or more co-expression networks."""
    from gcen2.io import read_network
    import gcen2

    networks = [read_network(p) for p in inputs]
    merged = gcen2.network.merge(networks, cor_cutoff=cor_cutoff)
    merged.to_edgelist(output)
    typer.echo(
        f"Merged network: {merged.node_count} nodes, {merged.edge_count} edges → {output}"
    )
```

### Step 4: 运行测试确认通过

```bash
python3 -m pytest tests/test_cli.py::TestMergeCommand -v --tb=short 2>&1 | tail -10
```

期望：`PASSED`

---

## Task 2：CLI `gcen2 extract` 命令

**Files:**
- Modify: `gcen2/cli.py`
- Test: `tests/test_cli.py`（追加 `TestExtractCommand` 类）

### Step 1: 写失败测试

```python
@requires_typer
class TestExtractCommand:
    def test_extract_subnetwork(self, runner, tmp_path):
        net_file = tmp_path / "net.network"
        genes_file = tmp_path / "genes.txt"
        out = tmp_path / "sub.network"
        net_file.write_text(
            "gene_a\tgene_b\tcorrelation\tp_value\n"
            "A\tB\t0.9\t0.001\nA\tC\t0.85\t0.001\nD\tE\t0.88\t0.001\n"
        )
        genes_file.write_text("A\nB\n")
        result = runner.invoke(app, [
            "extract", "-i", str(net_file), "-g", str(genes_file), "-o", str(out)
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        df = pd.read_csv(out, sep="\t")
        # D-E edge should be excluded (neither D nor E in gene list)
        all_genes = set(df["gene_a"]) | set(df["gene_b"])
        assert "D" not in all_genes or "E" not in all_genes  # not both
```

### Step 2: 运行测试确认失败

```bash
python3 -m pytest tests/test_cli.py::TestExtractCommand -v --tb=short 2>&1 | tail -10
```

### Step 3: 实现命令

在 `cmd_merge` 后插入：

```python
@app.command("extract")
def cmd_extract(
    input: Path = typer.Option(..., "-i", "--input", help="Input network TSV file."),
    genes: Path = typer.Option(..., "-g", "--genes",
                               help="Gene list file (one gene per line)."),
    output: Path = typer.Option(..., "-o", "--output", help="Output subnetwork TSV file."),
) -> None:
    """Extract a subnetwork for a list of genes."""
    from gcen2.io import read_network
    import gcen2

    net = read_network(input)
    gene_list = genes.read_text().strip().splitlines()
    sub = gcen2.network.extract(net, gene_list)
    sub.to_edgelist(output)
    typer.echo(
        f"Subnetwork: {sub.node_count} nodes, {sub.edge_count} edges → {output}"
    )
```

### Step 4: 运行测试确认通过

```bash
python3 -m pytest tests/test_cli.py::TestExtractCommand -v --tb=short 2>&1 | tail -10
```

---

## Task 3：CLI `gcen2 shuffle` 命令

**Files:**
- Modify: `gcen2/cli.py`
- Test: `tests/test_cli.py`（追加 `TestShuffleCommand` 类）

### Step 1: 写失败测试

```python
@requires_typer
class TestShuffleCommand:
    def test_shuffle_preserves_edge_count(self, runner, tmp_path):
        net_file = tmp_path / "net.network"
        out = tmp_path / "shuffled.network"
        # 写一个有10条边的网络
        lines = ["gene_a\tgene_b\tcorrelation\tp_value"]
        genes = [f"G{i}" for i in range(10)]
        import itertools
        for a, b in list(itertools.combinations(genes, 2))[:10]:
            lines.append(f"{a}\t{b}\t0.9\t0.001")
        net_file.write_text("\n".join(lines) + "\n")
        result = runner.invoke(app, [
            "shuffle", "-i", str(net_file), "-o", str(out), "--swap-factor", "10"
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        orig = pd.read_csv(net_file, sep="\t")
        shuf = pd.read_csv(out, sep="\t")
        assert len(shuf) == len(orig)
```

### Step 2: 运行测试确认失败

```bash
python3 -m pytest tests/test_cli.py::TestShuffleCommand -v --tb=short 2>&1 | tail -10
```

### Step 3: 实现命令

在 `cmd_extract` 后插入：

```python
@app.command("shuffle")
def cmd_shuffle(
    input: Path = typer.Option(..., "-i", "--input", help="Input network TSV file."),
    output: Path = typer.Option(..., "-o", "--output", help="Output shuffled network TSV file."),
    swap_factor: int = typer.Option(100, "--swap-factor",
                                    help="Number of swaps = swap_factor × edge_count."),
    seed: int = typer.Option(42, "--seed", help="Random seed."),
) -> None:
    """Degree-preserving edge shuffle (null model generation)."""
    from gcen2.io import read_network
    import gcen2

    net = read_network(input)
    shuffled = gcen2.network.shuffle(net, swap_factor=swap_factor, seed=seed)
    shuffled.to_edgelist(output)
    typer.echo(
        f"Shuffled network: {shuffled.node_count} nodes, {shuffled.edge_count} edges → {output}"
    )
```

### Step 4: 运行测试确认通过

```bash
python3 -m pytest tests/test_cli.py::TestShuffleCommand -v --tb=short 2>&1 | tail -10
```

---

## Task 4：CLI `gcen2 diff` 命令

**Files:**
- Modify: `gcen2/cli.py`
- Test: `tests/test_cli.py`（追加 `TestDiffCommand` 类）

### Step 1: 写失败测试

```python
@requires_typer
class TestDiffCommand:
    def test_diff_two_networks(self, runner, tmp_path):
        net_a = tmp_path / "netA.network"
        net_b = tmp_path / "netB.network"
        out = tmp_path / "diff.tsv"
        net_a.write_text(
            "gene_a\tgene_b\tcorrelation\tp_value\n"
            "A\tB\t0.9\t0.001\nA\tC\t0.85\t0.001\n"
        )
        net_b.write_text(
            "gene_a\tgene_b\tcorrelation\tp_value\n"
            "A\tB\t0.5\t0.001\nB\tC\t0.88\t0.001\n"
        )
        result = runner.invoke(app, [
            "diff",
            "--network-a", str(net_a),
            "--network-b", str(net_b),
            "-o", str(out),
            "--n-a", "30", "--n-b", "30",
        ])
        assert result.exit_code == 0, result.output
        assert out.exists()
        import pandas as pd
        df = pd.read_csv(out, sep="\t")
        assert "category" in df.columns
        assert set(df["category"]).issubset(
            {"A_only", "B_only", "shared_changed", "shared_stable"}
        )
```

### Step 2: 运行测试确认失败

```bash
python3 -m pytest tests/test_cli.py::TestDiffCommand -v --tb=short 2>&1 | tail -10
```

### Step 3: 实现命令

在 `cmd_shuffle` 后插入：

```python
@app.command("diff")
def cmd_diff(
    network_a: Path = typer.Option(..., "--network-a", "-A",
                                   help="Network from condition A."),
    network_b: Path = typer.Option(..., "--network-b", "-B",
                                   help="Network from condition B."),
    output: Path = typer.Option(..., "-o", "--output",
                                help="Output differential results TSV file."),
    pval: float = typer.Option(0.05, "-p", "--pval", help="P-value cutoff."),
    n_a: int = typer.Option(30, "--n-a", help="Number of samples in condition A."),
    n_b: int = typer.Option(30, "--n-b", help="Number of samples in condition B."),
) -> None:
    """Differential co-expression analysis between two conditions."""
    from gcen2.io import read_network
    import gcen2

    net_a = read_network(network_a)
    net_b = read_network(network_b)
    result = gcen2.differential.compare(net_a, net_b, pval=pval, n_a=n_a, n_b=n_b)
    result.results.to_csv(output, sep="\t", index=False)

    counts = result.results["category"].value_counts().to_dict()
    typer.echo(f"Differential analysis complete → {output}")
    for cat, n in sorted(counts.items()):
        typer.echo(f"  {cat}: {n}")
```

### Step 4: 运行测试确认通过

```bash
python3 -m pytest tests/test_cli.py::TestDiffCommand -v --tb=short 2>&1 | tail -10
```

---

## Task 5：更新文档状态

**Files:**
- Modify: `docs/design/GCEN2_Requirements_Specification.md`
- Modify: `docs/design/GCEN2_Implementation_Document.md`
- Modify: `docs/design/GCEN2_Code_Guide.md`

### Step 1: 更新需求规格说明书

把 §3.2.5 从 `🔲` 改为 `✅`，因为 `_chunk_gene_indices` / `_build_chunked` / `max_genes_per_chunk` 已在 `single_cell/core.py` 实现。

找到并替换：
```
#### 3.2.5 Scalability Optimization [P1] 🔲
```
改为：
```
#### 3.2.5 Scalability Optimization [P1] ✅
```

把 §6.3 的 Remaining 部分：
```
- Chunked computation for >20k gene datasets (scalability optimization). 🔲
```
改为：
```
- Chunked computation for >20k gene datasets (scalability optimization). ✅
```

### Step 2: 更新实现文档

在 §1 Implementation Status Overview 表格中，`gcen2.single_cell` 行：
把 `Chunked computation pending` 改为 `Chunked computation implemented`

在 §7 Remaining Work 中：
把 `#### 7.2.1 Chunked Computation for Large Datasets` 整段标记为已完成（在标题后加 `✅ Done`），或删除该节。

### Step 3: 更新 Code Guide CLI 表格

在 `docs/design/GCEN2_Code_Guide.md` 的"已实现子命令"表格末尾加入新命令：

已实现表格补充：
```
| `gcen2 merge`   | `network_merge`  | 合并多个网络 |
| `gcen2 extract` | `network_extract`| 提取子网络   |
| `gcen2 shuffle` | `network_shuffle`| 随机化网络   |
| `gcen2 diff`    | （新增）          | 差异共表达分析 |
```

"计划中"表格对应行删除（这 4 条命令已实现）。

---

## Task 6：更新一致性检查器并全量验证

**Files:**
- Modify: `tools/check_doc_consistency.py`

### Step 1: 更新 C2 期望命令集

把 `check_c2_cli()` 中的 `doc_commands` 集合从：
```python
doc_commands = {
    "normalize", "filter", "build", "stats",
    "module", "annotate", "rwr", "enrich", "report",
}
```
改为（加入 4 条新命令）：
```python
doc_commands = {
    "normalize", "filter", "build", "stats",
    "module", "annotate", "rwr", "enrich", "report",
    "merge", "extract", "shuffle", "diff",
}
```

### Step 2: 运行一致性检查器

```bash
python3 tools/check_doc_consistency.py 2>&1
```

期望末尾：`差异: 0 项` / `全部通过，文档与代码一致。`

### Step 3: 运行完整测试套件

```bash
python3 -m pytest tests/ --ignore=tests/test_integration.py -q --tb=short 2>&1 | tail -5
```

期望：`passed`，`0 failed`，`0 errors`

---

## 完成标准

全部任务完成后须满足：
1. `python3 -m pytest tests/test_cli.py -v` — `TestMergeCommand`、`TestExtractCommand`、`TestShuffleCommand`、`TestDiffCommand` 全部 PASSED（或因 typer 未安装而 SKIPPED，但不 FAIL）
2. `python3 tools/check_doc_consistency.py` — `差异: 0 项`
3. `python3 -m pytest tests/ --ignore=tests/test_integration.py -q` — 0 failed，0 errors
4. 需求文档 §3.2.5 标记为 ✅
5. Code Guide CLI 表格"计划中"栏目清空，4 条新命令出现在"已实现"栏目
