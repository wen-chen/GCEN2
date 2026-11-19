# 文档-代码一致性检查与修复计划

> **For Claude:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task.

**Goal:** 编写一个自动化检查器，系统性发现 `docs/design/GCEN2_Code_Guide.md` 与实际代码之间的全部差异，然后逐条修复。

**Architecture:** 分三阶段：① 写检查器脚本（解析文档表格 + AST/文件系统核对）→ ② 运行检查器生成差异报告 → ③ 按报告逐条修复文档（以代码为准）。

**Tech Stack:** Python 3.12, ast, importlib, pathlib, re（解析 Markdown 表格）；无需新依赖。

---

## 检查范围说明

检查器需核对以下 6 类一致性：

| 编号 | 检查类型 | 来源（文档） | 核对方式 |
|------|---------|------------|---------|
| C1 | Python 类/函数存在性 | Code Guide 各模块表格 | `ast.parse` 提取实际符号 |
| C2 | CLI 命令注册 | Code Guide CLI 命令表 | 解析 `cli.py` 的 typer 注册 |
| C3 | 目录/文件存在性 | Code Guide 目录结构图 | `Path.exists()` |
| C4 | 样本数据文件 | Code Guide 数据文件表格 | `Path.exists()` |
| C5 | C++ 源文件 | Code Guide C++ 文件列表 | `Path.exists()` |
| C6 | 测试文件 | Code Guide 测试文件表格 | `Path.exists()` |

---

## Task 1：编写一致性检查器脚本

**Files:**
- Create: `tools/check_doc_consistency.py`

**Step 1: 创建 tools/ 目录并写检查器**

```python
#!/usr/bin/env python3
"""
GCEN2 文档-代码一致性检查器
用法: python tools/check_doc_consistency.py
输出: 所有差异条目，格式 [FAIL] / [PASS] / [WARN]
退出码: 0=全部通过, 1=存在差异
"""
from __future__ import annotations

import ast
import re
import sys
from pathlib import Path

ROOT = Path(__file__).parent.parent
GUIDE = ROOT / "docs/design/GCEN2_Code_Guide.md"

issues: list[str] = []
passed: list[str] = []


def fail(msg: str) -> None:
    issues.append(f"[FAIL] {msg}")


def warn(msg: str) -> None:
    issues.append(f"[WARN] {msg}")


def ok(msg: str) -> None:
    passed.append(f"[PASS] {msg}")


# ── 辅助：从 Python 文件提取所有顶层函数和类名 ──────────────────────────────
def get_symbols(rel_path: str) -> set[str]:
    p = ROOT / rel_path
    if not p.exists():
        return set()
    tree = ast.parse(p.read_text(encoding="utf-8"))
    return {
        n.name
        for n in ast.walk(tree)
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef))
    }


# ── 辅助：从 cli.py 提取实际注册的 typer 命令名 ─────────────────────────────
def get_cli_commands() -> set[str]:
    """
    解析 cli.py，提取 @app.command(name="xxx") 或 @app.command() 的函数名。
    """
    p = ROOT / "gcen2/cli.py"
    if not p.exists():
        return set()
    src = p.read_text(encoding="utf-8")
    tree = ast.parse(src)
    cmds: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for deco in node.decorator_list:
            # @app.command() or @app.command(name="xxx")
            if not (isinstance(deco, ast.Call) and
                    isinstance(deco.func, ast.Attribute) and
                    deco.func.attr == "command"):
                continue
            # Check for explicit name= keyword
            explicit = None
            for kw in deco.keywords:
                if kw.arg == "name" and isinstance(kw.value, ast.Constant):
                    explicit = kw.value.value
            if explicit:
                cmds.add(explicit)
            else:
                # derive name: cmd_normalize → normalize, cmd_annotate_go → annotate-go
                fn = node.name
                name = fn.removeprefix("cmd_").replace("_", "-")
                cmds.add(name)
    return cmds


# ── C1：Python API 符号检查 ──────────────────────────────────────────────────
def check_c1_api() -> None:
    """从 Code Guide 表格提取文档声明的符号，对比 ast 实际符号。"""
    # (模块路径, [文档中声明的符号名])
    checks: list[tuple[str, list[str]]] = [
        ("gcen2/models.py", [
            "ExpressionMatrix", "Network", "NetworkStats",
            "ModuleSet", "EnrichmentResult", "RWRResult", "DiffCoExResult",
        ]),
        ("gcen2/preprocessing.py", ["normalize", "filter_genes"]),
        ("gcen2/network.py", ["build", "stats", "merge", "extract", "shuffle"]),
        ("gcen2/module.py", ["identify"]),
        ("gcen2/annotation.py", [
            "GOTerm", "GOOntology", "load_go", "load_association",
            "enrich", "annotate_go", "annotate_kegg",
        ]),
        ("gcen2/rwr.py", ["predict"]),
        ("gcen2/differential.py", ["compare"]),
        ("gcen2/plotting.py", [
            "plot_network", "plot_module_heatmap", "plot_enrichment_bar",
            "plot_rwr_scores", "plot_degree_distribution",
        ]),
        ("gcen2/report.py", ["generate"]),
        ("gcen2/io.py", [
            "read_tsv", "read_csv", "read_h5ad", "read_network",
            "write_tsv", "write_csv",
        ]),
        ("gcen2/config.py", ["GCENConfig"]),
        ("gcen2/exceptions.py", [
            "GCENError", "GCENInputError", "GCENComputeError", "GCENIOError",
        ]),
        ("gcen2/single_cell/core.py", [
            "pseudobulk", "build_by_celltype", "build_dropout_aware",
            "_chunk_gene_indices", "_build_chunked",
        ]),
    ]
    for rel, doc_symbols in checks:
        actual = get_symbols(rel)
        for sym in doc_symbols:
            if sym in actual:
                ok(f"C1 {rel}: {sym}")
            else:
                fail(f"C1 {rel}: 文档声明 `{sym}` 但代码中不存在（实际有: {sorted(actual)}）")


# ── C2：CLI 命令检查 ──────────────────────────────────────────────────────────
def check_c2_cli() -> None:
    """Code Guide CLI 表格声明的命令是否在 typer app 中注册。"""
    doc_commands = {
        "normalize", "filter", "build", "stats", "merge",
        "extract", "shuffle", "module", "annotate-go", "annotate-kegg",
        "rwr", "diff", "report",
    }
    actual_cmds = get_cli_commands()
    for cmd in sorted(doc_commands):
        if cmd in actual_cmds:
            ok(f"C2 CLI: `gcen2 {cmd}`")
        else:
            fail(f"C2 CLI: 文档声明命令 `gcen2 {cmd}` 但 cli.py 中未注册（实际: {sorted(actual_cmds)}）")


# ── C3：项目关键文件/目录存在性 ───────────────────────────────────────────────
def check_c3_files() -> None:
    """Code Guide 目录结构图中列出的关键路径是否存在。"""
    paths = [
        "gcen2/__init__.py",
        "gcen2/models.py",
        "gcen2/preprocessing.py",
        "gcen2/network.py",
        "gcen2/module.py",
        "gcen2/annotation.py",
        "gcen2/rwr.py",
        "gcen2/differential.py",
        "gcen2/plotting.py",
        "gcen2/report.py",
        "gcen2/io.py",
        "gcen2/cli.py",
        "gcen2/config.py",
        "gcen2/exceptions.py",
        "gcen2/core/__init__.py",
        "gcen2/single_cell/__init__.py",
        "gcen2/single_cell/core.py",
        "gcen2/templates/report_default.html.jinja2",
        "CMakeLists.txt",
        "pyproject.toml",
        "docs/conf.py",
        "docs/index.rst",
        "docs/quickstart.rst",
        "docs/api.rst",
        "docs/migration.rst",
        "docs/Makefile",
        "docs/design/GCEN2_Code_Guide.md",
        "docs/design/GCEN2_Design_Document.md",
        "docs/design/GCEN2_Implementation_Document.md",
        "docs/design/GCEN2_Requirements_Specification.md",
        "conda-recipe/meta.yaml",
        "conda-recipe/build.sh",
        "conda-recipe/bld.bat",
        ".github/workflows/ci.yml",
    ]
    for rel in paths:
        p = ROOT / rel
        if p.exists():
            ok(f"C3 文件存在: {rel}")
        else:
            fail(f"C3 文件不存在: {rel}")


# ── C4：样本数据文件 ───────────────────────────────────────────────────────────
def check_c4_sample_data() -> None:
    """Code Guide 数据文件表格中列出的样本数据文件是否存在。"""
    files = [
        "data/sample/gene_expr.tsv",
        "data/sample/gene_expr.csv",
        "data/sample/gene_expr_norm.tsv",
        "data/sample/gene_expr_norm_filter.tsv",
        "data/sample/gene_co_expr.network",
        "data/sample/test_1.network",
        "data/sample/test_2.network",
        "data/sample/gene_go.assoc",
        "data/sample/K2ko.tsv",
        "data/sample/gene_kegg.assoc",
        "data/sample/module.txt",
        "data/sample/rwr_seed_genes.list",
        "data/sample/enrichment_gene.list",
        "data/sample/background_gene.list",
        "data/sample/gene_list.txt",
        "data/sample/rsem/rsem_sample.txt",
        "data/sample/stringtie/stringtie_sample.txt",
    ]
    for rel in files:
        p = ROOT / rel
        if p.exists():
            ok(f"C4 数据文件存在: {rel}")
        else:
            fail(f"C4 数据文件不存在: {rel}")


# ── C5：C++ 源文件 ────────────────────────────────────────────────────────────
def check_c5_cpp() -> None:
    """Code Guide C++ 部分列出的源文件是否存在。"""
    files = [
        "src/correlation.cpp",
        "src/correlation.hpp",
        "src/enrichment.cpp",
        "src/enrichment.hpp",
        "src/module_identify.cpp",
        "src/module_identify.hpp",
        "src/rwr.cpp",
        "src/rwr.hpp",
        "src/util/thread_pool.hpp",
        "src/util/strim.hpp",
        "src/bindings/module.cpp",
        "src/bindings/correlation.cpp",
        "src/bindings/enrichment.cpp",
        "src/bindings/module_identify.cpp",
        "src/bindings/rwr.cpp",
    ]
    for rel in files:
        p = ROOT / rel
        if p.exists():
            ok(f"C5 C++文件存在: {rel}")
        else:
            fail(f"C5 C++文件不存在: {rel}")


# ── C6：测试文件 ───────────────────────────────────────────────────────────────
def check_c6_tests() -> None:
    """Code Guide 测试章节列出的测试文件是否存在。"""
    files = [
        "tests/conftest.py",
        "tests/test_preprocessing.py",
        "tests/test_network.py",
        "tests/test_module.py",
        "tests/test_annotation.py",
        "tests/test_rwr.py",
        "tests/test_differential.py",
        "tests/test_correlation.py",
        "tests/test_enrichment.py",
        "tests/test_io.py",
        "tests/test_cli.py",
        "tests/test_plotting.py",
        "tests/test_report.py",
        "tests/test_single_cell.py",
        "tests/test_integration.py",
        "tests/test_benchmark.py",
    ]
    for rel in files:
        p = ROOT / rel
        if p.exists():
            ok(f"C6 测试文件存在: {rel}")
        else:
            fail(f"C6 测试文件不存在: {rel}")


# ── 主入口 ────────────────────────────────────────────────────────────────────
def main() -> int:
    print("=" * 60)
    print("GCEN2 文档-代码一致性检查器")
    print("=" * 60)

    check_c1_api()
    check_c2_cli()
    check_c3_files()
    check_c4_sample_data()
    check_c5_cpp()
    check_c6_tests()

    print(f"\n通过: {len(passed)} 项")
    print(f"差异: {len(issues)} 项\n")

    for line in issues:
        print(line)

    if issues:
        print("\n需要修复上述差异。")
        return 1
    else:
        print("全部通过，文档与代码一致。")
        return 0


if __name__ == "__main__":
    sys.exit(main())
```

**Step 2: 验证脚本语法**

```bash
python3 -c "import ast; ast.parse(open('tools/check_doc_consistency.py').read()); print('语法 OK')"
```

期望输出：`语法 OK`

---

## Task 2：运行检查器，收集差异报告

**Files:**
- Read: `tools/check_doc_consistency.py`（运行）
- Read: 差异报告输出

**Step 1: 运行检查器**

```bash
cd /mnt/d/code_base/gcen2
python3 tools/check_doc_consistency.py 2>&1 | tee /tmp/consistency_report.txt
```

**Step 2: 统计差异类别**

```bash
grep "^\[FAIL\]" /tmp/consistency_report.txt | sed 's/\[FAIL\] //' | cut -d' ' -f1 | sort | uniq -c
```

期望看到每类 FAIL 的数量，作为后续修复的工作量估算。

---

## Task 3：修复 C1 — plotting 函数名不符

**问题（已知）：** 文档声明 `plot_network`、`plot_module_heatmap`、`plot_enrichment_bar`、`plot_rwr_scores`、`plot_degree_distribution`，但代码实际是 `network_interactive`、`network_static`、`dotplot`、`barplot`、`heatmap`、`degree_distribution`、`diff_network`、`expression_boxplot`、`expression_violin`。

**修复策略：以代码为准，更新文档。**

**Files:**
- Modify: `docs/design/GCEN2_Code_Guide.md`（plotting 模块表格）

**Step 1: 读取 plotting.py 的实际公开函数**

```bash
python3 -c "
import ast
tree = ast.parse(open('gcen2/plotting.py').read())
fns = [n.name for n in ast.walk(tree)
       if isinstance(n, ast.FunctionDef) and not n.name.startswith('_')]
print(fns)
"
```

**Step 2: 同时读取每个函数的 docstring（第一行）**

```bash
python3 -c "
import ast
tree = ast.parse(open('gcen2/plotting.py').read())
for n in ast.walk(tree):
    if isinstance(n, ast.FunctionDef) and not n.name.startswith('_'):
        doc = ast.get_docstring(n) or ''
        print(f'{n.name}: {doc.splitlines()[0] if doc else \"(无 docstring)\"}')
"
```

**Step 3: 将 Code Guide 中 `gcen2/plotting.py` 表格替换为实际函数列表**

用 Edit 工具把旧表格替换为以实际代码为准的新表格（表格内容根据 Step 2 输出填写）。

**Step 4: 运行检查器验证 C1 plotting 项全部通过**

```bash
python3 tools/check_doc_consistency.py 2>&1 | grep "C1.*plotting"
```

期望：所有 plotting 相关项显示 `[PASS]`。

---

## Task 4：修复 C2 — CLI 命令缺失

**问题（已知）：** 文档声明 13 条命令，但 cli.py 中只有 `cmd_normalize`、`cmd_filter`、`cmd_build`、`cmd_module`、`cmd_annotate`、`cmd_rwr`、`cmd_enrich`、`cmd_stats`、`cmd_report` 共 9 个函数，且命令名推导后可能与文档声明不符。

**修复策略：以代码为准，更新文档 CLI 命令表。**

**Files:**
- Modify: `docs/design/GCEN2_Code_Guide.md`（CLI 命令表格）

**Step 1: 提取 cli.py 实际注册的命令及其说明**

```bash
python3 -c "
import ast
src = open('gcen2/cli.py').read()
tree = ast.parse(src)
for n in ast.walk(tree):
    if not isinstance(n, ast.FunctionDef):
        continue
    for deco in n.decorator_list:
        if (isinstance(deco, ast.Call) and
                isinstance(deco.func, ast.Attribute) and
                deco.func.attr == 'command'):
            explicit = next((kw.value.value for kw in deco.keywords
                             if kw.arg == 'name' and isinstance(kw.value, ast.Constant)), None)
            name = explicit or n.name.removeprefix('cmd_').replace('_', '-')
            doc = ast.get_docstring(n) or ''
            print(f'{name}: {doc.splitlines()[0] if doc else \"(无 docstring)\"}')
"
```

**Step 2: 将 Code Guide 中 CLI 命令表格替换为实际命令列表**

用 Edit 工具把旧表格替换为以实际代码为准的新表格。

**Step 3: 验证**

```bash
python3 tools/check_doc_consistency.py 2>&1 | grep "C2"
```

期望：所有 C2 项显示 `[PASS]`。

---

## Task 5：修复 C3/C4/C5/C6 — 文件存在性差异

**Files:**
- Modify: `docs/design/GCEN2_Code_Guide.md`（涉及不存在文件的条目）

**Step 1: 运行检查器只显示文件不存在的 FAIL**

```bash
python3 tools/check_doc_consistency.py 2>&1 | grep "\[FAIL\].*C[3456]"
```

**Step 2: 对每个 FAIL 判断处理方式**

- 若文件在文档中多余记录（代码未实现）→ 在文档中加 `（计划中）` 标注或删除条目
- 若文件应存在但缺失（代码 bug）→ 在文档中加 `（待实现）` 标注

**Step 3: 修改 Code Guide 相应表格**

用 Edit 工具按 Step 2 的判断逐条更新。

**Step 4: 验证**

```bash
python3 tools/check_doc_consistency.py 2>&1 | grep "\[FAIL\].*C[3456]"
```

期望：无输出（所有文件类检查通过）。

---

## Task 6：全量验证，生成最终报告

**Step 1: 运行完整检查器**

```bash
python3 tools/check_doc_consistency.py 2>&1
```

期望末尾输出：`全部通过，文档与代码一致。`（退出码 0）

**Step 2: 若仍有 FAIL，逐条修复直到清零**

每修复一项后，重新运行检查器确认该项变为 `[PASS]`。

**Step 3: 在 Code Guide 顶部更新版本号**

将 `> 版本：2026-03-15` 改为当天日期，并注明"已通过一致性检查器验证"。

---

## 预期差异汇总（已知）

执行前预判的主要差异：

| 编号 | 模块 | 差异描述 | 修复方向 |
|------|------|---------|---------|
| C1-plotting | `plotting.py` | 文档函数名与代码完全不符 | 以代码为准更新文档 |
| C2-cli | `cli.py` | 文档 13 命令，代码仅 9 个函数 | 以代码为准更新文档 |
| C3/C4/C5/C6 | 各文件 | 待运行检查器后确认 | 按实际情况标注 |
