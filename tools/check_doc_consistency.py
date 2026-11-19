#!/usr/bin/env python3
"""
GCEN2 文档-代码一致性检查器
用法: python tools/check_doc_consistency.py
输出: 所有差异条目，格式 [FAIL] / [PASS] / [WARN]
退出码: 0=全部通过, 1=存在差异
"""
from __future__ import annotations

import ast
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


def get_cli_commands() -> set[str]:
    """Commands registered with ``@app.command(...)`` in gcen2/cli/*.py."""
    cmds: set[str] = set()
    for p in sorted((ROOT / "gcen2/cli").glob("*.py")):
        cmds |= _commands_in(ast.parse(p.read_text(encoding="utf-8")))
    return cmds


def _commands_in(tree) -> set[str]:
    cmds: set[str] = set()
    for node in ast.walk(tree):
        if not isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef)):
            continue
        for deco in node.decorator_list:
            if not (isinstance(deco, ast.Call) and
                    isinstance(deco.func, ast.Attribute) and
                    deco.func.attr == "command"):
                continue
            explicit = None
            for kw in deco.keywords:
                if kw.arg == "name" and isinstance(kw.value, ast.Constant):
                    explicit = kw.value.value
            if explicit:
                cmds.add(explicit)
            else:
                fn = node.name
                name = fn.removeprefix("cmd_").replace("_", "-")
                cmds.add(name)
    return cmds


def check_c1_api() -> None:
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
            "network_interactive", "network_static", "dotplot", "barplot",
            "expression_boxplot", "expression_violin", "heatmap",
            "degree_distribution", "diff_network",
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
        ]),
        ("gcen2/spatial/core.py", [
            "build_spatial_filtered", "build_spatial_weighted", "build_spatial_regions",
        ]),
        ("gcen2/spatial/neighbors.py", ["build_neighbor_graph"]),
    ]
    for rel, doc_symbols in checks:
        actual = get_symbols(rel)
        for sym in doc_symbols:
            if sym in actual:
                ok(f"C1 {rel}: {sym}")
            else:
                fail(f"C1 {rel}: 文档声明 `{sym}` 但代码中不存在（实际有: {sorted(s for s in actual if not s.startswith('__'))}）")


def check_c2_cli() -> None:
    # Only check commands documented as already implemented (not "计划中")
    doc_commands = {
        "normalize", "filter", "build", "stats",
        "module", "annotate", "rwr", "enrich", "report",
        "merge", "extract", "shuffle", "diff", "spatial-build",
    }
    actual_cmds = get_cli_commands()
    for cmd in sorted(doc_commands):
        if cmd in actual_cmds:
            ok(f"C2 CLI: `gcen2 {cmd}`")
        else:
            fail(f"C2 CLI: 文档声明命令 `gcen2 {cmd}` 但 gcen2/cli/ 中未注册（实际: {sorted(actual_cmds)}）")


def check_c3_files() -> None:
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
        "gcen2/cli/__init__.py",
        "gcen2/cli/_common.py",
        "gcen2/cli/network.py",
        "gcen2/cli/plot.py",
        "gcen2/config.py",
        "gcen2/exceptions.py",
        "gcen2/core/__init__.py",
        "gcen2/single_cell/__init__.py",
        "gcen2/single_cell/core.py",
        "gcen2/spatial/__init__.py",
        "gcen2/spatial/neighbors.py",
        "gcen2/spatial/core.py",
        "tests/test_spatial.py",
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


def check_c4_sample_data() -> None:
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
    # data/sample (the full GCEN v1 sample package, 277 MB) is not in git; a
    # clone has data/example only, so the check applies where data/sample is
    if not (ROOT / "data" / "sample").is_dir():
        ok("C4 跳过：data/sample 不在本仓库副本中（未纳入 git 的 v1 完整样例包）")
        return
    for rel in files:
        p = ROOT / rel
        if p.exists():
            ok(f"C4 数据文件存在: {rel}")
        else:
            fail(f"C4 数据文件不存在: {rel}")


def check_c5_cpp() -> None:
    files = [
        "src/correlation.cpp",
        "src/correlation.hpp",
        "src/enrichment.cpp",
        "src/enrichment.hpp",
        "src/module_identify.cpp",
        "src/module_identify.hpp",
        "src/rwr.cpp",
        "src/rwr.hpp",
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


def check_c6_tests() -> None:
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
