# GCEN2 documentation map

Each fact is kept in one place; the other documents link to it.

## For users (Sphinx site, `docs/*.rst`)

| Document | What it is for |
|---|---|
| `index.rst` | Start page of the site (Read the Docs) |
| `quickstart.rst` | Installation and a first analysis in Python and on the command line |
| `tutorials/` | One workflow per page on the example data: bulk, single cell, spatial, trajectory and multi-omics, figures and reports |
| `cli.rst` | Every command and option (generated: `python tools/make_cli_docs.py`) |
| `api.rst` | The Python API (autodoc from the docstrings) |
| `migration.rst` | From GCEN v1: command mapping, verified equivalence, deliberate differences |

Also for users: `../README.md` (overview), `../data/example/README.md` (example
data and their sources), `../benchmarks/README.md` (comparisons with GCEN v1,
WGCNA, hdWGCNA, DGCA and RandomWalkRestartMH).

## For developers (`docs/design/`, Markdown)

| Document | Question it answers |
|---|---|
| `GCEN2_Requirements_Specification.md` | What must GCEN2 do, and is each requirement met? |
| `GCEN2_Design_Document.md` | How is it designed: architecture, data models, modules, algorithms? |
| `GCEN2_Implementation_Document.md` | How is it built and kept correct: build, C++ rules, determinism, tests, pitfalls? |
| `GCEN2_Code_Guide.md` (Chinese) | Where is what: every file of the repository and what it holds |
| `GCEN2_Website_Design.md` | What is the website at biochen.org/gcen2 for, what does it contain, how does it look? |
| `GCEN2_Website_Implementation.md` | How is the website built (`website/`), tested and deployed? |

How to contribute (setup, tests, rules, commits): `../CONTRIBUTING.md`.

## Project record (Chinese)

| Document | Content |
|---|---|
| `status.md` | Current state and what is left |
| `development_log.md` | The development by phase: decisions and problems found |
| `test_report.md` | Functional test results, environments and coverage |
| `../CHANGELOG.md` | Every user-visible change (English) |
| `plans/archive/` | Earlier plans and status snapshots, kept unchanged |
