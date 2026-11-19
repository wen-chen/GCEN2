# GCEN2

**A Python Toolkit for Gene Co-Expression Network Analysis with Single-Cell Resolution and Interactive Visualization**

---

**Implementation Document**

| | |
|---|---|
| Document Version | 3.0 |
| Date | November 2026 (2.0: November 2026; 1.0: March 2026) |
| Author | Wen Harold Chen |
| Affiliation | Hunan University of Chinese Medicine |
| Contact | chenwen@biochen.org |

---

How GCEN2 is built and kept correct: the build, the rules the C++ and Python
code follow, the test system, the pitfalls met, and the steps for adding
code. What the software must do is in the Requirements Specification, how it
is designed (modules, algorithms) in the Design Document, where each piece of
code lives in the Code Guide; the current state in `docs/status.md`, test
results in `docs/test_report.md`, the history in `docs/development_log.md`.

## Contents

- [1. Build](#1-build)
- [2. C++ layer](#2-c-layer)
- [3. Python layer](#3-python-layer)
- [4. Tests and continuous integration](#4-tests-and-continuous-integration)
- [5. Release](#5-release)
- [6. Pitfalls](#6-pitfalls)
- [7. Adding code](#7-adding-code)

---

## 1. Build

**Toolchain.** Python 3.9 to 3.14, a C++17 compiler (built and tested with
GCC 15; GCC 7, Clang 5 and MSVC 2017 are the stated minimums, not yet
verified), CMake >= 3.18, scikit-build-core and pybind11. The extension
`_gcen_core` is one shared library built from `src/` (kernels) and
`src/bindings/` (one binding file per submodule: `correlation`,
`enrichment`, `module`, `network`, `rwr`), installed into `gcen2/core/`.

**Development install** (setup in `CONTRIBUTING.md`):

```bash
pip install --no-build-isolation -e ".[dev]"
```

scikit-build-core (`pyproject.toml`: `build-dir = "build"`,
`editable.mode = "redirect"`, `editable.rebuild = false`) compiles into
`build/`, installs the extension into `gcen2/core/` of site-packages and
redirects the Python files to the repository. Run the command again after a
change to `src/`, and after adding or moving Python modules or packages (the
redirect mapping is made at install time).

**Compiler settings** (`CMakeLists.txt`):

- `-O3 -Wall -Wextra -ffp-contract=off` (MSVC `/O2 /W3`): no `a*b+c` is fused
  into an FMA anywhere, so results do not depend on the target CPU or
  compiler defaults (AArch64 GCC and `-march=native` contract by default).
- `GCEN2_NATIVE=ON` adds `-march=native` for local builds only.
- On x86 targets a second build of the TOM product kernel with `-mavx`
  (`src/tom_gemm_avx.cpp`) is chosen at run time when the CPU has AVX; it
  gives the same results as the portable kernel (SSE2 / NEON / plain C++).
  `GCEN2_SIMD=baseline` forces the portable kernel (the CI runs the TOM tests
  both ways).
- No third-party C++ library: the kernels use the standard library only
  (threads are `std::thread`) besides pybind11.

**Dependency groups** (`pyproject.toml`; lower bounds are the oldest versions
the tests were run with):

| Group | Packages |
|---|---|
| core | numpy >= 1.24, pandas >= 2.0, scipy >= 1.10, threadpoolctl >= 3.0, typer >= 0.9 |
| `single-cell` | anndata >= 0.10, scanpy >= 1.9 |
| `visualization` | matplotlib >= 3.7, seaborn, pyvis, networkx |
| `report` | Jinja2, weasyprint |
| `all` | all of the above (`cli` is kept as an empty group) |
| `dev` | pytest, pytest-benchmark, pytest-cov, sphinx, sphinx-rtd-theme, ruff, pre-commit |

The documentation builds without the extension (`docs/requirements.txt`,
`.readthedocs.yaml`). The source distribution includes the example data
(`sdist.include`), so its tests run.

---

## 2. C++ layer

**Bindings.** `src/bindings/module.cpp` creates the submodules; each
`src/bindings/<name>.cpp` validates its arguments (shapes, CSR structure:
`binding_util.hpp`), copies what it needs out of Python objects, releases
the GIL around the computation and builds the results after re-acquiring it.
Never use `py::call_guard<py::gil_scoped_release>()` on a function that
reads a numpy buffer: the buffer would be read without the GIL.

**Threads.** Kernels run on `std::thread`s (`WorkerPool` in
`correlation.hpp`, started once per call and reused across blocks;
`util/parallel.hpp` for simple parallel loops). Work is handed out in small
units; each thread collects its results locally and they are merged and then
sorted, so the output never depends on the number of threads or their timing.

**Determinism.** Besides the threads and no FP contraction: fixed summation
orders in the matrix-product kernels; numpy's PCG64 stream reproduced in C++
(`pcg64.hpp`, for `shuffle`), R's Mersenne Twister in Python
(`gcen2/_rrng.py`, for preservation and consensus calibration); floats
written as Python's `repr` (`py_repr.hpp`) or GCEN v1's `%f` / `%e`
(`edge_format.hpp`).

**Files.** Edge lists and matrices are read in large chunks
(`edge_text::for_each_line`) with pandas' rules for quoting and missing
values (`edge_text.hpp`), so a C++ reader gives the same table as
`pandas.read_csv`; when a file is outside what the C++ reader is sure of, it
returns `ok = false` and Python falls back to pandas.

**Reference implementations.** Every kernel has a Python / numpy version
used when the extension is missing (`gcen2.core.extension()` returns None)
and compared with it in the tests, byte for byte where files are written.

---

## 3. Python layer

**Light imports.** `import gcen2` loads submodules on first use (PEP 562).
`import gcen2.cli` imports neither numpy nor pandas; the commands that stream
files (`build` on a plain matrix, `module -i`, `stats`, `merge`, `extract`,
`shuffle`, `rwr`, `annotate`) read and write them through the extension and
`gcen2/_matrix_io.py` without pandas (`tests` check this in fresh
interpreters).

**Extension access.** Through `gcen2.core.extension()` (None when not built)
or `require_extension()` (`GCENInputError` when not built); both read the
module attribute at call time, so tests disable the extension with
`monkeypatch.setattr(gcen2.core, "_gcen_core", None)`.

**Errors.** `GCENError` and its subclasses (`GCENInputError` for bad input or
parameters, `GCENIOError` for files, `GCENComputeError`). Library code raises
them; the CLI wrapper (`gcen2/cli/_common.wrap_command`) prints them as
`ERROR: ...` with exit code 1 and user warnings as `WARNING: ...`.

**Determinism.** Every sort has an explicit tie-break (e.g. RWR ranks ties by
gene name), so results do not depend on the numpy / pandas version; random
steps take a seed.

**CLI annotations.** typer evaluates the annotations of command functions:
they use `Optional[X]`, not `X | None` (fails on Python 3.9).

---

## 4. Tests and continuous integration

How the suite is organised (the results are in `docs/test_report.md`):

- **Golden outputs** (`tests/golden/`): every `gcen2` command runs on small
  fixed inputs and its outputs must equal the frozen files (text exactly,
  numbers to a relative 1e-9; the same with 1 and 4 threads). This is the
  baseline for refactoring and performance work;
  `python -m tests.golden.make_golden` regenerates it only for an intended
  change of results, stated in the commit.
- **References**: GCEN v1 outputs (`tests/data/v1_annotate`, and the full v1
  network when `data/sample` is present), WGCNA, hdWGCNA and
  RandomWalkRestartMH results made by the R scripts in `tests/data/*/make_*.R`.
- **Two implementations agree**: each C++ path against its Python reference.
- **Completeness**: a new command fails the suite until it has a golden step,
  a new figure command until it is in `test_plot_commands.py`; `docs/cli.rst`
  must match the commands; every subpackage must be reachable as
  `gcen2.<name>`; every command reports a missing input as `ERROR:`.
- Tests that need the extension, anndata or matplotlib skip themselves when
  it is missing; nothing else is skipped.

**CI** (`.github/workflows/ci.yml`; not run yet, the repository has not been
pushed):

| Job | Platform | Python | Content |
|---|---|---|---|
| test | Ubuntu 22.04 | 3.9 to 3.14 | all tests; TOM tests again with `GCEN2_SIMD=baseline`; coverage on 3.12 |
| test | macOS 14, macOS 15 (ARM64), Windows Server 2022 | 3.12 | all tests |
| minimum versions | Ubuntu 22.04 | 3.9 | numpy 1.24, pandas 2.0, scipy 1.10, anndata 0.10, matplotlib 3.7, typer 0.9 |
| docs | Ubuntu 22.04 | 3.12 | `sphinx-build -W` |
| lint | Ubuntu 22.04 | 3.12 | `ruff check` |
| benchmark | Ubuntu 22.04 | 3.12 | pytest-benchmark (pushes to main) |
| wheels / publish | Linux, macOS, Windows | cp39 to cp314 | cibuildwheel, tests on the wheels; PyPI on a `v*` tag |

---

## 5. Release

1. Feature freeze; re-time the benchmarks on an idle machine
   (`benchmarks/RETIMING_PLAN.md`; done on 2026-11-17).
2. Update `CHANGELOG.md` and the version (`pyproject.toml`, `docs/conf.py`,
   `CITATION.cff`); tag `v2.0.0`.
3. The CI builds the sdist and wheels and uploads them to PyPI.
4. Fill the sha256 of the bioconda recipe (`conda-recipe/meta.yaml`) and
   submit it.
5. GitHub release with the changelog; import the repository on Read the Docs.

---

## 6. Pitfalls

| Pitfall | What to do |
|---|---|
| Releasing the GIL before reading a numpy buffer (`py::call_guard`) | Read the buffer first, release the GIL in a scope |
| FMA contraction changes the last bits between machines | Keep `-ffp-contract=off`; check new kernels against the numpy version |
| Results depend on the thread count | Merge per-thread results, then sort |
| Ties sorted differently by numpy / pandas versions | Sort with an explicit tie-break |
| `X \| None` in a command's annotation breaks the CLI on Python 3.9 | `Optional[X]` |
| typer < 0.12 passes an empty tuple for a list option left out | The CLI wrapper turns it into None |
| pandas 2 and 3 parse and format the last digit differently | Compare with pandas to 1e-12, not byte for byte, under pandas 2 |
| matplotlib's backend probe can hang on a dead X display | `gcen2.plot` uses Agg outside interactive sessions |
| A new module or package is missing after an editable install | Run the install again |
| FDR on a subset of the edges | Adjust over all `n(n - 1)/2` pairs (the kernels keep every p-value <= cutoff) |

---

## 7. Adding code

**A C++ function.** Declare it in `src/<name>.hpp`, implement it in
`src/<name>.cpp` (add new files to `CMakeLists.txt`), bind it in
`src/bindings/<submodule>.cpp`, rebuild, and test it against a Python
reference.

**A Python module.** Add it to `gcen2/`, to `_SUBMODULES` and `__all__` in
`gcen2/__init__.py` if public, to `docs/api.rst`, and a test file.

**A command.** Add it to the matching `gcen2/cli/<group>.py` (raise
`GCENInputError` for bad input; import heavy libraries inside the function),
a golden step in `tests/golden/cases.py`, and run `python
tools/make_cli_docs.py`. A figure command: a `_plot(...)` declaration in
`gcen2/cli/plot.py` and a case in `tests/test_plot_commands.py`.

Then a line in `CHANGELOG.md`; the checks to run before a commit are in
`CONTRIBUTING.md`.

---

*End of Document*
