# Contributing to GCEN2

## Development setup

```bash
git clone https://github.com/wen-chen/GCEN2
cd GCEN2
python -m venv .venv && source .venv/bin/activate
pip install scikit-build-core pybind11 cmake ninja
pip install --no-build-isolation -e ".[dev]"
```

After a change to the C++ code (`src/`), run the last command again to rebuild
the extension (`gcen2.core._gcen_core`). A C++17 compiler is needed.

## Tests

```bash
python -m pytest -q                        # everything (1–2 min on a multi-core machine)
python -m pytest -q tests/test_golden.py   # golden outputs (seconds)
```

All tests must pass on every supported Python (3.9 to 3.14) and with the
oldest dependencies in `pyproject.toml` (the CI runs both). Tests that need
the compiled extension, anndata or matplotlib skip themselves when it is
missing; nothing else should be skipped.

Rules that the test suite enforces:

- **Golden outputs** (`tests/golden/`): every `gcen2` command runs on small
  fixed inputs and its outputs must equal the frozen ones in
  `tests/golden/expected/` (text exactly, numbers to a relative 1e-9). This
  is the baseline for performance work: an optimisation must leave it
  passing. Regenerate the expected outputs only when a change of results is
  intended, and say what changed and why in the commit message:

  ```bash
  python -m tests.golden.make_golden            # expected outputs
  python -m tests.golden.make_golden --inputs   # inputs too (rarely)
  ```

- **A new `gcen2` command** needs a step in `tests/golden/cases.py`; a new
  `gcen2 plot` command a case in `tests/test_plot_commands.py`. Both lists
  are checked against the registered commands.
- **CLI commands** live in `gcen2/cli/`, one module per group of commands
  (`network`, `modules`, `single_cell`, `spatial`, ...; figures in
  `plot.py`). A command raises `GCENInputError` / `GCENIOError` for bad
  input; the wrapper around every command prints it as `ERROR: ...`. Heavy
  libraries are imported inside the command, so `import gcen2.cli` stays
  free of numpy and pandas.
- **CLI annotations**: typer evaluates the annotations of command
  functions, so they use `Optional[X]`, not `X | None` (which fails on
  Python 3.9; checked by a test).
- **The compiled extension** is reached through `gcen2.core.extension()`
  (None when it is not built) or `require_extension()`; light file helpers
  (`require_file`, `tsv_field`, `write_modules`) are in `gcen2._matrix_io`.
- **Two implementations of one step** (the C++ kernel and the pandas
  reference path) must give the same output; tests compare them byte for
  byte (to 1e-12 with pandas 2, whose parser and float formatting differ in
  the last digit).
- **Results must not depend on the environment**: the same output with any
  number of threads and with every supported numpy / pandas version. Sort
  with an explicit tie-break (e.g. by gene name), never rely on how a sort
  algorithm orders ties.

## Code style

```bash
ruff check gcen2 tests tools benchmarks
```

The ruff configuration in `pyproject.toml` follows the code's conventions
(matrix names such as `X` and `W`, `typer.Option` defaults, `Optional` in CLI
annotations). Lines are kept to about 100 characters; the code is not
formatted with black.

## Documentation

- Every public function has a docstring (NumPy style); `docs/api.rst`
  documents the modules with autodoc.
- `docs/cli.rst` is generated from the CLI: run
  `python tools/make_cli_docs.py` after changing a command or option (a
  test checks that it is up to date).
- `python tools/check_doc_consistency.py` checks that the documents and the
  code agree (commands, modules, files); run it before committing.
- Build the documentation with `sphinx-build -W -b html docs docs/_build/html`
  (`pip install -r docs/requirements.txt` is enough; Read the Docs builds it
  that way, `.readthedocs.yaml`).
- Add a line to `CHANGELOG.md` for every user-visible change.

## Benchmarks

`benchmarks/` compares GCEN2 with GCEN v1, WGCNA, hdWGCNA, DGCA,
RandomWalkRestartMH and R implementations of its other analyses (scripts, raw
results and a README per comparison; `benchmarks/overview.py` builds the
summary table). The environments are in `benchmarks/environment-r.yml` and
`benchmarks/requirements-python.txt`, the downloads in each README; scripts
run from a data directory and never refer to a fixed path (`RSCRIPT`,
`GCEN2_STRING_GOLD`). Times are measured on an idle machine, after a feature
freeze, following `benchmarks/RETIMING_PLAN.md`.

## Example data

`data/example/` holds the example data (see its README for sources and
licences); `tools/make_example_data.py` rebuilds it. Files written for older
anndata versions use `portable()` in that script (object strings, not pandas
3's string dtype).

## Commits

One topic per commit, a subject such as `feat(plot): ...`, `fix: ...`,
`perf(build): ...`, `docs: ...`, `test: ...`, and a body that says what
changed and why, including any change of results.
