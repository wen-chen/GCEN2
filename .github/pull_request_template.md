## What and why

<!-- What this changes and why. Link the issue, if any. -->

## Checklist

See [CONTRIBUTING.md](../CONTRIBUTING.md).

- [ ] `python -m pytest -q` passes
- [ ] `ruff check gcen2 tests tools benchmarks` passes
- [ ] `python tools/check_doc_consistency.py` passes
- [ ] Results unchanged (golden outputs pass), or the change of results is intended and
      explained here and in the commit message (`python -m tests.golden.make_golden`)
- [ ] A new command has a golden step (`tests/golden/cases.py`); a new plot command a case
      in `tests/test_plot_commands.py`
- [ ] A changed command or option: `python tools/make_cli_docs.py` run
- [ ] Docstrings and documentation updated; a line in `CHANGELOG.md` for a user-visible change
