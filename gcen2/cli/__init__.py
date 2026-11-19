"""Command-line interface for GCEN2 (backward compatible with GCEN v1).

One module per group of commands, registered in this order (``gcen2
--help`` lists them so):

- :mod:`~gcen2.cli.network`: the GCEN v1 commands (normalize, filter, build,
  module, annotate, enrich, rwr, merge, extract, shuffle, stats, diff);
- :mod:`~gcen2.cli.modules`: WGCNA modules (consensus-modules,
  module-preservation, project-modules; ``module -e``);
- :mod:`~gcen2.cli.single_cell`: celltype, celltype-modules;
- :mod:`~gcen2.cli.spatial`: spatial-modules, spatial-build;
- :mod:`~gcen2.cli.trajectory`, :mod:`~gcen2.cli.multilayer`: trajectory,
  fuse, rwr-multiplex;
- :mod:`~gcen2.cli.report`: report;
- :mod:`~gcen2.cli.plot`: ``gcen2 plot ...``.

Helpers shared by the commands, and the error handling around each of
them, are in :mod:`gcen2.cli._common`.
"""

from __future__ import annotations

try:
    import typer  # noqa: F401
except ImportError as _err:
    raise ImportError(
        "typer is required for the CLI. Install with: pip install gcen2[cli]"
    ) from _err

import importlib  # noqa: E402

from gcen2.cli._common import app, plot_app, wrap_command  # noqa: E402

# each module registers its commands when imported, in this order
for _name in ("network", "modules", "single_cell", "spatial", "trajectory", "multilayer",
              "report", "plot"):
    importlib.import_module(f"gcen2.cli.{_name}")

app.add_typer(plot_app, name="plot")


@app.callback()
def main() -> None:
    """GCEN2 - Gene Co-Expression Network Analysis Toolkit."""


for _info in [*app.registered_commands, *plot_app.registered_commands]:
    _info.callback = wrap_command(_info.callback)

__all__ = ["app", "plot_app"]
