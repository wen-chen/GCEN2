"""Helpers shared by the gcen2 commands.

Every command is wrapped by :func:`wrap_command`: an error of GCEN2
(:class:`gcen2.GCENError`) is printed as ``ERROR: <message>`` (exit code 1)
and warnings as ``WARNING: <message>``, so the commands themselves only
raise. Heavy libraries are imported inside the commands: importing the CLI
must not import numpy or pandas (``build``, ``stats``, ``merge``, ... start
in a fraction of a second without them).
"""

from __future__ import annotations

from pathlib import Path
from typing import NoReturn, Optional

import typer

app = typer.Typer(
    name="gcen2",
    help="GCEN2: Gene Co-Expression Network Analysis Toolkit",
    no_args_is_help=True,
    add_completion=True,
)

plot_app = typer.Typer(
    help="Publication figures from GCEN2 result files (PDF, PNG and SVG; journal widths, "
         "editable text).", no_args_is_help=True)

THREADS_HELP = "Number of threads."
PRESETS = ("wgcna", "hdwgcna")


def fail(message: str) -> NoReturn:
    """Print ``ERROR: message`` and exit with status 1."""
    typer.echo(f"ERROR: {message}", err=True)
    raise typer.Exit(1)


def read_h5ad(path: Path, backed: bool = False):
    """An AnnData file (``backed``: kept on disk, read-only)."""
    try:
        import anndata as ad
    except ImportError:
        fail("anndata is required: pip install gcen2[single-cell]")
    if not Path(path).is_file():
        fail(f"file not found: {path}")
    return ad.read_h5ad(path, backed="r") if backed else ad.read_h5ad(path)


def check_groups(groupby, groups) -> None:
    """``--groupby`` and ``--group`` are given together or not at all."""
    if (groupby is None) != (not groups):
        fail("--groupby and --group go together.")


def h5ad_only(**options) -> None:
    """Fail when an option that applies to .h5ad input only was given."""
    if any(v is not None and v is not False and v not in ([], ()) for v in options.values()):
        names = [f"--{n}" for n in options]
        fail(f"{', '.join(names[:-1])} and {names[-1]} need an .h5ad input.")


def module_options(preset: str, defaults: dict, **given) -> dict:
    """Module-building options: the ``hdwgcna`` preset's ``defaults`` (none for
    ``wgcna``) with the options given (not None) on top; ``signed=True``
    adds a signed network to either preset."""
    if preset not in PRESETS:
        fail("--preset must be wgcna or hdwgcna.")
    opts = dict(defaults) if preset == "hdwgcna" else {}
    if given.pop("signed", False):
        opts["signed"] = True
    opts.update({k: v for k, v in given.items() if v is not None})
    return opts


def echo_tom(tom: Optional[dict]) -> None:
    """The soft power of a TOM and, when it was picked, the scale-free fit table."""
    if tom is None:
        return
    if not tom["fits"]:
        typer.echo(f"TOM power {tom['power']:g} (given)")
        return
    typer.echo("power  R^2     slope   mean(k)")
    for f in tom["fits"]:
        mark = " <-" if f["power"] == tom["power"] and tom["estimated"] else ""
        typer.echo(f"{f['power']:5g}  {f['r2']:.3f}  {f['slope']:6.2f}  {f['mean_k']:.4g}{mark}")
    if tom["estimated"] and tom.get("rule"):
        typer.echo(f"TOM power {tom['power']:g} ({tom['rule']})")
    elif tom["estimated"]:
        typer.echo(f"TOM power {tom['power']:g} (first with scale-free fit R^2 > 0.85)")
    else:
        typer.echo(f"TOM power {tom['power']:g} (no power met the R^2 cutoff: WGCNA FAQ "
                   f"default for this number of samples)")


def fmt_power(p) -> str:
    return p if isinstance(p, str) else f"{p:g}"


def read_module_table(path: Path):
    """gene -> module of a modules table: GCEN2's (module, gene; modules
    0..K-1, named ME0.. here) or any TSV with gene / gene_name and module
    columns (hdWGCNA's GetModules)."""
    from gcen2.plot._io import modules

    return modules(path).rename(None)


def wrap_command(callback):
    """Wrap a command callback:

    - errors of GCEN2 become ``ERROR: ...`` (exit 1), user warnings
      ``WARNING: ...`` on stderr;
    - the parent directories of its output files are created first: the
      ``output`` option and the options whose help starts with "write" /
      "also write" (``-e: write ...``);
    - list options left out get their default None: typer < 0.12 passes an
      empty tuple (``--group`` not given would mean "no group" instead of "all").
    """
    import functools
    import inspect
    import re

    from typer.models import OptionInfo

    writes = re.compile(r"^(?:-\w+: )?(?:also )?write\b", re.IGNORECASE)  # not "--write-x"
    params = inspect.signature(callback).parameters.items()
    outs = [n for n, p in params
            if isinstance(p.default, OptionInfo)
            and (n == "output" or writes.match(p.default.help or ""))]
    none_lists = [n for n, p in params
                  if isinstance(p.default, OptionInfo) and p.default.default is None]

    @functools.wraps(callback)
    def wrapper(*args, **kwargs):
        import warnings

        from gcen2.exceptions import GCENError

        for n in none_lists:
            if isinstance(kwargs.get(n), (tuple, list)) and not kwargs[n]:
                kwargs[n] = None
        for n in outs:
            v = kwargs.get(n)
            if isinstance(v, (str, Path)) and str(v):
                Path(v).parent.mkdir(parents=True, exist_ok=True)
        with warnings.catch_warnings(record=True) as caught:
            warnings.simplefilter("always")
            try:
                return callback(*args, **kwargs)
            except GCENError as exc:
                error = str(exc)
            finally:  # the warnings, also before an error
                for w in caught:
                    if issubclass(w.category, UserWarning):
                        typer.echo(f"WARNING: {w.message}", err=True)
        fail(error)

    return wrapper
