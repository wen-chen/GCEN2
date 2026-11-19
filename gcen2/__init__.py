"""
GCEN2: A Python Toolkit for Gene Co-Expression Network Analysis
with Single-Cell Resolution and Interactive Visualization.
"""

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version("gcen2")
except PackageNotFoundError:
    __version__ = "2.0.0.dev0"

__author__ = "Wen Harold Chen"
__email__ = "chenwen@biochen.org"
__license__ = "GPL-3.0-or-later"

# Submodules and the top-level helpers are loaded on first use (PEP 562), so
# ``import gcen2`` stays cheap: pandas, scipy and plotting libraries are only
# imported by the parts that need them. ``gcen2.network`` etc. work as before.
import importlib as _importlib

_SUBMODULES = {
    "preprocessing", "network", "module", "tom_modules", "consensus", "preservation",
    "annotation", "rwr", "differential", "single_cell", "harmony", "spatial", "trajectory",
    "multilayer", "plot", "plotting", "report", "io", "models", "config", "exceptions",
    "core",
}
_ATTRS = {
    "GCENConfig": "config",
    "GCENError": "exceptions",
    "GCENInputError": "exceptions",
    "GCENComputeError": "exceptions",
    "GCENIOError": "exceptions",
    "read_tsv": "io",
    "read_csv": "io",
    "read_h5ad": "io",
    "read_gene_list": "io",
    "read_modules": "io",
    "write_tsv": "io",
    "write_csv": "io",
}


def __getattr__(name):
    if name in _SUBMODULES:
        return _importlib.import_module(f"gcen2.{name}")
    if name in _ATTRS:
        value = getattr(_importlib.import_module(f"gcen2.{_ATTRS[name]}"), name)
        globals()[name] = value
        return value
    if name == "_HAS_CPP_EXTENSION":
        return _importlib.import_module("gcen2.core")._HAVE_CORE
    raise AttributeError(f"module 'gcen2' has no attribute {name!r}")


def __dir__():
    return sorted(set(globals()) | _SUBMODULES | set(_ATTRS))


__all__ = [
    "__version__",
    "preprocessing",
    "network",
    "module",
    "tom_modules",
    "consensus",
    "preservation",
    "annotation",
    "rwr",
    "differential",
    "single_cell",
    "harmony",
    "plot",
    "plotting",
    "report",
    "spatial",
    "trajectory",
    "multilayer",
    "GCENConfig",
    "GCENError",
    "GCENInputError",
    "GCENComputeError",
    "GCENIOError",
    "read_tsv",
    "read_gene_list",
    "read_modules",
    "read_csv",
    "read_h5ad",
    "write_tsv",
    "write_csv",
]
