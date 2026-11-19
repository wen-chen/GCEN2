"""The compiled C++ extension ``_gcen_core`` (pybind11).

Submodules: ``correlation`` (correlation kernels, network building, TOM),
``enrichment`` (hypergeometric tests), ``module`` (Jaccard modules,
hierarchical clustering, edge-file adjacency), ``network`` (statistics, edge
files: merge, extract, shuffle, diff, layers) and ``rwr``.

Code that can do without the extension calls :func:`extension` (``None``
when it is not built); code that needs it calls :func:`require_extension`.
Both read the module attribute ``_gcen_core`` at call time, so tests can
disable the extension with ``monkeypatch.setattr(gcen2.core, "_gcen_core",
None)``.
"""

from __future__ import annotations

try:
    from . import _gcen_core
except ImportError:
    _gcen_core = None  # type: ignore[assignment]

_HAVE_CORE = _gcen_core is not None


def extension():
    """The compiled extension, or ``None`` when it is not available."""
    return globals().get("_gcen_core")


def require_extension():
    """The compiled extension; :class:`gcen2.GCENInputError` when it is not built."""
    ext = extension()
    if ext is None:
        from gcen2.exceptions import GCENInputError

        raise GCENInputError(
            "C++ extension '_gcen_core' is not built. "
            "Please install gcen2 with: pip install gcen2 (which compiles C++ extensions).")
    return ext


__all__ = ["_gcen_core", "_HAVE_CORE", "extension", "require_extension"]
