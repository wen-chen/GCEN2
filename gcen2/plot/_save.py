"""Saving figures: vector PDF and SVG with editable text, 300 dpi PNG (and
600 dpi TIFF on request)."""

from __future__ import annotations

from pathlib import Path

from gcen2.plot._style import _STATE

_EXT = {"pdf", "png", "svg", "tiff", "tif", "eps"}
_DPI = {"png": 300, "tiff": 600, "tif": 600}


def save_figure(fig, output, formats=None, close: bool = False) -> list[Path]:
    """Save ``fig``. ``output`` with a known extension (.pdf, .png, .svg,
    .tiff, .eps) is written as given; without one it is a stem and the
    figure is written in every format of ``formats`` (default: the GCEN2
    defaults, PDF, PNG and SVG). Returns the paths written."""
    import matplotlib.pyplot as plt

    from gcen2.plot._style import rc_params

    out = Path(output)
    if out.suffix.lower().lstrip(".") in _EXT and formats is None:
        paths = [out]
    else:
        stem = out.with_suffix("") if out.suffix.lower().lstrip(".") in _EXT else out
        fmts = formats or _STATE["formats"]
        paths = [stem.with_name(f"{stem.name}.{f.lower().lstrip('.')}") for f in fmts]
    for p in paths:
        p.parent.mkdir(parents=True, exist_ok=True)
    import matplotlib as mpl

    with mpl.rc_context(rc_params()):
        for p in paths:
            ext = p.suffix.lower().lstrip(".")
            meta = _metadata(ext)
            # metadata only where there is some: TIFF / EPS take none before matplotlib 3.9
            fig.savefig(p, dpi=_DPI.get(ext, 300), bbox_inches="tight", pad_inches=0.02,
                        **({"metadata": meta} if meta is not None else {}))
    if close:
        plt.close(fig)
    return paths


def _metadata(ext: str) -> dict | None:
    from gcen2 import __version__

    # no creation dates: the same figure gives the same file
    if ext == "pdf":
        return {"Creator": f"GCEN2 {__version__}", "CreationDate": None}
    if ext == "svg":
        return {"Creator": f"GCEN2 {__version__}", "Date": None}
    if ext == "png":
        return {"Software": f"GCEN2 {__version__}"}
    return None
