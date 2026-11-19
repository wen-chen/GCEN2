"""One figure style for every GCEN2 plot: journal sizes (Nature's 89 / 120 /
183 mm widths), 5-7 pt sans-serif text kept as editable text in PDF and SVG,
thin lines, colour-blind-aware palettes, and one colour per module across all
figures (WGCNA's standardColors by default)."""

from __future__ import annotations

import contextlib
import re

from gcen2.plot._wgcna_colors import STANDARD_COLORS

MM = 1 / 25.4
WIDTHS = {"single": 89 * MM, "onehalf": 120 * MM, "double": 183 * MM}

# Okabe & Ito (2008), colour-blind safe; black last
OKABE_ITO = ["#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7",
             "#000000"]
GREY = "#BEBEBE"   # R's grey: unassigned genes
GOLD = "#FFD700"   # WGCNA's random "gold" module
DIVERGING = "RdBu_r"
SEQUENTIAL = "viridis"

_WGCNA = dict(STANDARD_COLORS)
_STATE = {"palette": "wgcna", "formats": ("pdf", "png", "svg")}
_FONTS = ["Arial", "Helvetica", "Liberation Sans", "Nimbus Sans", "DejaVu Sans"]


def set_defaults(palette: str | None = None, formats=None) -> None:
    """Change the module palette ("wgcna" or "okabe-ito") and the formats
    figures are saved in (e.g. ("pdf", "png", "svg"))."""
    from gcen2.exceptions import GCENInputError

    if palette is not None:
        if palette not in ("wgcna", "okabe-ito"):
            raise GCENInputError("palette must be 'wgcna' or 'okabe-ito'.")
        _STATE["palette"] = palette
    if formats is not None:
        formats = tuple(str(f).lower().lstrip(".") for f in formats)
        bad = [f for f in formats if f not in ("pdf", "png", "svg", "tiff", "tif", "eps")]
        if bad or not formats:
            raise GCENInputError(f"unknown figure format(s): {', '.join(bad) or 'none'}.")
        _STATE["formats"] = formats


def get_defaults() -> dict:
    return dict(_STATE)


def _font_family() -> list[str]:
    from matplotlib import font_manager

    have = {f.name for f in font_manager.fontManager.ttflist}
    found = [f for f in _FONTS if f in have]
    return found or ["sans-serif"]


def rc_params() -> dict:
    """The matplotlib settings of the GCEN2 style."""
    return {
        "font.family": "sans-serif",
        "font.sans-serif": _font_family(),
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "figure.titlesize": 7,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "legend.fontsize": 6,
        "legend.title_fontsize": 6,
        "axes.linewidth": 0.5, "lines.linewidth": 0.75, "patch.linewidth": 0.5,
        "xtick.major.width": 0.5, "ytick.major.width": 0.5,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5,
        "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
        "xtick.major.pad": 1.5, "ytick.major.pad": 1.5, "axes.labelpad": 2,
        "axes.titlepad": 3, "axes.spines.top": False, "axes.spines.right": False,
        "legend.frameon": False, "legend.handlelength": 1.2, "legend.borderaxespad": 0.3,
        "image.cmap": SEQUENTIAL,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "svg.hashsalt": "gcen2",  # reproducible SVG element ids
        "savefig.dpi": 300, "figure.dpi": 150, "savefig.transparent": False,
        "figure.facecolor": "white", "axes.facecolor": "white",
    }


def headless_backend() -> None:
    """Outside an interactive session, with no backend chosen (no MPLBACKEND,
    no matplotlib.use), draw with Agg: matplotlib's automatic choice probes
    for a display, which can hang for minutes when DISPLAY points to an X
    server that does not answer (e.g. a stale ssh session)."""
    import sys

    import matplotlib

    get = getattr(matplotlib.rcParams, "_get_backend_or_none", None)
    chosen = get() if get is not None else matplotlib.rcParams.get("backend")
    interactive = hasattr(sys, "ps1") or "ipykernel" in sys.modules or sys.flags.interactive
    if chosen is None and not interactive:
        matplotlib.use("Agg")


@contextlib.contextmanager
def style():
    """Use the GCEN2 style inside a with block (the user's global matplotlib
    settings are left as they were). Every GCEN2 figure is made inside it,
    so it also picks the headless backend when none was chosen."""
    import matplotlib as mpl

    headless_backend()
    with mpl.rc_context(rc_params()):
        yield


def figsize(size="single", height: float | None = None, aspect: float = 0.75):
    """(width, height) in inches: ``size`` "single" (89 mm), "onehalf"
    (120 mm), "double" (183 mm) or a width in mm; ``height`` in mm, else
    width x ``aspect`` (capped at Nature's 170 mm)."""
    from gcen2.exceptions import GCENInputError

    if isinstance(size, (int, float)):
        w = float(size) * MM
    elif size in WIDTHS:
        w = WIDTHS[size]
    else:
        raise GCENInputError("size must be single, onehalf, double or a width in mm.")
    h = height * MM if height is not None else min(w * aspect, 170 * MM)
    return w, h


def new_figure(size="single", height=None, aspect=0.75, **kwargs):
    import matplotlib.pyplot as plt

    return plt.figure(figsize=figsize(size, height, aspect), **kwargs)


_ME = re.compile(r"^(?:ME|M|module)?\s*(\d+)$", re.IGNORECASE)


def module_index(name) -> int | None:
    """0-based module number of a module name as GCEN2 writes them ("ME3",
    "3", 3), else None."""
    if isinstance(name, (int,)) and not isinstance(name, bool):
        return int(name)
    m = _ME.match(str(name).strip())
    return int(m.group(1)) if m else None


def is_unassigned(name) -> bool:
    return str(name).strip().lower() in ("grey", "gray", "me_unassigned", "unassigned", "megrey",
                                         "na", "nan")


def module_color(name, palette: str | None = None) -> str:
    """The colour of a module in every figure: hdWGCNA / WGCNA colour names
    are their own colour; GCEN2's module k (0-based, "ME<k>") is the
    (k + 1)-th WGCNA standard colour (turquoise, blue, brown, ...), or the
    k-th Okabe-Ito colour with palette "okabe-ito"; unassigned genes grey,
    WGCNA's random module gold."""
    palette = palette or _STATE["palette"]
    s = str(name).strip()
    low = s.lower()
    if low.startswith("me") and low[2:] in _WGCNA:
        low = low[2:]
    if is_unassigned(low):
        return GREY
    if low == "gold":
        return GOLD
    k = module_index(s)
    if k is not None:
        if palette == "okabe-ito":
            return OKABE_ITO[k % len(OKABE_ITO)]
        return STANDARD_COLORS[k % len(STANDARD_COLORS)][1]
    if low in _WGCNA:
        return _WGCNA[low]
    # any other name: a stable colour from its text
    return OKABE_ITO[sum(map(ord, s)) % len(OKABE_ITO)]


def module_color_name(name) -> str | None:
    """WGCNA's colour name of a GCEN2 module ("ME0" -> "turquoise")."""
    k = module_index(name)
    if k is not None:
        return STANDARD_COLORS[k % len(STANDARD_COLORS)][0]
    low = str(name).strip().lower()
    return low if low in _WGCNA or low in ("grey", "gold") else None


def module_label(name) -> str:
    """How a module is named in figures: "ME0 (turquoise)" for GCEN2's
    numbered modules (with the WGCNA palette), else the name itself, so a
    module never depends on colour alone."""
    s = str(name).strip()
    if is_unassigned(s):
        return "unassigned"
    k = module_index(s)
    if k is not None:
        base = f"ME{k}"
        return f"{base} ({module_color_name(k)})" if _STATE["palette"] == "wgcna" else base
    return s
