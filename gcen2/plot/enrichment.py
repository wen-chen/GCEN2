"""Functional enrichment figures (GO / KEGG of ``gcen2 enrich`` and
``gcen2 annotate``): bar chart of one result, dot plot and heatmap of the
modules' top terms (hdWGCNA EnrichrBarPlot / EnrichrDotPlot)."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.plot._common import canvas, figure, finish
from gcen2.plot._io import table as read_table
from gcen2.plot._style import OKABE_ITO, module_color, module_label, style

_NS = {"biological_process": ("BP", OKABE_ITO[2]), "molecular_function": ("MF", OKABE_ITO[0]),
       "cellular_component": ("CC", OKABE_ITO[1])}


def _short(text: str, width: int = 45) -> str:
    text = str(text)
    return text if len(text) <= width else text[:width - 1] + "…"


def _term_table(x) -> pd.DataFrame:
    t = read_table(x).copy()
    cols = {c.lower(): c for c in t.columns}
    pcol = next((cols[c] for c in ("p_adj", "fdr", "q_value", "p_value", "pvalue", "p")
                 if c in cols), None)
    if pcol is None:
        raise GCENInputError("an enrichment table needs a p-value column.")
    name = cols.get("name") or cols.get("description") or cols.get("term_id") or t.columns[0]
    t = t.assign(_p=pd.to_numeric(t[pcol], errors="coerce"), _name=t[name].astype(str),
                 _pcol=pcol)
    return t.dropna(subset=["_p"])


def _p_label(tabs) -> str:
    """Colour-bar label: the p-value column the tables share (e.g. fdr), else p."""
    cols = {t["_pcol"].iat[0] for t in tabs if not t.empty}
    return f"−log10 {cols.pop()}" if len(cols) == 1 else "−log10 p"


def enrich_bar(table, top: int = 15, size="single", height=None, output=None, formats=None):
    """The ``top`` terms of an enrichment result (``gcen2 enrich`` / one
    ``gcen2 annotate`` file) by p-value: bars of -log10 p, coloured by GO
    namespace, with the number of genes."""
    t = _term_table(table).sort_values("_p").head(top)
    if t.empty:
        raise GCENInputError("no enriched terms.")
    with canvas(None, size, height or max(30, 8 + 4.2 * len(t)), 0.7) as (fig, ax):
        y = np.arange(len(t))
        ns = t["namespace"] if "namespace" in t.columns else pd.Series([""] * len(t))
        cols = [_NS.get(n, ("", OKABE_ITO[4]))[1] for n in ns]
        v = -np.log10(np.maximum(t["_p"].to_numpy(float), 1e-300))
        ax.barh(y, v, color=cols, height=0.7)
        ax.set_yticks(y, [_short(n) for n in t["_name"]])
        ax.invert_yaxis()
        ax.set_xlabel(f"−log10 {t['_pcol'].iat[0]}")
        if "study_count" in t.columns:
            for yi, vi, c in zip(y, v, t["study_count"]):
                ax.text(vi, yi, f" {int(c)}", va="center", fontsize=4.5)
        present = [k for k in _NS if k in set(ns)]
        if present:
            from matplotlib.patches import Patch

            ax.legend([Patch(color=_NS[k][1]) for k in present], [_NS[k][0] for k in present],
                      loc="lower right")
    return finish(fig, output, formats)


def _module_tables(results) -> dict:
    """module -> term table, from a ``gcen2 annotate`` directory
    (module_<k>.go / .kegg) or a dict module -> table or EnrichmentResult
    (``gcen2.annotation.annotate_go`` of a ModuleSet: keys ``module_<k>``).
    Modules without enriched terms are left out, as in the directory."""
    if isinstance(results, dict):
        out = {}
        for k, v in results.items():
            v = getattr(v, "results", v)  # EnrichmentResult -> its table
            if isinstance(v, pd.DataFrame) and v.empty:
                continue
            t = _term_table(v)
            if not t.empty:
                key = str(k)
                out[key[len("module_"):] if key.startswith("module_") else key] = t
        return out
    d = Path(results)
    if not d.is_dir():
        raise GCENInputError(f"Not a directory: {d}")
    out = {}
    for f in sorted(d.glob("module_*.*")):
        if f.suffix in (".go", ".kegg", ".tsv") and f.stat().st_size > 0:
            k = f.stem.split("_", 1)[1]
            try:
                t = _term_table(f)
            except (GCENInputError, pd.errors.EmptyDataError):
                continue
            if not t.empty:
                out.setdefault(k, []).append(t)
    out = {k: pd.concat(v, ignore_index=True) for k, v in out.items()}
    return dict(sorted(out.items(), key=lambda kv: int(kv[0]) if kv[0].isdigit() else 10 ** 9))


def enrich_dot(results, modules=None, top: int = 3, max_modules: int = 12, size="double",
               height=None, output=None, formats=None):
    """The ``top`` terms of each module (hdWGCNA EnrichrDotPlot): a dot per
    module and term, size = genes in the term, colour = -log10 p (or of the
    table's adjusted p-value column, e.g. fdr); modules in
    their colours along the axis. ``results``: a ``gcen2 annotate`` output
    directory or a dict module -> table."""
    tabs = _module_tables(results)
    keys = list(tabs) if modules is None else [str(m).replace("ME", "") for m in modules]
    keys = [k for k in keys if k in tabs][:max_modules]
    if not keys:
        raise GCENInputError("no module with enriched terms.")
    rows = []
    for k in keys:
        t = tabs[k].sort_values("_p").head(top)
        for _, r in t.iterrows():
            rows.append((k, _short(r["_name"], 50), r["_p"], r.get("study_count", np.nan)))
    d = pd.DataFrame(rows, columns=["module", "term", "p", "n"])
    terms = list(dict.fromkeys(d["term"]))
    with style():
        fig = figure(size, height or max(40, 10 + 3.6 * len(terms)))
        from matplotlib.patches import Rectangle

        ax = fig.add_subplot(111)
        xi = d["module"].map({k: i for i, k in enumerate(keys)})
        yi = d["term"].map({t: i for i, t in enumerate(terms)})
        v = -np.log10(np.maximum(d["p"].to_numpy(float), 1e-300))
        nn = d["n"].fillna(5).to_numpy(float)
        sc = ax.scatter(xi, yi, s=8 + 50 * nn / max(nn.max(), 1), c=v, cmap="viridis",
                        edgecolors="black", linewidths=0.25)
        labels = [module_label(f"ME{k}" if k.isdigit() else k) for k in keys]
        ax.set_xticks(range(len(keys)), labels, rotation=45, ha="right")
        ax.set_yticks(range(len(terms)), terms)
        ax.set_xlim(-0.6, len(keys) - 0.4)
        ax.set_ylim(len(terms) - 0.5, -0.5)
        for i, k in enumerate(keys):
            ax.add_patch(Rectangle((i - 0.4, len(terms) - 0.45), 0.8, 0.35, clip_on=False,
                                       color=module_color(f"ME{k}" if k.isdigit() else k)))
        cb = fig.colorbar(sc, ax=ax, shrink=0.5, aspect=15)
        cb.set_label(_p_label(tabs[k] for k in keys))
        cb.outline.set_linewidth(0.3)
    return finish(fig, output, formats)


def enrich_heatmap(results, top: int = 2, max_modules: int = 12, size="double", height=None,
                   output=None, formats=None):
    """Modules x terms heatmap of -log10 p (the union of each module's
    ``top`` terms), to compare modules' functions at a glance; the colour
    scale stops at the 95th percentile so that a few extreme p-values do not
    hide the rest (marked by the arrow on the colour bar)."""
    tabs = _module_tables(results)
    keys = list(tabs)[:max_modules]
    if not keys:
        raise GCENInputError("no module with enriched terms.")
    terms = []
    for k in keys:
        terms.extend(tabs[k].sort_values("_p").head(top)["_name"].tolist())
    terms = list(dict.fromkeys(terms))
    mat = np.zeros((len(terms), len(keys)))
    for j, k in enumerate(keys):
        best = tabs[k].groupby("_name")["_p"].min()
        for i, t in enumerate(terms):
            if t in best.index:
                mat[i, j] = -np.log10(max(best[t], 1e-300))
    with canvas(None, size, height or max(40, 10 + 3.6 * len(terms)), 0.8) as (fig, ax):
        nz = mat[mat > 0]
        vmax = float(np.quantile(nz, 0.95)) if nz.size else 1.0
        im = ax.imshow(mat, cmap="viridis", aspect="auto", interpolation="nearest", vmin=0,
                       vmax=vmax)
        ax.set_xticks(range(len(keys)), [module_label(f"ME{k}" if k.isdigit() else k)
                                         for k in keys], rotation=45, ha="right")
        ax.set_yticks(range(len(terms)), [_short(t, 50) for t in terms])
        cb = fig.colorbar(im, ax=ax, shrink=0.6, extend="max" if mat.max() > vmax else "neither")
        cb.set_label(_p_label(tabs[k] for k in keys))
        cb.outline.set_linewidth(0.3)
    return finish(fig, output, formats)
