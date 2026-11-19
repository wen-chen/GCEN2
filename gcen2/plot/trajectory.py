"""Trajectory network figures (``gcen2 trajectory``): the networks of the
pseudotime windows, the edge patterns, and the genes whose connectivity
changes most."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.plot._common import figure, finish
from gcen2.plot._style import OKABE_ITO, style


def trajectory_overview(directory, top: int = 25, size="double", height=None, output=None,
                        formats=None):
    """Networks along pseudotime: genes and edges per window against the
    window centre, edges per pattern (stable, gained, lost, transient,
    intermittent), and the degree of the ``top`` most variable genes across
    the windows (scaled per gene)."""
    d = Path(directory)
    for f in ("windows.tsv", "edge_dynamics.tsv", "degree_dynamics.tsv"):
        if not (d / f).is_file():
            raise GCENInputError(f"{d / f} not found (a gcen2 trajectory output directory?).")
    w = pd.read_csv(d / "windows.tsv", sep="\t")
    e = pd.read_csv(d / "edge_dynamics.tsv", sep="\t", usecols=lambda c: c == "pattern")
    deg = pd.read_csv(d / "degree_dynamics.tsv", sep="\t", index_col=0)
    with style():
        fig = figure(size, height or 110)
        gs = fig.add_gridspec(2, 2, height_ratios=[1, 1.4])
        a1 = fig.add_subplot(gs[0, 0])
        a1.plot(w["center"], w["edge_count"], marker="o", ms=2.5, color=OKABE_ITO[4],
                label="edges")
        a1.set_xlabel("Pseudotime (window centre)")
        a1.set_ylabel("Edges", color=OKABE_ITO[4])
        a1b = a1.twinx()
        a1b.plot(w["center"], w["node_count"], marker="s", ms=2.5, color=OKABE_ITO[5],
                 label="genes")
        a1b.set_ylabel("Genes", color=OKABE_ITO[5])
        a1b.spines["right"].set_visible(True)
        a2 = fig.add_subplot(gs[0, 1])
        cnt = e["pattern"].value_counts()
        a2.bar(range(len(cnt)), cnt.to_numpy(), color=OKABE_ITO[2])
        a2.set_xticks(range(len(cnt)), cnt.index, rotation=30, ha="right")
        a2.set_ylabel("Edges")
        a2.set_title("Edge patterns")
        a3 = fig.add_subplot(gs[1, :])
        v = deg.to_numpy(float)
        sd = v.std(axis=1)
        idx = np.argsort(-sd)[:top]
        sub = v[idx]
        z = (sub - sub.mean(axis=1, keepdims=True)) / np.where(sub.std(axis=1, keepdims=True)
                                                               > 0,
                                                               sub.std(axis=1, keepdims=True),
                                                               1)
        im = a3.imshow(z, aspect="auto", cmap="RdBu_r", vmin=-2.5, vmax=2.5,
                       interpolation="nearest")
        a3.set_yticks(range(len(idx)), deg.index[idx], fontsize=4.5)
        a3.set_xticks(range(v.shape[1]), [f"{c:.1f}" for c in w["center"]][:v.shape[1]],
                      fontsize=5)
        a3.set_xlabel("Window (pseudotime centre)")
        a3.set_title(f"Degree of the {len(idx)} most variable genes (z-score per gene)")
        cb = fig.colorbar(im, ax=a3, shrink=0.6)
        cb.set_label("Degree (z)")
        cb.outline.set_linewidth(0.3)
    return finish(fig, output, formats)
