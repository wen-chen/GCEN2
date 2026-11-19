"""Multi-omics network fusion figures (``gcen2 fuse``)."""

from __future__ import annotations

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.plot._common import figure, finish
from gcen2.plot._io import table as read_table
from gcen2.plot._style import OKABE_ITO, style


def fusion_overview(fused, size="double", height=None, output=None, formats=None):
    """A fused network: edges supported by 1, 2, ... layers, edges per layer
    and the overlap of every pair of layers, and the fused score
    distribution."""
    f = read_table(fused)
    wcols = [c for c in f.columns if c.startswith("w_")]
    if "n_layers" not in f.columns or not wcols:
        raise GCENInputError("not a gcen2 fuse output (n_layers, w_<layer> columns).")
    layers = [c[2:] for c in wcols]
    present = f[wcols].notna() & (f[wcols].astype(float) != 0)
    n = len(layers)
    ov = np.zeros((n, n), dtype=np.int64)
    for i in range(n):
        for j in range(n):
            ov[i, j] = int((present.iloc[:, i] & present.iloc[:, j]).sum())
    with style():
        fig = figure(size, height or 55)
        a1, a2, a3 = fig.subplots(1, 3)
        cnt = f["n_layers"].value_counts().sort_index()
        a1.bar(cnt.index.astype(str), cnt.to_numpy(), color=OKABE_ITO[4])
        a1.set_xlabel("Layers supporting the edge")
        a1.set_ylabel("Edges")
        for xi, v in enumerate(cnt.to_numpy()):
            a1.text(xi, v, f"{v:,}", ha="center", va="bottom", fontsize=4.5)
        a1.set_ylim(0, cnt.max() * 1.15)
        a1.ticklabel_format(axis="y", style="sci", scilimits=(-3, 4), useMathText=True)
        im = a2.imshow(ov, cmap="viridis")
        a2.set_xticks(range(n), layers, rotation=45, ha="right")
        a2.set_yticks(range(n), layers)
        for i in range(n):
            for j in range(n):
                a2.text(j, i, f"{ov[i, j]:,}", ha="center", va="center", fontsize=4.5,
                        color="white" if ov[i, j] < ov.max() * 0.6 else "black")
        a2.set_title("Edges shared by layers")
        fig.colorbar(im, ax=a2, shrink=0.6).outline.set_linewidth(0.3)
        a3.hist(pd.to_numeric(f["score"], errors="coerce").dropna(), bins=60,
                color=OKABE_ITO[2], edgecolor="none")
        a3.set_xlabel("Fused score")
        a3.set_ylabel("Edges")
    return finish(fig, output, formats)
