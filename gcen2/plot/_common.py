"""Helpers shared by the plot functions."""

from __future__ import annotations

import contextlib

from gcen2.plot._style import new_figure, style


def _require_matplotlib():
    """Import matplotlib (a clear error when it is missing) with a headless
    backend outside interactive sessions; also used by :mod:`gcen2.plotting`."""
    try:
        import matplotlib  # noqa: F401
    except ImportError as exc:
        from gcen2.exceptions import GCENInputError

        raise GCENInputError("matplotlib is required for plotting. "
                             "Install with: pip install gcen2[visualization]") from exc
    from gcen2.plot._style import headless_backend

    headless_backend()


def figure(size="single", height=None, aspect: float = 0.75, layout="constrained"):
    """A new figure of the given width (``size``) and height in mm (else
    width x ``aspect``), constrained layout by default (``layout=None``:
    none)."""
    _require_matplotlib()
    return new_figure(size, height, aspect, layout=layout)


@contextlib.contextmanager
def canvas(ax=None, size="single", height=None, aspect=0.75):
    """Yield (fig, ax) in the GCEN2 style: the given axes, or a new figure."""
    _require_matplotlib()
    with style():
        if ax is None:
            fig = new_figure(size, height, aspect)
            ax = fig.add_subplot(111)
        else:
            fig = ax.figure
        yield fig, ax


def finish(fig, output, formats=None):
    """Save the figure if an output is given; return it."""
    if output is not None:
        from gcen2.plot._save import save_figure

        save_figure(fig, output, formats=formats)
    return fig


def r_linkage(merge, height):
    """SciPy linkage matrix of an R-convention dendrogram (merge: n - 1 x 2,
    negative = leaf -i, positive = step j; height)."""
    import numpy as np

    merge = np.asarray(merge, dtype=np.int64)
    height = np.asarray(height, dtype=np.float64)
    n = merge.shape[0] + 1
    z = np.zeros((n - 1, 4))
    size = np.zeros(n - 1, dtype=np.int64)
    for i, (a, b) in enumerate(merge):
        ia = -a - 1 if a < 0 else n + a - 1
        ib = -b - 1 if b < 0 else n + b - 1
        sa = 1 if a < 0 else size[a - 1]
        sb = 1 if b < 0 else size[b - 1]
        size[i] = sa + sb
        z[i] = (ia, ib, height[i], size[i])
    return z


def tree_layout(linkage, hang: float = 0.03):
    """Leaf order and the drawing segments of a SciPy linkage, without
    recursion (SciPy's dendrogram recurses once per level and fails on trees
    of thousands of genes). Leaves are at x = 0.5, 1.5, ...; each merge is
    drawn as the usual upside-down U."""
    import numpy as np

    z = np.asarray(linkage)
    n = z.shape[0] + 1
    # leaf order: depth-first, left child first, with an explicit stack
    order = []
    stack = [2 * n - 2]
    while stack:
        c = stack.pop()
        if c < n:
            order.append(c)
        else:
            a, b = int(z[c - n, 0]), int(z[c - n, 1])
            stack.append(b)
            stack.append(a)
    xpos = np.empty(2 * n - 1)
    xpos[np.asarray(order)] = np.arange(n) + 0.5
    ypos = np.zeros(2 * n - 1)
    hs = z[:, 2]
    span = (hs.max() - hs.min()) if n > 2 else 1.0
    drop = hang * (span if span > 0 else 1.0)
    segs = []
    bottom = np.inf
    for i in range(n - 1):
        a, b, h = int(z[i, 0]), int(z[i, 1]), z[i, 2]
        xa, xb = xpos[a], xpos[b]
        # leaves hang a little below their merge (R's hang), clusters end at
        # their own height
        ya = h - drop if a < n else ypos[a]
        yb = h - drop if b < n else ypos[b]
        bottom = min(bottom, ya, yb)
        segs.append([(xa, ya), (xa, h), (xb, h), (xb, yb)])
        xpos[n + i] = (xa + xb) / 2
        ypos[n + i] = h
    return np.asarray(order), segs, (float(bottom) if n > 1 else 0.0,
                                     float(hs.max()) if n > 1 else 1.0)


def draw_dendrogram(ax, linkage, color="black", lw=0.3, rasterize_above=2000):
    """Draw a SciPy linkage as a dendrogram on ``ax`` (leaves at x = 0.5, 1.5,
    ...); returns the leaf order."""
    from matplotlib.collections import LineCollection

    order, segs, (lo, hi) = tree_layout(linkage)
    lc = LineCollection(segs, colors=color, linewidths=lw)
    if len(order) > rasterize_above:
        lc.set_rasterized(True)
    ax.add_collection(lc)
    ax.set_xlim(0, len(order))
    pad = 0.02 * (hi - lo if hi > lo else 1)
    ax.set_ylim(lo - pad, hi + pad)
    return order


def color_rows(ax, rows, labels, rasterize_above=2000):
    """Colour bars under a dendrogram: ``rows`` lists of hex colours (one per
    leaf, in leaf order), named by ``labels``."""
    import numpy as np
    from matplotlib.colors import to_rgb

    img = np.array([[to_rgb(c) for c in r] for r in rows])
    im = ax.imshow(img, aspect="auto", interpolation="nearest",
                   extent=(0, img.shape[1], len(rows), 0))
    if img.shape[1] > rasterize_above:
        im.set_rasterized(True)
    ax.set_yticks(np.arange(len(rows)) + 0.5, labels)
    ax.set_xticks([])
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(axis="y", length=0)


def repel_labels(ax, xy, labels, fontsize: float = 4.5, color="black", point_px: float = 4.0,
                 iters: int = 300, leader: bool = True, bbox: bool = True, **text_kw):
    """Place text labels near their points without overlaps: each label
    starts just above its point; overlapping labels (and labels covering a
    point) are pushed apart in display space, with a gentle pull back to the
    anchor. Labels are kept as offsets in points from their anchors, so they
    stay apart when the axes are resized afterwards; labels that end up away
    from their point get a thin leader line. Returns the annotations."""
    import numpy as np

    xy = np.asarray(xy, dtype=float).reshape(-1, 2)
    if len(xy) == 0:
        return []
    fig = ax.figure
    ax.autoscale_view()
    ax.apply_aspect()
    box_kw = ({"boxstyle": "round,pad=0.12", "fc": "white", "ec": "none", "alpha": 0.75}
              if bbox else None)
    probe = [ax.text(0, 0, str(lab), fontsize=fontsize, **text_kw) for lab in labels]
    r = fig.canvas.get_renderer()
    pad = fontsize * fig.dpi / 72 * 0.5  # the white box around the text
    wh = np.array([[t.get_window_extent(r).width, t.get_window_extent(r).height]
                   for t in probe]) + pad
    for t in probe:
        t.remove()
    anchor = ax.transData.transform(xy)
    home = anchor + np.column_stack([np.zeros(len(xy)), wh[:, 1] / 2 + point_px])
    pos = home.copy()
    box = ax.get_window_extent(r)
    for _ in range(max(iters, 3 * len(xy))):
        moved = 0.0
        for i in range(len(pos)):  # label-label overlaps
            d = pos[i] - pos
            ox = (wh[i, 0] + wh[:, 0]) / 2 - np.abs(d[:, 0])
            oy = (wh[i, 1] + wh[:, 1]) / 2 - np.abs(d[:, 1])
            hit = (ox > 0) & (oy > 0)
            hit[i] = False
            for j in np.flatnonzero(hit):
                if oy[j] < ox[j]:
                    s = np.sign(d[j, 1]) or (1 if i < j else -1)
                    pos[i, 1] += s * oy[j] / 2
                    pos[j, 1] -= s * oy[j] / 2
                else:
                    s = np.sign(d[j, 0]) or (1 if i < j else -1)
                    pos[i, 0] += s * ox[j] / 2
                    pos[j, 0] -= s * ox[j] / 2
                moved += min(ox[j], oy[j])
        for i in range(len(pos)):  # labels off all points
            d = pos[i] - anchor
            ox = wh[i, 0] / 2 + point_px - np.abs(d[:, 0])
            oy = wh[i, 1] / 2 + point_px - np.abs(d[:, 1])
            for j in np.flatnonzero((ox > 0) & (oy > 0)):
                s = np.sign(d[j, 1]) or 1
                pos[i, 1] += s * oy[j] * 0.5
                moved += oy[j] * 0.5
        pos += (home - pos) * 0.01  # a little pull back, stay inside the axes
        pos[:, 0] = np.clip(pos[:, 0], box.x0 + wh[:, 0] / 2, box.x1 - wh[:, 0] / 2)
        pos[:, 1] = np.clip(pos[:, 1], box.y0 + wh[:, 1] / 2, box.y1 - wh[:, 1] * 1.5)
        if moved < 0.5:
            break
    # last pass: labels apart from each other (the point push above can stall
    # them in a dense cluster; a label over a point is the smaller harm)
    for _ in range(500):
        moved = 0.0
        for i in range(len(pos)):
            d = pos[i] - pos
            ox = (wh[i, 0] + wh[:, 0]) / 2 - np.abs(d[:, 0])
            oy = (wh[i, 1] + wh[:, 1]) / 2 - np.abs(d[:, 1])
            hit = (ox > 0) & (oy > 0)
            hit[i] = False
            for j in np.flatnonzero(hit):
                s = np.sign(d[j, 1]) or (1 if i < j else -1)
                pos[i, 1] += s * (oy[j] / 2 + 0.25)
                pos[j, 1] -= s * (oy[j] / 2 + 0.25)
                moved += oy[j]
        pos[:, 1] = np.clip(pos[:, 1], box.y0 + wh[:, 1] / 2, box.y1 - wh[:, 1] * 1.5)
        if moved < 0.5:
            break
    off = (pos - anchor) * 72 / fig.dpi
    far = np.hypot(*(pos - anchor).T) > wh[:, 1] + point_px + 2
    out = []
    for k, lab in enumerate(labels):
        arrow = ({"arrowstyle": "-", "color": "#707070", "lw": 0.3, "shrinkA": 0, "shrinkB": 2}
                 if leader and far[k] else None)
        out.append(ax.annotate(str(lab), xy[k], xytext=tuple(off[k]),
                               textcoords="offset points", fontsize=fontsize, color=color,
                               ha="center", va="center", zorder=6, bbox=box_kw,
                               arrowprops=arrow, annotation_clip=False, **text_kw))
    return out


def pack_components(graph, seed: int = 0, weight: str = "weight"):
    """Positions of a networkx graph with each connected component laid out
    on its own (spring layout) and the components packed in rows, largest
    first, so small components do not drift far from the rest."""
    import networkx as nx
    import numpy as np

    comps = sorted(nx.connected_components(graph), key=len, reverse=True)
    pos, boxes = {}, []
    for comp in comps:
        sub = graph.subgraph(comp)
        if len(comp) == 1:
            p = {next(iter(comp)): np.zeros(2)}
        else:
            p = nx.spring_layout(sub, seed=seed, weight=weight,
                                 k=1.0 / np.sqrt(len(comp)), iterations=100)
        a = np.array(list(p.values()))
        a = a - a.min(axis=0)
        span = a.max(axis=0)
        scale = np.sqrt(len(comp)) / max(span.max(), 1e-9) if len(comp) > 1 else 0
        a = a * scale
        boxes.append((list(p.keys()), a, a.max(axis=0) + 1.0))
    width = max(np.sqrt(sum(b[2][0] * b[2][1] for b in boxes)) * 1.3, boxes[0][2][0])
    x = y = row_h = 0.0
    for keys, a, (w, h) in boxes:
        if x + w > width and x > 0:
            x, y, row_h = 0.0, y - row_h, 0.0
        for k, v in zip(keys, a):
            pos[k] = np.array([x + v[0] + 0.5, y - v[1] - 0.5])
        x += w
        row_h = max(row_h, h)
    return pos


def plain_log_axis(axis, subs=(1.0,)) -> None:
    """Log-scale tick labels as plain numbers (300, 1,000) instead of 3×10²;
    ``subs`` the ticks within each decade (e.g. (1, 2, 5) for a short axis)."""
    from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter

    axis.set_major_locator(LogLocator(base=10, subs=subs))
    axis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:,.0f}" if v >= 1 else f"{v:g}"))
    axis.set_minor_formatter(NullFormatter())


_VENN2 = {"centers": [(-0.55, 0.0), (0.55, 0.0)],
          "regions": {1: (-1.02, 0.0), 2: (1.02, 0.0), 3: (0.0, 0.0)},
          "names": [(-0.25, 1.15, "right"), (0.25, 1.15, "left")]}
_VENN3 = {"centers": [(-0.58, 0.38), (0.58, 0.38), (0.0, -0.6)],
          "regions": {1: (-1.05, 0.62), 2: (1.05, 0.62), 4: (0.0, -1.2), 3: (0.0, 0.85),
                      5: (-0.68, -0.38), 6: (0.68, -0.38), 7: (0.0, 0.05)},
          "names": [(-0.95, 1.5, "center"), (0.95, 1.5, "center"), (0.0, -1.72, "center")]}


def draw_venn(ax, names, counts: dict, colors, extra: dict | None = None,
              fontsize: float = 5.5) -> None:
    """A 2- or 3-set Venn diagram with unit circles (areas not to scale):
    ``counts`` by region, the bit mask of the sets it belongs to (1 = first
    set, 2 = second, 4 = third); ``extra``: a second line per region (e.g.
    how many are kept downstream)."""
    from matplotlib.patches import Circle

    lay = _VENN2 if len(names) == 2 else _VENN3
    for (cx, cy), c in zip(lay["centers"], colors):
        ax.add_patch(Circle((cx, cy), 1.0, fc=c, ec="none", alpha=0.28))
        ax.add_patch(Circle((cx, cy), 1.0, fc="none", ec=c, lw=0.8))
    for mask, (x, y) in lay["regions"].items():
        n = counts.get(mask, 0)
        txt = f"{n:,}"
        if extra is not None and mask in extra:
            txt += f"\n{extra[mask]}"
        ax.text(x, y, txt, ha="center", va="center", fontsize=fontsize, linespacing=1.1)
    for (x, y, ha), nm in zip(lay["names"], names):
        ax.text(x, y, nm, ha=ha, va="center", fontsize=fontsize + 0.5, fontweight="bold")
    ax.set_aspect("equal")
    ax.set_xlim(-1.75, 1.75)
    ax.set_ylim(-1.9 if len(names) == 3 else -1.15, 1.75 if len(names) == 3 else 1.5)
    ax.axis("off")


def draw_upset(fig, names, counts: dict, colors, top: int = 15, extra: dict | None = None):
    """UpSet-style plot for any number of sets: bars of the ``top``
    largest intersections (exclusive regions, ``counts`` by bit mask) over a
    dot matrix of the sets they belong to; ``extra`` (by mask) adds a
    darker bar part (e.g. kept downstream). Returns (bar axes, matrix axes)."""
    import numpy as np

    items = sorted(((m, c) for m, c in counts.items() if c > 0), key=lambda t: -t[1])[:top]
    a1, a2 = fig.subplots(2, 1, height_ratios=[2.2, 0.35 * len(names) + 0.3], sharex=True)
    x = np.arange(len(items))
    v = [c for _, c in items]
    a1.bar(x, v, color="#BBBBBB", width=0.7, label="all")
    if extra is not None:
        a1.bar(x, [extra.get(m, 0) for m, _ in items], color="#555555", width=0.7,
               label="in merged network")
        a1.legend(loc="upper right")
    for xi, vi in zip(x, v):
        a1.text(xi, vi, f"{vi:,}", ha="center", va="bottom", fontsize=4.5, rotation=90)
    a1.set_ylabel("Edges")
    a1.set_ylim(0, max(v) * 1.4 if v else 1)
    from matplotlib.ticker import FuncFormatter

    a1.yaxis.set_major_formatter(FuncFormatter(lambda t, _: f"{t:,.0f}"))
    a1.tick_params(axis="x", length=0)
    for j in range(len(names)):
        on = [bool(m >> j & 1) for m, _ in items]
        a2.scatter(x, np.full(len(x), j), s=14, c=[colors[j] if o else "#E0E0E0" for o in on],
                   zorder=2)
    for xi, (m, _) in zip(x, items):
        js = [j for j in range(len(names)) if m >> j & 1]
        if len(js) > 1:
            a2.plot([xi, xi], [min(js), max(js)], color="#555555", lw=0.8, zorder=1)
    a2.set_yticks(range(len(names)), names)
    a2.set_ylim(-0.6, len(names) - 0.4)
    a2.invert_yaxis()
    a2.set_xticks([])
    for s in a2.spines.values():
        s.set_visible(False)
    return a1, a2
