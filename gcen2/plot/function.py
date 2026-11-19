"""Gene function prediction figures: the RWR ranking (``gcen2 rwr`` and
``gcen2 rwr-multiplex``) and the network around the seeds."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from gcen2.exceptions import GCENInputError
from gcen2.plot._common import canvas, finish
from gcen2.plot._io import table as read_table
from gcen2.plot._style import OKABE_ITO


def _rwr(table) -> pd.DataFrame:
    t = read_table(table).copy()
    t.columns = [str(c).lstrip("#") for c in t.columns]
    score = next((c for c in ("probability", "score", "p", "rwr") if c in t.columns),
                 t.columns[-1])
    return t.rename(columns={t.columns[0]: "gene", score: "score"})


def _seeds(seeds) -> set:
    if seeds is None:
        return set()
    if isinstance(seeds, (str, Path)) and Path(seeds).is_file():
        return {ln.split()[0] for ln in Path(seeds).read_text(encoding="utf-8").splitlines()
                if ln.strip() and not ln.startswith("#")}
    return set(map(str, seeds))


def rwr_rank(table, seeds=None, top: int = 30, size="single", height=None, output=None,
             formats=None):
    """The ``top`` genes by random-walk score, seeds marked: a lollipop chart
    (score on a log axis)."""
    t = _rwr(table).sort_values("score", ascending=False)
    if t.empty:
        raise GCENInputError("the RWR result is empty.")
    s = _seeds(seeds)
    t = t.head(top)
    with canvas(None, size, height or max(35, 8 + 3.4 * len(t)), 0.7) as (fig, ax):
        y = np.arange(len(t))
        v = t["score"].to_numpy(float)
        is_seed = t["gene"].astype(str).isin(s).to_numpy()
        ax.hlines(y, 0, v, color="#BBBBBB", lw=0.6)
        ax.scatter(v[~is_seed], y[~is_seed], s=16, color=OKABE_ITO[4], zorder=3,
                   label="candidate")
        if is_seed.any():
            ax.scatter(v[is_seed], y[is_seed], s=18, marker="D", color=OKABE_ITO[5], zorder=3,
                       label="seed")
        ax.set_yticks(y, t["gene"].astype(str))
        ax.invert_yaxis()
        ax.set_xlim(0, v.max() * 1.05)
        ax.ticklabel_format(axis="x", style="sci", scilimits=(-3, 3), useMathText=True)
        ax.set_xlabel("RWR score")
        ax.legend(loc="upper left", bbox_to_anchor=(1.0, 1.0))
    return finish(fig, output, formats)


def rwr_subnetwork(network, table, seeds, top: int = 20, seed: int = 0, size="single",
                   height=None, output=None, formats=None):
    """The seeds and the ``top`` candidates with their edges in ``network``:
    candidates coloured by score, seeds as diamonds."""
    try:
        import networkx as nx
    except ImportError as exc:
        raise GCENInputError("rwr_subnetwork needs networkx.") from exc
    from gcen2.plot.network import _edges

    t = _rwr(table).sort_values("score", ascending=False)
    s = _seeds(seeds)
    cand = [g for g in t["gene"].astype(str) if g not in s][:top]
    keep = s | set(cand)
    e = _edges(network)
    a, b = e.columns[:2]
    e = e[e[a].astype(str).isin(keep) & e[b].astype(str).isin(keep)]
    g = nx.Graph()
    g.add_nodes_from(keep)
    g.add_edges_from(zip(e[a].astype(str), e[b].astype(str)))
    g.remove_nodes_from([n for n in list(g.nodes) if g.degree(n) == 0 and n not in s])
    linked = [n for n in g.nodes if g.degree(n) > 0]
    lone = [n for n in g.nodes if g.degree(n) == 0]
    pos = nx.spring_layout(g.subgraph(linked), seed=seed,
                           k=1.5 / np.sqrt(max(len(linked), 1))) if linked else {}
    for i, n in enumerate(lone):  # unconnected seeds: a row below the network
        pos[n] = np.array([-1 + 2 * (i + 0.5) / len(lone), -1.35])
    score = dict(zip(t["gene"].astype(str), t["score"].astype(float)))
    named = set(s) | set(cand[:8])
    with canvas(None, size, height or 80, 1.0) as (fig, ax):
        for u, v in g.edges():
            ax.plot([pos[u][0], pos[v][0]], [pos[u][1], pos[v][1]], color="#909090", lw=0.4,
                    zorder=1)
        cn = [n for n in g.nodes if n not in s]
        sc = ax.scatter([pos[n][0] for n in cn], [pos[n][1] for n in cn], s=40,
                        c=[score.get(n, 0) for n in cn], cmap="viridis", edgecolors="black",
                        linewidths=0.3, zorder=2)
        sn = [n for n in g.nodes if n in s]
        ax.scatter([pos[n][0] for n in sn], [pos[n][1] for n in sn], s=50, marker="D",
                   color=OKABE_ITO[5], edgecolors="black", linewidths=0.3, zorder=2)
        for n in g.nodes:
            if n in named:
                ax.annotate(n, pos[n], xytext=(3, 3), textcoords="offset points",
                            fontsize=4.5)
        ax.axis("off")
        cb = fig.colorbar(sc, ax=ax, shrink=0.6, format="%.2e")
        cb.set_label("RWR score")
        cb.outline.set_linewidth(0.3)
        cb.outline.set_linewidth(0.3)
        ax.set_title(f"{len(sn)} seeds (diamonds), {len(cn)} candidates")
    return finish(fig, output, formats)
