"""Multi-omics network integration: heterogeneous layers over shared genes.

A *layer* is an undirected gene network of any kind (co-expression, PPI,
regulatory, ...), given as a :class:`~gcen2.models.Network`, a DataFrame or
an edge-list file. Its edge weights are ``|correlation|`` when the layer has
a ``correlation`` column, else ``|weight|`` when it has a ``weight`` column
(or the column named by ``weight_column``), else 1.

Two integration strategies (requirement §3.7.1):

* :func:`fuse` — weighted fusion into one network: per edge, the weighted
  mean of its (normalised) layer weights, an edge absent from a layer
  counting as 0.
* :func:`rwr_multiplex` — random walk with restart on the multiplex network
  (RWR-M, Valdeolivas et al., Bioinformatics 2019; by default the transition
  matrix of their RandomWalkRestartMH): one copy of every gene per layer,
  walks move within a layer or jump to the same gene in another layer, and
  the per-layer scores are aggregated per gene.
"""

from __future__ import annotations

from collections.abc import Mapping
from pathlib import Path
from typing import TYPE_CHECKING, Union

import numpy as np

from gcen2._matrix_io import require_file, tsv_header
from gcen2.core import extension
from gcen2.exceptions import GCENInputError, GCENIOError

if TYPE_CHECKING:
    import pandas as pd

    from gcen2.models import Network, RWRResult

LayerLike = Union["Network", "pd.DataFrame", str, Path]

NORMALIZE = ("max", "rank", "none")
AGGREGATE = ("geometric", "arithmetic", "sum")
TRANSITION = ("supra", "layer")


def read_layer(path: Union[str, Path], weight_column: Union[str, int, None] = None,
               sep: str = "\t") -> pd.DataFrame:
    """Read an edge-list file as a layer table ``gene_a, gene_b, weight``.

    The first two columns are the genes. A header line is recognised when the
    file has a third column whose first value is not a number (GCEN network
    files, STRING ``protein1 protein2 combined_score`` files, ...). The weight
    is the column ``weight_column`` (name, or 0-based position), else
    ``correlation`` or ``weight`` if present, else the third column if it is
    numeric, else 1.

    The C++ extension reads the file when its pandas reading is certain (see
    :func:`read_layer_codes`); the table is the same either way.
    """
    import pandas as pd

    names, a, b, w = read_layer_codes(path, weight_column, sep)
    return pd.DataFrame({"gene_a": pd.Series(names[a], dtype=object).astype(str),
                         "gene_b": pd.Series(names[b], dtype=object).astype(str),
                         "weight": w})


def read_layer_codes(path: Union[str, Path], weight_column: Union[str, int, None] = None,
                     sep: str = "\t"):
    """:func:`read_layer` as ``(names, a, b, w)``: gene names (object array),
    int64 codes of the two gene columns into ``names`` and float64 weights.

    The first record is read by pandas to decide the header and the weight
    column; the rows are then read by the C++ extension
    (``_gcen_core.network.read_layer_file``), which builds no string per
    row. Inputs it declines (quotes, a byte order mark, missing values,
    numbers pandas might parse differently, ...) and installs without the
    extension are read by pandas.
    """
    import pandas as pd

    path = require_file(path)
    core = extension()
    if core is not None:
        try:
            first = pd.read_csv(path, sep=sep, header=None, dtype=str, comment=None,
                                skip_blank_lines=True, nrows=1)
        except Exception as exc:
            raise GCENIOError(f"Failed to read layer file {path}: {exc}") from exc
        has_header, header, col = _layout(first, weight_column, path)
        if col is None or header.count(col) == 1:
            r = core.network.read_layer_file(str(path), sep, has_header, first.shape[1],
                                             -1 if col is None else header.index(col))
            if r["ok"]:
                return (np.asarray(r["names"], dtype=object), r["a"].astype(np.int64),
                        r["b"].astype(np.int64), r["w"])
    t = _read_layer_pandas(path, weight_column, sep)
    codes, names = pd.factorize(pd.concat([t["gene_a"], t["gene_b"]], ignore_index=True))
    m = len(t)
    return (np.asarray(names, dtype=object), codes[:m].astype(np.int64),
            codes[m:].astype(np.int64), t["weight"].to_numpy(dtype=np.float64))


def _layout(raw: pd.DataFrame, weight_column, path):
    """(has_header, column names, weight column name or None) of a layer
    file from its first records ``raw`` (read with header=None, dtype=str)."""
    import pandas as pd

    if raw.shape[1] < 2:
        raise GCENInputError(f"{path}: an edge list needs at least two columns.")
    has_header = bool(raw.shape[1] >= 3 and len(raw)
                      and pd.isna(pd.to_numeric(raw.iloc[0, 2], errors="coerce")))
    header = ([str(c) for c in raw.iloc[0]] if has_header
              else [f"col{i}" for i in range(raw.shape[1])])
    col = None
    if weight_column is not None:
        if isinstance(weight_column, int):
            if not 0 <= weight_column < raw.shape[1]:
                raise GCENInputError(f"{path}: no column {weight_column}.")
            col = header[weight_column]
        elif weight_column in header:
            col = weight_column
        else:
            raise GCENInputError(f"{path}: no column '{weight_column}' (has {header}).")
    else:
        for name in ("correlation", "weight"):
            if name in header:
                col = name
                break
        if col is None and raw.shape[1] >= 3:
            col = header[2]
    return has_header, header, col


def _read_layer_pandas(path: Path, weight_column, sep: str) -> pd.DataFrame:
    """:func:`read_layer` by pandas (the reference the C++ reader follows)."""
    import pandas as pd

    try:
        raw = pd.read_csv(path, sep=sep, header=None, dtype=str, comment=None,
                          skip_blank_lines=True)
    except Exception as exc:
        raise GCENIOError(f"Failed to read layer file {path}: {exc}") from exc
    has_header, header, col = _layout(raw, weight_column, path)
    if has_header:
        raw = raw.iloc[1:].reset_index(drop=True)
    raw.columns = header
    weight = (pd.to_numeric(raw[col], errors="coerce").abs() if col is not None
              else pd.Series(1.0, index=raw.index))
    if weight.isna().any():
        raise GCENInputError(f"{path}: column '{col}' has non-numeric weights.")
    return pd.DataFrame({"gene_a": raw.iloc[:, 0].astype(str),
                         "gene_b": raw.iloc[:, 1].astype(str),
                         "weight": weight.to_numpy(dtype=np.float64)})


def as_layer(layer: LayerLike) -> pd.DataFrame:
    """``gene_a, gene_b, weight`` table of a Network, DataFrame or file."""
    import pandas as pd

    from gcen2.models import Network

    if isinstance(layer, (str, Path)):
        return read_layer(layer)
    edges = layer.edges if isinstance(layer, Network) else layer
    if not isinstance(edges, pd.DataFrame):
        raise GCENInputError("A layer must be a Network, a DataFrame or an edge-list path.")
    missing = {"gene_a", "gene_b"} - set(edges.columns)
    if missing:
        raise GCENInputError(f"Layer table lacks columns {sorted(missing)}.")
    if "correlation" in edges.columns:
        w = edges["correlation"].abs()
    elif "weight" in edges.columns:
        w = edges["weight"].abs()
    else:
        w = pd.Series(1.0, index=edges.index)
    return pd.DataFrame({"gene_a": edges["gene_a"].astype(str).to_numpy(),
                         "gene_b": edges["gene_b"].astype(str).to_numpy(),
                         "weight": w.to_numpy(dtype=np.float64)})


def _layer_items(layers) -> list:
    if isinstance(layers, Mapping):
        items = list(layers.items())
    else:
        items = [(f"layer{i}", lay) for i, lay in enumerate(layers)]
    if not items:
        raise GCENInputError("At least one layer is required.")
    names = [str(k) for k, _ in items]
    if len(set(names)) != len(names):
        raise GCENInputError(f"Layer names must be unique, got {names}.")
    return [(n, lay) for n, (_, lay) in zip(names, items)]


def _layer_codes(layer: LayerLike):
    """``(names, a, b, w)`` of one layer (see :func:`read_layer_codes`)."""
    import pandas as pd

    if isinstance(layer, (str, Path)):
        return read_layer_codes(layer)
    t = as_layer(layer)
    codes, names = pd.factorize(pd.concat([t["gene_a"], t["gene_b"]], ignore_index=True))
    m = len(t)
    return (np.asarray(names, dtype=object), codes[:m].astype(np.int64),
            codes[m:].astype(np.int64), t["weight"].to_numpy(dtype=np.float64))


def _coded_layers(layers, threads: int = 4):
    """Layer names, all genes (sorted object array) and per layer the int64
    codes ``a``, ``b`` into those genes and the weights ``w``.

    Files are read in parallel (the C++ reader releases the GIL); only the
    small per-layer name lists are merged in Python.
    """
    from concurrent.futures import ThreadPoolExecutor

    items = _layer_items(layers)
    workers = max(1, min(int(threads), sum(isinstance(l, (str, Path)) for _, l in items)))
    if workers > 1:
        with ThreadPoolExecutor(workers) as pool:
            coded = list(pool.map(_layer_codes, [lay for _, lay in items]))
    else:
        coded = [_layer_codes(lay) for _, lay in items]
    local = [np.asarray(c[0], dtype=str) if len(c[0]) else np.empty(0, dtype=str)
             for c in coded]
    genes_str = np.unique(np.concatenate(local)) if local else np.empty(0, dtype=str)
    genes = genes_str.astype(object)
    out = []
    for names_l, (_, a, b, w) in zip(local, coded):
        remap = np.searchsorted(genes_str, names_l).astype(np.int64)
        out.append((remap[a], remap[b], np.asarray(w, dtype=np.float64)))
    return [n for n, _ in items], genes, out


def fuse(
    layers: Union[Mapping[str, LayerLike], list],
    weights: Union[Mapping[str, float], None] = None,
    normalize: str = "max",
    min_layers: int = 1,
    min_score: float = 0.0,
) -> Network:
    """Weighted fusion of heterogeneous networks into one network.

    Parameters
    ----------
    layers:
        ``{name: layer}`` (or a list, named ``layer0``, ``layer1``, ...).
        Each layer is a Network, a ``gene_a, gene_b[, weight]`` DataFrame or
        an edge-list path (see :func:`read_layer`). Edges are undirected;
        repeated edges within a layer keep their largest weight.
    weights:
        Layer importance ``{name: w >= 0}`` (default 1 for every layer).
    normalize:
        Per-layer edge-weight scaling before fusion: ``max`` (divide by the
        layer's largest weight, so all layers range up to 1), ``rank``
        (percentile rank in (0, 1], robust to differently distributed
        scores) or ``none``.
    min_layers:
        Keep edges present in at least this many layers.
    min_score:
        Keep edges whose fused score is at least this value.

    Returns
    -------
    Network
        Edges ``gene_a, gene_b, score, n_layers, layers`` and one column
        ``w_<name>`` per layer (its normalised weight, 0 when absent).
        ``score = sum_l weights[l] * w_l / sum_l weights[l]``; ``gene_a`` <=
        ``gene_b``; sorted by score (descending), then gene names.
    """
    import pandas as pd

    from gcen2.models import Network

    f = _fuse_arrays(layers, weights, normalize, min_layers, min_score)
    out = pd.DataFrame({
        "gene_a": f["genes"][f["a"]],
        "gene_b": f["genes"][f["b"]],
        "score": f["score"],
        "n_layers": f["n_layers"],
        "layers": np.asarray(f["labels"], dtype=object)[f["label_idx"]],
    })
    for k, n in enumerate(f["names"]):
        out[f"w_{n}"] = f["wide"][k]
    return Network(edges=out, method="fused", node_count=f["node_count"],
                   metadata={"layers": f["names"], "weights": f["alpha"],
                             "normalize": normalize})


def fuse_to_file(
    layers: Union[Mapping[str, LayerLike], list],
    output: Union[str, Path],
    weights: Union[Mapping[str, float], None] = None,
    normalize: str = "max",
    min_layers: int = 1,
    min_score: float = 0.0,
    threads: int = 4,
) -> dict:
    """:func:`fuse` written straight to ``output`` (used by ``gcen2 fuse``).

    Same file as ``fuse(...).to_edgelist(output)``, but the table is written
    by the C++ extension from gene codes, without building the gene-name
    columns in Python (falls back to that path without the extension).
    Returns ``node_count``, ``edge_count`` and ``layer_counts`` (edges per
    number of supporting layers).
    """
    core = extension()
    if core is None:
        net = fuse(layers, weights, normalize, min_layers, min_score)
        net.to_edgelist(output)
        counts = net.edges["n_layers"].value_counts().sort_index()
        return {"node_count": net.node_count, "edge_count": net.edge_count,
                "layer_counts": {int(k): int(v) for k, v in counts.items()}}

    f = _fuse_arrays(layers, weights, normalize, min_layers, min_score, threads)
    columns = ["gene_a", "gene_b", "score", "n_layers", "layers"] + \
        [f"w_{n}" for n in f["names"]]
    table = [("float", f["score"]), ("int", f["n_layers"]),
             ("label", f["label_idx"], f["labels"])]
    table += [("float", col) for col in f["wide"]]
    core.network.write_coded_table(str(output), tsv_header(columns),
                                   [str(g) for g in f["genes"]], f["a"], f["b"], table,
                                   max(1, int(threads)))
    counts = np.bincount(f["n_layers"], minlength=len(f["names"]) + 1)
    return {"node_count": f["node_count"], "edge_count": len(f["score"]),
            "layer_counts": {k: int(c) for k, c in enumerate(counts) if c}}


def _unique_ints(x: np.ndarray, return_inverse: bool = False):
    """``np.unique`` of an integer array by sorting (numpy >= 2.3 hashes
    integers in np.unique, several times slower on millions of values)."""
    if return_inverse:
        order = np.argsort(x)
        xs = x[order]
    else:
        xs = np.sort(x)
    first = np.empty(len(xs), dtype=bool)
    if len(xs):
        first[0] = True
        np.not_equal(xs[1:], xs[:-1], out=first[1:])
    uniq = xs[first]
    if not return_inverse:
        return uniq
    inv = np.empty(len(x), dtype=np.int64)
    inv[order] = np.cumsum(first) - 1
    return uniq, inv


def _fuse_arrays(layers, weights, normalize, min_layers, min_score, threads=4) -> dict:
    """The fused table as arrays: gene codes ``a`` / ``b`` into ``genes``
    (sorted names), ``score``, ``n_layers``, ``label_idx`` into ``labels``
    and ``wide`` (one weight array per layer), rows already filtered and
    ordered. Computed by the C++ extension (src/layer_fuse.cpp), else by
    :func:`_fuse_arrays_numpy`, with the same values."""
    if normalize not in NORMALIZE:
        raise GCENInputError(f"normalize must be one of {NORMALIZE}, got '{normalize}'.")
    names, genes, coded = _coded_layers(layers, threads)
    tables = dict.fromkeys(names)
    alpha = dict.fromkeys(tables, 1.0)
    if weights is not None:
        unknown = set(map(str, weights)) - set(tables)
        if unknown:
            raise GCENInputError(f"Weights given for unknown layers {sorted(unknown)}.")
        alpha.update({str(k): float(v) for k, v in weights.items()})
    if any(v < 0 for v in alpha.values()) or sum(alpha.values()) <= 0:
        raise GCENInputError("Layer weights must be >= 0 with a positive sum.")
    if min_layers < 1 or min_layers > len(tables):
        raise GCENInputError(f"min_layers must be in [1, {len(tables)}], got {min_layers}.")

    core = extension()
    f = None
    if core is not None:
        r = core.network.fuse_layer_codes(
            coded, len(genes), normalize, [alpha[n] for n in names], int(min_layers),
            float(min_score), max(1, int(threads)))
        if r["ok"]:
            f = r
            f["labels"] = [",".join(n for k, n in enumerate(names) if v >> k & 1)
                           for v in r["patterns"]]
    if f is None:
        f = _fuse_arrays_numpy(names, len(genes), coded, normalize, alpha, min_layers,
                               min_score)
    f.update(genes=genes, names=names, alpha=alpha)
    return f


def _fuse_arrays_numpy(names, n_genes, coded, normalize, alpha, min_layers,
                       min_score) -> dict:
    """numpy reference of :func:`_fuse_arrays` (and fallback of the C++
    ``fuse_layer_codes``)."""
    import pandas as pd

    # an undirected pair as one int64 (lo * n + hi); genes are codes in sorted
    # name order, so pair order is (gene_a, gene_b) name order
    per_layer = []
    for a, b, w_rows in coded:
        keep = a != b  # self-loops carry no pair information
        lo, hi = np.minimum(a, b)[keep], np.maximum(a, b)[keep]
        codes, inv = _unique_ints(lo.astype(np.int64) * n_genes + hi, return_inverse=True)
        w = np.full(len(codes), -np.inf)
        np.maximum.at(w, inv, w_rows[keep])
        if normalize == "max" and len(w) and w.max() > 0:
            w = w / w.max()
        elif normalize == "rank" and len(w):
            w = pd.Series(w).rank(method="average", pct=True).to_numpy()
        per_layer.append((codes, w))

    union = _unique_ints(np.concatenate([c for c, _ in per_layer]))
    wide = np.zeros((len(union), len(names)))
    pres = np.zeros((len(union), len(names)), dtype=bool)
    for k, (codes, w) in enumerate(per_layer):
        pos = np.searchsorted(union, codes)
        wide[pos, k] = w
        pres[pos, k] = True
    # summed term by term in layer order (as the C++ code; a BLAS matrix
    # product may round differently)
    total = sum(alpha.values())
    score = np.zeros(len(union))
    for k, n in enumerate(names):
        score = score + wide[:, k] * alpha[n]
    score = score / total
    n_layers = pres.sum(axis=1)
    keep = (n_layers >= min_layers) & (score >= min_score)
    # union is in (gene_a, gene_b) order, so a stable sort by score gives the
    # documented order
    order = np.flatnonzero(keep)[np.argsort(-score[keep], kind="stable")]
    kept = union[order]
    a, b = kept // n_genes, kept % n_genes
    # "a,b" label per presence pattern (few distinct patterns)
    bits = pres[order] @ (1 << np.arange(len(names), dtype=np.int64))
    uniq_bits, label_idx = _unique_ints(bits, return_inverse=True)
    labels = [",".join(n for k, n in enumerate(names) if int(v) >> k & 1) for v in uniq_bits]
    return {
        "a": a, "b": b, "score": score[order],
        "n_layers": n_layers[order].astype(np.int64), "label_idx": label_idx.astype(np.int64),
        "labels": labels, "wide": [np.ascontiguousarray(wide[order, k])
                                   for k in range(len(names))],
        "node_count": len(np.union1d(a, b)),
    }


def rwr_multiplex(
    layers: Union[Mapping[str, LayerLike], list],
    seed_genes: list[str],
    gamma: float = 0.5,
    delta: float = 0.5,
    weighted: bool = False,
    aggregate: str = "geometric",
    prob_cutoff: float = 0.0,
    transition: str = "supra",
    max_iter: int = 1_000_000,
    tol: float = 1e-10,
) -> RWRResult:
    """Random walk with restart on a multiplex network (RWR-M).

    Every gene of any layer has one node per layer. From a node, the walker
    restarts with probability ``gamma`` (to the seeds, spread equally over
    the layers); otherwise it moves along the supra-adjacency matrix: the
    layers' adjacency matrices scaled by ``1 - delta`` on the diagonal
    blocks, ``delta / (L - 1)`` between the copies of a gene in two layers.
    How that matrix becomes transition probabilities is ``transition``:

    * ``supra`` (default): the whole matrix is column-normalised, as in
      RandomWalkRestartMH (Valdeolivas et al., the authors' implementation
      of RWR-M), whose scores it reproduces. A gene of degree ``d`` jumps
      with probability ``delta / ((1 - delta) d + delta)``, so hubs stay in
      their layer.
    * ``layer``: each layer is column-normalised first, so every gene
      follows an edge with probability ``1 - delta`` and jumps with
      probability ``delta``, whatever its degree.

    A node without edges in its layer always jumps. With one layer both are
    :func:`gcen2.rwr.predict`.

    Parameters
    ----------
    layers:
        ``{name: layer}`` or a list, as in :func:`fuse`.
    seed_genes:
        Known genes; those absent from every layer are ignored.
    gamma:
        Restart probability, in (0, 1).
    delta:
        Weight of the jumps between layers, in [0, 1] (their probability with
        ``transition="layer"``).
    weighted:
        Use the layer edge weights (see the module docstring) for intra-layer
        transitions instead of 1 per edge.
    aggregate:
        Per-gene score from its per-layer probabilities: ``geometric`` mean
        (as RWR-M; favours genes supported in every layer), ``arithmetic``
        mean, or ``sum`` (the total probability, summing to 1 over genes).
    prob_cutoff:
        Minimum aggregated score of a reported gene.
    transition:
        ``supra`` or ``layer``, see above.

    Returns
    -------
    RWRResult
        ``ranked_genes`` holds ``gene, probability`` (the aggregated score)
        and ``p_<layer>`` per layer, seeds excluded, sorted by probability.
    """
    import pandas as pd

    from gcen2.models import RWRResult
    from gcen2.rwr import _check_inputs, _rwr_python

    _check_inputs(gamma, seed_genes)
    if not 0.0 <= delta <= 1.0:
        raise GCENInputError(f"delta must be in [0, 1], got {delta}.")
    if aggregate not in AGGREGATE:
        raise GCENInputError(f"aggregate must be one of {AGGREGATE}, got '{aggregate}'.")
    if transition not in TRANSITION:
        raise GCENInputError(f"transition must be one of {TRANSITION}, got '{transition}'.")
    layer_names, genes, coded = _coded_layers(layers)
    genes = genes.tolist()
    L = len(layer_names)
    n = len(genes)

    W = _transition_matrix(coded, n, delta, weighted, transition)

    idx = {g: i for i, g in enumerate(genes)}
    seeds = [g for g in dict.fromkeys(seed_genes) if g in idx]
    if not seeds:
        raise GCENInputError("None of the seed genes are present in the network.")
    P0 = np.zeros(L * n)
    for k in range(L):
        P0[[idx[g] + k * n for g in seeds]] = 1.0
    P0 /= P0.sum()

    core = extension()
    if core is not None:
        PT, steps = core.rwr.random_walk_restart(
            W.indptr.astype(np.int64), W.indices.astype(np.int32), W.data, P0,
            gamma, max_iter, tol)
        PT = np.asarray(PT)
    else:
        PT, steps = _rwr_python(W, P0, gamma, max_iter, tol)

    per_layer = PT.reshape(L, n)
    if aggregate == "geometric":
        with np.errstate(divide="ignore"):
            score = np.exp(np.log(per_layer).mean(axis=0))
    elif aggregate == "arithmetic":
        score = per_layer.mean(axis=0)
    else:
        score = per_layer.sum(axis=0)

    df = pd.DataFrame({"gene": genes, "probability": score})
    for k, name in enumerate(layer_names):
        df[f"p_{name}"] = per_layer[k]
    seed_set = set(seeds)
    df = df[~df["gene"].isin(seed_set) & (df["probability"] >= prob_cutoff)]
    df = df.sort_values(["probability", "gene"], ascending=[False, True]).reset_index(drop=True)
    return RWRResult(ranked_genes=df, seed_genes=list(seed_genes), gamma=gamma, steps=int(steps))


def _transition_matrix(coded, n: int, delta: float, weighted: bool, transition: str):
    """Column-stochastic CSR matrix of the walk on the multiplex network.

    Supra-adjacency: ``(1 - delta)`` times each layer on the diagonal blocks
    (weights divided by the layer's largest, as RandomWalkRestartMH rescales
    them), ``delta / (L - 1)`` between the copies of a gene (1 / (L - 1) for
    a gene without edges in its layer); then column-normalised, each layer
    block first for ``transition="layer"``.
    """
    from scipy.sparse import coo_matrix

    L = len(coded)
    stay = 1.0 - delta if L > 1 else 1.0
    rows, cols, vals = [], [], []
    has_edges = np.zeros((L, n), dtype=bool)
    for k, (i, j, w_rows) in enumerate(coded):
        w = w_rows if weighted else np.ones(len(i))
        if weighted and transition == "supra" and len(w):
            w = w / w.max()
        i, j, w = np.concatenate([i, j]), np.concatenate([j, i]), np.concatenate([w, w])
        A = coo_matrix((w, (i, j)), shape=(n, n)).tocsr()
        col_sums = np.asarray(A.sum(axis=0)).ravel()
        has_edges[k] = col_sums > 0
        A = A.tocoo()
        rows.append(A.row + k * n)
        cols.append(A.col + k * n)
        vals.append(A.data * stay / (col_sums[A.col] if transition == "layer" else 1.0))
    for k in range(L if L > 1 else 0):
        jump = np.where(has_edges[k], delta, 1.0) / (L - 1)
        for m in range(L):
            if m != k:
                rows.append(np.arange(n) + m * n)
                cols.append(np.arange(n) + k * n)
                vals.append(jump)
    W = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(L * n, L * n)).tocsr()
    if transition == "supra":
        W = W.multiply(1.0 / np.asarray(W.sum(axis=0)).ravel()).tocsr()
    return W
