"""Text files without pandas: expression matrices, edge-list columns, TSV fields.

Importing pandas costs about 50 MB and several seconds of start-up, more than
most commands need, so the commands that stream files (``build``, ``module``,
``stats``, ``merge``, ...) read and write them with these helpers.
:func:`read_matrix` matches :func:`gcen2.io.read_tsv` (same header
detection, values and gene order; gene names kept exactly as written), and
:func:`tsv_field` quotes a field as ``DataFrame.to_csv(sep="\t")`` does.
"""

from __future__ import annotations

import warnings
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    import numpy as np

# numpy is imported inside read_matrix: header detection and iter_lines are
# also used by the numpy-free build path (gcen2._build_file).

# pandas' default missing-value strings (pandas.read_csv ``na_values``)
_NA_STRINGS = frozenset({
    "", "#N/A", "#N/A N/A", "#NA", "-1.#IND", "-1.#QNAN", "-NaN", "-nan",
    "1.#IND", "1.#QNAN", "<NA>", "N/A", "NA", "NULL", "NaN", "None", "n/a",
    "nan", "null",
})


def require_file(path) -> Path:
    """``path`` as a :class:`~pathlib.Path`; :class:`gcen2.GCENIOError` when
    it does not exist or is a directory."""
    path = Path(path)
    if not path.exists() or path.is_dir():
        from gcen2.exceptions import GCENIOError

        raise GCENIOError(f"File not found: {path}")
    return path


def tsv_field(value) -> str:
    """A field quoted as pandas ``to_csv`` quotes it (QUOTE_MINIMAL, tab separator)."""
    value = str(value)
    if any(ch in value for ch in '\t"\n\r'):
        return '"' + value.replace('"', '""') + '"'
    return value


def tsv_header(names) -> str:
    """A header line as ``DataFrame.to_csv(sep="\t")`` writes it."""
    return "\t".join(tsv_field(c) for c in names) + "\n"


def write_modules(path, modules) -> None:
    """A module file (``gcen2 module`` format): a ``module\tgene`` header and
    one line per gene, modules numbered 0, 1, ... in the given order, genes in
    the given order."""
    with open(path, "w", encoding="utf-8", newline="") as fh:
        fh.write("module\tgene\n")
        for i, mod in enumerate(modules):
            fh.write("".join(f"{i}\t{tsv_field(g)}\n" for g in mod))


def iter_lines(path, encoding: str = "utf-8"):
    """Yield the lines of a text file (line endings kept, ``\r\n`` included).

    Reads in binary with a 1 MB buffer and decodes each line. Iterating a
    text-mode file is several times slower on network/9p file systems
    (e.g. WSL ``/mnt/*``: 0.9 s instead of 0.26 s for the 34 MB GO OBO file).
    Lines are split on ``\n`` only (a lone ``\r`` does not end a line).
    A gzip-compressed file (e.g. ``go-basic.obo.gz``) is decompressed on the
    fly, whatever its name.
    """
    with open(path, "rb", buffering=1 << 20) as fh:
        if fh.peek(2)[:2] == b"\x1f\x8b":
            import gzip

            fh = gzip.GzipFile(fileobj=fh)  # noqa: PLW2901
        for raw in fh:
            yield raw.decode(encoding)


def _is_number(token: str) -> bool:
    try:
        float(token)
    except ValueError:
        return False
    return True


def _detect_header(path: Path, sep: str) -> bool:
    """Decide whether the first line of an expression matrix is a header.

    Mirrors GCEN v1 (lines starting with ``#`` are headers/annotations; all
    other lines are data) while still accepting ordinary pandas-style headers:

    * first line starts with ``#``                   -> header
    * every field after the first parses as a number -> no header (v1 data row)
    * otherwise                                      -> header

    The line is split as the reader splits it (quoted fields may hold ``sep``).
    """
    import csv

    with open(path, encoding="utf-8") as fh:
        first = fh.readline().rstrip("\r\n")
    if first.startswith("#"):
        return True
    fields = [f for f in next(csv.reader([first], delimiter=sep), [])[1:] if f != ""]
    if fields and all(_is_number(f) for f in fields):
        if all(float(f).is_integer() for f in fields):
            warnings.warn(
                f"{path.name}: first line has only integer values after the first "
                "column; treating it as data (no header). If these are sample "
                "names, pass header=True.",
                stacklevel=3,
            )
        return False
    return True


def _unquote(field: str) -> str:
    if len(field) >= 2 and field[0] == '"' and field[-1] == '"':
        return field[1:-1]
    return field


def _to_float(field: str) -> float:
    field = _unquote(field.strip())
    if field in _NA_STRINGS:
        return float("nan")
    return float(field)


def read_matrix(
    path, sep: str = "\t", header: bool | None = None
) -> tuple[np.ndarray, list[str], list[str]]:
    """Read a genes x samples matrix; return ``(values, gene_names, sample_names)``.

    ``values`` is a C-contiguous float64 array. ``header`` works as in
    :func:`gcen2.io.read_tsv`. Blank lines are skipped, missing values
    (pandas' default NA strings) become NaN, and short rows are padded with
    NaN. Header-less files get sample names ``sample_1`` ... ``sample_k``.
    """
    import numpy as np

    from gcen2.exceptions import GCENIOError

    path = require_file(path)
    has_header = _detect_header(path, sep) if header is None else header

    # Rows are converted to float arrays as they are read; keeping the split
    # strings instead would cost several times the matrix itself.
    genes: list[str] = []
    rows: list[np.ndarray] = []
    head: list[str] | None = None
    lines = iter_lines(path)
    if has_header:
        head = [_unquote(f) for f in next(lines, "").rstrip("\r\n").split(sep)]
    for line in lines:
        line = line.rstrip("\r\n")
        if not line.strip():
            continue
        fields = line.split(sep)
        genes.append(_unquote(fields[0]))
        try:
            try:  # fast path: plain numbers convert directly
                rows.append(np.array(fields[1:], dtype=np.float64))
            except ValueError:
                rows.append(np.array([_to_float(f) for f in fields[1:]], dtype=np.float64))
        except ValueError as exc:
            raise GCENIOError(f"Failed to read matrix file {path}: line for gene "
                              f"'{genes[-1]}': {exc}") from exc

    width = max((len(r) for r in rows), default=0)
    if head is not None:
        # pandas: a header one field shorter than the rows names only the samples
        width = max(width, len(head) - 1)
    values = np.full((len(rows), width), np.nan, dtype=np.float64)
    for i, r in enumerate(rows):
        values[i, :len(r)] = r
    del rows

    if head is None:
        samples = [f"sample_{i + 1}" for i in range(width)]
    elif len(head) == width:
        samples = head
    else:
        samples = head[1:]
    if head is not None and len(samples) != width:
        raise GCENIOError(f"Failed to read matrix file {path}: header has "
                          f"{len(head)} fields but rows have up to {width + 1}")
    return values, genes, samples


# First-line column names that mark a network file header (after lower-casing
# and replacing spaces by "_"); shared with :func:`gcen2.io.read_network`.
NETWORK_HEADER_KEYS = frozenset({"gene_a", "#node1", "node1", "gene", "genea"})


def network_has_header(path) -> bool:
    """Whether an edge-list file starts with a header line (as ``read_network``)."""
    with open(path, encoding="utf-8") as fh:
        first = fh.readline().rstrip("\r\n")
    columns = {_unquote(c).lower().replace(" ", "_") for c in first.split("\t")}
    return bool(NETWORK_HEADER_KEYS & columns)


_HEADERLESS_COLUMNS = {
    4: ["gene_a", "gene_b", "correlation", "p_value"],
    5: ["gene_a", "gene_b", "correlation", "p_value", "fdr"],
    6: ["gene_a", "gene_b", "correlation", "p_value", "fdr", "direction"],
}


def network_columns(first_columns: list[str]) -> tuple[list[str], bool]:
    """Canonical column names of an edge-list file from its first line's
    fields, and whether that line is a header (``read_network`` rules).

    Header names are lower-cased with spaces as underscores and GCEN v1
    names mapped (``#node1``/``node1`` -> ``gene_a``, ``node2`` -> ``gene_b``,
    ``pvalue``/``p-value`` -> ``p_value``); a header-less file gets names by
    its column count.
    """
    columns = [c.lower().replace(" ", "_") for c in first_columns]
    col_set = set(columns)
    if not (NETWORK_HEADER_KEYS & col_set):
        n = len(columns)
        return _HEADERLESS_COLUMNS.get(n, [f"c{i}" for i in range(n)]), False
    rename = {}
    if "#node1" in col_set and "gene_a" not in col_set:
        rename["#node1"] = "gene_a"
    if "node1" in col_set and "gene_a" not in col_set:
        rename["node1"] = "gene_a"
    if "node2" in col_set and "gene_b" not in col_set:
        rename["node2"] = "gene_b"
    for alias in ("pvalue", "p-value"):  # v1 header: "#node1 node2 correlation p-value"
        if alias in col_set and "p_value" not in col_set:
            rename[alias] = "p_value"
    return [rename.get(c, c) for c in columns], True


def network_file_columns(path) -> tuple[list[str], bool]:
    """:func:`network_columns` of a file, read without pandas (the first
    line split at tabs, fields unquoted)."""
    with open(path, encoding="utf-8") as fh:
        first = fh.readline().rstrip("\r\n")
    return network_columns([_unquote(c) for c in first.split("\t")])


def correlation_specs(paths) -> list | None:
    """``(path, has_header, correlation column)`` of each edge-list file for
    the extension's edge-file readers (merge, diff), or None when a file's
    gene columns are not the first two or it has no correlation column (the
    pandas pipeline handles those)."""
    specs = []
    for p in paths:
        names, has_header = network_file_columns(p)
        if names[:2] != ["gene_a", "gene_b"] or "correlation" not in names:
            return None
        specs.append((str(p), has_header, names.index("correlation")))
    return specs
