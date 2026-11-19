"""Reading GCEN2 result files for plotting. Every reader also takes the
object itself (a DataFrame / Series), so plots work on files and in memory."""

from __future__ import annotations

import pandas as pd

from gcen2._matrix_io import require_file
from gcen2.exceptions import GCENInputError


def table(x, index_col=None, **kw) -> pd.DataFrame:
    """A DataFrame, or a TSV file read into one."""
    if isinstance(x, pd.DataFrame):
        return x
    return pd.read_csv(require_file(x), sep="\t", index_col=index_col, **kw)


def modules(x) -> pd.Series:
    """gene -> module name: "ME<k>" for GCEN2's numbered modules
    (modules.tsv: module, gene), the colour for hdWGCNA tables (gene /
    gene_name, module), or a Series / dict as given."""
    if isinstance(x, pd.Series):
        return x.astype(str)
    if isinstance(x, dict):
        return pd.Series(x).astype(str)
    t = table(x, dtype=str)
    cols = {c.lower(): c for c in t.columns}
    g = cols.get("gene") or cols.get("gene_name")
    if g is None or "module" not in cols:
        raise GCENInputError("a modules table needs gene and module columns.")
    m = t[cols["module"]].astype(str)
    if m.str.fullmatch(r"\d+").all():
        m = "ME" + m
    return pd.Series(m.to_numpy(), index=t[g].astype(str).to_numpy(), name="module")


def eigengenes(x) -> pd.DataFrame:
    """samples (or cells) x eigengenes (eigengenes.tsv, cell_eigengenes.tsv,
    MEs from hdWGCNA)."""
    return table(x, index_col=0)


def kme(x) -> pd.DataFrame:
    """kme.tsv / cell_kme.tsv: genes x (module, kME<k>, p.kME<k>)."""
    t = x.copy() if isinstance(x, pd.DataFrame) else table(x, index_col=0,
                                                           dtype={"module": str})
    t.index = t.index.map(str)
    if "module" in t.columns:
        t["module"] = t["module"].astype(object).where(t["module"].notna(), "NA").astype(str)
    return t
