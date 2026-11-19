"""File I/O utilities for GCEN2."""

from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Union

if TYPE_CHECKING:
    import pandas as pd

    from gcen2.models import ModuleSet, Network

from gcen2._matrix_io import _detect_header, require_file
from gcen2.exceptions import GCENInputError, GCENIOError
from gcen2.models import ExpressionMatrix


def read_tsv(
    path: Union[str, Path], index_col: int = 0, header: bool | None = None
) -> ExpressionMatrix:
    """Read a tab-separated expression matrix (genes × samples).

    Parameters
    ----------
    path:
        Path to the TSV file. First column is expected to be gene identifiers.
    index_col:
        Column to use as the row index (gene names). Default 0.
    header:
        Whether the first line holds sample names. ``None`` (default)
        auto-detects: header-less GCEN v1 matrices (first line is numeric data)
        are read without losing the first gene. Pass ``True`` explicitly if
        your sample names are purely numeric.

    Returns
    -------
    ExpressionMatrix
    """
    import pandas as pd

    path = require_file(path)
    try:
        has_header = _detect_header(path, "\t") if header is None else header
        df = pd.read_csv(
            path, sep="\t", index_col=index_col, header=0 if has_header else None
        )
    except Exception as exc:
        raise GCENIOError(f"Failed to read TSV file {path}: {exc}") from exc

    if not has_header:
        df.columns = [f"sample_{i + 1}" for i in range(df.shape[1])]
        df.index.name = None

    import numpy as np

    return ExpressionMatrix(
        data=df,
        gene_names=np.array(df.index, dtype=str),
        sample_names=np.array(df.columns, dtype=str),
        is_single_cell=False,
    )


def read_csv(
    path: Union[str, Path], index_col: int = 0, header: bool | None = None
) -> ExpressionMatrix:
    """Read a comma-separated expression matrix (genes × samples).

    Quoted fields are handled correctly. ``header`` behaves as in
    :func:`read_tsv` (auto-detects header-less GCEN v1 matrices).
    """
    import pandas as pd

    path = require_file(path)
    try:
        has_header = _detect_header(path, ",") if header is None else header
        df = pd.read_csv(
            path, sep=",", index_col=index_col, quoting=0, header=0 if has_header else None
        )
    except Exception as exc:
        raise GCENIOError(f"Failed to read CSV file {path}: {exc}") from exc

    if not has_header:
        df.columns = [f"sample_{i + 1}" for i in range(df.shape[1])]
        df.index.name = None

    import numpy as np

    return ExpressionMatrix(
        data=df,
        gene_names=np.array(df.index, dtype=str),
        sample_names=np.array(df.columns, dtype=str),
        is_single_cell=False,
    )


def read_h5ad(path: Union[str, Path]) -> ExpressionMatrix:
    """Read an AnnData (.h5ad) file for single-cell data."""
    try:
        import anndata as ad
    except ImportError as exc:
        raise GCENInputError(
            "anndata is required to read .h5ad files. "
            "Install with: pip install gcen2[single-cell]"
        ) from exc

    path = require_file(path)

    try:
        adata = ad.read_h5ad(path)
    except Exception as exc:
        raise GCENIOError(f"Failed to read h5ad file {path}: {exc}") from exc

    import numpy as np

    return ExpressionMatrix(
        data=adata,
        gene_names=np.array(adata.var_names, dtype=str),
        sample_names=np.array(adata.obs_names, dtype=str),
        is_single_cell=True,
    )


def read_network(path: Union[str, Path], genes_only: bool = False) -> Network:
    """Read a GCEN-format edge list (TSV: gene_a, gene_b, correlation, p_value).

    Compatible with GCEN v1 ``.network`` files and GCEN2 ``network.build()`` output.

    ``genes_only=True`` reads just the two gene columns, as categoricals
    (each gene name is stored once). Enough for annotation and needs much
    less memory for large networks.
    """
    import pandas as pd

    from gcen2.models import Network

    path = require_file(path)
    try:
        first = pd.read_csv(path, sep="\t", nrows=0)
    except Exception as exc:
        raise GCENIOError(f"Failed to read network file {path}: {exc}") from exc

    from gcen2._matrix_io import network_columns

    names, has_header = network_columns(list(first.columns))
    return Network(edges=_read_network_table(path, names, has_header, genes_only, path))


def _read_network_table(src, names: list[str], has_header: bool, genes_only: bool = False,
                        label=None) -> pd.DataFrame:
    """Parse an edge list (path or buffer) whose columns are ``names``
    (see :func:`gcen2._matrix_io.network_columns`).

    Gene columns are read as strings, so numeric gene IDs (e.g. Entrez)
    keep their text and match string gene lists.
    """
    import pandas as pd

    header = 0 if has_header else None
    kwargs = {}
    if genes_only:
        kwargs = dict(usecols=[0, 1], dtype={names[0]: "category", names[1]: "category"})
        names = names[:2]
    else:
        kwargs = dict(dtype={names[0]: str, names[1]: str})
    try:
        return pd.read_csv(src, sep="\t", header=header, names=names, **kwargs)
    except Exception as exc:
        raise GCENIOError(f"Failed to read network file {label}: {exc}") from exc


def read_gene_list(path: Union[str, Path]) -> list[str]:
    """Read a gene list file, one gene per line.

    GCEN v1 compatible: blank lines and lines starting with ``#`` (headers /
    comments) are skipped, and only the first tab/whitespace-separated column
    is used, so v1 seed files such as ``#gene<TAB>probability`` with
    ``GENE<TAB>1.0`` rows are read correctly.
    """
    path = require_file(path)
    from gcen2._matrix_io import iter_lines

    genes = []
    for line in iter_lines(path):
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        genes.append(line.split()[0])
    return genes


def read_modules(path: Union[str, Path]) -> ModuleSet:
    """Read a module file in GCEN v1 or GCEN2 format.

    * v1 (``module_identify``): one module per line, genes tab-separated.
    * GCEN2 (``gcen2 module``): header ``module<TAB>gene`` then one
      ``module_id<TAB>gene`` pair per line.

    Module order follows the file (v1 numbers annotation output ``module_<i>``
    by line).
    """
    from gcen2.models import ModuleSet

    path = require_file(path)
    from gcen2._matrix_io import iter_lines

    lines = [l.rstrip("\n\r") for l in iter_lines(path) if l.strip() and not l.startswith("#")]
    if lines and lines[0].split("\t")[:2] == ["module", "gene"]:
        groups: dict[str, set[str]] = {}
        for line in lines[1:]:
            mod, gene = line.split("\t")[:2]
            groups.setdefault(mod.strip(), set()).add(gene.strip())
        modules = list(groups.values())
    else:
        modules = [{g for g in line.split("\t") if g.strip()} for line in lines]
    return ModuleSet(modules=modules)


def write_tsv(df: pd.DataFrame, path: Union[str, Path]) -> None:
    """Write a DataFrame to a tab-separated file."""
    try:
        df.to_csv(path, sep="\t", lineterminator="\n")
    except Exception as exc:
        raise GCENIOError(f"Failed to write TSV file {path}: {exc}") from exc


def write_csv(df: pd.DataFrame, path: Union[str, Path]) -> None:
    """Write a DataFrame to a comma-separated file."""
    try:
        df.to_csv(path, sep=",", lineterminator="\n")
    except Exception as exc:
        raise GCENIOError(f"Failed to write CSV file {path}: {exc}") from exc
