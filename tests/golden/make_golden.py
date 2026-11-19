"""Write the golden inputs (``--inputs``: once; they are committed) and the
expected outputs of every step (``expected/``, gzipped, with manifest.json).

    python -m tests.golden.make_golden            # expected outputs only
    python -m tests.golden.make_golden --inputs   # inputs too (from data/example)

Regenerate the expected outputs only for an intended change of results, and
say so in the commit message: they are the baseline of every optimisation."""

from __future__ import annotations

import argparse
import gzip
import json
import shutil
import tempfile
from pathlib import Path

import numpy as np
import pandas as pd

from tests.golden.golden import EXPECTED, INPUTS, run_steps

ROOT = Path(__file__).resolve().parents[2]


def portable(a):
    """Strings of obs / var (index and columns) as object arrays, so the file
    is written with the string-array encoding that older anndata (0.8-0.10,
    e.g. on Python 3.9) reads; pandas 3's str dtype would be written as
    nullable-string-array, which they cannot read."""
    import pandas as pd

    for df in (a.obs, a.var):
        df.index = pd.Index(np.asarray(df.index, dtype=object), dtype=object,
                            name=df.index.name)
        for c in df.columns:
            if pd.api.types.is_string_dtype(df[c].dtype) and \
                    not isinstance(df[c].dtype, pd.CategoricalDtype):
                df[c] = df[c].astype(object)
    return a


def make_inputs(d: Path) -> None:
    d.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(20261116)
    # bulk: 300 well-expressed genes of the GCEN sample data (15 samples)
    expr = pd.read_csv(ROOT / "data" / "example" / "bulk" / "gene_expr.tsv", sep="\t", header=None,
                       index_col=0)  # GCEN v1 format: no header line
    expr.index.name = "gene"
    expr.columns = [f"S{i + 1:02d}" for i in range(expr.shape[1])]
    expr = expr.loc[expr.mean(axis=1).sort_values(ascending=False, kind="stable").index[:300]]
    expr.to_csv(d / "bulk.tsv", sep="\t", lineterminator="\n")
    expr.iloc[:60].to_csv(d / "bulk_small.tsv", sep="\t", lineterminator="\n")
    norm = expr.div(expr.sum(axis=0), axis=1) * 1e6
    norm.iloc[:, :8].to_csv(d / "set1.tsv", sep="\t", lineterminator="\n")
    norm.iloc[:, 8:].to_csv(d / "set2.tsv", sep="\t", lineterminator="\n")
    pd.DataFrame({"weight": rng.normal(70, 10, 15).round(1),
                  "group": np.tile([0, 1], 8)[:15]}, index=expr.columns).to_csv(
        d / "traits.tsv", sep="\t", lineterminator="\n")
    genes = list(expr.index)
    (d / "seeds.txt").write_text("\n".join(genes[:3]) + "\n", encoding="utf-8")
    # a small GO: 3 namespaces, 60 terms each in a random DAG, genes annotated
    # to random leaves (the enrichment propagates to the ancestors)
    lines = ["format-version: 1.2", ""]
    terms = []
    for ns_i, ns in enumerate(("biological_process", "molecular_function",
                               "cellular_component")):
        ids = [f"GO:{ns_i + 1}{i:06d}" for i in range(60)]
        for i, t in enumerate(ids):
            lines += ["[Term]", f"id: {t}", f"name: {ns.split('_')[0]} term {i}",
                      f"namespace: {ns}"]
            if i:
                for p in sorted(set(rng.choice(i, size=min(i, 1 + (i % 2)), replace=False))):
                    lines.append(f"is_a: {ids[p]} ! parent")
            lines.append("")
        terms += ids
    (d / "go.obo").write_text("\n".join(lines), encoding="utf-8")
    assoc = []
    for i, g in enumerate(genes):
        if i % 5 == 4:  # unannotated: annotate -n infers them from their neighbours
            continue
        k = rng.integers(1, 5)
        assoc.append("\t".join([g, *sorted(set(rng.choice(terms, size=k)))]))
    (d / "gene_go.assoc").write_text("\n".join(assoc) + "\n", encoding="utf-8")
    (d / "query.txt").write_text("\n".join(rng.choice(genes, 40, replace=False)) + "\n", encoding="utf-8")
    (d / "background.txt").write_text("\n".join(genes) + "\n", encoding="utf-8")
    mods = [genes[i:i + 25] for i in range(0, 100, 25)]
    (d / "annot_modules.txt").write_text("\n".join("\t".join(m) for m in mods) + "\n", encoding="utf-8")
    a, b = np.triu_indices(80, 1)
    keep = rng.random(a.size) < 0.08
    pd.DataFrame({"gene_a": np.array(genes)[a[keep]], "gene_b": np.array(genes)[b[keep]],
                  "correlation": rng.uniform(0.6, 1, keep.sum()).round(6),
                  "p_value": 1e-5}).to_csv(d / "annot.network", sep="\t", index=False, lineterminator="\n")
    _single_cell(d, rng)
    _spatial(d, rng)


def _single_cell(d: Path, rng) -> None:
    import anndata as ad
    import scipy.sparse as sp

    def cells(n, seed_shift, prefix):
        n_mod, per, noise, rare = 3, 25, 40, 10
        n_genes = n_mod * per + noise + rare
        x = np.zeros((n, n_genes), dtype=np.float32)
        types = np.array(["A", "B"])[np.arange(n) % 2]
        for t in ("A", "B"):
            c = np.flatnonzero(types == t)
            f = rng.normal(size=(len(c), n_mod))
            logmu = np.log(rng.uniform(0.5, 2.0, n_genes))[None, :].repeat(len(c), 0)
            for m in range(n_mod):
                logmu[:, m * per:(m + 1) * per] += 0.8 * f[:, [m]]
            logmu[:, -rare:] = np.log(0.01)
            x[c] = rng.poisson(np.exp(logmu))
        obs = pd.DataFrame({"ct": types,
                            "batch": np.where(np.arange(n) % 4 < 2, "b1", "b2"),
                            "age": rng.normal(50, 10, n).round(1),
                            "pseudotime": rng.random(n).round(6)},
                           index=[f"{prefix}{i}" for i in range(n)])
        var = pd.DataFrame({"highly_variable": [True] * (n_genes - rare) + [False] * rare},
                           index=[f"g{i}" for i in range(n_genes)])
        return ad.AnnData(X=sp.csr_matrix(x), obs=obs, var=var)

    portable(cells(1200, 0, "c")).write_h5ad(d / "sc.h5ad")
    portable(cells(600, 1, "q")).write_h5ad(d / "sc_query.h5ad")


def _spatial(d: Path, rng) -> None:
    import anndata as ad

    n = 30
    row, col = np.divmod(np.arange(n * n), n)
    coords = np.column_stack([col, row]).astype(float)
    n_genes = 80
    logmu = np.log(rng.uniform(0.5, 2.0, n_genes))[None, :].repeat(n * n, 0)
    for k in range(4):  # spatial patterns: gradients and blobs
        f = (np.sin(row / (3 + k)) + np.cos(col / (4 + k))) * 0.6
        logmu[:, k * 12:(k + 1) * 12] += f[:, None]
    x = rng.poisson(np.exp(logmu)).astype(np.float32)
    obs = pd.DataFrame({"array_row": row, "array_col": col,
                        "region": np.array(["top", "middle", "bottom"])[row * 3 // n]},
                       index=[f"s{i}" for i in range(n * n)])
    a = ad.AnnData(X=x, obs=obs, var=pd.DataFrame(index=[f"G{i}" for i in range(n_genes)]))
    a.obsm["spatial"] = coords
    portable(a).write_h5ad(d / "sp.h5ad")


def write_expected(manifest: dict, out: Path) -> None:
    if EXPECTED.exists():
        shutil.rmtree(EXPECTED)
    for files in manifest.values():
        for rel in files:
            dst = EXPECTED / (rel + ".gz")
            dst.parent.mkdir(parents=True, exist_ok=True)
            with open(out / rel, "rb") as src, gzip.GzipFile(dst, "wb", mtime=0) as fh:
                fh.write(src.read())
    (EXPECTED / "manifest.json").write_text(json.dumps(manifest, indent=1) + "\n", encoding="utf-8")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--inputs", action="store_true", help="also (re)write the inputs")
    args = ap.parse_args()
    if args.inputs:
        make_inputs(INPUTS)
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp)
        manifest = run_steps(out)
        write_expected(manifest, out)
    n = sum(len(v) for v in manifest.values())
    print(f"{len(manifest)} steps, {n} files → {EXPECTED}")


if __name__ == "__main__":
    main()
