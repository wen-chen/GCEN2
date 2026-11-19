#!/usr/bin/env python
"""Benchmark GCEN v1 against GCEN2: runtime, peak memory and result agreement.

Both tools are run through their command-line interfaces in separate
processes, so wall time and peak RSS (``/usr/bin/time``) are measured the same
way for both and include file I/O.

Steps (per gene-count subset of the expression matrix):
  build    network_build  vs  gcen2 build     -> identical edge set? max |dr|
  module   module_identify vs gcen2 module   -> identical module sets?
  rwr      rwr            vs  gcen2 rwr       -> probabilities after scaling
                                                 v1 by n_seeds (v1 does not
                                                 normalise P0)
  annotate annotate -n    vs  gcen2 annotate  -> runtime only: same queries and
                                                 counts, but v2 tests one-sided
                                                 and v1 two-sided
  stats    network_stat   vs  gcen2 stats     -> the statistics side by side
  extract  network_extract vs gcen2 extract   -> identical edge set?
  merge    network_merge  vs  gcen2 merge     -> runtime only: v1 merges unordered
                                                 pairs and sums signs apart, v2
                                                 sums ordered pairs
  shuffle  network_shuffle vs gcen2 shuffle   -> runtime only (different
                                                 algorithms; v1 counts successful
                                                 swaps, v2 attempts): degrees kept?

stats / extract / merge / shuffle read the same input for both tools: the v1
network of the build step (merge adds a second network, built by gcen2 at
cor - 0.05, header removed, so that v1 can read it).

Example:
  python benchmarks/compare_v1_v2.py --v1 /path/to/GCEN-0.6.3 \\
      --sizes 5000 10000 20000 0 --threads 4 --out benchmarks/results
(size 0 = all genes.)
"""

from __future__ import annotations

import argparse
import json
import os
import resource
import shutil
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]


# --------------------------------------------------------------------------
# process runner
# --------------------------------------------------------------------------
def run(cmd: list[str], log: Path, mem_limit_gb: float | None) -> dict:
    """Run ``cmd`` under /usr/bin/time; return wall seconds, peak MB, status."""
    tfile = log.with_suffix(".time")

    def limit():
        if mem_limit_gb:
            b = int(mem_limit_gb * 1024**3)
            resource.setrlimit(resource.RLIMIT_AS, (b, b))

    with open(log, "w") as fh:
        proc = subprocess.run(
            ["/usr/bin/time", "-f", "%e %M", "-o", str(tfile), *cmd],
            stdout=fh, stderr=subprocess.STDOUT, preexec_fn=limit,
        )
    wall, rss = (float(x) for x in tfile.read_text().split()[-2:])
    out = {"seconds": wall, "peak_mb": round(rss / 1024, 1), "ok": proc.returncode == 0}
    if not out["ok"]:
        tail = log.read_text(errors="replace").strip().splitlines()[-1:]
        out["error"] = tail[0][:200] if tail else f"exit {proc.returncode}"
    return out


# --------------------------------------------------------------------------
# readers
# --------------------------------------------------------------------------
def _pair_key(a, b) -> pd.Index:
    a = np.asarray(a, dtype=object)
    b = np.asarray(b, dtype=object)
    return pd.Index(np.where(a < b, a + "|" + b, b + "|" + a))


def read_edges(path: Path) -> pd.Series:
    """Correlation indexed by undirected gene pair (v1 or v2 network file)."""
    d = pd.read_csv(path, sep="\t", comment="#", header=None, usecols=[0, 1, 2], low_memory=False)
    if d.iloc[0, 0] == "gene_a":  # v2 header line
        d = d.iloc[1:]
    return pd.Series(d[2].astype(float).to_numpy(), index=_pair_key(d[0], d[1]))


def read_modules_v1(path: Path) -> set[frozenset]:
    return {frozenset(l.split()) for l in open(path) if l.strip() and not l.startswith("#")}


def read_modules_v2(path: Path) -> set[frozenset]:
    d = pd.read_csv(path, sep="\t")
    if d.empty:
        return set()
    return {frozenset(g) for _, g in d.groupby("module")["gene"]}


def read_rwr(path: Path) -> pd.Series:
    d = pd.read_csv(path, sep="\t", comment="#", header=None)
    if d.iloc[0, 0] == "gene":
        d = d.iloc[1:]
    return pd.Series(d[1].astype(float).to_numpy(), index=d[0].to_numpy())


def read_pairs(path: Path) -> list[tuple[str, str]]:
    """Gene pairs (first two columns) of a v1 or v2 network file."""
    out = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            f = line.rstrip("\r\n").split("\t")
            if f[0] == "gene_a":  # v2 header line
                continue
            out.append((f[0], f[1]))
    return out


def degrees(pairs) -> dict[str, int]:
    """Undirected simple-graph degrees (duplicate / reversed pairs once)."""
    deg: dict[str, int] = {}
    for a, b in {(min(a, b), max(a, b)) for a, b in pairs}:
        deg[a] = deg.get(a, 0) + 1
        deg[b] = deg.get(b, 0) + 1
    return deg


def read_gene_list(path: Path) -> list[str]:
    return [l.split()[0] for l in open(path) if l.strip() and not l.startswith("#")]


# --------------------------------------------------------------------------
# comparisons
# --------------------------------------------------------------------------
def compare_networks(f1: Path, f2: Path) -> dict:
    e1, e2 = read_edges(f1), read_edges(f2)
    common = e1.index.intersection(e2.index)
    return {
        "edges_v1": len(e1), "edges_v2": len(e2),
        "only_v1": int(len(e1.index.difference(e2.index))),
        "only_v2": int(len(e2.index.difference(e1.index))),
        "max_abs_dr": float(np.abs(e1[common].to_numpy() - e2[common].to_numpy()).max())
        if len(common) else None,
    }


def compare_modules(f1: Path, f2: Path) -> dict:
    m1, m2 = read_modules_v1(f1), read_modules_v2(f2)
    return {"modules_v1": len(m1), "modules_v2": len(m2), "common": len(m1 & m2),
            "identical": m1 == m2}


def compare_rwr(f1: Path, f2: Path, seeds: list[str]) -> dict:
    p1 = read_rwr(f1) / len(seeds)  # v1 P0 is not normalised
    p1 = p1.drop([s for s in seeds if s in p1.index])  # v1 also reports seeds
    p2 = read_rwr(f2)
    p2 = p2[p2 > 0]
    common = p1.index.intersection(p2.index)
    d = np.abs(p1[common] - p2[common])
    top1 = list(p1.sort_values(ascending=False).index[:100])
    top2 = list(p2.sort_values(ascending=False).index[:100])
    return {
        "genes_v1": len(p1), "genes_v2": len(p2), "common": len(common),
        "max_abs_dp": float(d.max()) if len(d) else None,
        "spearman": float(p1[common].rank().corr(p2[common].rank())) if len(common) > 2 else None,
        "top100_overlap": len(set(top1) & set(top2)),
    }


def compare_extract(f1: Path, f2: Path) -> dict:
    p1 = {frozenset(p) for p in read_pairs(f1)}
    p2 = {frozenset(p) for p in read_pairs(f2)}
    return {"edges_v1": len(p1), "edges_v2": len(p2), "identical": p1 == p2}


def compare_shuffle(src: Path, f1: Path, f2: Path) -> dict:
    d0 = degrees(read_pairs(src))
    e0 = {frozenset(p) for p in read_pairs(src)}
    out = {}
    for tool, f in (("v1", f1), ("v2", f2)):
        pairs = read_pairs(f)
        e = {frozenset(p) for p in pairs}
        out[f"degrees_kept_{tool}"] = degrees(pairs) == d0
        out[f"edges_changed_{tool}"] = round(1 - len(e & e0) / max(len(e0), 1), 4)
    return out


_V1_STAT = {
    "nodes": r"Number of nodes: ([^;]+)", "edges": r"Number of edges: ([^;]+)",
    "components": r"Connected components: ([^;]+)", "giant": r"The giant component: ([^;]+)",
    "density": r"Network density: ([^;]+)", "average_degree": r"Average node degree: ([^;]+)",
    "clustering": r"Average clustering coefficient: ([^;]+)", "transitivity": r"Transitivity: ([^;.]+(?:\.\d+)?)",
}
_V2_STAT = {
    "nodes": r"Nodes:\s+(\S+)", "edges": r"Edges:\s+(\S+)",
    "components": r"Connected components:\s+(\S+)", "giant": r"Giant component size:\s+(\S+)",
    "density": r"Density:\s+(\S+)", "average_degree": r"Average degree:\s+(\S+)",
    "clustering": r"Clustering coefficient:\s+(\S+)", "transitivity": r"Transitivity:\s+(\S+)",
}


def compare_stats(log1: Path, log2: Path) -> dict:
    """The statistics printed by v1 network_stat and gcen2 stats, side by side."""
    import re

    t1, t2 = log1.read_text(), log2.read_text()
    out = {}
    for k in _V1_STAT:
        m1, m2 = re.search(_V1_STAT[k], t1), re.search(_V2_STAT[k], t2)
        out[k] = [float(m1.group(1)) if m1 else None, float(m2.group(1)) if m2 else None]
    return out


# --------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------
def make_subset(expr: Path, n: int, keep: set[str], out: Path) -> int:
    """First ``n`` genes of a header-less v1 matrix, always keeping ``keep``."""
    with open(expr) as fh:
        lines = [l for l in fh if l.strip() and not l.startswith("#")]
    if n <= 0 or n >= len(lines):
        chosen = lines
    else:
        head = lines[:n]
        have = {l.split("\t", 1)[0] for l in head}
        extra = [l for l in lines[n:] if l.split("\t", 1)[0] in keep - have]
        chosen = head + extra
    out.write_text("".join(chosen))
    return len(chosen)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--v1", required=True, type=Path, help="GCEN v1 directory containing bin/")
    ap.add_argument("--gcen2", default=shutil.which("gcen2") or "gcen2", help="gcen2 executable")
    ap.add_argument("--data", type=Path, default=ROOT / "data" / "sample")
    ap.add_argument("--sizes", type=int, nargs="+", default=[5000, 10000, 20000, 0])
    ap.add_argument("--threads", type=int, default=4)
    ap.add_argument("--steps", nargs="+",
                    default=["build", "module", "rwr", "annotate", "stats", "extract", "merge",
                             "shuffle"])
    ap.add_argument("--swap-factor", type=int, default=10,
                    help="shuffle: v1 -s / gcen2 --swap-factor (v1 counts successful swaps)")
    ap.add_argument("--cor", type=float, default=0.9)
    ap.add_argument("--pval", type=float, default=0.001)
    ap.add_argument("--similarity", type=float, default=0.5)
    ap.add_argument("--mem-limit-gb", type=float, default=6.0,
                    help="address-space limit per process (protects the host)")
    ap.add_argument("--out", type=Path, default=ROOT / "benchmarks" / "results")
    args = ap.parse_args()

    v1bin = args.v1 / "bin"
    expr = args.data / "gene_expr.tsv"
    seeds_file = args.data / "rwr_seed_genes.list"
    seeds = read_gene_list(seeds_file)
    args.out.mkdir(parents=True, exist_ok=True)
    results = []

    # Interpreter + import start-up of the gcen2 CLI, included in every v2 row.
    # It is dominated by disk reads; on WSL with the venv under /mnt/* it can
    # be several seconds, on a native Linux file system ~1 s.
    r = run([args.gcen2, "--help"], args.out / "v2_startup.log", None)
    results.append({"genes": 0, "step": "startup", "tool": "v2", **r})
    print(f"v2 CLI start-up: {r['seconds']:.1f} s", flush=True)

    for size in args.sizes:
        tag = f"n{size or 'all'}"
        wd = args.out / tag
        wd.mkdir(exist_ok=True)
        sub = wd / "expr.tsv"
        n_genes = make_subset(expr, size, set(seeds), sub)
        print(f"== {tag}: {n_genes} genes", flush=True)
        net1, net2 = wd / "v1.network", wd / "v2.network"
        mod1, mod2 = wd / "v1.module", wd / "v2.module"
        lim = args.mem_limit_gb

        def record(step, tool, res, n_genes=n_genes, **extra):
            row = {"genes": n_genes, "step": step, "tool": tool, **res, **extra}
            results.append(row)
            print(f"   {step:8s} {tool}: {res['seconds']:8.1f} s {res['peak_mb']:9.1f} MB"
                  f"{'' if res['ok'] else '  FAILED: ' + res.get('error', '')}", flush=True)

        if "build" in args.steps:
            common = ["-m", "spearman", "-c", str(args.cor), "-p", str(args.pval), "-t", str(args.threads)]
            r1 = run([str(v1bin / "network_build"), "-i", str(sub), "-o", str(net1), *common],
                     wd / "v1_build.log", lim)
            record("build", "v1", r1)
            r2 = run([args.gcen2, "build", "-i", str(sub), "-o", str(net2), *common],
                     wd / "v2_build.log", lim)
            record("build", "v2", r2)
            if r1["ok"] and r2["ok"]:
                results[-1]["agreement"] = compare_networks(net1, net2)
                print(f"     agreement: {results[-1]['agreement']}", flush=True)

        if "module" in args.steps and net1.exists() and net2.exists():
            r1 = run([str(v1bin / "module_identify"), "-i", str(net1), "-o", str(mod1),
                      "-s", str(args.similarity), "-t", str(args.threads)], wd / "v1_module.log", lim)
            record("module", "v1", r1)
            r2 = run([args.gcen2, "module", "-i", str(net2), "-o", str(mod2),
                      "-s", str(args.similarity), "-t", str(args.threads)], wd / "v2_module.log", lim)
            record("module", "v2", r2)
            if r1["ok"] and r2["ok"]:
                results[-1]["agreement"] = compare_modules(mod1, mod2)
                print(f"     agreement: {results[-1]['agreement']}", flush=True)

        if "rwr" in args.steps and net1.exists() and net2.exists():
            o1, o2 = wd / "v1.rwr", wd / "v2.rwr"
            r1 = run([str(v1bin / "rwr"), "-n", str(net1), "-g", str(seeds_file), "-o", str(o1),
                      "-p", "0"], wd / "v1_rwr.log", lim)
            record("rwr", "v1", r1)
            r2 = run([args.gcen2, "rwr", "-n", str(net2), "-g", str(seeds_file), "-o", str(o2),
                      "--prob-cutoff", "0"], wd / "v2_rwr.log", lim)
            record("rwr", "v2", r2)
            if r1["ok"] and r2["ok"]:
                results[-1]["agreement"] = compare_rwr(o1, o2, seeds)
                print(f"     agreement: {results[-1]['agreement']}", flush=True)

        if "annotate" in args.steps and net1.exists() and net2.exists():
            go, assoc = args.data / "go-basic.obo", args.data / "gene_go.assoc"
            for tool, cmd, d in [
                ("v1", [str(v1bin / "annotate"), "-g", str(go), "-a", str(assoc), "-n", str(net1),
                        "-t", str(args.threads)], wd / "v1_go"),
                ("v2", [args.gcen2, "annotate", "-g", str(go), "-a", str(assoc), "-n", str(net2),
                        "-t", str(args.threads)],
                 wd / "v2_go"),
            ]:
                shutil.rmtree(d, ignore_errors=True)
                r = run([*cmd, "-o", str(d)], wd / f"{tool}_annotate.log", lim)
                files = len(os.listdir(d)) if d.exists() else 0
                record("annotate", tool, r, output_files=files)
                shutil.rmtree(d, ignore_errors=True)  # tens of thousands of small files

        # stats / extract / merge / shuffle: both tools read the v1 network
        util = args.v1 / "util"
        if "stats" in args.steps and net1.exists():
            l1, l2 = wd / "v1_stats.log", wd / "v2_stats.log"
            record("stats", "v1", run([str(util / "network_stat"), "-i", str(net1)], l1, lim))
            r2 = run([args.gcen2, "stats", "-i", str(net1), "-t", str(args.threads)], l2, lim)
            record("stats", "v2", r2)
            if results[-2]["ok"] and r2["ok"]:
                results[-1]["agreement"] = compare_stats(l1, l2)
                print(f"     agreement: {results[-1]['agreement']}", flush=True)

        if "extract" in args.steps and net1.exists():
            genes = args.data / "gene_list.txt"
            o1, o2 = wd / "v1_extract.network", wd / "v2_extract.network"
            r1 = run([str(util / "network_extract"), "-i", str(net1), "-g", str(genes),
                      "-o", str(o1)], wd / "v1_extract.log", lim)
            record("extract", "v1", r1)
            r2 = run([args.gcen2, "extract", "-i", str(net1), "-g", str(genes), "-o", str(o2)],
                     wd / "v2_extract.log", lim)
            record("extract", "v2", r2)
            if r1["ok"] and r2["ok"]:
                results[-1]["agreement"] = compare_extract(o1, o2)
                print(f"     agreement: {results[-1]['agreement']}", flush=True)

        if "merge" in args.steps and net1.exists():
            # second network: gcen2 build at cor - 0.05, header removed for v1
            net_b, tmp = wd / "second.network", wd / "second_v2.network"
            prep = run([args.gcen2, "build", "-i", str(sub), "-o", str(tmp), "-m", "spearman",
                        "-c", str(round(args.cor - 0.05, 4)), "-p", str(args.pval),
                        "-t", str(args.threads)], wd / "second_build.log", lim)
            if prep["ok"]:
                with open(tmp) as fi, open(net_b, "w") as fo:
                    next(fi)
                    fo.writelines(fi)
                tmp.unlink()
                o1, o2 = wd / "v1_merge.network", wd / "v2_merge.network"
                r1 = run([str(util / "network_merge"), "-i", f"{net1},{net_b}", "-o", str(o1),
                          "-c", "0.5"], wd / "v1_merge.log", lim)
                record("merge", "v1", r1)
                r2 = run([args.gcen2, "merge", "-i", str(net1), "-i", str(net_b), "-o", str(o2),
                          "-c", "0.5"], wd / "v2_merge.log", lim)
                record("merge", "v2", r2)
                if r1["ok"] and r2["ok"]:
                    results[-1]["agreement"] = {"edges_v1": len(read_pairs(o1)),
                                                "edges_v2": len(read_pairs(o2))}
                    print(f"     edges: {results[-1]['agreement']}", flush=True)

        if "shuffle" in args.steps and net1.exists():
            o1, o2 = wd / "v1_shuffle.network", wd / "v2_shuffle.network"
            sf = str(args.swap_factor)
            r1 = run([str(util / "network_shuffle"), "-i", str(net1), "-o", str(o1), "-s", sf],
                     wd / "v1_shuffle.log", lim)
            record("shuffle", "v1", r1)
            r2 = run([args.gcen2, "shuffle", "-i", str(net1), "-o", str(o2), "--swap-factor", sf],
                     wd / "v2_shuffle.log", lim)
            record("shuffle", "v2", r2)
            if r1["ok"] and r2["ok"]:
                results[-1]["agreement"] = compare_shuffle(net1, o1, o2)
                print(f"     agreement: {results[-1]['agreement']}", flush=True)

        (args.out / "results.json").write_text(json.dumps(results, indent=1))

    # summary table
    df = pd.DataFrame(results)
    tab = df.pivot_table(index=["genes", "step"], columns="tool",
                         values=["seconds", "peak_mb"], aggfunc="first")
    tab.columns = [f"{v}_{t}" for v, t in tab.columns]
    ok = df.pivot_table(index=["genes", "step"], columns="tool", values="ok", aggfunc="first")
    for t in ok.columns:
        tab[f"ok_{t}"] = ok[t]
    if {"seconds_v1", "seconds_v2"} <= set(tab.columns):
        tab["speedup"] = (tab["seconds_v1"] / tab["seconds_v2"]).round(2)
    tab.to_csv(args.out / "summary.tsv", sep="\t")
    print(tab.to_string())
    return 0


if __name__ == "__main__":
    sys.exit(main())
