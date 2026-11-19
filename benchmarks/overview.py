#!/usr/bin/env python
"""One table of every GCEN2 benchmark, built from the committed results.

Reads
  results.json, summary.tsv                    GCEN v1 vs GCEN2 (compare_v1_v2.py)
  wgcna_gtex/results/{tom_modules,time}.tsv    WGCNA, GTEx whole blood
  hdwgcna_pbmc/results/tom_modules.tsv,        hdWGCNA single-cell, PBMC 10k
    gcen2_time_gemm.txt, hdwgcna_time.txt
  hdwgcna_hlca/results/{compare,time}.tsv      hdWGCNA single-cell, lung atlas
  hdwgcna_consensus/results/compare.tsv        hdWGCNA consensus, harmonized MEs, DMEs
  hdwgcna_preservation/results/compare.tsv     hdWGCNA projection and preservation
  hdwgcna_spatial/results/compare.tsv,         hdWGCNA spatial, Visium,
    {slideseqv2,visium_hd}_compare.tsv         Slide-seqV2 and Visium HD
  dgca_gtex/results/compare.tsv                DGCA, differential co-expression
  multilayer_string/results/compare.tsv        RandomWalkRestartMH; fusion in R
  trajectory_paul15/results/compare.tsv        trajectory networks in R
and writes overview.tsv plus the Markdown table between the
``<!-- overview:begin -->`` / ``<!-- overview:end -->`` markers of README.md.

  python benchmarks/overview.py           # rewrite overview.tsv and README.md
  python benchmarks/overview.py --check   # exit 1 if either is out of date
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
BEGIN, END = "<!-- overview:begin -->", "<!-- overview:end -->"
COLUMNS = ["group", "data", "task", "gcen2_s", "gcen2_mb", "other_s", "other_mb",
           "faster", "less_memory", "agreement"]


# --------------------------------------------------------------------------
# readers and formatting
# --------------------------------------------------------------------------
def read_tsv(path: Path) -> list[dict]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def read_time(path: Path, name_cols: int) -> list[tuple[str, float, float]]:
    """(name, seconds, MB) of /usr/bin/time lines "<name...> <s> [s] <KB> [KB]"."""
    out = []
    for line in path.read_text(encoding="utf-8").splitlines():
        f = [x for x in line.split() if x not in ("s", "KB")]
        if len(f) >= name_cols + 2:
            out.append((" ".join(f[:name_cols]), float(f[-2]), float(f[-1]) / 1024))
    return out


def fnum(x: float) -> str:
    """0.77, 7.5, 23.0, 213: 2 significant digits below 1, 1 decimal below 100."""
    if x >= 99.5:
        return f"{x:,.0f}"
    return f"{x:.2g}" if x < 1 else f"{x:.1f}"


def fs(lo: float, hi: float | None = None) -> str:
    if hi is None or fnum(lo) == fnum(hi):
        return f"{fnum(lo)} s"
    return f"{fnum(lo)}–{fnum(hi)} s"


def fmb(lo: float | None, hi: float | None = None) -> str:
    """MB below 1 GB, else GB (1.43 GB, 10.6 GB); a range keeps both units."""
    if lo is None:
        return "–"

    def one(x):
        if x < 1000:
            return f"{x:.0f}", "MB"
        g = x / 1024
        return (f"{g:.2f}" if g < 9.995 else f"{g:.1f}"), "GB"
    (a, ua), (b, ub) = one(lo), one(lo if hi is None else hi)
    if (a, ua) == (b, ub):
        return f"{a} {ua}"
    return f"{a}–{b} {ua}" if ua == ub else f"{a} {ua}–{b} {ub}"


def ratio(lo: float, hi: float | None = None) -> str:
    def one(x):
        return f"{x:.1f}×" if x < 9.95 else f"{x:.0f}×"
    if hi is None or one(lo) == one(hi):
        return one(lo) if lo >= 0.95 else f"{one(1 / lo)} more"
    if lo < 0.95:
        return f"{one(1 / lo)} more – {one(hi)} less"
    return f"{one(lo)[:-1]}–{one(hi)}"


def row(group, data, task, v2, other, agreement, v2_mb=None, other_mb=None):
    """v2 / other: seconds (scalar) or a list over runs; same for the MB."""
    v2 = v2 if isinstance(v2, list) else [v2]
    other = other if isinstance(other, list) else [other]
    sp = [o / g for g, o in zip(v2, other)]
    out = {"group": group, "data": data, "task": task,
           "gcen2_s": fs(min(v2), max(v2)), "other_s": fs(min(other), max(other)),
           "faster": ratio(min(sp), max(sp))}
    if v2_mb is not None:
        v2_mb = v2_mb if isinstance(v2_mb, list) else [v2_mb]
        out["gcen2_mb"] = fmb(min(v2_mb), max(v2_mb))
    else:
        out["gcen2_mb"] = "–"
    if other_mb is not None and v2_mb is not None:
        other_mb = other_mb if isinstance(other_mb, list) else [other_mb]
        out["other_mb"] = fmb(min(other_mb), max(other_mb))
        m = [o / g for g, o in zip(v2_mb, other_mb)]
        out["less_memory"] = ratio(min(m), max(m))
    else:
        out["other_mb"] = "–"
        out["less_memory"] = "–"
    out["agreement"] = agreement
    return out


# --------------------------------------------------------------------------
# the benchmarks
# --------------------------------------------------------------------------
V1_AGREEMENT = {
    "build": lambda a: (f"identical edges ({a['edges_v2']:,}), max |Δr| = {a['max_abs_dr']:.2g}"
                        if a["only_v1"] == a["only_v2"] == 0 else
                        f"{a['only_v1']} / {a['only_v2']} edges only in v1 / v2"),
    "module": lambda a: (f"identical modules ({a['modules_v2']:,})" if a["identical"]
                         else f"{a['common']} of {a['modules_v1']} modules shared"),
    "rwr": lambda a: (f"max |Δp| = {a['max_abs_dp']:.1e}, top 100 "
                      + ("identical" if a["top100_overlap"] == 100
                         else f"{a['top100_overlap']} shared")),
    "stats": lambda a: ("same statistics" if all(abs(x - y) <= 5e-4 * max(abs(x), 1e-12)
                                               for x, y in a.values())
                        else "statistics differ"),
    "extract": lambda a: ("identical edges" if a["identical"] else "edges differ"),
    "merge": lambda a: f"same edge count ({a['edges_v2']:,}); sums differ by design",
    "shuffle": lambda a: ("both keep every degree; different algorithms"
                          if a["degrees_kept_v1"] and a["degrees_kept_v2"]
                          else "degrees changed"),
}
V1_STEPS = ["build", "module", "rwr", "annotate", "stats", "extract", "merge", "shuffle"]


def gcen_v1() -> tuple[list[dict], list[str]]:
    res = json.loads((HERE / "results.json").read_text(encoding="utf-8"))
    by = {(r["genes"], r["step"], r["tool"]): r for r in res}
    sizes = sorted({r["genes"] for r in res if r["genes"]})
    big = sizes[-1]
    edges = by[(big, "build", "v2")]["agreement"]["edges_v2"]
    rows = []
    for step in V1_STEPS:
        a, b = by.get((big, step, "v1")), by.get((big, step, "v2"))
        if not (a and b and a["ok"] and b["ok"]):
            continue
        agr = (V1_AGREEMENT[step](b["agreement"]) if "agreement" in b
               else "runtime only (one- vs two-sided test)")
        rows.append(row("GCEN v1 0.6.3", f"{big:,} genes × 15 samples, {edges / 1e6:.2g}M edges",
                        step, b["seconds"], a["seconds"], agr, b["peak_mb"], a["peak_mb"]))
    # where GCEN2 is not ahead, per step over all subsets
    notes = []
    for step in V1_STEPS:
        pairs = [(g, by[(g, step, "v1")], by[(g, step, "v2")]) for g in sizes
                 if (g, step, "v1") in by and (g, step, "v2") in by
                 and by[(g, step, "v1")]["ok"] and by[(g, step, "v2")]["ok"]]
        slow = [(g, a, b) for g, a, b in pairs if b["seconds"] > a["seconds"]]
        mem = [(g, a, b) for g, a, b in pairs if b["peak_mb"] > 1.05 * a["peak_mb"]]
        what = []
        if slow:
            g, a, b = slow[0]
            what.append(f"slower on {len(slow)} of {len(pairs)} subsets "
                        f"({b['seconds']:.2f} s vs {a['seconds']:.2f} s on {g:,} genes)")
        if mem:
            what.append(f"more memory on {len(mem)} of {len(pairs)} "
                        f"(GCEN2 {fmb(min(b['peak_mb'] for _, _, b in mem), max(b['peak_mb'] for _, _, b in mem))}, "
                        f"v1 {fmb(min(a['peak_mb'] for _, a, _ in mem), max(a['peak_mb'] for _, a, _ in mem))})")
        if what:
            notes.append(f"`{step}`: " + "; ".join(what))
    return rows, notes


def wgcna() -> list[dict]:
    d = HERE / "wgcna_gtex/results"
    tm = read_tsv(d / "tom_modules.tsv")
    # the largest run: pc5_<genes> with the most genes
    big = max((r for r in tm if r["run"].split("_")[-1].isdigit()),
              key=lambda r: int(r["run"].split("_")[-1]))
    n = big["run"].split("_")[-1]
    ari = [float(r["ARI"]) for r in tm]
    agr = (f"same modules in {sum(a == 1.0 for a in ari)} of {len(ari)} runs (ARI 1.0)")
    rows = [row("WGCNA 1.74", f"GTEx whole blood, {int(n):,} genes × 803 samples",
                "WGCNA modules (`module -e`)", float(big["gcen2_s"]), float(big["wgcna_s"]),
                agr, float(big["gcen2_GB"]) * 1024, float(big["wgcna_GB"]) * 1024)]
    t = {r["run"]: r for r in read_tsv(d / "time.tsv")}[big["run"]]
    rows.append(row("WGCNA 1.74", f"GTEx whole blood, {int(n):,} genes × 803 samples",
                    "correlation network + modules, FDR (`build` + `module`)",
                    float(t["gcen2_s"]), float(t["wgcna_s"]),
                    "different method (r and FDR, not TOM)",
                    float(t["gcen2_GB"]) * 1024, float(t["wgcna_GB"]) * 1024))
    return rows


def hdwgcna_pbmc() -> list[dict]:
    d = HERE / "hdwgcna_pbmc/results"
    tm = read_tsv(d / "tom_modules.tsv")
    # the six cell types, without the split halves (_A, _B) and metacell seeds
    main = [r for r in tm if not re.search(r"_(A|B|seed\d+)$", r["run"])]
    ari = [float(r["ARI"]) for r in tm]
    rows = [row("hdWGCNA 0.4.12", f"PBMC 10k, {len(main)} cell types, 44–1,000 metacells",
                "hdWGCNA modules (`module -e --preset hdwgcna`)",
                [float(r["gcen2_s"]) for r in main],
                [float(r["hdwgcna_network_s"]) for r in main],
                f"same modules in {sum(a == 1.0 for a in ari)} of {len(ari)} runs (ARI 1.0); "
                "hdWGCNA time: soft power + network only",
                [float(r["gcen2_GB"]) * 1024 for r in main])]
    g = {}
    for name, s, mb in read_time(d / "gcen2_time_gemm.txt", 2):
        g.setdefault(name.split()[0], []).append((s, mb))
    h = {name.split()[1]: (s, mb) for name, s, mb in read_time(d / "hdwgcna_time.txt", 2)}
    cts = [c for c in g if c in h]
    mean = {c: [sum(x[i] for x in g[c]) / len(g[c]) for i in (0, 1)] for c in cts}
    rows.append(row("hdWGCNA 0.4.12", f"PBMC 10k, {len(cts)} cell types, 355–3,391 cells",
                    "correlation network + modules per cell type (`03_gcen2.py`)",
                    [mean[c][0] for c in cts], [h[c][0] for c in cts],
                    "different method (r, no TOM); whole script vs hdWGCNA's",
                    [mean[c][1] for c in cts], [h[c][1] for c in cts]))
    return rows


def hdwgcna_hlca() -> list[dict]:
    """The whole atlas, the 10 largest cell types (both tools)."""
    d = HERE / "hdwgcna_hlca/results"
    cmp = [r for r in read_tsv(d / "compare.tsv") if r["size"] == "all"]
    t = {(r["run"], r["size"]): r for r in read_tsv(d / "time.tsv")}
    g, h = t[("gcen2_top10", "all")], t[("hdwgcna_top10", "all")]
    cells = sum(int(r["cells"]) for r in cmp)
    same = sum(r["same_genes"] == "True" for r in cmp)
    net = sum(float(r["network_step_ARI"]) == 1.0 for r in cmp)
    return [row("hdWGCNA 0.4.12", f"Human Lung Cell Atlas core (584,944 cells), the "
                f"{len(cmp)} largest cell types ({cells:,} cells)",
                "single-cell workflow (`celltype-modules`)",
                float(g["seconds"]), float(h["seconds"]),
                f"same genes in {same} of {len(cmp)} types; same modules on hdWGCNA's "
                f"metacells in {net} of {len(cmp)} (ARI 1.0)",
                float(g["peak_GB"]) * 1024, float(h["peak_GB"]) * 1024)]


def hdwgcna_consensus() -> list[dict]:
    """Consensus network over datasets with harmonized eigengenes (HLCA NK cells)."""
    c = {r["metric"]: r["value"] for r in read_tsv(HERE / "hdwgcna_consensus/results/compare.tsv")}
    ari = [float(x) for x in c["ARI_gcen2_vs_hdwgcna"].split(",")]
    return [row("hdWGCNA 0.4.12", "HLCA natural killer cells, 16,006 cells in 4 datasets, "
                f"{int(c['genes_gcen2']):,} genes",
                "consensus network, harmonized MEs, DMEs, traits (`celltype-modules "
                "--consensus-by --harmonize-by --dme-by --trait`)",
                float(c["gcen2_s"]), float(c["hdwgcna_s"]),
                f"same modules on hdWGCNA's metacells (ARI {float(c['network_step_ARI']):.1f}); "
                f"harmonized MEs r >= {float(c['hME_r_gcen2_vs_hdwgcna_min']):.3f}; "
                "same DMEs and trait correlations; "
                f"own metacells: ARI {min(ari):.2f}–{max(ari):.2f}, as seed to seed",
                float(c["gcen2_GB"]) * 1024, float(c["hdwgcna_GB"]) * 1024)]


def hdwgcna_preservation() -> list[dict]:
    """Module projection and preservation (HLCA NK cells, one dataset onto another)."""
    c = {r["metric"]: r["value"] for r in
         read_tsv(HERE / "hdwgcna_preservation/results/compare.tsv")}
    return [row("hdWGCNA 0.4.12", "HLCA natural killer cells, Meyer_2019 (6,732 cells) onto "
                f"Krasnow_2020 (4,303 cells), {int(c['genes']):,} genes",
                "reference modules, projection and preservation, 500 permutations "
                "(`celltype-modules`, `project-modules`, `module-preservation`)",
                float(c["gcen2_s"]), float(c["hdwgcna_s"]),
                f"projected MEs within {float(c['projected_ME_max_abs_diff']):.0e}; "
                "on hdWGCNA's metacells the same preservation statistics (Z within "
                f"{float(c['Z_max_rel_diff']):.0e} relative)",
                float(c["gcen2_GB"]) * 1024, float(c["hdwgcna_GB"]) * 1024)]


def hdwgcna_spatial() -> list[dict]:
    rows = read_tsv(HERE / "hdwgcna_spatial/results/compare.tsv")
    label = {"dlpfc": "Visium DLPFC 151673", "mouse_anterior": "Visium mouse brain anterior"}
    out = []
    for r in rows:
        spots = {"dlpfc": "3,611", "mouse_anterior": "2,695"}.get(r["sample"], "?")
        a = float(r["ARI_final"])
        agr = (f"same modules (ARI {a:.1f})" if a == 1.0 else f"modules ARI {a:.2f}")
        agr += f", spot MEs within {float(r['spot_ME_max_abs_diff']):.0e} (TSV precision)"
        out.append(row("hdWGCNA 0.4.12", f"{label.get(r['sample'], r['sample'])}, {spots} spots, "
                       f"{int(r['genes']):,} genes", "spatial workflow (`spatial-modules`)",
                       float(r["gcen2_s"]), float(r["hdwgcna_s"]), agr,
                       float(r["gcen2_GB"]) * 1024, float(r["hdwgcna_GB"]) * 1024))
    knn = [("slideseqv2", "Slide-seqV2 hippocampus (no grid), 41,786 beads", 100),
           ("visium_hd", "Visium HD mouse brain, 16 µm bins (square grid), 98,917 bins", 64)]
    for name, data, k in knn:
        ss = {r["metric"]: r[name]
              for r in read_tsv(HERE / f"hdwgcna_spatial/results/{name}_compare.tsv")}
        ari = [float(x) for x in ss["ARI_gcen2_vs_hdwgcna"].split(",")]
        out.append(row("hdWGCNA 0.4.12", f"{data}, {int(ss['genes']):,} genes",
                       f"spatial workflow, knn metaspots (`spatial-modules -k {k}`)",
                       float(ss["gcen2_s"]), float(ss["hdwgcna_s"]),
                       f"same modules on hdWGCNA's metacells "
                       f"(ARI {float(ss['network_step_ARI']):.1f}); "
                       f"own metacells: ARI {min(ari):.2f}–{max(ari):.2f}, as seed to seed",
                       float(ss["gcen2_GB"]) * 1024, float(ss["hdwgcna_GB"]) * 1024))
    return out


def metrics(path: str) -> dict:
    """{metric: value} of a two-column compare.tsv."""
    return {r["metric"]: r["value"] for r in read_tsv(HERE / path)}


def dgca() -> list[dict]:
    rows = read_tsv(HERE / "dgca_gtex/results/compare.tsv")
    out = []
    for m in ("pearson", "spearman"):
        c = {r["metric"]: r[f"{m}_5000"] for r in rows}
        out.append(row("DGCA 1.0.3", f"GTEx whole blood, {int(c['genes']):,} genes, "
                       f"{c['samples_A']} male vs {c['samples_B']} female",
                       f"differential co-expression of every pair, {m.capitalize()} "
                       "(`build` × 2 + `diff`)",
                       float(c["gcen2_s"]), float(c["dgca_s"]),
                       f"same test on {int(c['pairs_common']) / 1e6:.1f}M pairs: r within "
                       f"{float(c['max_abs_dr_A']):.0e} (6 decimals), z within "
                       f"{float(c['max_abs_dz_same_r']):.0e} on the same r; "
                       f"{c['significant_differ']} pairs differ at p = 0.05, "
                       f"{c['fdr_significant_differ']} at FDR = 0.05",
                       float(c["gcen2_GB"]) * 1024, float(c["dgca_GB"]) * 1024))
    return out


def multilayer() -> list[dict]:
    c = metrics("multilayer_string/results/compare.tsv")
    data = (f"GTEx co-expression + STRING + TRRUST, {int(c['genes']):,} genes, "
            "3 layers")
    return [row("RandomWalkRestartMH 1.13.2", data, "multiplex RWR (`rwr-multiplex`)",
                float(c["rwr_gcen2_s"]), float(c["rwr_rwrmh_s"]),
                f"scores within {float(c['rwr_max_abs_diff']):.0e}, "
                f"top 100 {'identical' if c['rwr_top100_shared'] == '100' else 'differ'}",
                float(c["rwr_gcen2_GB"]) * 1024, float(c["rwr_rwrmh_GB"]) * 1024),
            row("R (data.table)", data, "network fusion (`fuse`)",
                float(c["fuse_gcen2_s"]), float(c["fuse_r_s"]),
                f"same {int(c['fuse_edges_shared']):,} edges, scores within "
                f"{float(c['fuse_max_abs_score_diff']):.0e}",
                float(c["fuse_gcen2_GB"]) * 1024, float(c["fuse_r_GB"]) * 1024)]


def trajectory() -> list[dict]:
    c = metrics("trajectory_paul15/results/compare.tsv")
    same = c["same_edges_every_window"] == "True" and c["same_pattern_first_last"] == "True"
    return [row("R (WGCNA)", f"paul15 myeloid differentiation, 2,730 cells, 2,000 genes, "
                f"{c['windows']} windows", "networks along pseudotime (`trajectory`)",
                float(c["gcen2_s"]), float(c["r_reference_s"]),
                (f"same edges in every window, same patterns of "
                 f"{int(c['dynamics_edges_shared']):,} edges" if same else "edges differ"),
                float(c["gcen2_GB"]) * 1024, float(c["r_reference_GB"]) * 1024)]


# --------------------------------------------------------------------------
def markdown(rows: list[dict], notes: list[str]) -> str:
    head = ["vs", "data", "task", "GCEN2", "other tool", "faster", "less memory", "results"]
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for r in rows:
        cells = [r["group"], r["data"], r["task"],
                 r["gcen2_s"] + (f", {r['gcen2_mb']}" if r["gcen2_mb"] != "–" else ""),
                 r["other_s"] + (f", {r['other_mb']}" if r["other_mb"] != "–" else ""),
                 r["faster"], r["less_memory"], r["agreement"].replace("|", "\\|")]
        lines.append("| " + " | ".join(cells) + " |")
    if notes:
        lines += ["", "Where GCEN2 is behind GCEN v1 on the gene subsets "
                  "(the fixed cost is the Python + numpy start-up, about 0.1 s and 60 MB):", ""]
        lines += [f"- {n}" for n in notes]
    return "\n".join(lines)


def build() -> tuple[str, str]:
    v1, notes = gcen_v1()
    rows = (v1 + wgcna() + hdwgcna_pbmc() + hdwgcna_hlca() + hdwgcna_consensus()
            + hdwgcna_preservation() + hdwgcna_spatial() + dgca() + multilayer()
            + trajectory())
    buf = io.StringIO()
    w = csv.DictWriter(buf, COLUMNS, delimiter="\t", lineterminator="\n")
    w.writeheader()
    w.writerows(rows)
    return buf.getvalue(), markdown(rows, notes)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--check", action="store_true",
                    help="only check that overview.tsv and README.md are up to date")
    args = ap.parse_args()
    tsv, md = build()
    readme = (HERE / "README.md").read_text(encoding="utf-8")
    i, j = readme.index(BEGIN) + len(BEGIN), readme.index(END)
    new_readme = readme[:i] + "\n" + md + "\n" + readme[j:]
    tsv_path = HERE / "overview.tsv"
    if args.check:
        stale = [p.name for p, new in ((tsv_path, tsv), (HERE / "README.md", new_readme))
                 if not p.exists() or p.read_text(encoding="utf-8") != new]
        if stale:
            print("out of date (run benchmarks/overview.py): " + ", ".join(stale))
            return 1
        print("benchmark overview up to date")
        return 0
    tsv_path.write_text(tsv, encoding="utf-8")
    (HERE / "README.md").write_text(new_readme, encoding="utf-8")
    print(md)
    return 0


if __name__ == "__main__":
    sys.exit(main())
