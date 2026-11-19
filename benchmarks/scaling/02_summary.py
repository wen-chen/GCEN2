"""Summarise 01_threads.sh: mean time and peak memory per command and thread
count, speed-up over one thread and parallel efficiency.

Usage: python 02_summary.py OUT_DIR   ->  results/threads.tsv, results/fig_threads.*
"""

import sys
from pathlib import Path

import pandas as pd

from gcen2.plot._common import figure, finish

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
t = pd.read_csv(Path(sys.argv[1]) / "time.txt", sep=" ",
                names=["command", "threads", "rep", "seconds", "peak_kb"])
s = t.groupby(["command", "threads"]).agg(seconds=("seconds", "mean"),
                                          runs=("rep", "size"),
                                          peak_kb=("peak_kb", "max")).reset_index()
one = s[s["threads"] == 1].set_index("command")["seconds"]
s["speedup"] = (s["command"].map(one) / s["seconds"]).round(2)
s["efficiency"] = (s["speedup"] / s["threads"]).round(2)
s["peak_GB"] = (s["peak_kb"] / 1024 ** 2).round(2)
s["seconds"] = s["seconds"].round(2)
order = ["build", "module_e", "diff", "celltype"]
s["command"] = pd.Categorical(s["command"], order, ordered=True)
s = s.sort_values(["command", "threads"])[["command", "threads", "runs", "seconds",
                                           "speedup", "efficiency", "peak_GB"]]
s.to_csv(OUT / "threads.tsv", sep="\t", index=False)
print(s.to_string(index=False))

labels = {"build": "build (FDR), 11,944 genes", "module_e": "module -e, 11,944 genes",
          "diff": "diff, 12.5M pairs", "celltype": "celltype-modules, 78,816 cells"}
fig = figure("single", 62)
ax = fig.add_subplot(111)
threads = sorted(s["threads"].unique())
ax.plot(threads, threads, color="grey", lw=0.6, ls="--", label="linear")
for c in order:
    d = s[s["command"] == c]
    ax.plot(d["threads"], d["speedup"], marker="o", ms=3, lw=1, label=labels[c])
ax.set_xscale("log", base=2)
ax.set_yscale("log", base=2)
ax.set_xticks(threads, [str(x) for x in threads])
ax.set_yticks(threads, [str(x) for x in threads])
ax.set_xlabel("threads (-t)")
ax.set_ylabel("speed-up over 1 thread")
ax.legend(frameon=False, fontsize=6)
finish(fig, OUT / "fig_threads", ["pdf", "png"])
