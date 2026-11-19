"""Make the gallery images of the GCEN2 website from the example data.

Runs, in a temporary directory, the setup commands of each group and then the
command of each figure, exactly as listed in ``website/content.py`` (the page
shows the same commands), and writes ``website/static/img/gallery/<id>.png``
(downscaled to at most 1,200 px wide, 8-bit palette where lossless enough).

usage: python website/gallery/make_gallery.py [--keep DIR]
Needs gcen2 with the single-cell and visualization extras on PATH, and Pillow.
"""
from __future__ import annotations

import argparse
import os
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
SITE = HERE.parent
REPO = SITE.parent
sys.path.insert(0, str(SITE))
import content  # noqa: E402

OUT = SITE / "static" / "img" / "gallery"
MAX_W = 1200


def run(cmd: str, cwd: Path, env: dict) -> None:
    args = [os.path.expandvars(a) if "$" in a else a for a in shlex.split(cmd)]
    args = [a.replace("$EX", env["EX"]).replace("$SC", env["SC"]).replace("$SP", env["SP"])
            for a in args]
    t = time.time()
    r = subprocess.run(args, cwd=cwd, env=env, capture_output=True, text=True,
                       encoding="utf-8")
    if r.returncode:
        sys.stderr.write(r.stdout[-2000:] + r.stderr[-4000:])
        raise SystemExit(f"failed: {cmd}")
    print(f"  {time.time() - t:6.1f} s  {cmd[:100]}")


def shrink(src: Path, dst: Path) -> None:
    from PIL import Image

    im = Image.open(src).convert("RGB")
    if im.width > MAX_W:
        im = im.resize((MAX_W, round(im.height * MAX_W / im.width)), Image.LANCZOS)
    im.save(dst, optimize=True)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--keep", type=Path, help="work in this directory and keep it")
    a = ap.parse_args()
    ex = REPO / "data" / "example"
    env = dict(os.environ, EX=str(ex / "bulk"), SC=str(ex / "single_cell" / "pbmc_subset.h5ad"),
               SP=str(ex / "spatial" / "mouse_brain_anterior_subset.h5ad"),
               MPLBACKEND="Agg", PYTHONIOENCODING="utf-8")
    work = a.keep or Path(tempfile.mkdtemp(prefix="gcen2-gallery-"))
    work.mkdir(parents=True, exist_ok=True)
    (work / "fig").mkdir(exist_ok=True)
    (work / "wgcna").mkdir(exist_ok=True)
    OUT.mkdir(parents=True, exist_ok=True)
    try:
        groups = sorted({g[1] for g in content.GALLERY})
        for group in groups:
            print(f"[{group}] setup")
            for cmd in content.GALLERY_SETUP.get(group, []):
                run(cmd, work, env)
        for fid, group, _kind, title, _desc, cmd in content.GALLERY:
            print(f"[{group}] {title}")
            run(cmd, work, env)
            out = shlex.split(cmd)[shlex.split(cmd).index("-o") + 1]
            shrink(work / f"{out}.png", OUT / f"{fid}.png")
        stale = {p.stem for p in OUT.glob("*.png")} - {g[0] for g in content.GALLERY}
        for s in sorted(stale):
            (OUT / f"{s}.png").unlink()
            print(f"removed stale {s}.png")
        print(f"{len(content.GALLERY)} images in {OUT}")
    finally:
        if a.keep is None:
            shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
