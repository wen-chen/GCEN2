"""Build the GCEN2 website into website/_site/.

usage: python website/build.py [--no-docs] [--out DIR]

Pages are Jinja2 templates (website/templates/) filled with facts read from
the repository (version, benchmark and installation tables of README.md,
number of commands and figures, citation) and with website/content.py; the
documentation is the Sphinx project in docs/, built with the site's theme
(website/sphinx_theme/gcen2). See docs/design/GCEN2_Website_Implementation.md.
"""
from __future__ import annotations

import argparse
import datetime
import html
import re
import shutil
import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
REPO = HERE.parent
sys.path.insert(0, str(HERE))
import content  # noqa: E402

BOOTSTRAP = "5.3.9"
SITE_PATH = "/gcen2/"  # where the site is served; only the 404 page uses it

PAGES = [
    # template, output, nav key, title, description
    ("index.html", "index.html", "home", "GCEN2",
     "GCEN2: gene co-expression network analysis from bulk to single-cell and spatial "
     "transcriptomes. Command-line tool and Python package; successor of GCEN."),
    ("install.html", "install.html", "install", "Install",
     "Install GCEN2 from PyPI, Bioconda or the source code on Linux, Windows and macOS."),
    ("gallery.html", "gallery.html", "gallery", "Gallery",
     "Figures made by gcen2 plot from the example data, each with its command."),
    ("benchmarks.html", "benchmarks.html", "benchmarks", "Benchmarks",
     "GCEN2 compared with GCEN v1, WGCNA, hdWGCNA, DGCA and RandomWalkRestartMH: "
     "the same results, faster."),
    ("faq.html", "faq.html", "faq", "FAQ",
     "Frequently asked questions about GCEN2."),
    ("404.html", "404.html", "", "Page not found",
     "The page was not found on the GCEN2 website."),
]


# --- facts from the repository ------------------------------------------------------
def read_version() -> str:
    text = (REPO / "pyproject.toml").read_text(encoding="utf-8")
    m = re.search(r'^version\s*=\s*"([^"]+)"', text, re.M)
    if not m:
        raise SystemExit("no version in pyproject.toml")
    return m.group(1)


def read_extras() -> dict[str, list[str]]:
    try:
        import tomllib
    except ModuleNotFoundError:  # Python < 3.11
        import tomli as tomllib
    data = tomllib.loads((REPO / "pyproject.toml").read_text(encoding="utf-8"))
    return data["project"]["optional-dependencies"]


def markdown_table(section: str, first_header: str) -> tuple[list[str], list[list[str]]]:
    """The first Markdown table in README section `section` whose header starts with
    `first_header`, as (header, rows); cells keep inline Markdown converted to HTML."""
    text = (REPO / "README.md").read_text(encoding="utf-8")
    m = re.search(rf"^## {re.escape(section)}\n(.*?)(?=^## |\Z)", text, re.M | re.S)
    if not m:
        raise SystemExit(f"README.md has no section '## {section}'")
    lines = [ln for ln in m.group(1).splitlines() if ln.startswith("|")]
    for i, ln in enumerate(lines):
        cells = [c.strip() for c in ln.strip("|").split("|")]
        if cells[0] == first_header or (first_header == "" and cells[0] == ""):
            header = cells
            rows = []
            for ln2 in lines[i + 2:]:
                rows.append([inline_md(c.strip()) for c in ln2.strip().strip("|").split("|")])
            return header, rows
    raise SystemExit(f"README.md section '{section}' has no table starting with '{first_header}'")


def inline_md(s: str) -> str:
    s = html.escape(s, quote=False)
    s = re.sub(r"`([^`]+)`", r"<code>\1</code>", s)
    s = re.sub(r"\*\*([^*]+)\*\*", r"<strong>\1</strong>", s)
    s = re.sub(r"\[([^\]]+)\]\((https?://[^)]+)\)", r'<a href="\2">\1</a>', s)
    return s


def cli_counts() -> tuple[int, int]:
    """(number of gcen2 commands, number of `gcen2 plot` figures), from the CLI itself."""
    import typer.main

    from gcen2.cli import app

    cmd = typer.main.get_command(app)
    plot = cmd.commands["plot"]
    return len(cmd.commands), len(plot.commands)


def read_citation(version: str) -> str:
    text = (REPO / "CITATION.cff").read_text(encoding="utf-8")
    head = text.split("\nreferences:")[0]
    given = re.search(r"given-names:\s*(.+)", head).group(1).strip()
    family = re.search(r"family-names:\s*(.+)", head).group(1).strip()
    title = re.search(r'^title:\s*"(.+)"', head, re.M).group(1)
    url = re.search(r'^repository-code:\s*"(.+)"', head, re.M).group(1)
    initials = "".join(p[0] for p in given.split())
    return f"{family} {initials}. {title}. Version {version}. {url}"


V1_CITATION = ("Chen W, Li J, Huang S, Li X, Zhan X, Hu X, Xiang S, Liu C. GCEN: an easy-to-use "
               "toolkit for gene co-expression network analysis and lncRNAs annotation. "
               "Curr Issues Mol Biol. 2022;44(4):1479–1487.")


# --- assembling the pages ---------------------------------------------------------------
def gallery_items() -> list[dict]:
    from PIL import Image

    items = []
    for fid, _group, kind, title, desc, cmd in content.GALLERY:
        path = HERE / "static" / "img" / "gallery" / f"{fid}.png"
        if not path.exists():
            raise SystemExit(f"missing {path}: run website/gallery/make_gallery.py")
        w, h = Image.open(path).size
        items.append(dict(id=fid, kind=kind, title=title, description=desc, command=cmd,
                          alt=f"{title}: {desc}", w=w, h=h))
    return items


def headlines(bench_rows: list[list[str]]) -> list[dict]:
    out = []
    for key, value, label in content.HEADLINES:
        parts = [p.strip() for p in key.split("|")]
        rows = [r for r in bench_rows if all(any(c.startswith(p) or p in c for c in r)
                                             for p in parts)]
        if not rows:
            raise SystemExit(f"headline '{key}': no such row in the README benchmark table")
        if not any(value in re.sub("<[^>]+>", "", c) for c in rows[0]):
            raise SystemExit(f"headline '{value}' is not in the README row: {rows[0]}")
        out.append(dict(value=value, label=label))
    return out


EXTRA_TEXT = {
    "single-cell": "single-cell and spatial analyses, .h5ad files",
    "visualization": "figures (gcen2 plot) and interactive network views",
    "report": "HTML and PDF reports (gcen2 report)",
    "all": "everything above",
}


def context() -> dict:
    version = read_version()
    n_commands, n_plots = cli_counts()
    bench_header, bench_rows = markdown_table("Benchmarks", "Compared with")
    inst_header, inst_rows = markdown_table("Installation", "")
    extras = read_extras()
    extra_rows = []
    for name, adds in EXTRA_TEXT.items():
        pkgs = extras[name]
        if name == "all":
            shown = "the four extras"
        else:
            shown = ", ".join(re.split(r"[<>=!~\[ ;]", p)[0] for p in pkgs)
        extra_rows.append(dict(name=f"gcen2[{name}]", adds=adds, packages=shown))
    released = content.RELEASED
    site = dict(
        url=content.SITE_URL, github=content.GITHUB, email=content.CONTACT_EMAIL,
        version=version, released=released, install_default="pip" if released else "source",
        n_commands=n_commands, n_plots=n_plots, year=datetime.date.today().year,
    )
    caps = []
    for c in content.CAPABILITIES:
        c = dict(c)
        c["text"] = c["text"].format(n_plots=n_plots)
        caps.append(c)
    faq = [dict(question=q, answer=a.format(v1_site=content.GCEN_V1_SITE, github=content.GITHUB,
                                            email=content.CONTACT_EMAIL, n_plots=n_plots))
           for q, a in content.FAQ]
    gallery = gallery_items()
    kinds_used = [k for k in content.KINDS if any(g["kind"] == k for g in gallery)]
    return dict(
        site=site, kinds=content.KINDS, capabilities=caps,
        headlines=headlines(bench_rows),
        bench=dict(header=bench_header, rows=bench_rows),
        platforms=[h.replace(" / AMD64", "").replace("-", "\u2011") for h in inst_header[1:]],
        matrix=[dict(method=r[0], cells=[c.split(" (")[0] for c in r[1:]]) for r in inst_rows],
        matrix_notes=[c.split(" (", 1)[1].rstrip(")") for r in inst_rows for c in r[1:]
                      if " (" in c],
        extras=extra_rows, faq=faq, gallery=gallery, gallery_kinds=kinds_used,
        video=content.VIDEO,
        citation=dict(software=read_citation(version), v1=V1_CITATION),
    )


def render_pages(out: Path, ctx: dict) -> None:
    import jinja2

    env = jinja2.Environment(loader=jinja2.FileSystemLoader(HERE / "templates"),
                             autoescape=True, undefined=jinja2.StrictUndefined,
                             trim_blocks=False, lstrip_blocks=False, keep_trailing_newline=True)
    for tpl, name, active, title, desc in PAGES:
        root = SITE_PATH if name == "404.html" else ""
        canonical = "" if name == "index.html" else name
        page = env.get_template(tpl).render(root=root, active=active, title=title,
                                            description=desc, canonical=canonical, **ctx)
        (out / name).write_text(page, encoding="utf-8", newline="\n")
        print(f"  {name}")


def build_docs(out: Path, ctx: dict) -> None:
    from pygments.formatters import HtmlFormatter

    docs_out = out / "docs"
    cmd = [sys.executable, "-m", "sphinx", "-b", "html", "-q", "-E",
           "-D", "html_theme=gcen2",
           "-D", f"html_theme_path={HERE / 'sphinx_theme'}",
           "-D", f"templates_path={HERE / 'templates' / 'partials'}",
           "-A", f"github={content.GITHUB}",
           "-A", f"site_year={ctx['site']['year']}",
           "-d", str(HERE / "_doctrees"),
           str(REPO / "docs"), str(docs_out)]
    subprocess.run(cmd, check=True, cwd=REPO)
    css = HtmlFormatter(style="github-dark").get_style_defs('[data-bs-theme="dark"] .highlight')
    (docs_out / "_static" / "pygments-dark.css").write_text(css + "\n", encoding="utf-8")
    for junk in (".buildinfo", ".buildinfo.bak"):
        (docs_out / junk).unlink(missing_ok=True)
    print("  docs/")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--no-docs", action="store_true", help="skip the Sphinx documentation")
    ap.add_argument("--out", type=Path, default=HERE / "_site")
    a = ap.parse_args()
    ctx = context()
    out = a.out
    if out.exists():
        shutil.rmtree(out)
    out.mkdir(parents=True)
    shutil.copytree(HERE / "static", out / "static")
    print(f"building {out} (GCEN2 {ctx['site']['version']}, Bootstrap {BOOTSTRAP}, "
          f"{ctx['site']['n_commands']} commands, {ctx['site']['n_plots']} figures)")
    render_pages(out, ctx)
    if not a.no_docs:
        build_docs(out, ctx)


if __name__ == "__main__":
    main()
