"""Static checks of the built GCEN2 website (no browser).

usage: python website/tests/check_site.py [SITE_DIR]   (default website/_site)

Checks every HTML file: internal links and resources exist (and #fragments
exist as ids), images have alt, width and height, ids are unique; site pages
(not the Sphinx-generated ones) have lang="en", one <h1>, a title, a
description and the viewport meta; no external scripts or stylesheets; every
page loads Bootstrap 5.3.9. Exits 1 on any problem.
"""
from __future__ import annotations

import sys
from html.parser import HTMLParser
from pathlib import Path
from urllib.parse import unquote, urlsplit

BOOTSTRAP = "static/vendor/bootstrap-5.3.9/"
SITE_PATH = "/gcen2/"


class Page(HTMLParser):
    def __init__(self):
        super().__init__(convert_charrefs=True)
        self.links, self.ids, self.imgs, self.scripts, self.styles = [], [], [], [], []
        self.h1 = 0
        self.lang = None
        self.title = False
        self.meta = {}
        self._in_title = False

    def handle_starttag(self, tag, attrs):
        a = dict(attrs)
        if "id" in a:
            self.ids.append(a["id"])
        if tag == "html":
            self.lang = a.get("lang")
        elif tag == "a" and "href" in a:
            self.links.append(a["href"])
        elif tag == "img":
            self.imgs.append(a)
            self.links.append(a.get("src", ""))
        elif tag == "script" and "src" in a:
            self.scripts.append(a["src"])
            self.links.append(a["src"])
        elif tag == "link" and a.get("rel") in ("stylesheet", "icon"):
            if a.get("rel") == "stylesheet":
                self.styles.append(a["href"])
            self.links.append(a["href"])
        elif tag == "h1":
            self.h1 += 1
        elif tag == "title":
            self._in_title = True
        elif tag == "meta" and "name" in a:
            self.meta[a["name"]] = a.get("content", "")
        elif tag == "form" and "action" in a:
            self.links.append(a["action"])
        elif tag == "video":
            # posters are in the site; the MP4 files are served from media/, outside it
            for k in ("poster", "data-gc-poster-light", "data-gc-poster-dark"):
                if k in a:
                    self.links.append(a[k])

    def handle_endtag(self, tag):
        if tag == "title":
            self._in_title = False

    def handle_data(self, data):
        if self._in_title and data.strip():
            self.title = True


def main() -> int:
    site = Path(sys.argv[1] if len(sys.argv) > 1 else Path(__file__).resolve().parents[1] / "_site")
    if not (site / "index.html").exists():
        print(f"no site at {site}; run website/build.py")
        return 1
    pages = {p: Page() for p in sorted(site.rglob("*.html"))}
    for p, pg in pages.items():
        pg.feed(p.read_text(encoding="utf-8"))
    ids = {p.resolve(): set(pg.ids) for p, pg in pages.items()}
    problems = []

    def bad(p, msg):
        problems.append(f"{p.relative_to(site)}: {msg}")

    for p, pg in pages.items():
        is_docs = p.relative_to(site).parts[0] == "docs"
        seen = set()
        for i in pg.ids:
            if i in seen:
                bad(p, f"duplicate id '{i}'")
            seen.add(i)
        for href in pg.links:
            u = urlsplit(href)
            if u.scheme in ("http", "https", "mailto") or href.startswith("//"):
                continue
            if u.scheme or (not u.path and not u.fragment):
                if href not in ("#", ""):
                    bad(p, f"odd link '{href}'")
                continue
            if u.path.startswith(SITE_PATH):
                target = site / unquote(u.path[len(SITE_PATH):])
            elif u.path.startswith("/"):
                bad(p, f"absolute link outside the site '{href}'")
                continue
            elif u.path:
                target = (p.parent / unquote(u.path))
            else:
                target = p
            if target.is_dir():
                target = target / "index.html"
            if not target.exists():
                bad(p, f"broken link '{href}'")
                continue
            if (u.fragment and target.suffix == ".html"
                    and u.fragment not in ids.get(target.resolve(), set())):
                bad(p, f"missing anchor '{href}'")
        for img in pg.imgs:
            if "alt" not in img:
                bad(p, f"<img {img.get('src')}> without alt")
            if not is_docs and not ("width" in img and "height" in img):
                bad(p, f"<img {img.get('src')}> without width/height")
        for s in pg.scripts + pg.styles:
            if urlsplit(s).scheme or s.startswith("//"):
                bad(p, f"external resource '{s}'")
        if not any(BOOTSTRAP in s for s in pg.styles) or not any(BOOTSTRAP in s for s in pg.scripts):
            bad(p, "Bootstrap 5.3.9 not loaded")
        if pg.lang != "en":
            bad(p, f'lang="{pg.lang}"')
        if not pg.title:
            bad(p, "no <title>")
        if "viewport" not in pg.meta:
            bad(p, "no viewport meta")
        if not is_docs:
            if pg.h1 != 1:
                bad(p, f"{pg.h1} <h1> elements")
            if not pg.meta.get("description"):
                bad(p, "no description")
    for msg in problems:
        print(msg)
    print(f"{len(pages)} pages checked, {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())
