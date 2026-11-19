"""Visual and in-browser tests of the built GCEN2 website (Playwright, Chromium).

usage: python website/tests/visual_test.py [--site DIR] [--out DIR] [--pages a,b]

Serves the site at http://127.0.0.1:<port>/gcen2/ (as on biochen.org), then for
every page, at 1440, 768 and 390 px wide, in light and dark mode:
- saves a full-page screenshot (website/tests/screenshots/<page>-<width>-<mode>.png);
- checks: no horizontal overflow, no console errors or failed requests, every
  image decoded, WCAG contrast >= 4.5 for text (3 for large text) of headings,
  paragraphs, list items, links, buttons, table cells and labels.
Then the interactions: the navbar toggler (390 px), the colour-mode switch, the
gallery filter and dialog, the copy button. Exits 1 on any failure.
Needs: pip install playwright && python -m playwright install chromium.
"""
from __future__ import annotations

import argparse
import functools
import http.server
import socketserver
import sys
import tempfile
import threading
from pathlib import Path

HERE = Path(__file__).resolve().parent
SITE = HERE.parent / "_site"

PAGES = {
    "home": "index.html",
    "install": "install.html",
    "gallery": "gallery.html",
    "benchmarks": "benchmarks.html",
    "faq": "faq.html",
    "404": "404.html",
    "docs": "docs/index.html",
    "quickstart": "docs/quickstart.html",
    "tutorial-sc": "docs/tutorials/single_cell.html",
    "cli": "docs/cli.html",
    "api": "docs/api.html",
    "search": "docs/search.html?q=module",
}
WIDTHS = (1440, 768, 390)
MODES = ("light", "dark")

CONTRAST_JS = r"""
() => {
  function parse(c) {
    const m = c.match(/rgba?\(([^)]+)\)/); if (!m) return null;
    const p = m[1].split(/[ ,\/]+/).filter(Boolean).map(Number);
    return {r: p[0], g: p[1], b: p[2], a: p.length > 3 ? p[3] : 1};
  }
  function lum(c) {
    const f = v => { v /= 255; return v <= 0.03928 ? v / 12.92 : Math.pow((v + 0.055) / 1.055, 2.4); };
    return 0.2126 * f(c.r) + 0.7152 * f(c.g) + 0.0722 * f(c.b);
  }
  function blend(top, bottom) {
    const a = top.a;
    return {r: top.r * a + bottom.r * (1 - a), g: top.g * a + bottom.g * (1 - a), b: top.b * a + bottom.b * (1 - a), a: 1};
  }
  function bg(el) {
    const stack = [];
    for (let e = el; e; e = e.parentElement) {
      const s = getComputedStyle(e);
      if (s.backgroundImage && s.backgroundImage !== "none" && !s.backgroundImage.startsWith("url")) return null;
      const c = parse(s.backgroundColor);
      if (c && c.a > 0) { stack.push(c); if (c.a >= 1) break; }
    }
    let out = {r: 255, g: 255, b: 255, a: 1};
    if (document.documentElement.getAttribute("data-bs-theme") === "dark" && !stack.some(c => c.a >= 1)) {
      out = parse(getComputedStyle(document.body).backgroundColor) || out;
    }
    for (let i = stack.length - 1; i >= 0; i--) out = blend(stack[i], out);
    return out;
  }
  const sel = "h1,h2,h3,h4,h5,h6,p,li,a,button,th,td,label,dt,dd,figcaption,.gc-stat-value,.gc-stat-label,.gc-badge,code";
  const res = [];
  for (const el of document.querySelectorAll(sel)) {
    if (!el.offsetParent && getComputedStyle(el).position !== "fixed") continue;
    const r = el.getBoundingClientRect(); if (r.width === 0 || r.height === 0) continue;
    const own = Array.from(el.childNodes).some(n => n.nodeType === 3 && n.textContent.trim());
    if (!own) continue;
    if (el.closest(".visually-hidden,.visually-hidden-focusable,[aria-hidden=true],.highlight,.gc-figure,.modal:not(.show)")) continue;
    const s = getComputedStyle(el);
    if (s.visibility === "hidden" || parseFloat(s.opacity) < 0.5) continue;
    const fg = parse(s.color); const b = bg(el); if (!fg || !b) continue;
    const f = blend(fg, b);
    const L1 = lum(f), L2 = lum(b);
    const ratio = (Math.max(L1, L2) + 0.05) / (Math.min(L1, L2) + 0.05);
    const size = parseFloat(s.fontSize), bold = parseInt(s.fontWeight) >= 700;
    const large = size >= 24 || (bold && size >= 18.66);
    const need = large ? 3 : 4.5;
    if (ratio < need) res.push({tag: el.tagName.toLowerCase(), cls: el.className.toString().slice(0, 40), text: el.textContent.trim().slice(0, 40), ratio: Math.round(ratio * 100) / 100, need});
  }
  return res;
}
"""


def serve(site: Path):
    root = Path(tempfile.mkdtemp(prefix="gcen2-site-"))
    (root / "gcen2").symlink_to(site.resolve(), target_is_directory=True)
    class Quiet(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *args):
            pass

    handler = functools.partial(Quiet, directory=str(root))
    httpd = socketserver.TCPServer(("127.0.0.1", 0), handler)
    threading.Thread(target=httpd.serve_forever, daemon=True).start()
    return httpd, f"http://127.0.0.1:{httpd.server_address[1]}/gcen2/"


def shoot(page, path: Path, **kw) -> None:
    """A screenshot, retried: headless Chromium on some hosts now and then never
    delivers the frame a screenshot waits for."""
    for attempt in range(3):
        try:
            page.screenshot(path=str(path), timeout=20000, **kw)
            return
        except Exception as e:  # playwright TimeoutError
            if "Timeout" not in str(e) or attempt == 2:
                raise
            page.wait_for_timeout(1000)


def main() -> int:
    from playwright.sync_api import sync_playwright

    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--site", type=Path, default=SITE)
    ap.add_argument("--out", type=Path, default=HERE / "screenshots")
    ap.add_argument("--pages", help="comma-separated page keys (default: all)")
    a = ap.parse_args()
    pages = {k: v for k, v in PAGES.items() if not a.pages or k in a.pages.split(",")}
    a.out.mkdir(parents=True, exist_ok=True)
    httpd, base = serve(a.site)
    failures: list[str] = []

    def fail(msg):
        failures.append(msg)
        print("FAIL", msg)

    with sync_playwright() as pw:
        # Chromium's software rasterizer hangs screenshots on some headless Linux hosts.
        browser = pw.chromium.launch(args=["--disable-gpu", "--disable-software-rasterizer"])
        for mode in MODES:
            for width in WIDTHS:
                ctx = browser.new_context(viewport={"width": width, "height": 900},
                                          color_scheme=mode, device_scale_factor=1)
                page = ctx.new_page()
                errors: list[str] = []
                page.on("console", lambda m, errs=errors: errs.append(
                    f"console {m.type}: {m.text}") if m.type == "error" else None)
                page.on("pageerror", lambda e, errs=errors: errs.append(f"page error: {e}"))
                page.on("requestfailed", lambda r, errs=errors: errs.append(
                    f"request failed: {r.url}"))
                page.on("response", lambda r, errs=errors: errs.append(
                    f"HTTP {r.status}: {r.url}") if r.status >= 400 else None)
                for key, path in pages.items():
                    errors.clear()
                    tag = f"{key}-{width}-{mode}"
                    page.goto(base + path, wait_until="networkidle")
                    page.evaluate("document.querySelectorAll('img').forEach(i => { i.loading = 'eager'; i.decoding = 'sync'; })")
                    page.wait_for_load_state("networkidle")
                    page.wait_for_function("Array.from(document.images).every(i => i.complete)")
                    theme = page.evaluate("document.documentElement.getAttribute('data-bs-theme')")
                    if theme != mode:
                        fail(f"{tag}: data-bs-theme={theme}")
                    broken = page.evaluate("Array.from(document.images).filter(i => i.offsetParent && i.naturalWidth === 0).map(i => i.src)")
                    if broken:
                        fail(f"{tag}: images not decoded {broken}")
                    over = page.evaluate("[document.documentElement.scrollWidth, window.innerWidth]")
                    if over[0] > over[1]:
                        wide = page.evaluate("""() => Array.from(document.querySelectorAll('body *')).filter(e => e.getBoundingClientRect().right > window.innerWidth + 1 && !e.closest('.table-responsive,.gc-figure-scroll,pre,.highlight,table.docutils')).slice(0, 5).map(e => e.tagName + '.' + e.className.toString().slice(0, 30))""")
                        fail(f"{tag}: horizontal overflow {over[0]} > {over[1]} {wide}")
                    low = page.evaluate(CONTRAST_JS)
                    if low:
                        fail(f"{tag}: low contrast {low[:4]}")
                    for e in errors:
                        if key == "404" or "favicon" not in e:
                            fail(f"{tag}: {e}")
                    # Without a GPU, headless Chromium paints only what has been in the
                    # viewport: make the viewport as tall as the page for the screenshot.
                    h = page.evaluate("document.documentElement.scrollHeight")
                    page.set_viewport_size({"width": width, "height": min(h, 16000)})
                    page.wait_for_timeout(300)
                    shoot(page, a.out / f"{tag}.png")
                    page.set_viewport_size({"width": width, "height": 900})
                ctx.close()

        # Interactions
        ctx = browser.new_context(viewport={"width": 390, "height": 844}, color_scheme="light")
        ctx.grant_permissions(["clipboard-read", "clipboard-write"])
        page = ctx.new_page()
        page.goto(base + "index.html")
        if page.is_visible("#gcNav .nav-link[href='install.html']"):
            fail("navbar: links visible before toggling at 390 px")
        page.click(".navbar-toggler")
        page.wait_for_selector("#gcNav.show")
        page.wait_for_timeout(400)
        if not page.is_visible("#gcNav .nav-link[href='install.html']"):
            fail("navbar: toggler does not show the links")
        shoot(page, a.out / "interaction-navbar-open-390.png")
        page.click("[data-gc-theme-toggle]")
        if page.evaluate("document.documentElement.getAttribute('data-bs-theme')") != "dark":
            fail("colour mode: switch does not change to dark")
        page.reload()
        if page.evaluate("document.documentElement.getAttribute('data-bs-theme')") != "dark":
            fail("colour mode: choice not kept after reload")
        page.evaluate("localStorage.clear()")
        page.reload()
        page.click("#pane-source .gc-copy") if page.is_visible("#pane-source .gc-copy") else page.click(".gc-copy >> nth=0")
        page.wait_for_timeout(200)
        copied = page.evaluate("navigator.clipboard.readText()")
        if "gcen2" not in copied.lower() or "$ " in copied:
            fail(f"copy button: clipboard has {copied!r}")
        ctx.close()

        ctx = browser.new_context(viewport={"width": 1440, "height": 900}, color_scheme="light")
        page = ctx.new_page()
        page.goto(base + "gallery.html")
        total = page.locator("[data-gc-kind]").count()
        page.click("[data-gc-show=spatial]")
        shown = page.locator("[data-gc-kind]:not(.d-none)").count()
        spatial = page.locator("[data-gc-kind=spatial]").count()
        if shown != spatial or shown == total:
            fail(f"gallery filter: {shown} shown, {spatial} spatial of {total}")
        if page.inner_text("#galleryStatus") != f"{spatial} figures":
            fail("gallery filter: status text not updated")
        page.click("[data-gc-show=all]")
        page.click(".gc-thumb-img >> nth=2")
        page.wait_for_selector("#figureModal.show")
        page.wait_for_timeout(400)
        cmd = page.inner_text("[data-gc-modal-cmd]")
        if not cmd.startswith("gcen2 plot"):
            fail(f"gallery dialog: command {cmd!r}")
        shoot(page, a.out / "interaction-gallery-dialog-1440.png")
        page.keyboard.press("Escape")
        page.wait_for_selector("#figureModal", state="hidden")
        ctx.close()
        browser.close()
    httpd.shutdown()
    n = len(pages) * len(WIDTHS) * len(MODES)
    print(f"{n} page views and the interactions checked; {len(failures)} failures; "
          f"screenshots in {a.out}")
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
