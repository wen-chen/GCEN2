# GCEN2 website — implementation

How the site designed in `GCEN2_Website_Design.md` is built, tested and
deployed.

## 1. Layout of `website/`

```
website/
  README.md                 build, test and deploy in short
  build.py                  builds the site into website/_site/
  content.py                page content that is not in the repository elsewhere
  templates/                Jinja2 templates
    base.html               <head>, navbar, footer, scripts
    partials/_chrome.html   head assets, navbar, footer, icons (shared with the Sphinx theme)
    partials/_blocks.html   code block with copy button, page header, pre-release notice
    index.html install.html gallery.html benchmarks.html faq.html 404.html
  static/
    css/gcen2.css           the site's styles on top of Bootstrap
    js/gcen2.js             colour mode, copy buttons, gallery filter and dialog
    img/                    logo.svg, favicon.svg, workflow.png, gallery/*.png
    vendor/bootstrap-5.3.9/ bootstrap.min.css, bootstrap.bundle.min.js, LICENSE
  sphinx_theme/gcen2/       Sphinx theme: Bootstrap layout around the documentation
  gallery/make_gallery.py   runs the tutorial commands on the example data, writes img/gallery/
  workflow/make_workflow.py the workflow figure (PowerPoint), rendered to img/workflow.png
  tests/
    check_site.py           static checks: links, images, HTML ids, page structure
    visual_test.py          Playwright: screenshots and in-browser checks
  _site/, _doctrees/        build output (ignored by git)
```

## 2. Build (`python website/build.py`)

1. Read facts from the repository (R4):
   - version from `pyproject.toml`;
   - the benchmark table from the `## Benchmarks` section of `README.md`
     (a Markdown table parsed into rows; the build fails if it is missing);
   - the number of commands and of figure commands from the CLI
     (`gcen2.cli.app` through Typer's Click command), so "24 commands" and
     "43 figures" cannot drift;
   - the citation from `CITATION.cff`.
2. Render each page template with Jinja2 into `_site/` (strict undefined, so a
   missing value fails the build). Links are relative (`./install.html`,
   `../index.html` from the documentation), so the site works under `/gcen2/`,
   at the root of any host and from disk.
3. Copy `static/` to `_site/static/`.
4. Build the documentation: `sphinx-build -b html docs _site/docs` with
   `-D html_theme=gcen2 -D html_theme_path=website/sphinx_theme` (the
   repository's `docs/conf.py` is otherwise unchanged; Read the Docs keeps its
   theme). `--no-docs` skips this step for quick work on the pages.

Release switch: `content.RELEASED` (False until GCEN2 is on PyPI and
Bioconda). While False, the pages show the pre-release notice and give the
source installation first.

### Sphinx theme

`sphinx_theme/gcen2/theme.toml` inherits Sphinx's `basic` theme.
`layout.html` replaces the page frame with the site's: Bootstrap 5.3.9, the
site navbar and footer (the macros of `templates/partials/_chrome.html`,
which `build.py` puts on Sphinx's `templates_path`; the GitHub URL and the
year come in with `-A`), a left sidebar with the global table of contents
(`toctree(collapse=False, maxdepth=2)`) and the search box, the document body,
and previous/next links. `static/gcen2-docs.css` styles Sphinx's generated
markup (admonitions, code blocks, tables, API signatures) to match the site.
`basic`'s own JavaScript (search, highlighting) is kept. Code is coloured
with Pygments' `default` style, and in dark mode with `github-dark`, whose
CSS `build.py` writes scoped to `[data-bs-theme="dark"]`
(`_static/pygments-dark.css`). The theme accepts, and ignores, the two
options `docs/conf.py` sets for `sphinx_rtd_theme`.

## 3. Content and assets

- **Workflow figure.** `workflow/make_workflow.py` (python-pptx) draws the
  slide; LibreOffice renders it to PDF and `pdftoppm` to a 2,400 px PNG
  (`img/workflow.png`). Both the PNG and the script are committed, so the
  site builds without LibreOffice.
- **Gallery.** `gallery/make_gallery.py` runs, in a temporary directory, the
  setup commands of the tutorials on `data/example/` and then each figure's
  `gcen2 plot` command exactly as `content.GALLERY` lists it (the gallery page
  shows the same text), and writes the PNGs, downscaled to at most 1,200 px
  wide, into `img/gallery/`. The images are committed (16 figures, 3.0 MB).
  Figures whose example data make them unreadable (32 bulk modules with
  Ensembl IDs) were replaced by the same figure on a smaller result.
- **Bootstrap.** The npm package `bootstrap@5.3.9`, checked against the
  registry's SHA-512, `dist/css/bootstrap.min.css` and
  `dist/js/bootstrap.bundle.min.js` (with Popper) and the MIT licence; the
  `sourceMappingURL` comments are removed (the maps are not shipped, and a
  browser would ask for them).
- **Colour mode.** A small script in `<head>` sets `data-bs-theme` from the
  saved choice or `prefers-color-scheme` before the page paints; the navbar
  switch toggles light and dark and saves the choice in `localStorage` (in
  `try`, so a blocked storage only loses the choice); without a saved choice
  the page follows the system setting, also when it changes.
- **Colours.** Brand teal `#27756b` for links and buttons (contrast 5.5 on
  white, 5.0 on the tinted page headers; `#6cc5b8` in dark mode); every text
  colour was chosen with `visual_test.py`'s contrast check.
- **Copy buttons.** `navigator.clipboard.writeText` with a fallback that
  selects the text; the button says "Copied" for 1.5 s.

## 4. Tests

`website/tests/check_site.py` (no browser; runs after every build):
- every internal `href` and `src` points to an existing file (fragments to an
  existing `id`);
- every `<img>` has `alt`, and `width` and `height`;
- every page has one `<h1>`, `lang="en"`, a `<title>`, a description and the
  viewport meta; ids are unique;
- no external `<script>` or stylesheet (R5);
- Bootstrap is 5.3.9 on every page.

`website/tests/visual_test.py` (Playwright, Chromium headless, started with
`--disable-gpu --disable-software-rasterizer`: on the development machine the
software rasterizer now and then never delivers a frame, and screenshots
hang; screenshots are also retried. Without a GPU Chromium paints only what
has been in the viewport, so each screenshot is taken with the viewport as
tall as the page):
- every page at 1440, 768 and 390 px wide, in light and dark mode, saved as
  full-page screenshots in `website/tests/screenshots/` (ignored by git) for
  inspection;
- in the browser: no horizontal overflow (`scrollWidth <= innerWidth`), no
  console errors or failed requests, all images decoded, the navbar
  collapses below 768 px and expands on the toggler, the colour-mode switch
  changes `data-bs-theme`, the gallery filter shows only the chosen analysis
  and the dialog opens, the copy button copies;
- WCAG contrast of body text, links and buttons in both modes, computed from
  the rendered colours.

The screenshots are then looked at, page by page (overflow, overlap,
alignment, spacing, contrast) before a release of the site.

## 5. Deployment

`_site/` is copied as is to the directory that the web server serves at
`/gcen2/` (for the Django site at biochen.org: a static location or a
`django.views.static.serve` route for `gcen2/`; `404.html` for its 404
handler). Nothing else is needed on the server.

## 6. Maintenance

- New benchmark rows, a new version or new commands: rebuild; the pages follow.
- New figure types: add them to `content.GALLERY` and run
  `make_gallery.py`.
- At release: set `content.RELEASED = True`, rebuild, deploy.
- Bootstrap updates: replace `static/vendor/bootstrap-<version>/` and the
  version in `base.html` and the theme; `check_site.py` checks the version.
