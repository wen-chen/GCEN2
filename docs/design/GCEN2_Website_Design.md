# GCEN2 website — design

The project website at <https://www.biochen.org/gcen2/>. This document says what
the site is for, what it contains and how it looks; the implementation document
(`GCEN2_Website_Implementation.md`) says how it is built and tested.

## 1. Purpose and audience

The site is where a biologist who heard of GCEN2 decides whether to use it and
gets it running, and where a user comes back for reference.

| Reader | Comes for | Finds it on |
|---|---|---|
| New user, biologist | what GCEN2 does, whether it fits the data | Home |
| New user | installing it on Linux, Windows or macOS | Install |
| GCEN v1 user | what changes, whether old pipelines still work | Home, FAQ, migration guide |
| User at work | commands, options, functions, tutorials | Documentation, Tutorials |
| User choosing figures | what each figure looks like and the command | Gallery |
| Method developer, evaluator | evidence of correctness and speed | Benchmarks |
| Everyone | citation, licence, contact, source code | Home, FAQ, footer |

## 2. Requirements

- **R1 Static.** Plain HTML, CSS and JavaScript; no server code, no database.
  It can be served by the existing biochen.org server (or any static host)
  under `/gcen2/`, and opened from disk.
- **R2 Bootstrap 5.3.9.** The front-end framework is Bootstrap v5.3.9 (CSS
  and the JavaScript bundle), on every page including the documentation.
- **R3 English.**
- **R4 One source of truth.** Facts that live in the repository (the
  benchmark table, the installation matrix, the version, the number of
  commands and figures, the documentation) are read from it at build time,
  not copied by hand.
- **R5 Self-contained.** Bootstrap and every asset are served from the site;
  no CDN, web fonts or third-party scripts (many readers are in mainland
  China, where common CDNs and Google services are slow or blocked). No
  analytics by default.
- **R6 Honest status.** Until GCEN2 is on PyPI and Bioconda, the site says so
  and gives the source installation; one switch changes this at release.
- **R7 Public content only.** The site shows only what the public repository
  holds.
- **R8 Accessible and responsive.** Readable from 360 px phones to wide
  screens without horizontal scrolling; light and dark colour modes; WCAG AA
  contrast; keyboard navigation; text alternatives for images.
- **R9 Continuity with GCEN.** A GCEN user recognises the site: the same
  BioChen colours, the same kinds of pages (home, documentation, tutorial,
  download, visualization, FAQ) and the same workflow figure language.

## 3. Information architecture

URLs are relative to `https://www.biochen.org/gcen2/`.

| Page | URL | GCEN v1 counterpart | Content |
|---|---|---|---|
| Home | `index.html` (`/gcen2/`) | `/gcen` | what GCEN2 is, install command, workflow figure, six capabilities, headline benchmarks, quick start, citation |
| Install | `install.html` | `/gcen/download` | status, platform matrix, PyPI, Bioconda, source, optional parts, check, example data, GO / KEGG data |
| Documentation | `docs/index.html` | `/gcen/documentation` | the Sphinx documentation: quickstart, CLI reference, Python API, migration guide |
| Tutorials | `docs/tutorials/*.html` | `/gcen/tutorial` | bulk, single cell, spatial, trajectories and multi-omics, figures and reports |
| Gallery | `gallery.html` | `/gcen/visualization` | figures made from the example data, each with its command, filterable by analysis |
| Benchmarks | `benchmarks.html` | — | the comparison with GCEN v1, WGCNA, hdWGCNA, DGCA, RandomWalkRestartMH; conditions |
| FAQ | `faq.html` | `/gcen/faq` | questions, including the GCEN v1 ones that still apply |
| Not found | `404.html` | — | for the web server's 404 handler |

Navigation bar (all pages, the documentation included): brand "GCEN2", then
Home · Install · Documentation · Tutorials (menu) · Gallery · Benchmarks · FAQ,
and at the right GitHub and the colour-mode switch. The current page is marked
(`aria-current="page"`). On narrow screens the links collapse behind the
toggler.

Footer (all pages): GCEN2 and the BioChen home page, links to GitHub, issues,
GCEN v1 and the licence; "GPL-3.0-or-later"; copyright.

## 4. Content

### Home

1. **Hero.** "GCEN2", the tagline *Gene co-expression network analysis from
   bulk to single-cell and spatial transcriptomes*, two sentences on what it
   is, buttons *Install*, *Get started* (quickstart) and *GitHub*, and the
   install command with a copy button (tabs pip / conda / source).
2. **Workflow.** The GCEN2 workflow figure (the same drawing as the
   `GCEN2 workflow.pptx` slide), with a caption.
3. **What it does.** Six cards, one per kind of analysis: GCEN v1 networks;
   WGCNA modules; single cells (hdWGCNA); spatial transcriptomes;
   trajectories and multi-omics; figures and reports. Each names its commands
   and links to its tutorial.
4. **Same results, faster.** Four figures from the benchmark table (GCEN v1,
   WGCNA, the lung atlas, DGCA) and a link to Benchmarks.
5. **Quick start.** The GCEN v1 pipeline in seven commands, and the same in
   Python.
6. **Cite and contact.** Citation (GCEN2 per `CITATION.cff`, GCEN v1 paper),
   licence, contact.

### Install

A notice while unreleased (R6); the matrix PyPI / Bioconda / source × Linux
x86-64, Windows x86-64, macOS Apple Silicon; the three installation methods;
the optional parts (`single-cell`, `visualization`, `report`, `all`) and what
each adds; how to check the installation; the example data in the repository;
GO and KEGG annotation files (as on the GCEN v1 download page).

### Documentation and tutorials

The Sphinx documentation in `docs/` built with a site theme: the site's navbar
and footer, a left sidebar with the documentation tree, the page in the
centre, and Sphinx's search. The text is not duplicated on the site.

### Gallery

Figure cards in a grid: image, title, the analysis, and the `gcen2 plot`
command; filter buttons by analysis (All, Networks, WGCNA modules, Single
cell, Spatial, Function); click opens the figure large in a dialog. Every
figure is produced by `website/gallery/make_gallery.py` from the example data
with the commands of the tutorials, so each image matches its command.

### Benchmarks

The README's benchmark table (tool, data, results, speed), the conditions
(one 16-core machine, idle, the versions compared) and links to the scripts
and results in `benchmarks/` on GitHub.

### FAQ

What is GCEN2; how it differs from GCEN v1; whether GCEN v1 pipelines still
work; what data it takes (TSV, CSV, h5ad); which analyses; Windows and conda;
threads and reproducibility; large atlases; GO / KEGG annotation; figures;
how to cite; where to ask.

## 5. Visual design

- **Character.** A scientific software site in the manner of current
  bioinformatics tools: white page, generous space, one brand colour, quiet
  greys, real figures as the illustrations; no stock images, gradients or
  decorative stripes.
- **Colour.** BioChen's teal as the brand colour (navbar `#2d6760`, footer
  `#073839`, links `#27756b` — a darker step of the v1 link colour `#43aea0`
  for AA contrast on white); Bootstrap's greys. The workflow figure's lane
  colours (Okabe–Ito blue, vermilion, green, purple) mark the analysis kinds
  in the capability cards and the gallery filters. Dark mode uses Bootstrap
  5.3's `data-bs-theme="dark"` with lighter brand steps.
- **Type.** The system sans-serif stack (Bootstrap's default) for text, the
  system monospace stack for commands; no web fonts (R5).
- **Components.** Bootstrap navbar, buttons, cards, nav tabs, tables,
  accordion (FAQ), modal (gallery), badges; one custom element, the code
  block with a copy button.
- **Icons.** A few inline SVG icons (GitHub, copy, sun/moon); no icon font.
- **Logo and favicon.** A small SVG mark: three connected nodes in the brand
  teal, with the word GCEN2.

## 6. Out of scope

Blog, news, user accounts, comments, online analysis, Chinese translation,
analytics (can be added in one template if wanted).
