# GCEN2 website

The static website at <https://www.biochen.org/gcen2/> (Bootstrap 5.3.9, English).
Design: `docs/design/GCEN2_Website_Design.md`; implementation:
`docs/design/GCEN2_Website_Implementation.md`.

## Build

```bash
pip install -e ".[all]" sphinx jinja2 pillow    # gcen2 itself, for the command counts
python website/build.py                        # -> website/_site/ (pages + docs/)
python website/build.py --no-docs              # pages only, quick
```

## Test

```bash
python website/tests/check_site.py             # links, images, page structure; no browser
pip install playwright && python -m playwright install chromium
python website/tests/visual_test.py            # screenshots in website/tests/screenshots/
```

`visual_test.py` checks 12 pages at 1440, 768 and 390 px in light and dark
mode (overflow, console errors, failed requests, images, WCAG contrast) and
the interactions; then look at the screenshots.

## Update

| When | Do |
|---|---|
| benchmark table, version, commands change | rebuild (they are read from README.md, pyproject.toml and the CLI) |
| gallery figures change | edit `content.GALLERY`, run `python website/gallery/make_gallery.py` |
| workflow figure changes | edit `workflow/make_workflow.py`, render it (see its docstring), replace `static/img/workflow.png` |
| GCEN2 is released on PyPI and Bioconda | `content.RELEASED = True`, rebuild |

## Deploy

Copy `website/_site/` to the directory served at `/gcen2/`; serve `404.html`
for missing pages under `/gcen2/`. Everything (Bootstrap included) is in
`_site/`; nothing is loaded from other hosts.

The introduction video on the home page is not in the repository: its two MP4 files
(`GCEN2_intro.mp4`, light; `GCEN2_intro_dark.mp4`, dark) are served from
`/gcen2/media/` (paths in `content.VIDEO`). Keep that directory when deploying, e.g.
`rsync -a --delete --exclude media/ website/_site/ SERVER:/path/to/gcen2/`.
Serve them as `video/mp4` with range requests (the default of Apache and nginx), so
browsers can seek.
