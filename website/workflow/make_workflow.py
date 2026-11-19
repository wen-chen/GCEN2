"""GCEN2 workflow, one 16:9 slide of editable shapes (python-pptx).

Swimlanes: four kinds of data (rows) through five steps (columns), a shared
column for gene function and outputs, and the engine underneath. The shape
grammar follows the GCEN v1 workflow: data = snipped-corner box, program =
rounded box with the command in a monospace font, step = chevron header.
"""
import sys

from lxml import etree
from pptx import Presentation
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_CONNECTOR, MSO_SHAPE
from pptx.enum.text import MSO_ANCHOR, PP_ALIGN
from pptx.oxml.ns import qn
from pptx.util import Inches, Pt

OUT = sys.argv[1] if len(sys.argv) > 1 else "GCEN2_workflow.pptx"
SANS, MONO = "Arial", "Courier New"
INK, MUTED, LINE, HEAD = "222222", "5F6368", "9AA0A6", "E3E6EA"
ENGINE = "F1F3F4"

LANES = [  # name, colour, tint
    ("Bulk\nRNA-seq", "0072B2", "E8F1F8"),
    ("Single cells", "D55E00", "FBEEE6"),
    ("Spatial", "009E73", "E5F4EF"),
    ("Trajectories\n& multi-omics", "B4527F", "F7ECF2"),
]
STEPS = ["Input", "Preprocessing", "Network construction", "Module identification",
         "Comparison & dynamics"]

# cell = (kind, payload): data -> text; method -> text; cmds -> (commands, caption, row?)
CELLS = [
    [("data", "Expression matrix\ngenes × samples"),
     ("cmds", (["normalize", "filter"], "TMM · DESeq · UQ ·\nmedian · HKG", True)),
     ("cmds", (["build"], "Pearson / Spearman\np-value · FDR", False)),
     ("cmds", (["module", "module -e"], "GCEN · WGCNA", False)),
     ("cmds", (["diff", "consensus-modules"], "DGCA test · data sets", False))],
    [("data", "AnnData (.h5ad)\ncell types · backed"),
     ("method", "kNN metacells\neffective samples"),
     ("cmds", (["celltype"], "per cell type\nFDR · --diff", False)),
     ("cmds", (["celltype-modules"], "hdWGCNA: cell MEs,\nkME, Harmony, DMEs", False)),
     ("cmds", (["project-modules", "module-preservation"], "projection · Zsummary", False))],
    [("data", "Visium · Visium HD\nSlide-seq (.h5ad)"),
     ("method", "metaspots (grid, kNN)\nMoran's I test"),
     ("cmds", (["spatial-build"], "filtered · weighted ·\nregions", False)),
     ("cmds", (["spatial-modules"], "spot eigengenes, kME", False)),
     ("method", "module–region\ncorrelations")],
    [("data", "Pseudotime ·\nprior networks"),
     ("method", "sliding windows\nlayer weights"),
     ("cmds", (["trajectory", "fuse"], "", False)),
     ("method", "edge dynamics: gained\n· lost · transient"),
     ("cmds", (["rwr-multiplex"], "multiplex RWR", False))],
]
ENGINE_ITEMS = [
    ("C++17 kernels", "multithreaded, GIL released"),
    ("Same results", "as GCEN v1, WGCNA, hdWGCNA,\nDGCA and RandomWalkRestartMH"),
    ("Deterministic", "any thread count; Linux,\nWindows, macOS"),
    ("Atlas scale", "584,944 cells, 49 cell types\nin about 8 min"),
    ("CLI + Python API", "24 commands; every step\nalso in Python"),
]

# --- geometry (inches) --------------------------------------------------------
W, H = 13.333, 7.5
M = 0.45
LABEL_W = 1.42
X0 = M + LABEL_W + 0.12
SHARED_W = 1.90
SHARED_X = W - M - SHARED_W
GAP = 0.20
LANES_END = SHARED_X - 0.28
COL_W = (LANES_END - X0 - 4 * GAP) / 5
HEAD_Y, HEAD_H = 1.28, 0.40
LANE_Y0, LANE_H, LANE_GAP = 1.86, 1.02, 0.13
ENGINE_Y, ENGINE_H = LANE_Y0 + 4 * LANE_H + 3 * LANE_GAP + 0.20, 0.72
CHIP_H, CHIP_PT = 0.27, 9


def col_x(j):
    return X0 + j * (COL_W + GAP)


def lane_y(i):
    return LANE_Y0 + i * (LANE_H + LANE_GAP)


prs = Presentation()
prs.slide_width, prs.slide_height = Inches(W), Inches(H)
slide = prs.slides.add_slide(prs.slide_layouts[6])
NAMES = {}


def rgb(h):
    return RGBColor.from_string(h)


def _plain(s):
    st = s._element.find(qn("p:style"))
    if st is not None:
        s._element.remove(st)


def box(kind, x, y, w, h, fill=None, line=None, line_pt=0.75, dash=False, adj=None,
        name=None):
    s = slide.shapes.add_shape(kind, Inches(x), Inches(y), Inches(w), Inches(h))
    if fill:
        s.fill.solid()
        s.fill.fore_color.rgb = rgb(fill)
    else:
        s.fill.background()
    if line:
        s.line.color.rgb = rgb(line)
        s.line.width = Pt(line_pt)
        if dash:
            s.line.dash_style = 4  # dash
    else:
        s.line.fill.background()
    if adj is not None:
        s.adjustments[0] = adj
    _plain(s)
    if name:
        s.name = name
    tf = s.text_frame
    tf.margin_left = tf.margin_right = Inches(0.04)
    tf.margin_top = tf.margin_bottom = Inches(0.02)
    tf.word_wrap = True
    return s


def text(shape_or_frame, lines, size, color=INK, bold=False, font=SANS, align=PP_ALIGN.CENTER,
         anchor=MSO_ANCHOR.MIDDLE, italic=False, spacing=None):
    tf = shape_or_frame.text_frame if hasattr(shape_or_frame, "text_frame") else shape_or_frame
    tf.vertical_anchor = anchor
    if isinstance(lines, str):
        lines = lines.split("\n")
    for k, ln in enumerate(lines):
        p = tf.paragraphs[0] if k == 0 else tf.add_paragraph()
        p.alignment = align
        if spacing:
            p.line_spacing = spacing
        runs = ln if isinstance(ln, list) else [(ln, {})]
        for t, o in runs:
            r = p.add_run()
            r.text = t
            f = r.font
            f.name = o.get("font", font)
            f.size = Pt(o.get("size", size))
            f.bold = o.get("bold", bold)
            f.italic = o.get("italic", italic)
            f.color.rgb = rgb(o.get("color", color))


def textbox(x, y, w, h, lines, size, **kw):
    tb = slide.shapes.add_textbox(Inches(x), Inches(y), Inches(w), Inches(h))
    tf = tb.text_frame
    tf.word_wrap = True
    tf.margin_left = tf.margin_right = tf.margin_top = tf.margin_bottom = 0
    text(tb, lines, size, **kw)
    return tb


def arrow(x1, y1, x2, y2, color=LINE, pt=1.0):
    c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x1), Inches(y1),
                                   Inches(x2), Inches(y2))
    _plain(c)
    c.line.color.rgb = rgb(color)
    c.line.width = Pt(pt)
    ln = c.line._get_or_add_ln()
    tail = etree.SubElement(ln, qn("a:tailEnd"))
    tail.set("type", "triangle")
    tail.set("w", "med")
    tail.set("len", "med")
    return c


def chip(x, y, w, cmd, color):
    s = box(MSO_SHAPE.ROUNDED_RECTANGLE, x, y, w, CHIP_H, fill="FFFFFF", line=color,
            line_pt=1.0, adj=0.35)
    s.text_frame.margin_left = s.text_frame.margin_right = Inches(0)
    s.text_frame.word_wrap = False
    text(s, cmd, CHIP_PT, color=INK, bold=True, font=MONO)
    return s


def chip_w(cmd):
    return min(COL_W, 0.0752 * len(cmd) + 0.14)


# --- title and legend -----------------------------------------------------------
textbox(M, 0.34, 7.5, 0.5, "GCEN2 workflow", 26, bold=True, align=PP_ALIGN.LEFT,
        anchor=MSO_ANCHOR.TOP)
textbox(M, 0.84, 8.6, 0.3,
        "Gene co-expression network analysis from bulk to single-cell and spatial transcriptomes",
        13, color=MUTED, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP)

legend = [("data", "data", 0.36), ("program", "gcen2 command", 1.20),
          ("method", "method", 0.58), ("step", "step", 0.36)]
ly = 0.48
lx = W - M - sum(0.42 + lw for _, _, lw in legend) - 0.25 * (len(legend) - 1)
for kind, label, lw in legend:
    if kind == "data":
        s = box(MSO_SHAPE.SNIP_2_DIAG_RECTANGLE, lx, ly, 0.34, 0.22, fill="FFFFFF", line=INK,
                line_pt=0.75, adj=0.3)
    elif kind == "program":
        s = box(MSO_SHAPE.ROUNDED_RECTANGLE, lx, ly, 0.34, 0.22, fill="FFFFFF", line=INK,
                line_pt=1.0, adj=0.35)
    elif kind == "method":
        s = box(MSO_SHAPE.ROUNDED_RECTANGLE, lx, ly, 0.34, 0.22, fill="F1F3F4", line=MUTED,
                line_pt=0.75, dash=True, adj=0.2)
    else:
        s = box(MSO_SHAPE.CHEVRON, lx, ly, 0.34, 0.22, fill=HEAD, adj=0.35)
    textbox(lx + 0.42, ly - 0.01, lw, 0.24, label, 10, color=MUTED, align=PP_ALIGN.LEFT,
            font=MONO if kind == "program" else SANS)
    lx += 0.42 + lw + 0.25

# --- step headers -----------------------------------------------------------------
for j, step in enumerate(STEPS):
    s = box(MSO_SHAPE.PENTAGON if j == 0 else MSO_SHAPE.CHEVRON, col_x(j) - (0.10 if j else 0),
            HEAD_Y, COL_W + (0.22 if j else 0.12), HEAD_H, fill=HEAD, adj=0.28)
    text(s, step, 11, color=INK, bold=True)
s = box(MSO_SHAPE.CHEVRON, SHARED_X - 0.10, HEAD_Y, SHARED_W + 0.10, HEAD_H, fill=HEAD,
        adj=0.28)
text(s, "Function & output", 11, color=INK, bold=True)

# --- lanes ------------------------------------------------------------------------
for i, (lname, color, tint) in enumerate(LANES):
    y = lane_y(i)
    box(MSO_SHAPE.ROUNDED_RECTANGLE, M, y, LANES_END - M, LANE_H, fill=tint, adj=0.08,
        name=f"lane {i + 1}")
    dot = box(MSO_SHAPE.OVAL, M + 0.12, y + LANE_H / 2 - 0.08, 0.16, 0.16, fill=color)
    textbox(M + 0.34, y + 0.08, LABEL_W - 0.30, LANE_H - 0.16, lname, 11, bold=True,
            color=INK, align=PP_ALIGN.LEFT)
    cy = y + LANE_H / 2
    for j, (kind, payload) in enumerate(CELLS[i]):
        x = col_x(j)
        if kind == "data":
            s = box(MSO_SHAPE.SNIP_2_DIAG_RECTANGLE, x, cy - 0.31, COL_W, 0.62, fill="FFFFFF",
                    line=color, line_pt=1.0, adj=0.18)
            text(s, payload, 10, color=INK)
        elif kind == "method":
            s = box(MSO_SHAPE.ROUNDED_RECTANGLE, x, cy - 0.29, COL_W, 0.58, fill="FFFFFF",
                    adj=0.18)
            s.fill.fore_color.rgb = rgb("FFFFFF")
            s.fill.transparency = 0
            text(s, payload, 10, color=INK)
            s.fill.solid()
            s.fill.fore_color.rgb = rgb(tint)
            # a slightly deeper tint than the lane: the lane colour at low weight
            s.line.color.rgb = rgb(color)
            s.line.width = Pt(0.75)
            s.line.dash_style = 4
        else:
            cmds, caption, row = payload
            n_cap = len(caption.split("\n")) if caption else 0
            cap_h = 0.17 * n_cap
            if row:
                ws = [chip_w(c) for c in cmds]
                total = sum(ws) + 0.06 * (len(cmds) - 1)
                block_h = CHIP_H + (0.05 + cap_h if caption else 0)
                yy = cy - block_h / 2
                xx = x + (COL_W - total) / 2
                for c, w_ in zip(cmds, ws):
                    chip(xx, yy, w_, c, color)
                    xx += w_ + 0.06
                cap_y = yy + CHIP_H + 0.05
            else:
                block_h = len(cmds) * CHIP_H + 0.05 * (len(cmds) - 1) + (0.05 + cap_h if caption else 0)
                yy = cy - block_h / 2
                for c in cmds:
                    w_ = chip_w(c)
                    chip(x + (COL_W - w_) / 2, yy, w_, c, color)
                    yy += CHIP_H + 0.05
                cap_y = yy
            if caption:
                textbox(x, cap_y, COL_W, cap_h, caption, 9.5, color=MUTED)
        if j < 4:
            arrow(x + COL_W + 0.03, cy, x + COL_W + GAP - 0.03, cy)
    arrow(LANES_END + 0.03, cy, SHARED_X - 0.03, cy)

# --- shared column ----------------------------------------------------------------
top_y, bot_end = LANE_Y0, lane_y(3) + LANE_H
mid = top_y + 0.56 * (bot_end - top_y)
for (py, ph, title, body) in [
        (top_y, mid - top_y - 0.07, "Gene function", "top"),
        (mid + 0.07, bot_end - mid - 0.07, "Outputs", "bottom")]:
    p = box(MSO_SHAPE.ROUNDED_RECTANGLE, SHARED_X, py, SHARED_W, ph, fill="FFFFFF",
            line="C9CDD2", line_pt=0.75, adj=0.05)
    textbox(SHARED_X + 0.14, py + 0.10, SHARED_W - 0.28, 0.24, title, 11.5, bold=True,
            align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.TOP)
    yy = py + 0.42
    if body == "top":
        rows = [(["annotate", "enrich"], "GO / KEGG enrichment"),
                (["rwr"], "seed-gene ranking"),
                (["stats", "merge"], None), (["extract", "shuffle"], "network tools")]
    else:
        rows = [(None, "modules · eigengenes ·\nkME · hub genes · traits\n(GCEN v1 table formats)"),
                (["plot"], "43 figure types"),
                (["report"], "HTML / PDF")]
    for cmds, cap in rows:
        if cmds:
            ws = [chip_w(c) for c in cmds]
            xx = SHARED_X + 0.14
            for c, w_ in zip(cmds, ws):
                chip(xx, yy, w_, c, "5F6368")
                xx += w_ + 0.06
            if cap:
                if len(cmds) == 1:
                    textbox(xx + 0.04, yy + 0.02, SHARED_X + SHARED_W - xx - 0.16, CHIP_H, cap,
                            9.5, color=MUTED, align=PP_ALIGN.LEFT)
                    yy += CHIP_H + 0.12
                else:
                    textbox(SHARED_X + 0.14, yy + CHIP_H + 0.03, SHARED_W - 0.28, 0.18, cap,
                            9.5, color=MUTED, align=PP_ALIGN.LEFT)
                    yy += CHIP_H + 0.27
            else:
                yy += CHIP_H + 0.06
        else:
            nl = len(cap.split("\n"))
            textbox(SHARED_X + 0.14, yy, SHARED_W - 0.28, 0.18 * nl, cap, 10, color=INK,
                    align=PP_ALIGN.LEFT)
            yy += 0.18 * nl + 0.12

# --- engine -----------------------------------------------------------------------
e = box(MSO_SHAPE.ROUNDED_RECTANGLE, M, ENGINE_Y, W - 2 * M, ENGINE_H, fill=ENGINE, adj=0.12,
        name="engine")
textbox(M + 0.18, ENGINE_Y, 1.0, ENGINE_H, "Engine", 11.5, bold=True, align=PP_ALIGN.LEFT)
ex0 = M + 1.05
EW = [1.0, 1.45, 1.05, 1.15, 1.1]
EW = [v * (W - M - 0.15 - ex0) / sum(EW) for v in EW]
for k, (key, detail) in enumerate(ENGINE_ITEMS):
    x = ex0 + sum(EW[:k])
    ew = EW[k]
    if k:
        c = slide.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, Inches(x - 0.10),
                                       Inches(ENGINE_Y + 0.14), Inches(x - 0.10),
                                       Inches(ENGINE_Y + ENGINE_H - 0.14))
        _plain(c)
        c.line.color.rgb = rgb("C9CDD2")
        c.line.width = Pt(0.75)
    textbox(x, ENGINE_Y + 0.10, ew - 0.22, ENGINE_H - 0.20,
            [[(key, {"bold": True, "color": INK, "size": 10.5})]] +
            [[(ln, {"color": MUTED, "size": 9.5})] for ln in detail.split("\n")],
            10, align=PP_ALIGN.LEFT, anchor=MSO_ANCHOR.MIDDLE)

prs.save(OUT)
print(OUT)
