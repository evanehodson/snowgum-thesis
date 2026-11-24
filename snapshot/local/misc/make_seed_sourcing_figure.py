"""
Six-stratified seed-sourcing schematic.

Flat, 2D, native-MS-Office look: solid fills, smooth/striped gradients,
hairline strokes, no shadows / textures / transparency. All geometry is
deterministic.

The strips are colour-coded, because they are three different spaces:

    rows 1-4   climate space              blue
    row  5     Genome-wide Turnover       green
    row  6     Adaptive Allele Turnover   violet

Run:  python make_seed_sourcing_figure.py
Out:  seed_sourcing_strategies.png  (300 dpi)
"""

import numpy as np
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Circle, Polygon, Rectangle

# ---------------------------------------------------------------- settings --

SEED = 20260928
OUT = "seed_sourcing_strategies"

plt.rcParams.update(
    {
        "font.family": ["Segoe UI", "Calibri", "Arial", "DejaVu Sans"],
        "font.size": 7,
        "pdf.fonttype": 42,
        "svg.fonttype": "none",
    }
)

W, H = 10.00, 12.90                   # inches == data units, so 1 unit == 1 in
X0, X1 = 1.58, 9.45                   # gradient strip extent
SW = X1 - X0
STAR_FRAC = 0.196                     # the star keeps its place along a wider bar
STAR_X = X0 + STAR_FRAC * SW          # identical on every strip
DOT_SCALE = SW / 4.85                 # keep dot density as the strip widens
LABEL_X, LABEL_W = 0.12, 1.30         # narrow left gutter for the row labels
R_OUT, R_IN = 0.125, 0.125 * 0.382    # five-point star
DOT_R = 0.036                         # donor dot radius
CLIMATE_COLOUR = 0.62                 # 62% along the climate axis

INK = "#000000"
GRID = "#9E9E9E"
LABEL_BG, LABEL_ED = "#F4F5F7", "#C9CDD2"
STAR_FILL = "#FFC000"
MUTED, BODY = "#5A5F66", "#40444A"

# Four hand-specified Office ramps, one per quantity. All run dark -> white,
# so the strips share a shape and differ only in hue.
CLIMATE_BLUE = [
    "#0B1F51", "#123E8C", "#1E4FA8", "#2E6BC4",
    "#6C9BDE", "#B6CFEE", "#E6EEF9", "#FFFFFF",
]
GDM_GREEN = [
    "#08361F", "#0E5734", "#167A4A", "#259C60",
    "#57BB85", "#93D3AE", "#C8E7D5", "#FFFFFF",
]
LFMM_VIOLET = [
    "#2B0B4E", "#431A73", "#5E2B9B", "#7C41BE",
    "#9C6AD4", "#BFA2E4", "#DED0F2", "#FFFFFF",
]
CMAP_CLIMATE = LinearSegmentedColormap.from_list("office_blue", CLIMATE_BLUE, N=1024)
CMAP_GDM = LinearSegmentedColormap.from_list("office_green", GDM_GREEN, N=1024)
CMAP_LFMM = LinearSegmentedColormap.from_list("office_violet", LFMM_VIOLET, N=1024)

# Ink for the coloured annotations, so text and strip share a hue.
C_CLIMATE, C_GDM, C_LFMM = "#1E4FA8", "#167A4A", "#5E2B9B"

# What each strip is *not*: rows 1-4 share one climate space, rows 5 and 6
# are transformed spaces, so each gets its own name.
SPACE = {
    0: ("climate space  \u2192", C_CLIMATE),
    1: ("climate space  \u2192", C_CLIMATE),
    2: ("climate space  \u2192", C_CLIMATE),
    3: ("climate space  \u2192", C_CLIMATE),
    4: ("Neutral Genome-Wide Turnover  \u2192", C_GDM),
    5: ("Adaptive Polygenic Turnover  \u2192", C_LFMM),
}

# ------------------------------------------------------------- row layout ---

# (row_label, sublabel, y_top)
ROW_H = 0.75
ROWS = [
    ("Local", None, 0.35),
    ("Composite", None, 1.42),
    ("Admixture", None, 2.44),
    ("Climate-Adjusted", None, 4.10),
    ("GDM", None, 7.85),
    ("LFMM", None, 11.65),
]
ROW = [(a, b, t, t + ROW_H) for (a, b, t) in ROWS]

# --- GDM matrix, as a two-sheet paper stack ---------------------------------
# formatsitepair(bioFormat = 3) builds one dissimilarity matrix with the
# genomic block on top of the environmental block; here the two blocks are the
# two sheets, the environmental one behind, the genomic one in front and
# offset down-right so the depth reads.
N_SAMP, DM_CELL, DM_DEPTH = 5, 0.21, 0.30
DM_X0, DM_Y0 = 2.40, 5.50             # back sheet, top-left
DM_S = N_SAMP * DM_CELL              # sheet side
DM_FX, DM_FY = DM_X0 + DM_DEPTH, DM_Y0 + DM_DEPTH   # front sheet, top-left
DM_CAP_Y = DM_Y0 - 0.06               # "pairwise dissimilarity" caption
DM_BAR_X = DM_X0 - 0.06               # "environmental", on the left edge
# The blocks are abstract: no sample IDs are labelled.

# --- GDM turnover-curve inset ----------------------------------------------
# Both inset graphs share this one rectangle, so they always sit in the same
# place at the same size.
# Panel right edge stops short of the gradient's, so the inset never reads as a
# continuation of the strip: the two dotted axis_map connectors stay slanted
# and the shared X1 can shorten without moving the stack.
PANEL_CLEAR = 0.40                 # gap between panel right edge and X1
PANEL_W = 0.381 * SW               # inset widens with the strip, height does not
PANEL_X1 = X1 - PANEL_CLEAR
PANEL_X0 = PANEL_X1 - PANEL_W
PANEL_H = 1.75
GX0, GX1 = PANEL_X0, PANEL_X1
GY0, GY1 = 5.60, 5.60 + PANEL_H

# --- LFMM row --------------------------------------------------------------
# Five cells, and the brace spans exactly those five and nothing more.
BETA = np.array([0.07, -0.33, 0.19, 0.79, -0.21])
CELL_W, CELL_H = 0.46, 0.22
BX0, BY0 = 2.90, 9.415                          # beta column
BETA_H = len(BETA) * CELL_H                     # brace height == the 5 cells
BRACE_X, BRACE_W = BX0 + CELL_W + 0.08, 0.15    # curly brace feeding the graph
LX0, LX1 = PANEL_X0, PANEL_X1
LY0, LY1 = 9.20, 9.20 + PANEL_H

# ------------------------------------------------- GDM sigmoid (space B) ----
SIG_K = 10.0


def gd_spline(u):
    """Monotone logistic GDM transform. g(0) = 0, g(1) = 1, one inflection."""
    u = np.asarray(u, dtype=float)
    s0 = 1.0 / (1.0 + np.exp(SIG_K * 0.5))          # raw logistic at u = 0
    g = 1.0 / (1.0 + np.exp(-SIG_K * (u - 0.5)))
    return (g - s0) / (1.0 - s0)


_UG = np.linspace(0.0, 1.0, 200001)
_GG = gd_spline(_UG)


def gd_inverse(c):
    """Lower inverse: the smallest d at which the model has consumed c of its
    colour axis. Returns 0 for c <= 0 rather than the right edge of the
    flat-zero plateau, which is what a bare np.interp would give.
    """
    idx = np.searchsorted(_GG, np.clip(c, 0.0, 1.0), side="left")
    return _UG[np.minimum(idx, _UG.size - 1)]


D5 = float(gd_inverse(CLIMATE_COLOUR))
TIP4 = X0 + CLIMATE_COLOUR * SW
TIP5 = X0 + D5 * SW

# ------------------------------------------ LFMM adaptive-space strip (C) --
#
# LFMM defines a continuous quadratic distance (Mahalanobis-like metric) using
# the inner product of the effect matrix (B^T B). Unlike the GDM turnover
# curve, the response is LINEAR in the climate difference, so the strip is a
# straight ramp and the only thing the SNP column sets is its slope.

STRETCH = 1.9                                  # turnover reached 1/1.9 sooner


def lfmm_transform(u):
    return np.clip(np.asarray(u, dtype=float), 0.0, 1.0)


TIP6 = X0 + (CLIMATE_COLOUR / STRETCH) * SW

# d(d) = sqrt(1/L sum (b_i d)^2) factors as d * m, so the whole SNP column
# collapses to one coefficient m on the climate axis -- the slope of the
# straight line in row 6.
SLOPE = float(np.sqrt(np.mean(BETA ** 2)))

# ------------------------------------------------------------- primitives --


def star(ax, cx, cy, z=9):
    pts = [
        (cx + (R_OUT if i % 2 == 0 else R_IN) * np.cos(-np.pi / 2 + i * np.pi / 5),
         cy + (R_OUT if i % 2 == 0 else R_IN) * np.sin(-np.pi / 2 + i * np.pi / 5))
        for i in range(10)
    ]
    ax.add_patch(Polygon(pts, closed=True, fc=STAR_FILL, ec="#FFFFFF", lw=0.7, zorder=z))


def dot(ax, x, y, z=8):
    ax.add_patch(Circle((x, y), DOT_R, fc=INK, ec="none", zorder=z))


def arrow(ax, x0, y0, x1, y1, lw=0.6, head=5.0, z=7):
    ax.annotate(
        "", xy=(x1, y1), xytext=(x0, y0), zorder=z,
        arrowprops=dict(arrowstyle="-|>", color=INK, lw=lw, mutation_scale=head,
                        shrinkA=0, shrinkB=0),
    )


def gradient(ax, y0, y1, colour_of_u, cmap=CMAP_CLIMATE, n=2048):
    """Continuous strip drawn as n rectangles, so the raster stays smooth at
    300 dpi. Each rect is widened by a hair to close antialiasing seams; the
    colour step across one rect is 1/n of the ramp, so the overlap is
    invisible.
    """
    u = (np.arange(n) + 0.5) / n
    rgb = cmap(colour_of_u(u))[:, :3]
    w = SW / n
    bleed = w * 0.06
    for i in range(n):
        ax.add_patch(
            Rectangle((X0 + i * w, y0), w + bleed, y1 - y0, fc=rgb[i],
                      ec="none", antialiased=False, zorder=2)
        )


def border(ax, y0, y1):
    ax.add_patch(Rectangle((X0, y0), SW, y1 - y0, fc="none", ec=GRID,
                           lw=0.6, zorder=10))


def label_box(ax, title, sub, cy):
    """Box is centred on cy, i.e. on the vertical middle of its strip, and is
    the same height on every row so the left gutter reads as one column.
    """
    ax.add_patch(Rectangle((LABEL_X, cy - 0.35), LABEL_W, 0.70, fc=LABEL_BG,
                           ec=LABEL_ED, lw=0.6, zorder=3))
    if sub:
        ax.text(LABEL_X + 0.09, cy - 0.185, title, ha="left", va="center", fontsize=10.0,
                weight="bold", color=INK, zorder=4)
        ax.text(LABEL_X + 0.09, cy + 0.16, sub, ha="left", va="center", fontsize=6.5,
                color=MUTED, linespacing=1.25, zorder=4)
    else:
        ax.text(LABEL_X + 0.09, cy, title, ha="left", va="center", fontsize=10.0,
                weight="bold", color=INK, linespacing=1.3, zorder=4)


def panel(ax, x0, y0, x1, y1, xtitle, ytitle, caption=None):
    """White inset with the flat tick style, shared by the GDM turnover curve
    and the LFMM slope so the two read as the same kind of object.
    """
    ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fc="white", ec=GRID, lw=0.6, zorder=3))
    for t in (0.0, 0.5, 1.0):
        ax.plot([x0, x0 - 0.05], [y1 - t * (y1 - y0)] * 2, color=INK, lw=0.5, zorder=4)
        ax.plot([x0 + t * (x1 - x0)] * 2, [y1, y1 + 0.05], color=INK, lw=0.5, zorder=4)
    ax.text(x0 - 0.07, y1, "0", ha="right", va="center", fontsize=9.5, color=MUTED)
    ax.text(x0 - 0.07, y0, "1", ha="right", va="center", fontsize=9.5, color=MUTED)
    ax.text((x0 + x1) / 2, y1 + 0.09, xtitle, ha="center", va="top",
            fontsize=11.0, color=BODY)
    ax.text(x0 - 0.21, (y0 + y1) / 2, ytitle, ha="center", va="center",
            fontsize=11.0, color=BODY, rotation=90)
    if caption:
        ax.text(x0, y0 - 0.06, caption, ha="left", va="bottom", fontsize=9.5,
                color=BODY, linespacing=1.15)


def axis_map(ax, x0, x1, y_from, y_to):
    """Dotted pair showing how a panel's horizontal axis maps onto the strip's.

    The y axis is inverted, so a panel's y1 is its visual bottom and the x axis
    is drawn there. Pass y_from as that same y1 and the connectors leave from
    the two x-axis ends exactly, rather than floating clear of them.
    """
    for gx, sx in ((x0, X0), (x1, X1)):
        ax.plot([gx, sx], [y_from, y_to], color=BODY, lw=0.6, ls=(0, (1.2, 1.8)), zorder=5)


def tint(colour, amount):
    """Blend a hex colour toward white. 0 = unchanged, 1 = white."""
    r, g, b = (int(colour[i:i + 2], 16) for i in (1, 3, 5))
    mix = lambda v: int(round(v + (255 - v) * amount))
    return f"#{mix(r):02X}{mix(g):02X}{mix(b):02X}"


# ---- fitted-curve insets ----------------------------------------------------
RING_N = 14                       # dots clustered at each strip tip
FIT_CLOSE_N = RING_N              # unfainted points match the ring count
FIT_FAINT_N = round(30 * PANEL_W / 1.85)  # keeps cloud density as the inset widens
PANEL_T = 0.62                    # marked point, same fraction on both curves
PANEL_PT_R = 0.028
# Two tiers, binary: unfainted (the curve's own colour) or one faint tint.
FIT_FAINT_R, FIT_CLOSE_R = 0.014, 0.020
FIT_FAINT_TINT = 0.86
FIT_SIG_FAINT = 0.16              # faint spread, as a fraction of panel height
FIT_CLEAR = 0.055                 # keep a point from hiding under the marked dot
FIT_LABEL_DY = 0.22               # how far above the inset the drop line runs
FIT_SEP = 0.055                   # separation within the pool
FIT_MARGIN = 0.06                 # keep points off the panel border
FIT_Z_FAINT, FIT_Z_CLOSE = 4, 5   # both under the fit, which stays on top


def lf_curve(t):
    """The LFMM response is linear, so the fitted curve is the identity on the
    panel's normalised axes."""
    return np.asarray(t, dtype=float)


def curve_point(ax, x0, y0, x1, y1, f, colour, t=PANEL_T, label=None):
    """The marked point on a fitted curve, with its dashed drop line to the
    horizontal axis. Identical treatment on both insets.

    The y axis is inverted, so y1 is the inset's visual bottom and y0 its top.
    Given a label, the same dashed line carries on past the point and out of
    the top of the inset to that label, so the drop line and its caption read
    as one continuous mark.
    """
    lw_, lh_ = x1 - x0, y1 - y0
    px = x0 + t * lw_
    py = y1 - float(f(t)) * lh_
    y_end = y0 - FIT_LABEL_DY if label else py
    ax.plot([px, px], [y_end, y1], color=MUTED, lw=0.6, ls=(0, (1.2, 1.6)), zorder=6)
    if label:
        ax.text(px, y_end - 0.08, label, ha="center", va="bottom",
            fontsize=11.0, color=BODY)
    ax.add_patch(Circle((px, py), PANEL_PT_R, fc=colour, ec="none", zorder=7))
    return px, py


def curve_scatter(ax, x0, y0, x1, y1, f, colour, rng, px, py):
    """Scatter about a fitted curve, so the inset reads as a fit to observed
    data rather than a bare drawn line.

    One uniformly random pool, with no designed cluster around the marked point:
    every point is placed exactly the same way, scattered perpendicular to the
    curve, and the tiers are only decided afterwards by sorting the finished
    pool on distance to the marked point. The closest FIT_CLOSE_N come out
    unfainted in the curve's own colour, the remainder take a single faint tint.

    So the two tiers are binary, and which tier a point lands in depends purely
    on how near the marked point it happened to fall, not on any separate rule
    for how the unfainted set is generated.
    """
    lw_, lh_ = x1 - x0, y1 - y0
    n_total = FIT_CLOSE_N + FIT_FAINT_N
    pts, tries = [], 0
    while len(pts) < n_total and tries < 60000:
        tries += 1
        t = rng.uniform(0.02, 0.98)
        g0 = float(f(t))
        dt = 1e-3
        slope = (float(f(t + dt)) - g0) / dt          # dy_norm / dx_norm
        tx, ty = lw_, -slope * lh_                    # tangent in figure units
        norm = np.hypot(tx, ty) or 1.0
        off = rng.normal(0.0, FIT_SIG_FAINT) * lh_
        qx = x0 + t * lw_ - ty / norm * off
        qy = y1 - g0 * lh_ + tx / norm * off
        if not (x0 + FIT_MARGIN < qx < x1 - FIT_MARGIN
                and y0 + FIT_MARGIN < qy < y1 - FIT_MARGIN):
            continue
        if np.hypot(qx - px, qy - py) < FIT_CLEAR:    # would hide under the marked dot
            continue
        if not _clear(qx, qy, pts, FIT_SEP):
            continue
        pts.append((qx, qy))
    if len(pts) < n_total:
        raise RuntimeError(f"curve_scatter placed {len(pts)}/{n_total} points")

    pts.sort(key=lambda q: np.hypot(q[0] - px, q[1] - py))
    close_pts, faint_pts = pts[:FIT_CLOSE_N], pts[FIT_CLOSE_N:]
    fc_faint = tint(colour, FIT_FAINT_TINT)
    for qx, qy in faint_pts:
        ax.add_patch(Circle((qx, qy), FIT_FAINT_R, fc=fc_faint, ec="none",
                            zorder=FIT_Z_FAINT))
    for qx, qy in close_pts:
        ax.add_patch(Circle((qx, qy), FIT_CLOSE_R, fc=colour, ec="none",
                            zorder=FIT_Z_CLOSE))
    return faint_pts, close_pts


# ---- jittered donor dots ---------------------------------------------------
MIN_SEP = 0.115
MAX_TRIES = 60000

def _clear(x, y, pts, min_sep):
    return all(np.hypot(x - px, y - py) >= min_sep for px, py in pts)

def _in_strip(x, y, y0, y1, x_lo, x_hi):
    return x_lo < x < x_hi and y0 + 0.09 < y < y1 - 0.09

def ring_jitter(ax, cx, cy, n, r_in, r_out, rng, y0, y1, min_sep=MIN_SEP):
    """n dots, area-uniform in an annulus around (cx, cy), min-separation
    enforced. Deliberately NOT polar-spaced, so the cloud has no radial spokes.
    """
    pts, tries = [], 0
    while len(pts) < n and tries < MAX_TRIES:
        tries += 1
        a = rng.uniform(0, 2 * np.pi)
        r = np.sqrt(rng.uniform(r_in ** 2, r_out ** 2))
        x, y = cx + r * np.cos(a), cy + r * np.sin(a)
        if _in_strip(x, y, y0, y1, X0 + 0.09, X1 - 0.09) and _clear(x, y, pts, min_sep):
            pts.append((x, y))
    if len(pts) < n:
        raise RuntimeError(f"ring_jitter placed {len(pts)}/{n} dots at ({cx}, {cy})")
    for x, y in pts:
        dot(ax, x, y)
    return pts

def scatter_dots(ax, n, rng, y0, y1, min_sep=MIN_SEP, exclude=None, seed_pts=()):
    """n dots spread over the whole strip, honouring min_sep, `exclude`, and
    any already-placed points in `seed_pts`.
    """
    pts, tries = list(seed_pts), 0
    while len(pts) < n and tries < MAX_TRIES:
        tries += 1
        x, y = rng.uniform(X0 + 0.09, X1 - 0.09), rng.uniform(y0 + 0.09, y1 - 0.09)
        if exclude and exclude(x, y):
            continue
        if _clear(x, y, pts, min_sep):
            pts.append((x, y))
    if len(pts) < n:
        raise RuntimeError(f"scatter_dots placed {len(pts)}/{n} dots")
    for x, y in pts:
        dot(ax, x, y)
    return pts

def space_label(ax, i, y):
    text, colour = SPACE[i]
    ax.text(X0, y, text, ha="left", va="top", fontsize=10.5, color=colour)

def group_label(ax, y0, y1, text):
    ax.plot([X1 + 0.13, X1 + 0.13], [y0, y1], color=GRID, lw=0.6, zorder=3)
    ax.text(X1 + 0.23, (y0 + y1) / 2, text, rotation=90, ha="center", va="center",
            fontsize=10.5, color=BODY)

def climate_vector(ax, y0, y1, tip_x, label=None):
    cy = (y0 + y1) / 2
    arrow(ax, STAR_X + R_OUT + 0.02, cy, tip_x, cy, lw=1.1, head=9.0, z=9)
    ax.add_patch(Circle((tip_x, cy), 0.055, fc="none", ec=INK, lw=0.8, zorder=9))
    if label:
        ax.text(tip_x, y0 - 0.06, label, ha="center", va="bottom", fontsize=9.5,
                color=BODY)


def curly_brace(ax, x, y0, y1, w=0.15, color=INK, lw=0.9, z=6):
    """A right curly brace `}` spanning y0..y1: the arms start at x, and the
    cusp at mid-height sits at x + w, pointing at whatever it feeds.

    Two mirrored cubic Beziers meeting at a cusp. Each arm leaves its end
    horizontally, sweeps out to the right and back, then arrives at the cusp
    travelling down-and-right while its mirror leaves down-and-left, so the
    two tangents differ sharply and the junction reads as a point.
    """
    ym = 0.5 * (y0 + y1)
    half = 0.5 * (y1 - y0)
    flare, reach, drop = 0.80, 0.75, 0.22      # terminal flare, cusp reach, cusp drop

    def bezier(p0, p1, p2, p3, n=160):
        t = np.linspace(0.0, 1.0, n)[:, None]
        m = 1.0 - t
        return m ** 3 * p0 + 3 * m ** 2 * t * p1 + 3 * m * t ** 2 * p2 + t ** 3 * p3

    top = bezier((x, y0), (x + flare * w, y0),
                 (x + (1.0 - reach) * w, ym - drop * half), (x + w, ym))
    bot = bezier((x + w, ym), (x + (1.0 - reach) * w, ym + drop * half),
                 (x + flare * w, y1), (x, y1))
    path = np.vstack([top, bot[1:]])
    ax.plot(path[:, 0], path[:, 1], color=color, lw=lw, zorder=z,
            solid_capstyle="round")


def sheet(ax, x0, y0, vals, cmap, ec, z=3):
    """One n x n dissimilarity block, as a sheet in the paper stack.

    Cells and the sheet's own outline share a depth band, with the outline just
    above the cells, so the front sheet's border can never paint over the back
    sheet.
    """
    for i in range(N_SAMP):
        for j in range(N_SAMP):
            ax.add_patch(Rectangle((x0 + j * DM_CELL, y0 + i * DM_CELL), DM_CELL,
                                   DM_CELL, fc=cmap(0.98 - 0.72 * vals[i, j]),
                                   ec="white", lw=0.8, zorder=z))
    ax.add_patch(Rectangle((x0, y0), DM_S, DM_S, fc="none", ec=ec, lw=0.8,
                           zorder=z + 0.5))


# ------------------------------------------------------------------ build --
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, W)
ax.set_ylim(H, 0)
ax.set_aspect("auto")
ax.axis("off")
rng = np.random.default_rng(SEED)
# Separate stream for the inset scatters, so adding points to a panel never
# reshuffles the strip dots that the rest of the figure depends on.
rng_fit = np.random.default_rng(SEED + 101)

# --- Row 1: local -----------------------------------------------------------
_, _, y0, y1 = ROW[0]
cy = (y0 + y1) / 2
label_box(ax, *ROW[0][0:2], cy)
gradient(ax, y0, y1, lambda u: u)
border(ax, y0, y1)
ring_jitter(ax, STAR_X, cy, 8, 0.16, 0.32, rng, y0, y1)
star(ax, STAR_X, cy)
space_label(ax, 0, y1 + 0.05)

# --- Row 2: composite: heavy bias toward the restoration site ---------------
_, _, y0, y1 = ROW[1]
cy = (y0 + y1) / 2
label_box(ax, *ROW[1][0:2], cy)
gradient(ax, y0, y1, lambda u: u)
border(ax, y0, y1)
_near = ring_jitter(ax, STAR_X, cy, RING_N, 0.16, 0.42, rng, y0, y1, min_sep=0.105)
scatter_dots(ax, round((14 + 6) * DOT_SCALE), rng, y0, y1, min_sep=0.12,
             exclude=lambda x, y: np.hypot(x - STAR_X, y - cy) < 0.50,
             seed_pts=_near)
star(ax, STAR_X, cy)
space_label(ax, 1, y1 + 0.05)

# --- Row 3: admixture: broad even scatter, almost nothing near the star -----
_, _, y0, y1 = ROW[2]
cy = (y0 + y1) / 2
label_box(ax, *ROW[2][0:2], cy)
gradient(ax, y0, y1, lambda u: u)
border(ax, y0, y1)
scatter_dots(ax, round(30 * DOT_SCALE), rng, y0, y1, min_sep=0.16)
star(ax, STAR_X, cy)
space_label(ax, 2, y1 + 0.05)

# --- Row 4: climate-adjusted (baseline) -------------------------------------
_, _, y0, y1 = ROW[3]
cy = (y0 + y1) / 2
label_box(ax, *ROW[3][0:2], cy)
gradient(ax, y0, y1, lambda u: u)
border(ax, y0, y1)
climate_vector(ax, y0, y1, TIP4)
ring_jitter(ax, TIP4, cy, RING_N, 0.15, 0.30, rng, y0, y1)
star(ax, STAR_X, cy)
space_label(ax, 3, y1 + 0.05)
# The predicted-future point gets the same treatment as the GDM curve's: an
# unboxed label above, tied to its target by one dotted leader.  The leader
# runs all the way down to the circle rather than stopping at the strip's edge,
# so the label reads as naming the point instead of the bar.
_by = y0 - FIT_LABEL_DY
ax.plot([TIP4, TIP4], [_by, cy], color=MUTED, lw=0.6, ls=(0, (1.2, 1.6)), zorder=6)
ax.text(TIP4, _by - 0.08, "Predicted future climate", ha="center", va="bottom",
        fontsize=10.5, color=BODY)

# --- Row 5: GDM (dissimilarity stack -> turnover curve -> strip) -------------
_, _, y0, y1 = ROW[4]
cy = (y0 + y1) / 2
label_box(ax, *ROW[4][0:2], cy)
gradient(ax, y0, y1, gd_spline, cmap=CMAP_GDM)    # colour = g(d)
border(ax, y0, y1)
space_label(ax, 4, y1 + 0.05)

# Two blocks of the same dissimilarity matrix, drawn as a paper stack: the
# environmental block is genuinely climate-space data, so it stays blue; the
# genomic block on top is what GDM actually models, so it is green.
rng_dm = np.random.default_rng(SEED + 5)
_pos = np.sort(rng_dm.uniform(0.0, 1.0, N_SAMP))
D_GEN = np.abs(_pos[:, None] - _pos[None, :]) ** 0.6
D_ENV = np.abs(_pos[:, None] - _pos[None, :]) ** 0.42

# Both sheets share one colour, and sit in separate depth bands so the back
# sheet's outline stays behind the front sheet.
sheet(ax, DM_X0, DM_Y0, D_ENV, CMAP_GDM, C_GDM, z=2)             # behind
sheet(ax, DM_FX, DM_FY, D_GEN, CMAP_GDM, C_GDM, z=5)             # in front

ax.text(DM_X0 + DM_S / 2, DM_CAP_Y, "pairwise dissimilarity",
        ha="center", va="bottom", fontsize=12.0, color=BODY)
ax.text(DM_BAR_X, DM_Y0 + DM_S / 2, "environmental", ha="right",
        va="center", fontsize=12.0, color=BODY, rotation=90)
ax.text(DM_FX + DM_S / 2, DM_FY + DM_S + 0.02, "genomic", ha="center",
        va="top", fontsize=12.0, color=BODY, zorder=5)

# Turnover curve inset, to the right of the stack.
panel(ax, GX0, GY0, GX1, GY1, "environmental axis", "genomic distance")
_t = np.linspace(0, 1, 400)
ax.plot(GX0 + _t * (GX1 - GX0), GY1 - gd_spline(_t) * (GY1 - GY0),
        color=C_GDM, lw=1.2, solid_capstyle="round", zorder=6)
_gd_pt = curve_point(ax, GX0, GY0, GX1, GY1, gd_spline, C_GDM,
                    label="Predicted future climate")
curve_scatter(ax, GX0, GY0, GX1, GY1, gd_spline, C_GDM, rng_fit, *_gd_pt)

# Stack -> curve arrow, then curve -> strip connectors.
arrow(ax, DM_FX + DM_S + 0.10, DM_FY + DM_S / 2, GX0 - 0.42, DM_FY + DM_S / 2,
      lw=0.9, head=7.0, z=6)
axis_map(ax, GX0, GX1, GY1, y0)
climate_vector(ax, y0, y1, TIP5)
ring_jitter(ax, TIP5, cy, RING_N, 0.15, 0.30, rng, y0, y1)
star(ax, STAR_X, cy)

# --- Row 6: LFMM (one column of B -> brace -> slope -> straight strip) ------
_, _, y0, y1 = ROW[5]
cy = (y0 + y1) / 2
label_box(ax, *ROW[5][0:2], cy)
gradient(ax, y0, y1, lfmm_transform, cmap=CMAP_LFMM)
border(ax, y0, y1)
space_label(ax, 5, y1 + 0.05)

ax.text(BX0 + CELL_W / 2, BY0 - 0.06, "SNP effects", ha="center",
        va="bottom", fontsize=12.0, color=BODY)
for i, v in enumerate(BETA):
    ax.add_patch(Rectangle((BX0, BY0 + i * CELL_H), CELL_W, CELL_H,
                           fc=CMAP_LFMM(0.30 + 0.68 * min(1.0, abs(v))), ec="white",
                           lw=0.8, zorder=3))
    ax.text(BX0 + CELL_W / 2, BY0 + (i + 0.5) * CELL_H, f"{v:+.2f}".replace("+", " "),
            ha="center", va="center", fontsize=12.0, color=INK, zorder=4)
    ax.text(BX0 - 0.06, BY0 + (i + 0.5) * CELL_H, f"SNP {i + 1}", ha="right",
            va="center", fontsize=10.5, color=MUTED)
ax.text(BX0 + CELL_W / 2, BY0 + BETA_H + 0.02, "⋮", ha="center",
        va="top", fontsize=13.5, color=MUTED)

# Curly brace over the whole column: every SNP feeds the sum. The sum is the
# coefficient, so the brace points straight at the slope of the graph.
curly_brace(ax, BRACE_X, BY0, BY0 + BETA_H, w=BRACE_W)
arrow(ax, BRACE_X + BRACE_W + 0.05, BY0 + BETA_H / 2, LX0 - 0.42, BY0 + BETA_H / 2,
      lw=0.9, head=7.0, z=6)

# The response is linear, so this is a straight line and the summed SNP
# effects are its slope.
panel(ax, LX0, LY0, LX1, LY1, "environmental axis",
      "polygenic distance",
      "sum of SNP effects  =  slope")
lw_, lh_ = LX1 - LX0, LY1 - LY0
ax.plot([LX0, LX0 + 0.88 * lw_], [LY1, LY1 - 0.88 * lh_], color=C_LFMM, lw=1.2,
        solid_capstyle="round", zorder=6)
ax.plot([LX0 + 0.88 * lw_, LX1], [LY1 - 0.88 * lh_, LY0], color=C_LFMM, lw=0.7,
        ls=(0, (1.6, 1.8)), zorder=6)
_lf_pt = curve_point(ax, LX0, LY0, LX1, LY1, lf_curve, C_LFMM)
curve_scatter(ax, LX0, LY0, LX1, LY1, lf_curve, C_LFMM, rng_fit, *_lf_pt)

axis_map(ax, LX0, LX1, LY1, y0)
climate_vector(ax, y0, y1, TIP6)
ring_jitter(ax, TIP6, cy, RING_N, 0.15, 0.30, rng, y0, y1)
star(ax, STAR_X, cy)

# ------------------------------------------------------------------ output --

# The panel's right edge must never coincide with the gradient's, or the inset
# reads as a continuation of the strip instead of a separate object.
assert PANEL_X1 <= X1 - PANEL_CLEAR + 1e-9, (
    f"panel right edge {PANEL_X1} is within {PANEL_CLEAR} in of the gradient "
    f"right edge {X1}; the two would read as one object"
)
assert GX1 == LX1, "both insets must share one rectangle"

fig.savefig(f"{OUT}.png", dpi=300, facecolor="white")
print(f"wrote {OUT}.png")
print(f"panel right edge {PANEL_X1} vs gradient right edge {X1} "
      f"(gap {X1 - PANEL_X1:.2f} in)")
print(f"climate vector length  row4={TIP4 - STAR_X:.3f}  "
      f"row5={TIP5 - STAR_X:.3f}  row6={TIP6 - STAR_X:.3f} in")
print(f"sigmoid g^(-1)({CLIMATE_COLOUR}) = {D5:.4f}  (uniform would be {CLIMATE_COLOUR})")
print(f"LFMM slope m = sqrt(1/L sum beta^2) = {SLOPE:.4f}")
print(f"matrix stack footprint {DM_S + DM_DEPTH:.2f} x {DM_S + DM_DEPTH:.2f} in")
