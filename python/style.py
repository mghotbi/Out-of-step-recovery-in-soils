"""Shared Nature-style settings: Arial (falls back to the metrically identical Liberation Sans),
7-pt text everywhere, 9-pt bold lower-case panel letters, 183-mm double-column width."""
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import font_manager
import re, os

FS = 7          # all figure text (pt)
FS_LETTER = 9   # panel letters (pt)
W_DOUBLE = 183 / 25.4

_avail = {f.name for f in font_manager.fontManager.ttflist}
FONT = "Arial" if "Arial" in _avail else "Liberation Sans"
plt.rcParams.update({
    "font.family": "sans-serif", "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica"],
    "font.size": FS, "axes.titlesize": FS, "axes.labelsize": FS, "xtick.labelsize": FS,
    "ytick.labelsize": FS, "legend.fontsize": FS, "legend.title_fontsize": FS,
    "mathtext.fontset": "custom", "mathtext.rm": FONT, "mathtext.it": FONT + ":italic",
    "mathtext.bf": FONT + ":bold", "mathtext.sf": FONT, "mathtext.default": "regular",
    "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
    "xtick.minor.width": 0.4, "ytick.minor.width": 0.4,
    "xtick.major.size": 2.5, "ytick.major.size": 2.5, "xtick.minor.size": 1.5, "ytick.minor.size": 1.5,
    "xtick.major.pad": 2, "ytick.major.pad": 2, "axes.labelpad": 2.5,
    "axes.spines.top": False, "axes.spines.right": False, "legend.frameon": False,
    "legend.handlelength": 1.6, "legend.handletextpad": 0.5, "legend.borderaxespad": 0.2,
    "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none", "savefig.dpi": 600,
    "lines.linewidth": 1.2,
})

COL = {"Capacity": "#D08C2E", "Connectivity": "#2E8F95", "Kinetics": "#C2463B", "Memory": "#7E5CB0",
       "Debt": "#3F8A4E", "Bulk": "#8C8C8C", "ink": "#222222"}
HIST = {"Ox": "#1F78B4", "SA": "#E6A100", "LA": "#B2182B"}

def letter(ax, s, x=-0.02, y=1.0, fig=None):
    """Bold panel letter placed in figure coordinates at the axes' upper-left corner."""
    fig = fig or ax.figure
    bb = ax.get_position()
    fig.text(bb.x0 + x, bb.y1 + y * 0, s, fontsize=FS_LETTER, fontweight="bold",
             ha="right", va="bottom", transform=fig.transFigure)

def save(fig, outdir, name):
    os.makedirs(outdir, exist_ok=True)
    for ext in ("pdf", "png", "svg"):
        p = os.path.join(outdir, f"{name}.{ext}")
        fig.savefig(p, dpi=600)
        if ext == "svg":   # declare Arial first so the editable SVG opens in Arial
            s = open(p, encoding="utf-8").read()
            s = re.sub(r"font-family:\s*'?[^;'\"]*'?;", "font-family: Arial, 'Liberation Sans', Helvetica, sans-serif;", s)
            s = re.sub(r'font-family="[^"]*"', 'font-family="Arial, Liberation Sans, Helvetica, sans-serif"', s)
            open(p, "w", encoding="utf-8").write(s)
    plt.close(fig)
