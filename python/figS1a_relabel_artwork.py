"""Relabel the authors' artwork for Supplementary Fig. 1a (label changes only; drawing unchanged)."""
import numpy as np, cv2
from PIL import Image, ImageDraw, ImageFont
SRC = "../artwork/SuppFig1a_source_author_artwork.png"; OUT = "../figures/Supplementary_Figure_1a_artwork.png"
F = "/usr/share/fonts/truetype/liberation2/LiberationSans-"
def font(style, size):
    import os
    for base in (F, "/usr/share/fonts/truetype/liberation/LiberationSans-", "/Library/Fonts/Arial "):
        p = base + style + ".ttf"
        if os.path.exists(p): return ImageFont.truetype(p, size)
    return ImageFont.truetype("Arial.ttf", size)
im = Image.open(SRC).convert("RGBA"); A = np.array(im)

def inpaint_white(x0, y0, x1, y1, thr=200):
    rgb = cv2.cvtColor(A[y0:y1, x0:x1, :3].copy(), cv2.COLOR_RGB2BGR)
    m = ((A[y0:y1, x0:x1, :3].min(axis=2) > thr) & (A[y0:y1, x0:x1, 3] > 0)).astype(np.uint8) * 255
    m = cv2.dilate(m, np.ones((9, 9), np.uint8))
    out = cv2.inpaint(rgb, m, 9, cv2.INPAINT_TELEA)
    A[y0:y1, x0:x1, :3] = cv2.cvtColor(out, cv2.COLOR_BGR2RGB)

def fill_box(x0, y0, x1, y1, inset=22):
    c = np.median(A[y0 + inset:y0 + inset + 30, x0 + inset:x0 + inset + 200, :3].reshape(-1, 3), axis=0)
    A[y0 + inset:y1 - inset, x0 + inset:x1 - inset, :3] = c; A[y0 + inset:y1 - inset, x0 + inset:x1 - inset, 3] = 255

def clear_text(x0, y0, x1, y1):           # transparent background region with dark text
    A[y0:y1, x0:x1, 3] = 0

# 1 header boxes
fill_box(500, 410, 2000, 990); fill_box(4314, 410, 5740, 990)
# 2 CH4 bubble -> N2O
cx, cy = 622, 1118
red = A[cy - 45, cx, :3].copy()
yy, xx = np.ogrid[:A.shape[0], :A.shape[1]]
disc = (xx - cx) ** 2 + (yy - cy) ** 2 <= 50 ** 2
A[disc, :3] = red
# 3 gene and reservoir labels over soil
inpaint_white(780, 2490, 1050, 2610)
inpaint_white(1390, 2350, 1860, 2850)
inpaint_white(4490, 2360, 4830, 2750)
inpaint_white(5200, 2320, 5690, 2820)
# 4 misspelled line
clear_text(3110, 1274, 3895, 1326)

# 5 'ETC = EAC + EDC' and 'Total redox-active charge reservoir' -> assay exchange capacity
clear_text(2755, 1738, 3425, 1868)
# 6 DDQ structure -> generic quinone (drawn below)
clear_text(2045, 1668, 2432, 1998)
# 7 time axis labels -> moisture axis (hysteresis needs a forcing axis, as in Fig. 2c)
clear_text(2700, 2408, 3492, 2537)
# 8 downward arrow over core B (implied uptake) -> small upward arrow (smaller escape)
box = (slice(1000, 1155), slice(5488, 5578))
sub = A[box]; dark = (sub[..., :3].max(axis=2) < 70) & (sub[..., 3] > 0)
m = cv2.dilate(dark.astype(np.uint8) * 255, np.ones((7, 7), np.uint8))
rgb = cv2.inpaint(np.ascontiguousarray(sub[..., :3]), m, 7, cv2.INPAINT_TELEA)
alp = cv2.inpaint(np.ascontiguousarray(sub[..., 3]), m, 7, cv2.INPAINT_TELEA)
sub[..., :3] = rgb; sub[..., 3] = alp; A[box] = sub

im = Image.fromarray(A); d = ImageDraw.Draw(im)
INK = (20, 20, 20); W = (255, 255, 255)
def ctext(x, y, s, f, fill, anchor="mm"): d.text((x, y), s, font=f, fill=fill, anchor=anchor)
def eh_line(x, y, f, fs):
    fsub = font("Regular", int(fs * 0.62))
    wE = d.textlength("E", font=f); wh = d.textlength("h", font=fsub); wr = d.textlength(" = 200 mV", font=f)
    x0 = x - (wE + wh + wr) / 2
    d.text((x0, y), "E", font=f, fill=INK, anchor="lm")
    d.text((x0 + wE, y + fs * 0.22), "h", font=fsub, fill=INK, anchor="lm")
    d.text((x0 + wE + wh, y), " = 200 mV", font=f, fill=INK, anchor="lm")
for x, title, line3 in ((1250, "Core A", "N₂O sink activates late"), (5027, "Core B", "N₂O sink active")):
    ctext(x, 530, title, font("Bold", 120), INK)
    fs = 66; f1 = font("Regular", fs); fsub = font("Regular", int(fs * 0.62))
    a, b = "Matched bulk moisture, E", " and N₂O flux"
    w = d.textlength(a, font=f1) + d.textlength("h", font=fsub) + d.textlength(b, font=f1); x0 = x - w / 2
    d.text((x0, 700), a, font=f1, fill=INK, anchor="lm"); x0 += d.textlength(a, font=f1)
    d.text((x0, 700 + fs * 0.22), "h", font=fsub, fill=INK, anchor="lm"); x0 += d.textlength("h", font=fsub)
    d.text((x0, 700), b, font=f1, fill=INK, anchor="lm")
    ctext(x, 850, line3, font("Bold", 96), INK)
ctext(cx, cy, "N₂O", font("Regular", 50), W)
ctext(915, 2550, "nosZ", font("BoldItalic", 118), W)
ctext(1600, 2560, "N₂O\nreduction\ndelayed", font("Bold", 100), W)
ctext(4660, 2550, "nosZ", font("BoldItalic", 118), W)
ctext(5445, 2560, "N₂O\nreduction\nresumes", font("Bold", 100), W)
d.text((3120, 1294), "Spatially structured access", font=font("Italic", 47), fill=(52, 72, 120), anchor="lm")
# assay exchange capacity line
fb = font("BoldItalic", 48); parts = [("Assay exchange capacity: ", INK), ("EAC", (84, 130, 53)), (" + ", INK), ("EDC", (230, 130, 40))]
tw = sum(d.textlength(t, font=fb) for t, _ in parts); x = 3050 - tw / 2
for t, c in parts:
    d.text((x, 1800), t, font=fb, fill=c, anchor="lm"); x += d.textlength(t, font=fb)
# x axis of the hysteresis panel
d.text((3080, 2455), "Water-filled pore space", font=font("BoldItalic", 56), fill=(52, 72, 120), anchor="mm")
# generic 1,4-benzoquinone (replaces DDQ)
import math
cx, cy, r = 2235, 1835, 52; lw = 7; K = (25, 25, 25)
pts = [(cx + r * math.cos(math.radians(90 + 60 * i)), cy - r * math.sin(math.radians(90 + 60 * i))) for i in range(6)]
for i in range(6): d.line([pts[i], pts[(i + 1) % 6]], fill=K, width=lw)
for i, j in ((1, 2), (4, 5)):           # ring C=C double bonds (inner lines)
    (x1, y1), (x2, y2) = pts[i], pts[j]
    f = 0.80; mx, my = cx + (x1 - cx) * f, cy + (y1 - cy) * f; nx, ny = cx + (x2 - cx) * f, cy + (y2 - cy) * f
    d.line([(mx, my), (nx, ny)], fill=K, width=lw)
for (px, py), sgn in ((pts[0], -1), (pts[3], 1)):  # C=O at top and bottom
    d.line([(px - 9, py), (px - 9, py + sgn * 38)], fill=K, width=lw); d.line([(px + 9, py), (px + 9, py + sgn * 38)], fill=K, width=lw)
    d.text((px, py + sgn * 64), "O", font=font("Regular", 46), fill=K, anchor="mm")
d.text((cx, 2010), "Quinone", font=font("Bold", 50), fill=K, anchor="mm")
# small upward arrow over core B: escape continues, but smaller than from core A
d.line([(5532, 1100), (5532, 1045)], fill=(34, 34, 34), width=14)
d.polygon([(5510, 1050), (5554, 1050), (5532, 1012)], fill=(34, 34, 34))
# N2 released from core B (complete reduction)
d.ellipse([(5638, 1048), (5742, 1152)], fill=(43, 33, 24))
d.text((5690, 1100), "N₂", font=font("Bold", 50), fill=(255, 255, 255), anchor="mm")
d.line([(5690, 1046), (5690, 1012)], fill=(34, 34, 34), width=10)
d.polygon([(5674, 1016), (5706, 1016), (5690, 990)], fill=(34, 34, 34))
im.save(OUT, dpi=(600, 600)); print("saved", OUT)
