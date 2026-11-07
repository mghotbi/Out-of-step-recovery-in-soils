"""Compose Supplementary Fig. 1: a, relabelled author artwork; b, illustrative-example calculation (figS1_worked_example.py b)."""
from PIL import Image, ImageDraw, ImageFont
import numpy as np, os
DPI = 600; mm = DPI / 25.4; Wpx = int(183 * mm)
art = Image.open("../figures/Supplementary_Figure_1a_artwork.png")
bg = Image.new("RGB", art.size, "white"); bg.paste(art, mask=art.split()[3])
a = np.array(bg.convert("L")); ys, xs = np.where(a < 245)
bg = bg.crop((max(xs.min() - 20, 0), max(ys.min() - 20, 0), xs.max() + 20, ys.max() + 20))
bg = bg.resize((Wpx, int(bg.height * Wpx / bg.width)), Image.LANCZOS)
pb = Image.open("../figures/Supplementary_Figure_1_worked_example_panel_b.png").convert("RGB")
pb = pb.resize((int(95 * mm), int(pb.height * 95 * mm / pb.width)), Image.LANCZOS)
gap = int(4 * mm); H = bg.height + gap + pb.height
fig = Image.new("RGB", (Wpx, H), "white"); fig.paste(bg, (0, 0)); fig.paste(pb, ((Wpx - pb.width) // 2, bg.height + gap))
d = ImageDraw.Draw(fig)
fp = "/usr/share/fonts/truetype/liberation2/LiberationSans-Bold.ttf"
f = ImageFont.truetype(fp if os.path.exists(fp) else "Arial Bold.ttf", int(9 / 72 * DPI))
d.text((int(1 * mm), int(1 * mm)), "a", font=f, fill="black")
for ext in ("png", "pdf"):
    fig.save(f"../figures/Supplementary_Figure_1_worked_example.{ext}", dpi=(DPI, DPI))
print(fig.size)
