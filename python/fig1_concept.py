"""Figure 1 | Out-of-step recovery, cumulative functional deviation and hidden redox legacies.

Panels
  a  Delayed functional recovery: bulk indicators recover before connectivity, capacity and kinetics.
  b  Cumulative functional deviation: two trajectories that converge on the same endpoint
     accumulate different signed deviations (redox debt, D_W).
  c  Latent functional legacy: two soils (restored, legacy) with identical bulk conditions and function before a
     common second disturbance, but a displaced internal property in soil B, respond differently.
  d  Internal redox state (capacity, connectivity, kinetics) and the memory-update rule.

All trajectories are schematic: they are constructed to illustrate the concepts, not fitted to data
or predicted from equations (1)-(3). Panel titles are given in the figure legend, not on the figure.

Run from python/:  python fig1_concept.py  ->  ../figures/Figure_1_out_of_step.{pdf,svg,png}
Requires numpy and matplotlib; shared typography comes from style.py (Arial or Liberation Sans, 7 pt).
"""
import os

import numpy as np
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch, Rectangle

from style import COL, FS, FS_LETTER, plt

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "figures")
NAME = "Figure_1_out_of_step"

SOIL_A = "#285E87"      # internal state restored
SOIL_B = "#7B4F97"      # persistent legacy
EVENT_SHADE = "#ECEDEF"
WINDOW_SHADE = "#F3EBD9"
FW, FH = 183 / 25.4, 5.98  # figure width and height (inches); 183 mm = double column


def add_axes(fig, x, y_top, w, h):
    """Axes placed by left edge and top edge in inches from the figure's top-left corner."""
    return fig.add_axes([x / FW, (FH - y_top - h) / FH, w / FW, h / FH])


def panel_letter(fig, x, y_top, letter):
    fig.text(x / FW, 1 - y_top / FH, letter, ha="left", va="center",
             fontsize=FS_LETTER, fontweight="bold")


def smoothstep(x):
    x = np.clip(x, 0, 1)
    return x * x * (3 - 2 * x)


def dip(t, start, end, recover, depth):
    """Smooth decline between start and end, then recovery to baseline by `recover`."""
    down = smoothstep((t - start) / (end - start))
    up = 1 - smoothstep((t - end) / (recover - end))
    return depth * np.where(t < end, down, up)


# --------------------------------------------------------------------------------------------
# Panel a: recovery order of bulk indicators and internal-state components
# --------------------------------------------------------------------------------------------
def panel_a(fig):
    t = np.linspace(0, 12, 2401)
    bulk = 1 - np.exp(-t / 0.3)
    conn = 1 - np.exp(-t / 1.2)
    cap = 1 - 0.6 * np.exp(-t / 3)                      # reserve only partly depleted
    kin = np.where(t < 1.5, 0, 1 - np.exp(-np.maximum(t - 1.5, 0) / 0.8))
    func = conn * cap * kin                              # illustrative product form
    t_bulk = t[np.argmax(bulk >= 0.95)]
    t_func = t[np.argmax(func >= 0.95)]

    ax = add_axes(fig, 0.55, 0.22, 2.75, 1.55)
    ax.axvspan(t_bulk, t_func, color=WINDOW_SHADE, lw=0, zorder=0)
    for x, lab in ((t_bulk, "$t_{bulk}$"), (t_func, "$t_F$")):
        ax.axvline(x, color="0.42", ls=":", lw=0.7, zorder=1)
        ax.text(x, 1.16, lab, ha="center", va="center")
    curves = [(bulk, COL["Bulk"], 1.2, "bulk (water, E$_h$)"),
              (conn, COL["Connectivity"], 1.3, "connectivity, α"),
              (cap, COL["Capacity"], 1.3, "capacity, Q"),
              (kin, COL["Kinetics"], 1.3, "kinetics, k"),
              (func, COL["ink"], 1.9, "function")]
    for y, c, lw, lab in curves:
        ax.plot(t, y, color=c, lw=lw, label=lab, zorder=3)
    ax.text((t_bulk + t_func) / 2 - 0.4, 0.055, "out-of-step\nwindow", ha="center", va="bottom",
            fontweight="bold")
    ax.set(xlim=(0, 12), ylim=(0, 1.24), yticks=[0, 0.5, 1], xticks=[],
           xlabel="Time after disturbance (arbitrary units)", ylabel="Recovery score")
    ax.legend(loc="lower right", bbox_to_anchor=(1.01, 0.015), handlelength=1.35,
              labelspacing=0.18, borderaxespad=0.05)
    panel_letter(fig, 0.05, 0.12, "a")
    return t_bulk, t_func


# --------------------------------------------------------------------------------------------
# Panel b: cumulative deviation (redox debt) of two trajectories with the same endpoint
# --------------------------------------------------------------------------------------------
def panel_b(fig):
    t = np.linspace(0, 10, 2001)
    dt = t[1] - t[0]

    def trajectory(lag, tau, overshoot):
        y = np.ones_like(t)
        drop = (t >= 1) & (t < 2.2)
        y[drop] = 1 - 0.75 * (t[drop] - 1) / 1.2
        after = t >= 2.2
        s = t[after] - 2.2
        s_lag = np.maximum(s - lag, 0)
        y[after] = np.where(s < lag, 0.25,
                            1 - 0.75 * np.exp(-s_lag / tau) + overshoot * s_lag * np.exp(-s_lag / 0.9))
        return y

    fast, slow = trajectory(0, 0.45, 0), trajectory(1.6, 0.9, 1.1)
    debt_fast, debt_slow = np.sum(1 - fast) * dt, np.sum(1 - slow) * dt

    ax = add_axes(fig, 4.15, 0.22, 2.88, 1.55)
    ax.axvspan(1, 2.2, color="0.9", lw=0)
    ax.text(1.6, 1.39, "event", ha="center", va="center", color="0.3")
    ax.fill_between(t, np.minimum(fast, 1), 1, color="0.5", alpha=0.20, lw=0)
    ax.fill_between(t, np.minimum(slow, 1), 1, color=COL["Debt"], alpha=0.28, lw=0)
    ax.fill_between(t, 1, np.maximum(slow, 1), facecolor="none", edgecolor=COL["Debt"],
                    hatch="/////", lw=0)
    ax.plot(t, np.ones_like(t), color="0.5", ls="--", lw=0.75)
    ax.plot(t, fast, color="0.25", lw=1.4)
    ax.plot(t, slow, color=COL["Debt"], lw=1.5)
    arrow = dict(arrowstyle="-", lw=0.5, color="0.25")
    ax.annotate("fast path", xy=(2.62, 0.78), xytext=(3.25, 1.34), color="0.25", va="center",
                arrowprops=arrow)
    ax.text(4.55, 0.42, "slow path", color=COL["Debt"], va="center")
    ax.annotate("lag", xy=(3.1, 0.25), xytext=(3.1, 0.09), ha="center", va="center", arrowprops=arrow)
    ax.annotate("overshoot", xy=(5.3, 1.12), xytext=(6.55, 1.37), va="center", arrowprops=arrow)
    ax.text(9.92, 0.67, r"$D_W=s\int_W(F-F_{ref})\,dt$", ha="right")
    ax.text(9.92, 0.48, f"fast: $|D_W|$ = {debt_fast:.2f}", ha="right", color="0.25")
    ax.text(9.92, 0.31, f"slow: $|D_W|$ = {debt_slow:.2f}", ha="right", color=COL["Debt"])
    ax.set(xlim=(0, 10), ylim=(0, 1.5), xticks=[], yticks=[1], yticklabels=["$F_{ref}$"],
           xlabel="Time within window W (arbitrary units)", ylabel="Function F")
    panel_letter(fig, 3.62, 0.12, "b")
    return debt_fast, debt_slow


# --------------------------------------------------------------------------------------------
# Panel c: latent functional legacy across a first and a common second disturbance
# --------------------------------------------------------------------------------------------
def panel_c(fig):
    t = np.linspace(0, 11.5, 2301)
    event1, event2, check = (0.95, 1.85), (7.35, 8.05), 7.13   # check = pre-challenge assessment

    bulk = 1 - dip(t, *event1, 3.45, 0.51) - dip(t, *event2, 9.25, 0.50)   # identical in A and B
    func_first = 1 - dip(t, *event1, 4.65, 0.57)                           # identical before challenge
    func_a = func_first - dip(t, *event2, 9.55, 0.31)
    func_b = func_first - dip(t, *event2, 11.15, 0.70)
    state_a = 1 - dip(t, *event1, 4.45, 0.58)                               # state restored
    state_b = (1 - 0.58 * smoothstep((t - event1[0]) / (event1[1] - event1[0]))
               + 0.23 * smoothstep((t - event1[1]) / (5.8 - event1[1])))     # persistent shift
    before = t <= check                                                      # carrier not drawn after

    left, width, height = 1.18, 5.84, 0.53
    tops = [2.36, 2.99, 3.62]
    axes = [add_axes(fig, left, y, width, height) for y in tops]
    for k, ax in enumerate(axes):
        for ev in (event1, event2):
            ax.axvspan(*ev, color=EVENT_SHADE, lw=0, zorder=0)
        ax.axvline(check, color="#777777", ls=(0, (2, 2)), lw=0.8, zorder=1)
        ax.axhline(1, color="#B8B8B8", ls=(0, (3, 2)), lw=0.65, zorder=1)
        ax.set(xlim=(0, 11.5), ylim=(0.16, 1.17), xticks=[], yticks=[0.25, 1])
        ax.set_yticklabels(["0.25", "1"])
        ax.tick_params(axis="y", pad=1.5)
        ax.spines["left"].set_color("#777777")
        ax.spines["bottom"].set_visible(k == 2)

    top = axes[0]
    top.text(sum(event1) / 2, 1.30, "First disturbance", ha="center", va="bottom", clip_on=False)
    top.text(check - 0.12, 1.30, "Pre-challenge assessment", ha="right", va="bottom",
             fontweight="bold", clip_on=False)
    top.text(event2[0] + 0.05, 1.30, "Common second disturbance", ha="left", va="bottom", clip_on=False)
    top.plot(t, bulk, color="#555555", lw=1.55, zorder=3)
    top.plot(check, 1, "o", color="#555555", ms=3, zorder=4)
    top.text(4.95, 0.51, "Same bulk conditions", ha="center", va="center", color="#4C4C4C")

    mid = axes[1]
    mid.plot(t[before], state_a[before], color=SOIL_A, lw=1.7, zorder=3)
    mid.plot(t[before], state_b[before], color=SOIL_B, lw=1.7, ls=(0, (4, 2)), zorder=3)
    for y, c in ((state_a, SOIL_A), (state_b, SOIL_B)):
        mid.plot(check, np.interp(check, t, y), "o", ms=3.4, color=c, zorder=4)
    mid.text(4.52, 1.08, "State restored", color=SOIL_A, ha="center", va="bottom")
    mid.text(4.5, 0.38, "Persistent legacy", color=SOIL_B, ha="center", va="center")
    mid.annotate("", xy=(check - 0.16, 0.98), xytext=(check - 0.16, 0.66),
                 arrowprops=dict(arrowstyle="<->", lw=0.8, color=SOIL_B))
    mid.text(9.65, 0.73, "Candidate carrier\nmeasured before challenge", color="#666666",
             ha="center", va="center")

    bottom = axes[2]
    bottom.plot(t, func_a, color=SOIL_A, lw=1.7, zorder=3)
    bottom.plot(t, func_b, color=SOIL_B, lw=1.65, ls=(0, (4, 2)), zorder=4)
    bottom.plot(check, 1, "o", color=SOIL_A, ms=3, zorder=5)
    bottom.text(5.35, 0.53, "Same present function", color="#414141", ha="center", va="center")
    for lab, y_lab, y, c in (("restored", 0.72, func_a, SOIL_A), ("legacy", 0.23, func_b, SOIL_B)):
        bottom.annotate(lab, xy=(8.6, np.interp(8.6, t, y)), xytext=(9.0, y_lab), color=c,
                        fontweight="bold",
                        arrowprops=dict(arrowstyle="-", lw=0.65, color=c))
    bottom.set_xlabel("Time across disturbances (arbitrary units)", labelpad=3)

    for y, lab in zip(tops, ["Bulk recovery", "Internal state\n(candidate\ncarrier)", "Target function"]):
        fig.text(0.88 / FW, 1 - (y + height / 2) / FH, lab, ha="right", va="center", fontweight="bold")
    panel_letter(fig, 0.05, 2.05, "c")

    # The schematic must show what the legend states.
    i = np.searchsorted(t, check)
    assert np.isclose(bulk[i], 1.0), "bulk conditions not recovered at the assessment"
    assert np.isclose(func_a[i], func_b[i]) and np.isclose(func_a[i], 1.0), "functions differ before challenge"
    assert abs(state_a[i] - state_b[i]) > 0.1, "no internal-state difference at the assessment"
    assert np.max(np.abs(func_a[t >= event2[0]] - func_b[t >= event2[0]])) > 0.3, "no divergence after challenge"


# --------------------------------------------------------------------------------------------
# Panel d: internal redox state and memory-update rule
# --------------------------------------------------------------------------------------------
def panel_d(fig):
    ax = add_axes(fig, 0.10, 4.66, 7.0, 1.25)
    ax.set_axis_off()
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 34)

    def box(x, y, w, h, colour, fill=None):
        style = "round,pad=0,rounding_size=1.2"
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, fc=fill or colour, ec="none",
                                    alpha=0.13 if fill is None else 1))
        ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle=style, fc="none", ec=colour, lw=1))

    def tile(x, colour, name, symbol):
        box(x, 19.5, 20.5, 9.5, colour)
        ax.text(x + 10.25, 25.6, name, ha="center", va="center", fontweight="bold", color=colour)
        ax.text(x + 10.25, 21.9, symbol, ha="center", va="center", style="italic", color=colour)

    ax.add_patch(Rectangle((1.3, 18.2), 66.4, 12.1, fc="none", ec="0.6", ls=(0, (2, 2)), lw=0.6))
    ax.text(34.5, 31.2, "present state $z_n$", ha="center", va="bottom", color="0.3")
    tile(2.5, COL["Capacity"], "CAPACITY", "accepting inventory, Q")
    tile(24.5, COL["Connectivity"], "CONNECTIVITY", "connected fraction, α")
    tile(46.5, COL["Kinetics"], "KINETICS", "exchange coefficient, k")

    box(11, 1.5, 47, 8.0, "#5B7A99", fill="#EEF3F7")
    ax.text(34.5, 7.1, "Restricted accepting-pool example", ha="center", va="center", color="#425B72")
    ax.text(34.5, 4.0, r"$C_{acc}(\tau)=\Sigma_i Q_i\alpha_i(1-e^{-k_i\tau})$  (equation 1)",
            ha="center", va="center")
    for x in (12.75, 34.75, 56.75):
        ax.add_patch(FancyArrowPatch((x, 18.0), (34.5 + (x - 34.75) * 0.55, 9.8), arrowstyle="-|>",
                                     mutation_scale=6, color="0.35", lw=0.7))

    mem = COL["Memory"]
    box(73, 4, 25, 26.3, mem)
    ax.text(85.5, 26.6, "MEMORY", ha="center", va="center", fontweight="bold", color=mem)
    ax.text(85.5, 22.9, "update rule, M", ha="center", va="center", style="italic", color=mem)
    ax.text(85.5, 16.2, r"$z_{n+1}=M(z_n,H_n)$", ha="center", va="center")
    ax.text(85.5, 9, r"persistence ratio $R=T_m/\Delta$", ha="center", va="center", color=mem)
    ax.add_patch(FancyArrowPatch((67.9, 26), (72.8, 26), arrowstyle="-|>", mutation_scale=7,
                                 color=mem, lw=0.9))
    ax.text(70.4, 27.4, "$H_n$", ha="center", va="bottom", color=mem)
    ax.add_patch(FancyArrowPatch((72.8, 9), (61.5, 18), connectionstyle="arc3,rad=-.5",
                                 arrowstyle="-|>", mutation_scale=7, color=mem, lw=0.9))
    ax.text(65, 3.2, "update to $z_{n+1}$", ha="center", va="center", color=mem)
    panel_letter(fig, 0.05, 4.50, "d")


def main():
    fig = plt.figure(figsize=(FW, FH), facecolor="white")
    t_bulk, t_func = panel_a(fig)
    debt_fast, debt_slow = panel_b(fig)
    panel_c(fig)
    panel_d(fig)
    os.makedirs(OUT, exist_ok=True)
    for ext in ("pdf", "svg", "png"):
        fig.savefig(os.path.join(OUT, f"{NAME}.{ext}"), dpi=600, facecolor="white")
    plt.close(fig)
    print(f"Figure 1 written to figures/: t_bulk={t_bulk:.2f}, t_F={t_func:.2f}, "
          f"|D_W| fast={debt_fast:.2f}, slow={debt_slow:.2f}; panel c checks passed.")


if __name__ == "__main__":
    main()
