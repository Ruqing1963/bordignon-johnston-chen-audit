"""Figure for the note on Bordignon-Johnston-Starichkova, arXiv:2207.09452 (IJNT 21 (2025) 2497-2572).

Run (after verify_bjs_chen.py, which writes the two CSV files read here):
    python code/make_figures.py [--png preview.png]
Writes figures/fig_bjs_chen_audit.pdf (vector, 6.3 in wide, two panels).

(a) Case (b) margin of the final inequality, pi_2(N)/(U_N N/log^2 N) - 0, at N = X_2, against loglog X_2,
    with X_3 = X_2^{1/8} and the best delta at each point: the paper's computation, the corrected computation
    with the paper's lemmas, and the corrected computation with the sharpened Lemma 6.3.
(b) C_1(eps), C_2(eps): Table 1 against the values implied by Lemmas 2.4-2.8 (corrected gamma_{s0}).
"""
import csv, os, sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from mpmath import mpf, exp
from verify_bjs_chen import chen

DATA = os.path.join(HERE, "..", "data")
BLUE, VERM, GREEN, GREY, BLACK, ORANGE = "#0072B2", "#D55E00", "#009E73", "#8C8C8C", "#000000", "#E69F00"


def read_csv(name):
    with open(os.path.join(DATA, name)) as fh:
        rows = list(csv.DictReader(fh))
    return {k: np.array([float(r[k]) for r in rows]) for k in rows[0]}


def main():
    png = sys.argv[sys.argv.index("--png") + 1] if "--png" in sys.argv else None
    M = read_csv("margin_scan.csv")
    S = read_csv("sieve_constants.csv")
    # the paper's own point (delta = 1.478, X_3 = exp(exp(30.62))) and the same point corrected
    pr = float(chen(LX3=exp(mpf('30.62')))['res_b'])
    co = float(chen(LX3=exp(mpf('30.62')), ufix=True, eps2_arg="X3", C1=mpf('132.6'), C2=mpf('136.1'))['res_b'])
    print("paper's point: printed %.6g, corrected %.6g" % (pr, co))

    plt.rcParams.update({
        "font.family": "cmr10", "mathtext.fontset": "cm", "axes.formatter.use_mathtext": True,
        "font.size": 9, "axes.labelsize": 9, "axes.titlesize": 9, "xtick.labelsize": 8, "ytick.labelsize": 8,
        "legend.fontsize": 8, "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.minor.width": 0.4, "ytick.minor.width": 0.4, "xtick.direction": "in", "ytick.direction": "in",
        "pdf.fonttype": 42, "axes.unicode_minus": False,
    })
    fig, (ax, bx) = plt.subplots(1, 2, figsize=(6.3, 3.0))
    fig.subplots_adjust(left=0.085, right=0.985, bottom=0.15, top=0.92, wspace=0.30)

    # ------------------------------------------------------------ (a)
    x = M["loglogX2"]
    sel = (x >= 32.5 - 1e-9) & (x <= 33.0 + 1e-9)
    x = x[sel]
    sc = 1e3
    ax.axhline(0, color=GREY, lw=0.5, zorder=1)
    ax.axhline(0.2, color=BLACK, lw=0.7, ls=(0, (4, 2)), zorder=1)
    ax.axvline(32.7, color=GREY, lw=0.6, ls=(0, (1, 1.5)), zorder=1)
    ax.plot(x, M["printed_margin"][sel] * sc, color=GREY, lw=1.0, zorder=2)
    ax.plot(x, M["corrected_margin"][sel] * sc, color=VERM, lw=1.1, zorder=3)
    ax.plot(x, M["corrected_sharp_margin"][sel] * sc, color=BLUE, lw=1.1, ls=(0, (5, 1.5)), zorder=4)
    ax.plot([32.7], [pr * sc], ls="", marker="o", ms=4.2, mfc="white", mec=BLACK, mew=0.8, zorder=5)
    ax.plot([32.7], [co * sc], ls="", marker="D", ms=4.0, mfc=VERM, mec=BLACK, mew=0.5, zorder=5)
    ax.set_xlim(32.5, 33.0)
    ax.set_ylim(-16, 9)
    ax.set_xticks([32.5, 32.6, 32.7, 32.8, 32.9, 33.0])
    ax.set_xticklabels([r"$32.5$", r"$32.6$", r"$32.7$", r"$32.8$", r"$32.9$", r"$33.0$"])
    ax.set_yticks([-15, -10, -5, 0, 5])
    ax.set_xlabel(r"$\log\log X_2$")
    ax.set_ylabel(r"case (b) lower bound $\times 10^{3}$")
    ax.set_title(r"(a) Final inequality, case (b), $N=X_2$")
    ax.annotate(r"target $2\cdot10^{-4}$", xy=(32.57, 0.2), xytext=(32.53, -3.6), fontsize=8, ha="left", va="top",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=BLACK, shrinkA=1, shrinkB=1))
    ax.annotate(r"paper, $\delta=1.478$", xy=(32.7, pr * sc), xytext=(32.76, -2.2), fontsize=8, ha="left",
                va="center", arrowprops=dict(arrowstyle="-", lw=0.5, color=BLACK, shrinkA=1, shrinkB=3))
    ax.annotate(r"corrected", xy=(32.7, co * sc), xytext=(32.66, -8.3), fontsize=8, ha="right", va="center",
                arrowprops=dict(arrowstyle="-", lw=0.5, color=BLACK, shrinkA=1, shrinkB=3))
    handles = [
        Line2D([], [], color=GREY, lw=1.0, label=r"paper's conventions"),
        Line2D([], [], color=VERM, lw=1.1, label=r"corrected"),
        Line2D([], [], color=BLUE, lw=1.1, ls=(0, (5, 1.5)), label="corrected,\nLemma 41 sharpened"),
    ]
    ax.legend(handles=handles, loc="lower right", frameon=False, borderpad=0.3,
              handletextpad=0.4, labelspacing=0.3, handlelength=2.2)

    # ------------------------------------------------------------ (b)
    inv = S["inv_eps"]
    tab = inv <= 1e5 + 1
    bx.axvspan(219.38, 2e7, color="#F0F0F0", zorder=0, lw=0)
    bx.axvline(219.38, color=GREY, lw=0.6, ls=(0, (1, 1.5)), zorder=1)
    bx.plot(inv, S["C1_direct"], color=BLUE, lw=1.0, zorder=3)
    bx.plot(inv, S["C2_direct"], color=VERM, lw=1.0, ls=(0, (4, 1.5)), zorder=3)
    bx.plot(inv[tab], S["C1_printed"][tab], ls="", marker="o", ms=5.0, mfc="white", mec=BLUE, mew=0.8, zorder=4)
    bx.plot(inv[tab], S["C2_printed"][tab], ls="", marker="s", ms=2.2, mfc=VERM, mec=VERM, mew=0.0, zorder=5)
    fin = np.isfinite(S["C1_lemma28"])
    bx.plot(inv[fin], S["C1_lemma28"][fin], ls="", marker="+", ms=4.5, mec=BLACK, mew=0.7, zorder=5)
    bx.set_xscale("log"); bx.set_yscale("log")
    bx.set_xlim(60, 1.2e7)
    bx.set_ylim(90, 1000)
    bx.set_yticks([100, 150, 200, 300, 500, 800])
    bx.set_yticklabels([r"$%d$" % t for t in [100, 150, 200, 300, 500, 800]])
    bx.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())
    bx.set_xticks([1e2, 1e3, 1e4, 1e5, 1e6, 1e7])
    bx.set_xlabel(r"$1/\varepsilon$")
    bx.set_ylabel(r"$C_1(\varepsilon),\ C_2(\varepsilon)$")
    bx.set_title(r"(b) Linear-sieve constants (Table 1)")
    bx.text(290, 360, "Lemma 14 gives\nno bound here\n" r"($J(499)\geq 1$)", fontsize=8, ha="left", va="top")
    bx.annotate(r"$\varepsilon=1.452\cdot10^{-7}$", xy=(6.887e6, 132.6), xytext=(6.887e6, 205), fontsize=8,
                ha="right", va="bottom", arrowprops=dict(arrowstyle="-", lw=0.5, color=BLACK, shrinkA=1, shrinkB=2))
    handles = [
        Line2D([], [], color=BLUE, lw=1.0, label="$C_1$, Lemmas 8, 11, 12"),
        Line2D([], [], color=VERM, lw=1.0, ls=(0, (4, 1.5)), label="$C_2$, Lemmas 8, 11, 12"),
        Line2D([], [], ls="", marker="o", ms=5.0, mfc="white", mec=BLUE, mew=0.8, label=r"$C_1$, Table 1"),
        Line2D([], [], ls="", marker="s", ms=2.2, mfc=VERM, mec=VERM, mew=0.0, label=r"$C_2$, Table 1"),
        Line2D([], [], ls="", marker="+", ms=4.5, mec=BLACK, mew=0.7, label=r"$C_1$, Lemma 14"),
    ]
    bx.legend(handles=handles, loc="upper right", frameon=False, borderpad=0.3, handletextpad=0.4,
              labelspacing=0.25, handlelength=2.2)

    outdir = os.path.join(os.path.dirname(HERE), "figures")
    os.makedirs(outdir, exist_ok=True)
    pdf = os.path.join(outdir, "fig_bjs_chen_audit.pdf")
    fig.savefig(pdf)
    print("wrote", pdf)
    if png:
        fig.savefig(png, dpi=250)
        print("wrote", png)


if __name__ == "__main__":
    main()


