"""Figure: diploid graph size vs naive pairing of observed haplotypes, BRCA1 and MHC.

    python experiments/plot_real.py experiments/figures/real_loci.json --out experiments/figures/real_loci.png
"""

from __future__ import annotations

import argparse
import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.ticker import FixedLocator, NullLocator  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8a8984"
GRID = "#e6e5e0"
INPUT = "#b9b8b2"   # neutral: the input graph is context, not a series
OURS = "#2a78d6"    # categorical slot 1
NAIVE = "#eb6834"   # categorical slot 2


def human(n: float) -> str:
    if n < 1e6:
        return f"{n:,.0f}"
    for div, unit in ((1e9, " billion"), (1e6, " million")):
        if n >= div:
            v = n / div
            return (f"{v:.0f}" if v >= 10 else f"{v:.1f}") + unit
    return f"{n:,.0f}"


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("json")
    ap.add_argument("--out", default="experiments/figures/real_loci.png")
    args = ap.parse_args()
    data = json.load(open(args.json))
    loci = ["BRCA1", "MHC"]

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 11})
    fig, ax = plt.subplots(figsize=(9, 5.6), facecolor=SURFACE)
    ax.set_facecolor(SURFACE)

    width, gap = 0.26, 0.03
    series = [("Input pangenome graph", "input_nodes", INPUT),
              ("Our diploid graph", "diploid_nodes", OURS),
              ("Naive pairing of observed haplotypes", "naive_observed_nodes", NAIVE)]
    for k, (label, key, color) in enumerate(series):
        xs = [i + (k - 1) * (width + gap) for i in range(len(loci))]
        ys = [data[l][key] for l in loci]
        ax.bar(xs, ys, width=width, color=color, label=label, zorder=3)
        for x, y in zip(xs, ys):
            ax.text(x, y * 1.35, human(y), ha="center", va="bottom", fontsize=10, color=INK2, zorder=4)

    for i, l in enumerate(loci):
        d = data[l]
        ratio = d["naive_observed_nodes"] / d["diploid_nodes"]
        x = i  # centre of "our diploid graph" bar
        lo, hi = d["diploid_nodes"] * 4.5, d["naive_observed_nodes"]
        ax.annotate("", xy=(x, hi), xytext=(x, lo), zorder=5,
                    arrowprops=dict(arrowstyle="<->", color=OURS, lw=1.6, shrinkA=0, shrinkB=0))
        ax.plot([x, x + width / 2 + gap], [hi, hi], color=OURS, lw=1, ls=(0, (2, 2)), zorder=5)
        mid = (lo * hi) ** 0.5
        ax.text(x - 0.03, mid, f"~{round(ratio, -3):,.0f}\u00d7\nfewer nodes\nthan naive", ha="right", va="center",
                fontsize=12, fontweight="bold", color=OURS, zorder=5, linespacing=1.15)

    ax.set_yscale("log")
    ax.set_ylim(1e3, 1e13)
    ax.yaxis.set_major_locator(FixedLocator([1e3, 1e5, 1e7, 1e9, 1e11, 1e13]))
    ax.yaxis.set_minor_locator(NullLocator())
    ax.set_yticklabels(["1K", "100K", "10M", "1B", "100B", "10T"], color=MUTED)
    ax.set_ylabel("Nodes (log scale)", color=INK2)
    ax.set_xticks(range(len(loci)))
    ax.set_xticklabels([f"{l}\n{data[l]['input_nodes']:,} nodes, {data[l]['haplotypes']} haplotypes" for l in loci],
                       color=INK)
    ax.grid(axis="y", color=GRID, lw=0.8, zorder=0)
    for s in ("top", "right", "left"):
        ax.spines[s].set_visible(False)
    ax.spines["bottom"].set_color(MUTED)
    ax.tick_params(length=0)

    fig.suptitle("Diploid graph size vs. naive haplotype pairing", x=0.075, ha="left", fontsize=15,
                 fontweight="bold", color=INK, y=0.97)
    fig.text(0.075, 0.885, "HPRC v2.0 Minigraph-Cactus loci. Our graph represents every pair of haplotypes, "
             "including all recombinants;\nnaive pairing covers only the observed pairs.",
             ha="left", fontsize=10, color=INK2)
    ax.legend(loc="upper left", frameon=False, fontsize=10, bbox_to_anchor=(0.0, 1.0))
    b, m = data["BRCA1"], data["MHC"]
    fig.text(0.075, 0.02,
             f"Build time: {b['build_s'] * 1000:.0f} ms (BRCA1), {m['build_s']:.1f} s (MHC). "
             "Every added node is used by real haplotypes at both loci.", ha="left", fontsize=9.5, color=INK2)
    fig.subplots_adjust(left=0.075, right=0.98, top=0.8, bottom=0.17)
    fig.savefig(args.out, dpi=200, facecolor=SURFACE)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
