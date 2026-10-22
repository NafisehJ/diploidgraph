"""Plot synthetic_scaling.py results: diploid size vs naive, and runtime, against input size.

    python experiments/plot_synthetic.py experiments/results/synthetic.csv --out experiments/results/synthetic.png
"""

from __future__ import annotations

import argparse
import csv

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("csv")
    ap.add_argument("--out", default="experiments/results/synthetic.png")
    args = ap.parse_args()

    with open(args.csv) as fh:
        rows = list(csv.DictReader(fh))
    n = [int(r["n_vertices"]) for r in rows]
    vd = [int(r["diploid_vertices"]) for r in rows]
    naive = [float(r["log10_naive_ordered_all_paths"]) for r in rows]
    secs = [float(r["time_decompose_s"]) + float(r["time_build_s"]) for r in rows]

    fig, (a1, a2) = plt.subplots(1, 2, figsize=(10, 4))
    a1.scatter(n, [__import__("math").log10(v) for v in vd], s=14, label="this method: |V(G_D)|")
    a1.scatter(n, naive, s=14, marker="x", label="naive: N(N-1) path pairs")
    a1.set_xscale("log")
    a1.set_xlabel("input vertices |V|")
    a1.set_ylabel("log10(vertices)")
    a1.legend(frameon=False)
    a2.scatter(n, secs, s=14)
    a2.set_xscale("log")
    a2.set_yscale("log")
    a2.set_xlabel("input vertices |V|")
    a2.set_ylabel("decompose + build time (s)")
    for a in (a1, a2):
        a.spines[["top", "right"]].set_visible(False)
    fig.tight_layout()
    fig.savefig(args.out, dpi=200)
    print(f"wrote {args.out}")


if __name__ == "__main__":
    main()
