"""Runtime, memory and diploid graph size vs. input size on random series-parallel DAGs.

    python experiments/synthetic_scaling.py --out experiments/results/synthetic.csv

Each row is one random graph (several seeds per size). Columns match stats.run_stats.
"""

from __future__ import annotations

import argparse
import csv
import pathlib

from diploidgraph.generators import random_sp_dag
from diploidgraph.stats import run_stats

DEFAULT_SIZES = [100, 300, 1000, 3000, 10000, 30000]


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sizes", type=int, nargs="+", default=DEFAULT_SIZES, help="number of generator operations")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--max-branches", type=int, default=3)
    ap.add_argument("--out", default="experiments/results/synthetic.csv")
    args = ap.parse_args()

    out = pathlib.Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    rows = []
    for n_ops in args.sizes:
        for seed in range(args.seeds):
            g = random_sp_dag(n_ops, max_branches=args.max_branches, seed=seed)
            r = {"n_ops": n_ops, "seed": seed}
            r.update(run_stats(g))
            rows.append(r)
            print(f"ops={n_ops:>6} seed={seed} |V|={r['n_vertices']:>7} |V_D|={r['diploid_vertices']:>8} "
                  f"build={r['time_decompose_s'] + r['time_build_s']:.3f}s "
                  f"log10(naive)={r['log10_naive_ordered_all_paths']:.1f}")
    with out.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
