#!/usr/bin/env bash
# Reproduce the results table in README.md (BRCA1 and MHC, HPRC v2.0 Minigraph-Cactus).
#
# Usage:
#   bash reproduce.sh                 # both loci (MHC takes about 30 minutes, 1.5 GB of RAM)
#   bash reproduce.sh --quick         # BRCA1 only (under a minute)
#   DATA=~/brca1 bash reproduce.sh    # use a folder that already holds BRCA1.gfa / MHC.gfa
#
# Needs: python >= 3.10 with this package installed (python -m pip install -e .),
#        curl, and vg (https://github.com/vgteam/vg) only if the .gfa files are not there yet.
#        Set VG to the vg command if it is not on PATH, e.g. VG="conda run -n vgenv vg".

set -euo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
DATA="${DATA:-$ROOT/data}"
OUT="${OUT:-$ROOT/experiments/results}"
VG="${VG:-vg}"
URL="https://public.gi.ucsc.edu/cgl/ci/vg/vg-data/bakeoff"
LOCI="${LOCI:-BRCA1 MHC}"
[[ "${1:-}" == "--quick" ]] && LOCI="BRCA1"

mkdir -p "$DATA" "$OUT"
python -c "import diploidgraph" 2>/dev/null || {
  echo "diploidgraph is not installed: run  python -m pip install -e .  in $ROOT" >&2; exit 1; }

step() { printf '\n== %s\n' "$*"; }

for L in $LOCI; do
  GFA="$DATA/$L.gfa"
  if [[ ! -s "$GFA" ]]; then
    GBZ="$DATA/hprc-v2.0-mc-grch38.ec1M-$L.gbz"
    if [[ ! -s "$GBZ" ]]; then
      step "$L: downloading $(basename "$GBZ")"
      curl -fL -o "$GBZ" "$URL/$(basename "$GBZ")"
    fi
    step "$L: converting to GFA with vg"
    $VG version >/dev/null 2>&1 || {
      echo "vg not found. Install it (e.g. conda install -c bioconda vg) or set VG=..." >&2; exit 1; }
    $VG convert -f "$GBZ" > "$GFA"
  fi

  step "$L: build the diploid graph, verify it, measure time and memory"
  python -m diploidgraph.cli "$GFA" --json "$OUT/${L}_stats.json" > /dev/null

  step "$L: ablation without decomposition (same pairs)"
  python "$ROOT/experiments/no_decomposition.py" "$GFA" --json "$OUT/${L}_nodecomp.json" > /dev/null

  step "$L: check sampled haplotype pairs one by one (each must be decoded exactly once)"
  python "$ROOT/experiments/pair_check.py" "$GFA" --pairs "${PAIRS:-100}" --json "$OUT/${L}_pairs.json" > /dev/null

  step "$L: diploid graph and ablation sizes over ${ORDERS:-20} random vertex orders"
  python "$ROOT/experiments/order_sensitivity.py" "$GFA" --orders "${ORDERS:-20}" --json "$OUT/${L}_order.json" > /dev/null

  step "$L: state-space size of the PanGenie model on the same graph (size only)"
  python "$ROOT/experiments/external_state_space.py" "$GFA" --out "$OUT/${L}_external.csv" > /dev/null

  step "$L: which added vertices are used by observed haplotypes"
  python "$ROOT/experiments/explain_unsupported.py" "$GFA" --json "$OUT/${L}_unsupported.json" > /dev/null
done

step "Results (compare with README.md; times and memory depend on the machine)"
python "$ROOT/experiments/summarize.py" "$OUT"
if python -c "import matplotlib" 2>/dev/null; then
  python "$ROOT/experiments/plot_real.py" "$OUT/real_loci.json" --out "$OUT/real_loci.png" > /dev/null
  echo "figure: $OUT/real_loci.png"
fi
