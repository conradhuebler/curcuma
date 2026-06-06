#!/usr/bin/env bash
# sqm_bench.sh - Reproducible single-core benchmark for the native GFN1/GFN2 SCF.
#
# Compares curcuma (release/) against the Fortran references tblite and xtb on a
# molecule ladder, all pinned to ONE core with serial BLAS, min-of-N runs.
# For curcuma it parses the -verbosity 2 timing breakdown so the SCF step, the
# post-SCF energies (incl. D4 dispersion) and the gradient are separated from
# process startup. Energy+gradient workload throughout (curcuma -sp always does
# the gradient; xtb/tblite are asked for --grad to match).
#
# Usage:  scripts/sqm_bench.sh [N_REPEATS] [molecule ...]
#   N_REPEATS  number of timed repeats, min is reported (default 3)
#   molecule   subset of {caffeine triose complex} (default: all)
#
# Claude Generated. Not human production-tested.
set -u

REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
CURCUMA="$REPO/release/curcuma"
TBLITE="$REPO/release_tblite/_deps/tblite-build/app/tblite"
XTB="$(command -v xtb || true)"
MOLDIR="$REPO/test_cases/sqm_reference/molecules"

N=${1:-3}; shift 2>/dev/null || true
MOLS=("$@"); [ ${#MOLS[@]} -eq 0 ] && MOLS=(caffeine triose complex)
CORE=0
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OMP_MAX_ACTIVE_LEVELS=1
PIN=(taskset -c $CORE)
SCRATCH="$(mktemp -d /tmp/sqmbench.XXXXXX)"
strip() { sed -r 's/\x1b\[[0-9;]*m//g'; }

# min wall (ms) over N runs of "$@", stdout/stderr discarded, run in $SCRATCH
min_wall_ms() {
  local best=99999999 i s e w
  for ((i=0;i<N;i++)); do
    s=$(date +%s%N); ( cd "$SCRATCH" && "${PIN[@]}" "$@" >/dev/null 2>&1 ); e=$(date +%s%N)
    w=$(( (e-s)/1000000 )); (( w<best )) && best=$w
  done
  echo "$best"
}

printf "# sqm_bench  N=%d  core=%d  OMP/MKL=1\n" "$N" "$CORE"
printf "# curcuma=%s\n# tblite =%s\n# xtb    =%s\n\n" \
  "$([ -x "$CURCUMA" ] && echo yes || echo MISSING)" \
  "$([ -x "$TBLITE" ] && echo yes || echo MISSING)" \
  "$([ -n "$XTB" ] && echo "$XTB" || echo MISSING)"

for m in "${MOLS[@]}"; do
  xyz="$MOLDIR/$m.xyz"
  [ -f "$xyz" ] || { echo "skip $m (no $xyz)"; continue; }
  nat=$(head -1 "$xyz" | tr -d ' ')
  echo "================ $m ($nat atoms) ================"
  cp "$xyz" "$SCRATCH/$m.xyz"
  for meth in gfn1 gfn2; do
    echo "---- $meth ----"
    # curcuma: one verbose run for the breakdown, then min wall over N
    vb="$("${PIN[@]}" "$CURCUMA" -sp "$xyz" -method "$meth" -verbosity 2 2>&1 | strip)"
    iters=$(echo "$vb" | grep -oE "converged in [0-9]+ iterations" | grep -oE "[0-9]+" | head -1)
    setup=$(echo "$vb" | grep -E "setup *:" | grep -oE "[0-9]+\.[0-9]+" | head -1)
    scf=$(echo "$vb"   | grep -E "SCF \( *[0-9]+ it\)" | grep -oE "[0-9]+\.[0-9]+" | head -1)
    post=$(echo "$vb"  | grep -E "post-SCF E *:" | grep -oE "[0-9]+\.[0-9]+" | head -1)
    grd=$(echo "$vb"   | grep -E "gradient *:" | grep -oE "[0-9]+\.[0-9]+" | head -1)
    cwall=$([ -x "$CURCUMA" ] && min_wall_ms "$CURCUMA" -sp "$xyz" -method "$meth" -verbosity 0 || echo NA)
    printf "  curcuma : wall %6s ms | setup %7s | SCF %8s (%s it) | postSCF %7s | grad %7s\n" \
      "$cwall" "${setup:-?}" "${scf:-?}" "${iters:-?}" "${post:-?}" "${grd:-?}"
    # tblite
    if [ -x "$TBLITE" ]; then
      "${PIN[@]}" "$TBLITE" run --method "$meth" --grad "$SCRATCH/g.txt" "$xyz" >"$SCRATCH/tb.out" 2>&1
      tcyc=$(grep -E "^ *[0-9]+ +-?[0-9]" "$SCRATCH/tb.out" | tail -1 | awk '{print $1}')
      tscc=$(grep -E "^ *- scc" "$SCRATCH/tb.out" | grep -oE "[0-9]+\.[0-9]+" | head -1)
      tdsp=$(grep -E "^ *- dispersion" "$SCRATCH/tb.out" | grep -oE "[0-9]+\.[0-9]+" | head -1)
      twall=$(min_wall_ms "$TBLITE" run --method "$meth" --grad g.txt "$SCRATCH/$m.xyz")
      printf "  tblite  : wall %6s ms | scc %6ss (%s cyc) | dispersion %6ss\n" \
        "$twall" "${tscc:-?}" "${tcyc:-?}" "${tdsp:-?}"
    fi
    # xtb
    if [ -n "$XTB" ]; then
      gfn=${meth#gfn}
      # --norestart: xtb caches the converged wavefunction in xtbrestart, which
      # would make repeated runs in the same dir read it back and converge in ~1
      # iteration (warm). Force a cold single point for a fair comparison.
      ( cd "$SCRATCH" && rm -f xtbrestart && "${PIN[@]}" "$XTB" "$m.xyz" --gfn "$gfn" --grad --norestart >xtb.out 2>&1 )
      xscf=$(grep -A3 "^ SCF:" "$SCRATCH/xtb.out" | grep "wall-time" | grep -oE "[0-9]+\.[0-9]+ sec" | head -1)
      rm -f "$SCRATCH/xtbrestart"
      xwall=$(min_wall_ms "$XTB" "$m.xyz" --gfn "$gfn" --grad --norestart)
      printf "  xtb     : wall %6s ms | SCF %s\n" "$xwall" "${xscf:-?}"
    fi
  done
done
rm -rf "$SCRATCH"
