#!/usr/bin/env bash
# regression_bench.sh - compare two curcuma binaries numerically and by wall time.
#
# For every molecule x method x thread count it runs a single point with gradient
# (-dump_gradient) in a fresh directory (topology-cache miss), gfnff additionally a
# second time (cache hit), then reports for each run the energy difference (12 digits),
# the max |dgradient| and the wall-time ratio. Optional short MD runs (deterministic
# seed) catch anything a single point misses: an MD trajectory amplifies any 1-ulp force
# difference, so identical final energies are strong evidence of bit-identical forces.
#
# Usage:
#   scripts/regression_bench.sh REF_BINARY NEW_BINARY [quick] [OUTDIR]
#     quick  : small molecule subset, no MD
#     OUTDIR : where to keep logs/gradients (default: /tmp/curcuma_regbench)
#
# Molecules come from test_cases/ (sqm_reference set + polymer). Add your own by dropping
# NAME.xyz into $OUTDIR/mols before running (charge via CHG[NAME] below).
#
# Claude Generated (Sep 2026). Not human production-tested.
set -u
REF=$(readlink -f "$1"); NEW=$(readlink -f "$2"); QUICK=${3:-}; OUT=${4:-/tmp/curcuma_regbench}
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
mkdir -p "$OUT/mols"
for f in H2O caffeine H2S PH3 acetic_acid_dimer triose complex; do
  [ -f "$OUT/mols/$f.xyz" ] || cp "$REPO/test_cases/sqm_reference/molecules/$f.xyz" "$OUT/mols/"
done
[ -f "$OUT/mols/polymer.xyz" ] || cp "$REPO/test_cases/molecules/larger/polymer.xyz" "$OUT/mols/"
declare -A CHG=()
LIST="H2O caffeine H2S PH3 acetic_acid_dimer triose complex polymer"
[ "$QUICK" = quick ] && LIST="caffeine triose complex"
for extra in "$OUT"/mols/*.xyz; do b=$(basename "$extra" .xyz); [[ " $LIST " == *" $b "* ]] || LIST="$LIST $b"; done

strip() { sed -r 's/\x1b\[[0-9;]*m//g'; }
run_sp() { # bin tag mol method threads mode -> "energy ms"
  local bin=$1 tag=$2 mol=$3 m=$4 th=$5 mode=$6 c=${CHG[$3]:-0}
  local wd=$OUT/$tag/work.$mol.$m.t$th; mkdir -p "$wd"; cp "$OUT/mols/$mol.xyz" "$wd/"
  local s e; s=$(date +%s%N)
  ( cd "$wd" && timeout 3600 "$bin" -sp $mol.xyz -method $m -threads $th -charge $c -no_bmt \
      -dump_gradient "$OUT/$tag/$mol.$m.t$th.$mode.grad" > "$OUT/$tag/$mol.$m.t$th.$mode.log" 2>&1 )
  e=$(date +%s%N)
  echo "$(grep -m1 '^# energy' "$OUT/$tag/$mol.$m.t$th.$mode.grad" 2>/dev/null | awk '{print $3}') $(( (e-s)/1000000 ))"
}
run_md() { # bin tag mol method threads maxtime -> "lastEpot lastEtot ms"
  local bin=$1 tag=$2 mol=$3 m=$4 th=$5 mt=$6 c=${CHG[$3]:-0}
  local wd=$OUT/$tag/md.$mol.$m; mkdir -p "$wd"; cp "$OUT/mols/$mol.xyz" "$wd/"
  local s e; s=$(date +%s%N)
  ( cd "$wd" && timeout 3600 "$bin" -md $mol.xyz -method $m -maxtime $mt -threads $th -charge $c \
      -md.seed 42 -md.no_restart -md.rattle_12 false -md.print_frequency 50 -no_bmt > run.log 2>&1 )
  e=$(date +%s%N)
  echo "$(strip < "$wd/run.log" | grep -E '^\s+[0-9]+\.[0-9]+\s+-?[0-9]+\.[0-9]+\s+' | tail -1 | awk '{print $2, $6}') $(( (e-s)/1000000 ))"
}
cmp_grad() { # gradA gradB -> max|dg|
  python3 - "$1" "$2" <<'EOF'
import sys
def load(p):
    g=[]
    for l in open(p):
        if l.startswith('#'): continue
        t=l.split()
        if len(t)>=3: g.append([float(x) for x in t[:3]])
    return g
a,b=load(sys.argv[1]),load(sys.argv[2])
print("%.2e" % (max(abs(x-y) for ra,rb in zip(a,b) for x,y in zip(ra,rb)) if a and b and len(a)==len(b) else float('nan')))
EOF
}
rm -rf "$OUT/ref" "$OUT/new"; mkdir -p "$OUT/ref" "$OUT/new"
printf "%-34s %18s %12s %10s %8s %8s %7s\n" run energy_ref dE max_dgrad ms_ref ms_new ratio | tee "$OUT/summary.txt"
worst_e=0; worst_g=0
for mol in $LIST; do
  for m in gfnff gfn1 gfn2; do
    for th in 1 8; do
      [ "$mol" = polymer ] && [ "$m" != gfnff ] && [ "$th" = 1 ] && continue
      for mode in cold warm; do
        [ "$mode" = warm ] && [ "$m" != gfnff ] && continue
        read ea ta <<<"$(run_sp "$REF" ref $mol $m $th $mode)"
        read eb tb <<<"$(run_sp "$NEW" new $mol $m $th $mode)"
        de=$(python3 -c "print('%.3e' % (float('$eb')-float('$ea')))" 2>/dev/null || echo nan)
        dg=$(cmp_grad "$OUT/ref/$mol.$m.t$th.$mode.grad" "$OUT/new/$mol.$m.t$th.$mode.grad" 2>/dev/null || echo nan)
        ratio=$(python3 -c "print('%.2f' % ($ta/max($tb,1)))")
        printf "%-34s %18s %12s %10s %8s %8s %7s\n" "$mol.$m.t$th.$mode" "$ea" "$de" "$dg" "$ta" "$tb" "$ratio" | tee -a "$OUT/summary.txt"
      done
    done
  done
done
if [ "$QUICK" != quick ]; then
  for spec in "polymer gfnff 8 100" "caffeine gfn2 8 50" "triose gfnff 1 200"; do
    read mol m th mt <<<"$spec"
    read pa ta ma <<<"$(run_md "$REF" ref $mol $m $th $mt)"
    read pb tb mb <<<"$(run_md "$NEW" new $mol $m $th $mt)"
    same=$([ "$pa $ta" = "$pb $tb" ] && echo identical || echo "DIFFERENT ($pa/$ta vs $pb/$tb)")
    printf "%-34s %18s %12s %10s %8s %8s %7s\n" "$mol.$m.t$th.md${mt}fs" "$pa" "$same" "-" "$ma" "$mb" "$(python3 -c "print('%.2f' % ($ma/max($mb,1)))")" | tee -a "$OUT/summary.txt"
  done
fi
echo "results: $OUT/summary.txt"
