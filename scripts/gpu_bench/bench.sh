#!/bin/bash
# usage: bench.sh xyz method gpu threads tag
C=${CURCUMA:-/home/conrad/src/curcuma/release/curcuma}
log=$5_$(basename $1 .xyz)_$2_$3_t$4.log
rm -f $(basename $1 .xyz).topo.json
$(dirname $0)/tm.py $C -sp $1 -method $2 -gpu $3 -threads $4 -no_bmt -verbosity 2 > $log 2>&1
echo "== $log"; grep -a -E "Total  |TOTAL|SCF \(|setup  |WALL|fallback|GPU context|device-resident loop|gap =|not converged|Final Energy" $log | sed 's/\x1b\[[0-9;]*m//g' | head -12
