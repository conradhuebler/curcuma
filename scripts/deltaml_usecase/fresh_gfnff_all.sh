#!/bin/bash
# fresh GFN-FF energies (no cache) for every structure of the dataset, 8 jobs in parallel
cd /home/conrad/src/curcuma/build_deltaml/usecase/chk
python3 -I - <<'PY'
import json,os
d=json.load(open('/home/conrad/src/curcuma/build_deltaml/ablation_run/conf.deltaml.dataset.json'))
sym={1:'H',6:'C',7:'N',8:'O'}
os.makedirs('all',exist_ok=True)
for k,s in enumerate(d['structures']):
    el=s['elements']; g=s['geometry']; flat=[x for r in g for x in (r if isinstance(r,list) else [r])]
    os.makedirs(f'all/{k:03d}',exist_ok=True)
    with open(f'all/{k:03d}/s.xyz','w') as f:
        f.write(f"{len(el)}\n\n")
        for i,e in enumerate(el): f.write(f"{sym.get(e,e)} {flat[3*i]:.10f} {flat[3*i+1]:.10f} {flat[3*i+2]:.10f}\n")
PY
run_one() { k=$1; cd all/$k; ../../../curcuma_uc -sp s.xyz -method gfnff -gfnff.cache_topology false -threads 1 -no_bmt 2>&1 | sed 's/\x1b\[[0-9;]*m//g' | grep -m1 'Single Point Energy' | awk -v k=$k '{print k, $5}' > e.txt; }
export -f run_one
ls all | xargs -P 8 -I{} bash -c 'run_one {}'
cat all/*/e.txt > fresh_gfnff.txt
echo done > fresh_done.flag
