#!/bin/bash
# usage: label_one.sh <Mx/fr/NNN>
cd $1 || exit 1
B=/home/conrad/src/curcuma/build_deltaml/usecase3/curcuma_uc3
mkdir -p base ref
cp s.xyz base/x.xyz; cp s.xyz ref/x.xyz
(cd base && $B -sp x.xyz -method gfnff -gfnff.cache_topology false -dump_gradient g_base.txt -threads 1 -no_bmt > log.txt 2>&1; echo $? > rc.txt)
(cd ref && $B -sp x.xyz -method gfn2 -dump_gradient g_ref.txt -threads 1 -no_bmt > log.txt 2>&1; echo $? > rc.txt)
