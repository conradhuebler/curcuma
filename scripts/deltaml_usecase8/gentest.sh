cd /home/conrad/src/curcuma/build_deltaml/usecase8/gentest
U=/home/conrad/src/curcuma/build_deltaml/usecase8
for m in V1/A V1/B2 V2/A V2/B2; do
  t=$(echo $m|tr / _)
  mkdir -p $t; cp heldout_gfn2_s2.xyz $t/in.xyz; cd $t
  $U/curcuma_uc8 -activelearn in.xyz -al_mode test -model $U/models/$m/start.deltaml.json -method gfnff -ref_method gfn2 -threads 8 -no_bmt > log.txt 2>&1; echo "$t rc=$?"
  cd ..
done
