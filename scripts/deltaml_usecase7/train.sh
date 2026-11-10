cd /home/conrad/src/curcuma/build_deltaml/usecase7
DS=/home/conrad/src/curcuma/build_deltaml/usecase5/peo/md0/frames.deltaml.dataset.json
for t in A B2; do
 (cd models/$t; X=""; [ $t = B2 ] && X="-activelearn.n_centers 4000"
  s=$(date +%s.%N)
  ../../curcuma_uc7 -activelearn start.xyz -method gfnff -ref_method gfn2 -activelearn.dataset $DS -test_fraction 0 -activelearn.seed 1 -threads 8 -no_bmt $X > log.txt 2>&1
  echo "rc=$? wall=$(echo "$(date +%s.%N)-$s"|bc)" > time.txt) &
done
wait
