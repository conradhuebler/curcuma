cd /home/conrad/src/curcuma/build_deltaml/usecase/t2
(cd A && ../../curcuma_uc -activelearn train.xyz -method gfnff -ref_method gfn2 -activelearn.dataset $PWD/train.dataset.json -test_fraction 0 -activelearn.seed 1 -threads 8 -no_bmt > log.txt 2>&1) &
(cd B && ../../curcuma_uc -activelearn train.xyz -method gfnff -ref_method gfn2 -activelearn.dataset $PWD/train.dataset.json -test_fraction 0 -activelearn.seed 1 -threads 8 -no_bmt -activelearn.n_centers 4000 -activelearn.lambda 1e-7 -activelearn.acsf_rcut_radial 7 -activelearn.acsf_zeta 1,2,4,8,16 > log.txt 2>&1) &
wait
echo done > t2_train_done.flag
