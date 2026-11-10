# label.sh dir ref_method [frames.xyz]  : internal labelling, input order
cd $1
/home/conrad/src/curcuma/build_deltaml/usecase5/curcuma_uc5 -activelearn ${3:-frames.xyz} -method gfnff -ref_method $2 -test_fraction 0 -threads 8 -no_bmt > label.log 2>&1
echo "label done rc=$?" > label.done
