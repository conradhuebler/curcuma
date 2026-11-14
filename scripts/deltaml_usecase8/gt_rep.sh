cd /home/conrad/src/curcuma/build_deltaml/usecase8/gentest
for t in V1_A V1_B2 V2_A V2_B2; do echo "== $t"; sed 's/\x1b\[[0-9;]*m//g' $t/log.txt | grep -v CITE | grep -i "MAE\|RMSE\|rank\|structures\|warn\|error" | cut -c1-160 | tail -7; done
