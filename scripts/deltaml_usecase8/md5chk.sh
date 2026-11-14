cd /home/conrad/src/curcuma/build_deltaml/usecase8/runs
for t in $(ls | sed 's/_s[0-9]$//' | sort -u); do echo "$t $(md5sum ${t}_s?/md.snapshots/md.trj.xyz | awk '{print $1}' | sort -u | wc -l)"; done | awk '{printf "%s:%s ",$1,$2}'; echo
grep -h "rc=" ../STATUS.txt | grep -c "MD done.*rc=0"
