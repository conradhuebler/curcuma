#!/bin/bash
cd /home/conrad/src/curcuma/build_deltaml/usecase9
while pgrep -f "xargs -P 3 -L 1 bash scripts/al.sh" >/dev/null; do sleep 20; done
xargs -P 3 -L 1 bash scripts/al.sh < jobs_al2.txt
rm -f jobs_md.txt
for d in N0_s1 N0_s2 N0_s3 N1_s1 N1_s2 N1_s3 N2_s1 N2_s2 N3_s1 N3_s2 N3_s3; do for s in 1 2 3 4; do echo "$d $PWD/al/$d/start.deltaml.json $s" >> jobs_md.txt; done; done
xargs -P 3 -L 1 bash scripts/run_md9.sh < jobs_md.txt
echo "$(date +%H:%M:%S) ALL MD done" >> STATUS.txt
