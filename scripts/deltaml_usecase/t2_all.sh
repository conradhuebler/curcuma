cd /home/conrad/src/curcuma/build_deltaml/usecase
xargs -P 3 -L 1 sh -c 'python3 -I scripts/t2_job.py $0 $1 > /dev/null 2>t2/err_$0_$1.txt' < t2/jobs.txt
echo done > t2/t2_done.flag
