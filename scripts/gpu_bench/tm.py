#!/usr/bin/env python3
# Run a command, report wall time, peak host RSS and peak GPU memory (all GPUs, sampled 0.5 s).
import subprocess, sys, time, resource, threading
peak = {}
stop = False
def poll():
    while not stop:
        try:
            out = subprocess.run(["nvidia-smi","--query-gpu=index,memory.used","--format=csv,noheader,nounits"],capture_output=True,text=True).stdout
            for l in out.strip().splitlines():
                i,m = [int(x) for x in l.split(",")]
                peak[i] = max(peak.get(i,0), m)
        except Exception: pass
        time.sleep(0.5)
t=threading.Thread(target=poll); t.start()
t0=time.time()
rc=subprocess.call(sys.argv[1:])
stop=True; t.join()
r=resource.getrusage(resource.RUSAGE_CHILDREN)
print(f"WALL {time.time()-t0:.2f} s MAXRSS {r.ru_maxrss/1024/1024:.2f} GB GPU_PEAK_MiB {peak} RC {rc}")
