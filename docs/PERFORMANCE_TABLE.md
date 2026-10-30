# Performance table per machine, backend and method

One block per machine and per state of the code, made by `scripts/perf_table.py` and pasted below. Purpose: make gaps between backends visible (CPU, CUDA, ROCm) and see what a code change did on a given machine. Repeat the block after a hardware change or after a change to the GPU or SCF paths.

## What a row means

- Wall time of the whole process `curcuma -sp <system> -method <m> [-gpu <backend>] -threads <N> -no_bmt`. It includes start, plugin load and device initialisation, so small systems are dominated by the fixed cost of the GPU backends. It says nothing about the time per MD or optimisation step.
- n repeats per cell (default 5), backends alternating inside every repeat, one process at a time, a fresh directory per run (GFN-FF reuses `<name>.topo.json`), memory cap 24 GB through `systemd-run` where available. GPU runs use `-gpu_strict true`, so a silent CPU fallback fails the run instead of being timed as GPU.
- The columns are median, minimum and maximum of the n wall times; the speed-up is the CPU median divided by the backend median. A cell whose energy differs from the CPU value by more than 1e-6 Eh is marked.
- Each block names the binary, the commit the binary reports, the repository HEAD, the number of uncommitted files and the number of source files newer than the binary. Read those first: a block is only as current as its binary.

## Fixed system list

| system | atoms | methods |
|---|---|---|
| caffeine | 24 | gfnff, gfn2 |
| trisaccharide | 66 | gfnff, gfn2 |
| aaa-bglc.conf76 | 117 | gfnff, gfn2 |
| macrocycle-bgal | 231 | gfnff, gfn2 |
| peo201-chain | 1410 | gfnff, gfn2 |
| peo201-water500 | 2910 | gfnff |
| urea400-water1000 | 6200 | gfnff |
| peo201x2-water1500 | 7320 | gfnff |

The inputs are in `test_cases/structures/`. Changing the list means changing the table.

## How to add a block

```
python3 scripts/perf_table.py --curcuma <build>/curcuma --backends cpu,rocm --repeats 5 > block.md
```

`--backends cpu,cuda` on a CUDA machine, `--quick` runs three systems as a check of the setup, `--json file` keeps the single timings. Run it on an idle machine and with the canonical build (`release/`, `CLAUDE.md`). Paste the block under "Blocks".

## Blocks

No block recorded yet. The script ran once on the machine with an Intel Core Ultra 5 225 and an AMD Radeon RX 9060 XT (gfx1200), on caffeine only and with an older binary, as a check of the setup; that run is not an entry.

## Not covered

MD and optimisation step times, memory use, power, clock behaviour, multi-GPU, the first call separately from later calls.
