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
| jl22-bgal.conf6 | 173 | gfnff, gfn2 |
| peo201-chain | 1410 | gfnff, gfn2 |
| peo201-water500 | 2910 | gfnff, gfn2 (n = 3, timeout 1800 s) |
| peo201x2-water1500 | 7320 | gfnff, gfn2 (n = 1, timeout 3600 s, a trial: memory and time are not known beforehand) |

The inputs are in `test_cases/structures/`. Changing the list means changing the table.

## How to add a block

```
python3 scripts/perf_table.py --curcuma <build>/curcuma --backends cpu,rocm --repeats 5 > block.md
```

`--backends cpu,cuda` on a CUDA machine, `--quick` runs two systems as a check of the setup, `--json file` keeps the single timings. Run it on an idle machine and with the canonical build (`release/`, `CLAUDE.md`). Paste the block under "Blocks".

## Blocks

#### 2026-10-30: Intel(R) Core(TM) Ultra 5 225, AMD Radeon RX 9060 XT (gfx1200)

- Binary `/home/conrad/src/curcuma/release/curcuma` (built 2026-10-30 16:43, reports commit `77f033ca`); repository HEAD `77f033ca`, 2 uncommitted files (this script and this document, not the program), 0 source files newer than the binary.
- 5 repeats per cell unless the n column says otherwise, backends alternating, `-threads 10`, RAM 62 GB, load average at start 0.39 3.06 2.15, at end 1.00 1.50 2.37.
- Wall time of the whole process (start, plugin load, device initialisation included), `-sp`, fresh directory per run, GPU runs with `-gpu_strict true`. No GPU energy differed from the CPU value by more than 1e-6 Eh (no cell is marked).
- peo201x2-water1500 gfn2 was a trial with one repeat and a 24 GB memory cap; both backends were killed at the cap, so there is no time for that cell.

| system | atoms | method | backend | median s | min s | max s | n | speed-up vs CPU | note |
|---|---|---|---|---|---|---|---|---|---|
| caffeine | 24 | gfnff | cpu | 0.03 | 0.03 | 0.03 | 5 |  |  |
| caffeine | 24 | gfnff | rocm | 1.13 | 1.08 | 1.15 | 5 | 0.02 |  |
| caffeine | 24 | gfn2 | cpu | 0.03 | 0.03 | 0.04 | 5 |  |  |
| caffeine | 24 | gfn2 | rocm | 1.37 | 1.33 | 1.40 | 5 | 0.03 |  |
| jl22-bgal.conf6 | 173 | gfnff | cpu | 0.04 | 0.04 | 0.05 | 5 |  |  |
| jl22-bgal.conf6 | 173 | gfnff | rocm | 1.12 | 1.07 | 1.16 | 5 | 0.04 |  |
| jl22-bgal.conf6 | 173 | gfn2 | cpu | 0.22 | 0.21 | 0.23 | 5 |  |  |
| jl22-bgal.conf6 | 173 | gfn2 | rocm | 1.69 | 1.47 | 1.92 | 5 | 0.13 |  |
| peo201-chain | 1410 | gfnff | cpu | 0.14 | 0.13 | 0.15 | 5 |  |  |
| peo201-chain | 1410 | gfnff | rocm | 1.09 | 0.93 | 1.22 | 5 | 0.13 |  |
| peo201-chain | 1410 | gfn2 | cpu | 18.27 | 17.86 | 19.94 | 5 |  |  |
| peo201-chain | 1410 | gfn2 | rocm | 12.52 | 12.37 | 13.42 | 5 | 1.46 |  |
| peo201-water500 | 2910 | gfnff | cpu | 0.36 | 0.35 | 0.36 | 5 |  |  |
| peo201-water500 | 2910 | gfnff | rocm | 1.07 | 0.70 | 1.47 | 5 | 0.33 |  |
| peo201-water500 | 2910 | gfn2 | cpu | 135.81 | 135.51 | 136.33 | 3 |  |  |
| peo201-water500 | 2910 | gfn2 | rocm | 68.24 | 68.04 | 69.08 | 3 | 1.99 |  |
| peo201x2-water1500 | 7320 | gfnff | cpu | 1.59 | 1.58 | 1.59 | 5 |  |  |
| peo201x2-water1500 | 7320 | gfnff | rocm | 1.59 | 1.58 | 2.60 | 5 | 1.00 |  |
| peo201x2-water1500 | 7320 | gfn2 | cpu | - | - | - | 0 | - | killed at the 24 GB memory cap (kernel OOM killer in the journal) |
| peo201x2-water1500 | 7320 | gfn2 | rocm | - | - | - | 0 | - | killed at the 24 GB memory cap (kernel OOM killer in the journal) |

Reading this block (what the numbers show, not why):

- Every ROCm cell up to 2910 atoms with GFN-FF sits at about 1.1 s whatever the size (1.07 to 1.13 s median), so the time there is the fixed cost of the process with the GPU backend (plugin load, device initialisation), not the calculation. A single point on a small system is the wrong measure of GPU speed.
- GFN2: ROCm is slower at 24 and 173 atoms (0.03 and 0.13 times the CPU speed), 1.46 times faster at 1410 atoms and 1.99 times faster at 2910 atoms (n = 3 there).
- GFN-FF at 7320 atoms: CPU and ROCm both 1.59 s median. With `-verbosity 2` the ROCm run reports "using GPU acceleration (rocm)" and solves the EEQ on the CPU (ppcg, N = 7320, 1502 fragments), so the single point includes CPU steps; how much of the 1.59 s these take was not measured.
- One ROCm GFN-FF run at 7320 atoms took 2.60 s against a median of 1.59 s (max column); the other four were close to the median. The cause was not examined.

## Not covered

MD and optimisation step times, memory use, power, clock behaviour, multi-GPU, the first call separately from later calls.
