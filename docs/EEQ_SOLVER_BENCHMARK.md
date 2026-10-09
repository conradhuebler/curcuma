# EEQ solver benchmark (GFN-FF)

`scripts/eeq_solver_benchmark.py` — Claude Generated (Oct 2026), 🤖 AI-generated, ⚙️ machine-tested on caffeine and polymer_2x only.

## What it does

For each structure, EEQ solver variant and device (cpu / cuda), in a fresh directory each (no `.topo.json`, topology rebuilt):

- **sp**: curcuma single point, wall time and energy; with `--xtb` also `xtb --gfnff --sp` once per structure (the script sets
  `ulimit -s unlimited`, `OMP_NUM_THREADS=1`) and the energy difference.
- **md**: curcuma MD, per-step phase times from `-md_diagnostics_timing` (median/min/max over steps 1..n, step 0 dropped) and the final Epot.

Output: `results.json` (raw, with commit, threads, ids, sha256) and Markdown tables on stdout.

## Usage

```sh
scripts/eeq_solver_benchmark.py peo201x2-water1500 caffeine --curcuma release/curcuma --xtb /opt/bin/xtb \
    --solvers default,cholesky,ppcg,lu,pcg --devices cpu,cuda --md-fs 50 --threads 1 --out bench_out
```

- Structures: library id or alias from `test_cases/structures` (`polymer_2x` = `peo201x2-water1500`), or a path to an `.xyz`. Charge and spin
  from the manifest are passed on when recorded and non-zero.
- `--modes sp,md`, `--slow-md-fs` (MD length for `pcg`, default 3 fs), `--timeout` per run (default 3600 s), `--extra` for further curcuma flags.

## Solver variants (flag spelling matters)

| variant | flags | note |
|---|---|---|
| `default` | none | = `auto` |
| `auto`, `cholesky`, `ppcg`, `pcg`, `lu`, `ldlt`, `batched` | `-eeq_solver.solve_method <name>` | `batched` drops cross-fragment Coulomb: inexact; `lu`/`ldlt` are errors with `--devices cuda` |

Solver control and naming: [EEQ_SOLVER_CONTROL.md](EEQ_SOLVER_CONTROL.md). The measurement below was taken BEFORE that change, when `cholesky`
and `default` were switched to projected PCG at 500 atoms; its `cholesky_dense` row is today's `cholesky`.

## Result so far (Oct 9, 2026, commit `4da2b26d`, RTX 5080, 1 thread, n = 1 run per configuration)

polymer_2x (7320 atoms), single point: xtb 28.7 s; curcuma default 6.7 s, dense Cholesky 21 s, ldlt/lu 29 s, pcg 240 s; energies equal
to 1.3e-6 Eh of xtb. MD per step (CPU / CUDA, ms): default 2441 / 65, cholesky_dense (= today's cholesky) 5745 / 593, ldlt 7718 / 602, lu 7828 / 596,
pcg 61 786 / 36 309 (2 steps). Raw data and open points: lab journal `Labor/curcuma EEQ-Solver Benchmark.md` in the vault.

Not tested: other molecules, repeats (no spread), multi-GPU, thread scaling. Known open point: `batched` on CUDA reports "WP7-B batched
failed, falling through to WP7-A" and the MD ends after 3 steps with exit code 0 even with `-gpu_strict true`; cause not investigated.
