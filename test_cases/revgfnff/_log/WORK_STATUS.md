# WORK_STATUS — rev-gfnff work packages 1-5 (2026-09-18)
Packages done: 2/5

AI-generated, machine-tested. Repository `/home/conrad/src/curcuma_branches/curcuma`, branch
`reactff2-llm`, start HEAD `265a18b0`. Every measurement was taken with a FROZEN copy of the
binary; the copy's md5 stands next to each number and in every run directory's `wall.txt`.
Harness in `.../scratchpad/wk/h/` (`run_one.sh`, `seq.sh`, `grid_an.py`, `fals.sh`, `pathB.py`,
`fdchk3.py`, `adduct.py`, `joblist20.txt`).

**The grid has 20 live cells, not 22.** `ch4_H.xyz` holds 15 frames, so `ch4_H/T1000_f16` and
`ch4_H/T2000_f16` exit with rc = 1 and an empty log (confirmed again here). Both are dropped from
`joblist20.txt`; they contributed nothing, which the control arm below proves (it reproduces the
22-cell records exactly).

**Rebuild counting**: `max(REACT rebuild #N)` = the number of distinct rebuilds. `grep -c` gives
1.5x that (the line is printed twice for some rebuilds from calculator verbosity 2 on) — measured
here: c2h6/T2000_f16 default arm, max index 136, `grep -c` 204. The "204 / 126" fingerprints in
`HBUDGET_STATUS` / `FABLE_REVIEW_2` are `grep -c` numbers of the same runs as the 136 / 84 below.

---

## Package 1 — `rev_budget_fix_h` is the default, and the new baseline

### 1.1 The change (commit "rev-gfnff 3a(ii): make rev_budget_fix_h the default")

Four places write the default; all four flipped to `true`:
`gfnff.h` PARAM (line 363, help text rewritten), `gfnff.h` member `m_rev_budget_fix_h`,
`gfnff_method.cpp` `setupRevSettings()` fallback `m_parameters.value("rev_budget_fix_h", …)`,
`ff_workspace.h` `RevSettings::budget_fix_h`. `make GenerateParams` clean (exit 0, no warning),
`make -C build_rev -j8 curcuma` MAKE_EXIT=0.

Verified with `CURCUMA_REVDUMP=1`: no flag -> `budget_fix_h true`; `-gfnff.rev_budget_fix_h false`
-> `false`. **Binary md5 `5bfccf51c67366a8818192eb2ba8d984`** ("P1"), frozen as `scratchpad/wk/cur_p1`;
every number below is from that copy.

### 1.2 The 20-cell grid — new baseline (default) and control (`-gfnff.rev_budget_fix_h false`)

Protocol = the `hb/runs_{OFF,ON}` protocol: `-md -method revgfnff -gfnff.topology_mode react
-temperature {1000,2000} -maxtime 5000 -md.time_step 0.25 -md.thermostat csvr -md.coupling 10
-md.rattle_12 false -md.no_restart -md.seed 42 -threads 1 -verbosity 3 -md.print_frequency 1`,
fresh directory per cell, `*.topo.json` removed. Aggregate over the 20 cells:

| arm | rebuilds | (a) step max / kJ | (a) n >= 50 | (b) hard swaps | (c) jump max / kJ | (c) n >= 50 | jump median / p99 | T_max / K | wall / s |
|---|---:|---:|---:|---:|---:|---:|---|---:|---:|
| **default (fix_h on) = NEW BASELINE** | **1186** | **391.44** | **487** | **0 of 591** | **48.6** | **0** | 0.0 / 1.7 | **8 306** | 17.3 |
| control (fix_h off = pre-Sep-18) | 960 | 2593.60 | 324 | 5 of 478 | 471.1 | 3 | 0.0 / 1.8 | 1.594e8 | 17.0 |

(a) = max |Epot(t+dt) - Epot(t)| over adjacent printed steps, intervals containing a
`REACT rebuild` line excluded. Excluding in addition the t = 0.00025 ps initial step (the
`HBUDGET_STATUS` section 3 convention): default 222.53 / 485 events, control 2593.60 / 322.
n intervals: 398 817 (default) / 399 047 (control).

**The control reproduces the record exactly**: 960 rebuilds, max dE_jump 471.1, 3 events >= 50 kJ,
5 of 478 hard swaps, T_max 1.594e8 K, per-step max 2593.60 @3.20825 ps, 322 events >= 50 outside
the init step — i.e. `QP_STATUS` A.2 / `HBUDGET_STATUS` section 3 to the last digit, restricted to
20 cells. The two dropped cells therefore contributed nothing, as expected.

Per cell (def = new default, ctrl = control; `=` marks a cell that is bit-identical in both arms):

| cell | reb def/ctrl | step max def/ctrl | n>=50 def/ctrl | jump max def/ctrl | hard def/ctrl | T_max def/ctrl |
|---|---|---|---|---|---|---|
| c2h6/T1000_f0 = | 0/0 | 22.9/22.9 | 0/0 | 0.0/0.0 | 0-0/0-0 | 2444/2444 |
| c2h6/T1000_f8 = | 0/0 | 26.0/26.0 | 0/0 | 0.0/0.0 | 0-0/0-0 | 2469/2469 |
| c2h6/T1000_f16 = | 0/0 | 31.8/31.8 | 0/0 | 0.0/0.0 | 0-0/0-0 | 2622/2622 |
| c2h6/T2000_f0 | 86/36 | 57.4/61.3 | 8/6 | 0.3/0.2 | 0-43/0-18 | 5016/5254 |
| c2h6/T2000_f8 | 66/98 | 60.2/75.9 | 8/39 | 0.2/0.6 | 0-33/0-49 | 5346/5142 |
| **c2h6/T2000_f16** | **136/90** | **59.3/2593.6** | 17/22 | **0.6/471.1** | **0-68/1-45** | **5086/6.21e4** |
| ch3nh2/T1000_f0 | 10/26 | 31.2/29.4 | 0/0 | 0.2/0.3 | 0-5/0-13 | 3183/3059 |
| ch3nh2/T1000_f8 | 16/20 | 31.1/31.1 | 0/0 | 0.4/0.1 | 0-8/0-10 | 2752/2882 |
| ch3nh2/T1000_f16 | 20/30 | 26.7/33.6 | 0/0 | 0.2/0.4 | 0-10/0-15 | 2922/2744 |
| ch3nh2/T2000_f0 | 288/174 | 53.1/49.0 | 3/0 | 3.6/3.1 | 0-144/0-87 | 6181/5396 |
| ch3nh2/T2000_f8 | 230/222 | 71.4/71.3 | 7/15 | 16.8/2.1 | 0-115/0-111 | 5785/5632 |
| ch3nh2/T2000_f16 | 244/240 | 67.5/70.0 | 6/22 | 18.5/4.0 | 0-122/0-120 | 5384/5613 |
| ch4_H/T1000_f0 = | 0/0 | 10.6/10.6 | 0/0 | 0.0/0.0 | 0-0/0-0 | 2835/2835 |
| ch4_H/T1000_f5 = | 0/0 | 19.2/19.2 | 0/0 | 0.0/0.0 | 0-0/0-0 | 3018/3018 |
| ch4_H/T1000_f8 = | 2/2 | 15.7/15.7 | 0/0 | 0.0/0.0 | 0-0/0-0 | 3005/3005 |
| **ch4_H/T1000_f10 =** | 0/0 | **388.4/388.4** | **40/40** | 0.0/0.0 | 0-0/0-0 | 6812/6812 |
| ch4_H/T2000_f0 = | 0/0 | 33.1/33.1 | 0/0 | 0.0/0.0 | 0-0/0-0 | 5923/5923 |
| ch4_H/T2000_f5 = | 2/2 | 39.2/39.2 | 0/0 | 0.0/0.0 | 0-1/0-1 | 5850/5850 |
| ch4_H/T2000_f8 = | 2/2 | 45.3/45.3 | 0/0 | 0.0/0.0 | 0-0/0-0 | 6437/6437 |
| **ch4_H/T2000_f10** | **84/18** | **391.4/1089.0** | 398/180 | **48.6/310.4** | **0-42/4-9** | **8306/1.594e8** |

10 of the 20 cells are bit-identical in both arms. The two runaway cells are exactly the two the
review named. **The default arm's remaining (a) events are the CARBON budget, not the hydrogen**:
438 of the 487 sit in the two `ch4_H/*_f10` cells, and `ch4_H/T1000_f10` is bit-identical in both
arms (388.4 kJ, 40 events) — i.e. untouched by this package and the target of package 3.

### 1.3 Falsifiers — all reproduce, nothing moves

| falsifier | result |
|---|---|
| NH4+ / H3O+ / CH5+ / ClO4- / BF4- 1.394 A / BF4- 1.143 A | default = control to 12 digits (0.827303406347 / 1.207219902415 / 0.782692121320 / 0.039597768947 / -1.470184046880 / 0.118269033480 Eh), dE = 0.00000000 kcal, dGmax diff 0.00e+00, n = 6 |
| equilibrium toggle set (share x fix_h) x 5 systems | 20/20 dE = +0.000000000 kcal, dGmax 0.00e+00; caffeine -4.673521653477, benzene -2.363224128930, 2h2 -0.324422067632, n2_3h2 -0.928844746252, ch4_H -0.631216645034 |
| rkt06_h_h2, 11 points | rms(vs reactant) **2.7140** (record 2.714), barrier cur +3.41 @pt4 vs ref +2.57 @pt10; max abs(E_default - E_control) = **0.000e+00 Eh** |
| FD gradient, rkt06 point 10 (central, dx 1e-4 A, fresh dir per displacement) | worst abs(analytic - FD) = **1.194e-08 Eh/A** (record 1.194e-08) |
| `gfnff` identity | caffeine **-4.672737068614**, benzene **-2.362725526194**, `-dump_params` md5 **d297bc3b91ad75224f224b0cb4c3d189** / **77c134bbfd37284059c15d89c6125efc** = the records |

### 1.4 NEW falsifier — the four class-C radical-adduct approach scans

`test_cases/revgfnff/ref/C/{ch4_H,nh3_H,h2o_H,n2h4_H}`, 15 rigid points each, fresh SP per point,
reported as model MINUS reference in kcal/mol, both curves referenced to their own d = 3.00 A point
(`h/adduct.py`). This is the package-3 target; the numbers below are the state to be improved.

| scan | dev at 1.0 / 1.2 / 1.3 A | dev min | dev max | rms over 15 pts |
|---|---|---:|---:|---:|
| CH4 + H | -87.0 / -78.0 / -68.6 | **-87.0** | +32.2 | 41.8 |
| NH3 + H | -107.0 / -91.9 / -58.2 | **-107.0** | +1.0 | 47.1 |
| H2O + H | -54.1 / -53.7 / +5.4 | **-54.6** | +5.4 | 24.3 |
| N2H4 + H | -90.4 / -82.3 / -52.5 | **-90.4** | +2.1 | 41.0 |

This reproduces `FABLE_REVIEW_2` A.4 ("54-107 kcal/mol below the reference", n = 4) exactly and
independently of the review's offline arithmetic. **The hydrogen fix does not touch it**: the
control arm gives the identical four rows to the printed digit, as it must — the budget that
grows here is the carbon's / nitrogen's / oxygen's.

### 1.5 Tests

`ctest -R gfnff` **62/65**, `ctest -R cli_simplemd_` **19/22** — the same three failures in both,
and they are the known set:

| test | failure |
|---|---|
| `cli_simplemd_16_gfnff_rev_h2_smooth` | 578 events, median 0.0, max 0.5 kJ/mol (both inside their thresholds); fails only on `T_max 7946 K > 6400` |
| `cli_simplemd_18_gfnff_rev_nve_vs_gfnff` | rev NVE slope 3.42e-3 (dt 0.25) / 6.18e-3 (dt 0.125) vs threshold 1.6e-3; gfnff arm fine (-9.5e-6 / 2.7e-5); T_max 11924 K, 1030 rebuilds |
| `cli_simplemd_19_gfnff_rev_form_refuses_hbond` | mean Epot off by -0.3995 kcal/mol from the static run (tolerance 0.01) |

**Not caused by this flip**: `HBUDGET_STATUS` section 5 recorded the same three failing on a binary
whose default was OFF. They are the stage-3a calibration debt and are recalibrated in package 5.
`cli_simplemd_08/09` (which failed in that worktree) pass here.

### 1.6 Open after package 1

- The carbon/nitrogen/oxygen budget: 438 of the 487 remaining >= 50 kJ steps, and the whole class-C
  adduct error of 1.4. Package 3.
- The 49 events of the other 18 cells have no budget change and sit inside the thermal amplitude
  (`FABLE_REVIEW_2` A.3); a 50 kJ threshold at 5000-7000 K and dt = 0.25 fs is at the thermal
  ceiling. Scaling the threshold with 3N k T was suggested and is NOT done here.

---

## Package 2 — merge `origin/feature/multi-gpu`

Two commits: the merge itself with `coulomb_implicit` pinned **false**, then the flip to the
remote's default **true**. Binaries: P1 `5bfccf51` -> merge `e88c25235dbee736c2ca3d6bd420d0ae`
-> flip `6c0b0db0f1b422a7327163eaee689cbf` (three distinct md5s, so every "did not change"
below is on a proven rebuild).

### 2.1 The five conflicts

| file | resolution |
|---|---|
| `gfnff_method.cpp` | kept our `perceiveGeometricBonds()` call, moved the remote's row-parallel loop INTO that function. One implementation; the list stays (i, j) ascending because react and the topology-reuse check compare bond graphs. |
| `gfnff_gpu_method_impl.h` | both accessors kept (`getGFNFF`/`gfnffInstance` and `gpuDevice`). |
| `src/main.cpp` | our explicit `-batch true` block first (it returns), the remote's automatic multi-frame batch after it. |
| `.gitignore`, `AIChangelog.md` | both sides kept. |

**The `-sp` semantic change, checked rather than assumed**: a multi-frame `-sp` WITHOUT
`-batch` used to evaluate the first structure only and now evaluates all frames on `-threads`
workers. Every rev harness that passes a multi-frame file uses `-batch true`
(`revgfnff_classa/curves/contact/fit`); the ones on the plain path
(`revgfnff_barrier_terms.py`, `pathB.py`, `adduct.py`, `scan.py`, the ctest scripts) all write a
single structure. No script changes meaning.

**One build break in the remote, fixed here**: `XTB::reduceToStandardForm()` (065e2a29) called
`dsygst_` unconditionally, but that symbol is declared only inside the
`EIGEN_USE_BLAS || USE_BLAS || USE_MKL` guard — so a BLAS-less configuration did not compile,
and `build_rev` is one. The call now sits inside the same guard with the documented
triangular-solve route as the fallback (the two agree to 8e-15 elementwise per the function's
own comment, and both call sites are themselves inside the LAPACK guard, so the BLAS build is
untouched and the BLAS-less build never reaches the LAPACK branch).

### 2.2 Identity of the merge (`coulomb_implicit false`), all at `-threads 1`

| check | result |
|---|---|
| 20-cell grid, per cell | **identical in every field** (n_rebuild, step_max, n>=50, jump_max, n_hard/n_begin, T_max): 0 of 160 compared fields differ |
| 20-cell grid, aggregate | 1186 rebuilds / step max 391.44 / 487 events / 0 of 591 hard swaps / jump max 48.6 / T_max 8306 K = package 1 |
| grid logs, byte level | 143 799 status rows differ, and **only in columns 13 and 14** — `remaining` (wall-clock estimate) and accumulated wall ms/1000. All 13 physical columns bit-identical. |
| caffeine/benzene 12 digits, `revgfnff` + `gfnff` | bit-identical, `dump_params` md5 d297bc3b / 77c134bb unchanged |
| 6 hypervalent ions + BF4-, equilibrium toggle set | bit-identical |
| rkt06 | rms **2.7140**, max abs(dE) over 11 points **0.000e+00 Eh** |
| FD gradient at rkt06 pt 10 | **1.194e-08 Eh/A** |
| 4 class-C adduct scans | max difference **0.000e+00 kcal/mol** |
| class-A `ch4_C-H`, `--mode kept` + `-gfnff.topology_mode react` | rms **5.87**, D_e dev **-13.32**, r90 dev **-0.351** = the `BASELINE_HEAD` record (this is the `-batch` path through the merged `main.cpp`) |
| `ctest -R "gfnff|sqm_val|react"` | 95/98 |
| `ctest -R cli_simplemd_` | 19/22 — the same three known failures |

**Harness note**: `--mode kept` alone is NOT the `BASELINE_HEAD` protocol. `revgfnff_classa.py`
needs `--extra "-gfnff.topology_mode react"` as well; without it `ch4_C-H` reads 5.23 / -11.43 /
-0.335 instead of 5.87 / -13.32 / -0.351 (consistent with `BASELINE_HEAD` section 1c, which
records that the literal command without react differs on 30/32 bonds).

### 2.3 `-threads 4` — the merge does not change the thread dependence

| cell | T1 | T4 before merge | T4 after merge | T4 with `coulomb_implicit true` |
|---|---|---|---|---|
| c2h6/T2000_f16 | 136 reb / 59.26 kJ / 17 | 114 / 65.77 / 20 | **114 / 65.77 / 20** | 138 / 62.19 / 5 |
| ch4_H/T2000_f10 | 84 / 391.44 / 398 | 76 / 391.44 / 458 | **76 / 391.44 / 458** | 80 / 391.44 / 365 |

So: `-threads 4` already differed from `-threads 1` before the merge (pre-existing, Known Issue
#33 class), the merge leaves the T4 trajectories bit-identical to the pre-merge ones (again only
the two wall-clock columns differ), and the `coulomb_implicit` flip is what moves them.

### 2.4 The `coulomb_implicit true` flip, measured alone

At `-threads 1` **nothing** moves: grid identical cell by cell, hypervalent ions, equilibrium
set, 12-digit energies, rkt06 (2.7140, max dE 0.000e+00), FD gradient 1.194e-08, the four adduct
scans (0.000e+00) and class-A ch4_C-H (5.87 / -13.32 / -0.351) all unchanged; `ctest` unchanged.

Two things move, both predicted by `FABLE_REVIEW_2` section C:
1. **`-gfnff.dump_params` md5, by construction** — the Coulomb pair list is no longer built, so
   it is no longer dumped: caffeine **d297bc3b -> 4013d6fcd3398eb317af0543a0ea90d5**, benzene
   **77c134bb -> 6c3a87c8efa39c6c7b32718dd70065b2**. The energies behind them are bit-identical.
   **These are the yardstick md5s from here on.**
2. **The last ulp at `-threads` > 1** — the partition is by atom range instead of pair range, so
   the reduction order changes; a react MD trajectory amplifies it (table above).

### 2.5 GPU: merged, NOT verified

`build_rev/CMakeCache.txt` has `USE_CUDA:BOOL=OFF`, `USE_ROCM:BOOL=OFF`, `USE_VULKAN:BOOL=OFF`.
None of the merged GPU code is compiled in this build directory, so every GPU claim of the
remote branch is carried over unverified.
