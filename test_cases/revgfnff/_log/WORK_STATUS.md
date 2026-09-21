# WORK_STATUS — rev-gfnff work packages 1-11 (2026-09-18 / 21)
Packages done: 11/11 (9 = 9a + 9b)

> **Package 10 re-scales the time axis of everything before it.** The MD time step was integrated
> in the wrong unit until Sep 2026, so every "fs" and every duration written by packages 1-9 is
> **1.9516144x** larger than stated. Multiply before quoting any of them. Package 10 re-derives the
> react-MD dt characterization and the mg/mg2/mg3 smoothness comparison in true femtoseconds.

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

---

## Package 3 — the valence-conserving share with a charge-granted excess budget

New mode `-gfnff.rev_share_form delivered|conserving` (+ `-gfnff.rev_share_min_width`, default
0.1). **Default unchanged (`delivered`) and inert**: the 20-cell grid is identical to the
package-2 binary cell by cell (only the two wall-clock columns of the status rows differ), the
equilibrium set and the gfnff yardsticks are bit-identical. Binary `01e9b4e6862f8558a28a0f51e67f5e42`.

### 3.1 What `conserving` is

    S_i   = sum_k w_ik g_ik                         (the atom's whole claim, as before)
    X_i   = 0 (H, F) | 1 (group 13) | 6 - Val_Z (period >= 3, groups 15-17) | clip(Q_i)
    Val_i = Val_Z + min(G(S_i - Val_Z), X_i)
    f_i   = min(1, Val_i / S_i)                     PER ATOM
    c_ij  = 1 - g_ij (1 - f_i f_j)                  = f_i f_j for an ordinary pair

`Q_i` = the topological (phase-1 EEQ) charge of atom i plus that of its H partners in this
corner — a per-corner constant, so no chain rule runs through it. The delivered rule is a
LEFT-OVER rule and forfeits ALL of an atom's valence at one partner too many; this one spreads
it (`sum_j f_i w_ij = min(Val_i, S_i)` exactly). The PRODUCT, not the mean: A.5 measured that the
mean gives the rkt06 TS pair c = 0.75 and breaks the path.

**Two implementation points that matter**, both measured rather than assumed:
- The budget is built from the same **wide** `S_i` the share divides by, not from the settled
  count. That is what makes NH4+ come out at `Val = S` exactly: its residual is **0.0013
  kcal/mol**, not the +0.17 the review predicted from the settled-count budget.
- The smooth min (`shareMinOne`, a cubic join) is **exact on both sides**: literally 1 for
  `Val >= S` (so an equilibrium atom's term stays bit-identical) and literally `Val/S` below
  `1 - width` (so the conservation identity holds where the share bites). A softplus min
  satisfies neither — it is off by ln(2)/beta at x = 1.
- **`GFNFFParameters::periodic_group` uses MAIN-GROUP numbering 1-8, not IUPAC 1-18.** The
  review's "group 13 / groups 15-17" are 3 / 5-7 in that table. With the IUPAC numbers no
  element rule ever fired and ClO4- cost **+164 kcal/mol**; caught by the new `shareA` line of
  `CURCUMA_SHAREDUMP` (it prints Z, S, Val_Z, cap, Val, f, df/dS per atom).

### 3.2 Falsifiers, both modes side by side (same binary)

| falsifier | delivered | conserving | difference |
|---|---|---|---|
| NH4+ | 0.827303406347 | 0.827301258361 | **-0.00135 kcal/mol** |
| H3O+ | 1.207219902415 | 1.207219745650 | -0.0001 |
| CH5+ | 0.782692121320 | 0.781425328593 | **-0.795** |
| ClO4- | 0.039597768947 | 0.039597767376 | -0.000001 |
| BF4- 1.394 A | -1.470184046880 | -1.470184050200 | -0.000002 |
| BF4- 1.143 A (compressed) | 0.118269033480 | 0.147985923798 | **+18.6** |
| caffeine / benzene / 2h2 / n2_3h2 / ch4_H (equilibria) | — | — | **bit-identical, 0.000000000** |
| `gfnff` caffeine / benzene + dump_params md5 | — | — | untouched |
| rkt06_h_h2, 11 points | rms **2.7140** | rms **2.7614** | barrier +3.41 @pt4 in both |
| rkt03 (H + CH4 -> CH3 + H2, class P, 14 points) | rms 19.55 | rms 19.55 | max abs(dE) 3.6e-5 Eh over the path |

The BF4--compressed **+18.6** reproduces the review's offline prediction (+588.4 vs +569.7 =
+18.7) to 0.1 kcal/mol — an independent cross-check of both the review's arithmetic and this
implementation. CH5+ was "not evaluable offline" in A.5 (two corners); it is -0.795 kcal/mol.

### 3.3 The class-C adduct falsifier — this is what the mode is for

model minus r2SCAN-3c reference, kcal/mol, both curves referenced to their own d = 3.00 A point:

| scan | delivered at 1.0/1.2/1.3 A | conserving at 1.0/1.2/1.3 A | dev min del -> con | rms del -> con |
|---|---|---|---|---|
| CH4 + H | -87.0 / -78.0 / -68.6 | **+16.5 / +24.0 / +23.3** | -87.0 -> -1.5 | 41.8 -> **12.1** |
| NH3 + H | -107.0 / -91.9 / -58.2 | **+2.0 / +12.7 / +12.9** | -107.0 -> -1.4 | 47.1 -> **5.1** |
| H2O + H | -54.1 / -53.7 / +5.4 | **+32.4 / +23.9 / +5.4** | -54.6 -> -3.0 | 24.3 -> **13.3** |
| N2H4 + H | -90.4 / -82.3 / -52.5 | **+11.6 / +15.2 / +12.6** | -90.4 -> +0.0 | 41.0 -> **7.0** |

The artificial adduct is gone: the model goes from 54-107 kcal/mol BELOW the reference to
2-32 kcal/mol above it — exactly the "+2..+32" A.5 predicted offline, now measured from a build.

### 3.4 Gradient

Analytic vs central FD (dx 1e-4 A, fresh directory per displacement), worst component:

| geometry | delivered | conserving |
|---|---|---|
| rkt06 point 10 | 1.194e-08 | **1.360e-08** |
| c2h6 runaway-window frame (t = 3207.0 fs) | 1.210e-07 | **1.218e-07** |
| ch4_H frame 10 (the class-C MD start) | 2.689e-07 | **1.421e-08** |
| class-C adduct point CH4 + H, d = 1.20 A | 8.681e-09 | **3.815e-09** |

All far below the 1e-6 Eh/A acceptance. In `conserving` the whole geometry dependence rides the
existing Lambda pass over the term weights (the budget included), so `m_rev_share_dval` is 0 and
there is no `dc/dw` term — a simpler chain rule than the delivered rule's, which the FD confirms.

### 3.5 The 20-cell grid

| arm | rebuilds | step max / kJ | n >= 50 | hard swaps | jump max / kJ | T_max / K |
|---|---:|---:|---:|---:|---:|---:|
| delivered | 1186 | 391.44 | **487** | 0 of 591 | 48.6 | 8 306 |
| conserving | 902 | **216.62** | **72** | 0 of 449 | **21.3** | 11 139 |

**The review's INFERENCE about the two `f10` cells is confirmed by measurement**:

| cell | delivered | conserving |
|---|---|---|
| ch4_H/T1000_f10 | step max **388.4** kJ, **40** events, T_max 6812 | **33.5** kJ, **0** events, T_max 3196 |
| ch4_H/T2000_f10 | 391.4 kJ, 398 events, 84 rebuilds, jump max 48.6 | **216.6** kJ, **6** events, 12 rebuilds, jump max 21.3 |

The carbon-budget snap of the first MD step disappears. Cost, honestly: the three `ch3nh2`
T = 2000 K cells get *worse* on the per-step metric (3 -> 8, 6 -> 35, 7 -> 18 events; step max
53 -> 97 and 67 -> 87 kJ), and `ch4_H/T2000_f10`'s T_max rises 8306 -> 11139 K — the radical H
is no longer held in the artificial adduct well, so it leaves with its kinetic energy instead.
Net over the grid the metric improves by a factor of 6.8.

### 3.6 The systems A.5 names as NOT covered — measured, and the cost is real

Geometries built here and optimised with gfn2 (N-B 1.658 A, N-O 1.345, N-C 1.450, symmetric
O-H-O 1.225/1.225 and N-H-N 1.288/1.288 — all chemically sensible). Energies in Eh, differences
in kcal/mol; "share OFF" is `-gfnff.rev_valence_share false`, which isolates the share itself.

| system | chg | `gfnff` (pinned) | rev share OFF | rev delivered | rev conserving | del-off | **con-off** | con-del |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| H3N-BH3 | 0 | -0.897038980 | -0.765756400 | -0.765730080 | -0.615248850 | 0.017 | **+94.4** | +94.4 |
| H3N-O (amine oxide) | 0 | -0.667945910 | -0.600342710 | -0.600339100 | -0.483440260 | 0.002 | **+73.4** | +73.4 |
| H3N-CH2 (N-ylide) | 0 | -1.008339380 | -0.941139510 | -0.941134510 | -0.766636920 | 0.003 | **+109.5** | +109.5 |
| H5O2+ (Zundel) | 1 | 0.708382040 | 0.528414250 | 0.747824070 | 0.720470080 | +137.7 | +120.5 | **-17.2** |
| N2H7+ | 1 | 0.266793700 | 0.044595910 | 0.321976850 | 0.296027390 | +174.1 | +157.8 | **-16.3** |

- **The dative/ylide neutrals are the price.** Their donor N has four partners at a group charge
  well below 1, so `Val_N ~ 3.2` against `S_N ~ 4` and all four of its wells are scaled by ~0.8.
  The delivered share is inert there (0.002-0.017 kcal/mol), so this is a pure regression of the
  new mode: 73-110 kcal/mol on a 7-8 atom molecule.
- **The proton-shared dimers go the other way**: both modes deviate strongly from the pinned
  `gfnff` value (the bridging H genuinely is a 3c-2e case, which is what the share exists for),
  and conserving is 16-17 kcal/mol LESS repulsive, i.e. closer to `gfnff`.

### 3.7 Tests

`ctest -R "gfnff|sqm_val|react"` 95/98, `ctest -R cli_simplemd_` 19/22 — the same three known
failures, unchanged. The default path is inert, so this is expected rather than reassuring.

### 3.8 Recommendation

**Do not flip the default yet.** The mode does exactly what A.5 said on everything A.5 measured
— the adducts, the f10 cells, the grid metric, rkt06, the hypervalent ions — and the two
implementation refinements above make its hypervalent residuals smaller than predicted. But the
dative/ylide family is a 73-110 kcal/mol regression on molecules that are neither exotic nor
rare (amine boranes, amine oxides, ylides, and by extension sulfoxides/phosphine oxides, which
are NOT covered by the period >= 3 rule either because their donor is the period-2 partner).

What would close it, in order of cost:
1. **Extend the charge rule to a donor rule.** The missing physics is that a dative bond puts a
   full valence into the acceptor's empty orbital; the donor's *formal* charge is +1 even when
   its EEQ charge is +0.2. A per-corner test "this atom has a partner of group 13 or an
   sp3 partner with a formal octet deficit" would grant X = 1 there and cost nothing elsewhere.
   Cheap, but it is a new rule, so it needs its own falsifier set.
2. **Raise X_i by the charge of the WHOLE connected group, not just the H partners.** For
   H3N-BH3 the N-B-H6 group is neutral, so this does not help by itself.
3. Accept the regression and restrict the mode to the reactive path (it is already opt-in).

My recommendation to the operator: keep `delivered` as the default, and treat `conserving` as
the candidate for stage 3b, to be adopted together with (1). The adduct falsifier of 1.4 is now
the acceptance criterion for that work — nothing else in the reference sets measures it.

---

## Package 4 — stage 3a(iii) wells: MG and erf-Morse, both, switchable

`-gfnff.rev_well_form gauss|mg|erfmorse`. **Default `gauss`, bit-identical** (20-cell grid, the
falsifier set to 12 digits and the gfnff yardsticks unchanged on a proven rebuild). Final binary
**`7fef81bd9ef52d337d631a65c446f7e7`**; every number in this section is from that one binary.

### 4.1 The forms and where the parameters come from

Both new forms are `E = -D (2y - y^2)` with `y = exp(-(a x + beta x^2))` (MG) resp.
`y = erfc((x - u)/sigma) / erfc(-u/sigma)` (erf-Morse), `x = r - r0`. Both have `E(r0) = -D` and
`E'(r0) = 0` identically, and both are **curvature-pinned**: `a` (closed form) resp. `u`
(bisection) is set so that `E''(r0) = 2 alpha |k_b|`, the delivered Gaussian's own force
constant. So r_min and the force constant are reproduced by construction and only the depth scale
`s = D/|k_b|` and the tail are fitted — the review's step (1), two parameters per bond type.

`scripts/revgfnff_wellfit.py` (new) does the fit: a kept-topology scan per class-A bond type with
`CURCUMA_SHAREDUMP=1`, `E_rest = E_total - E_pair` frozen, `r0/alpha/k_b` recovered from
`ln(D/w)`, objective = break-side RMS against the r2SCAN-3c curve, **on the charge-frozen rest**
(the Coulomb drift's sign flips between the static and the react protocol, so a depth fitted on
the raw rest absorbs a charge-model error that is not the well's — `FABLE_REVIEW_2` B.2). It
aggregates to element pairs (median) and writes `rev_well_table.h`.

**The fit reproduces the review's offline B.2 table independently**, which is the best available
check on both:

| n = 32, medians | this fit | FABLE_REVIEW_2 B.2 |
|---|---:|---:|
| delivered Gaussian, break rms | **19.24** | 19.24 |
| MG, curvature pinned | **2.07** | 2.10 |
| erf-Morse, curvature pinned | **2.13** | 2.17 |
| n(rms < 3) | **19 / 19** | 19 / 19 |
| depth scale s, median (range) | **1.016 (0.459-2.357)** | 1.02 (0.46-2.36) |
| join, median / max abs(E_pair) at the last grid point | **0.000 / 0.69, 0.57** | 0.000 / 0.57-0.69 |

The well therefore ends by itself, so **the term weight `w` is taken OFF the well** in the two new
forms (`TOPO_REUSE_STATUS` part B's truncation); `w` stays on the angle/torsion terms and on the
share's own sums.

**Unit check, done rather than assumed**: the table is in Angstrom and the workspace in Bohr. The
implemented well reproduces the offline form to 0.36 kcal/mol on `ch3cl_C-Cl`, and the delivered
Gaussian reproduces its own offline form to 0.30 on the same scan — i.e. the residual is the
`ln(D/w)` recovery (0.91 % on D), not the conversion. The wrong-unit control (beta left in
Bohr^-2) is 23.5 kcal/mol off, 66x worse.

### 4.2 MG vs erf-Morse vs gauss, every acceptance row

| row | gauss | MG | erf-Morse |
|---|---:|---:|---:|
| class-A median rms (`--mode kept` + react, 32 bonds) | 24.49 | **19.50** | **19.21** |
| class-A median dev D_e | -25.53 | **-12.74** | **-13.13** |
| class-A median dev r90 | -0.330 | **-0.058** | **-0.068** |
| guard, pooled MAD (167 reactions) | **1.0341** | 1.0439 | 1.0427 |
| guard, worst set movement | — | ICONF +0.10 | ICONF +0.09 |
| max equilibrium bond-length shift (4 molecules, opt) | — | **0.0063 A** | **0.0067 A** |
| class D dE_MAD | **4.980** | 5.152 | 5.184 |
| class D grad_RMS | 16.599 | **16.376** | 16.494 |
| rkt06 rms | 2.7140 | 2.67 | 2.71 |
| adducts, dev min (CH4/NH3/H2O/N2H4 + H) | -87 / -107 / -55 / -90 | -89 / -110 / -89 / -94 | -89 / -110 / -88 / -93 |
| 20-cell grid: rebuilds / step max / n >= 50 / jump max | 1186 / 391.4 / 487 / 48.6 | 1073 / 396.9 / 458 / 49.7 | 1012 / 395.5 / 567 / **58.3** |
| grid wall time (20 cells, cached setup) | 17.39 s | 17.03 s | 17.33 s |
| FD gradient, worst of 4 geometries | 2.69e-07 | 2.77e-07 | 2.76e-07 |

Class-A reference: `BASELINE_HEAD` records 24.68 / -25.27 / -0.318 for the delivered form. Our
gauss arm gives 24.49 / -25.53 / -0.330, i.e. within 1 %. **Not the code**: the same run with
`-gfnff.rev_budget_fix_h false` is identical to the digit, so the package-1 flip is inert here;
the residual is the class-A reference set itself, which gained the rks/uks branches of four bonds
(untracked additions in this working tree) since that baseline.

**The guard does not open** (+0.9 % on a 1.03 kcal/mol MAD = 0.01 kcal/mol), but the review's
inferred r_eq shift is real and **twice its estimate**: measured by optimising four molecules,
c2h6 C-C 1.51301 -> 1.51930 / 1.51967 A, ch4 C-H 1.08854 -> 1.09269 / 1.09314, h2o O-H 0.97274 ->
**0.96690** / 0.96852 (the only one that shortens).

### 4.3 The adducts get worse before the share compensates — measured

As the package brief predicted: a deeper, wider well makes the partial-well sum at a radical
approach WORSE. `h2o_H` goes -54.6 -> -88.5 kcal/mol under the **delivered** share. With the
package-3 **conserving** share the well form is fully compensated and one scan improves further:

| scan, dev min / rms | gauss + delivered | mg + delivered | gauss + conserving | **mg + conserving** |
|---|---|---|---|---|
| CH4 + H | -87.0 / 41.8 | -89.4 / 43.1 | -1.5 / 12.1 | **-1.5 / 11.8** |
| NH3 + H | -107.0 / 47.1 | -110.1 / 48.7 | -1.4 / 5.1 | **-1.4 / 5.1** |
| H2O + H | -54.6 / 24.3 | -88.5 / 39.0 | -3.0 / 13.3 | **-3.0 / 13.3** |
| N2H4 + H | -90.4 / 41.0 | -93.6 / 42.8 | +0.0 / 7.0 | **+0.0 / 2.9** |

rkt06 is 2.7140 / 2.67 / 2.7614 / 2.72 across the same four combinations — unmoved either way.

### 4.4 MG or erf-Morse?

**They are indistinguishable on the data and MG is the cheaper one.** On the 32 class-A curves the
fitted rms differs by at most 0.85 and the medians by 0.06; every acceptance row above agrees to
within its own noise except two, and both favour MG: erf-Morse has the one 20-cell jump above
50 kJ/mol (58.3 vs 49.7) and its class-D dE_MAD is marginally worse.

**Cost**: after caching the per-bond setup, the two forms and the Gaussian are within noise of
each other (17.03 / 17.33 / 17.39 s over the grid). **Before** the cache, erf-Morse cost
**1.42x** the whole react-MD wall time and MG 1.02x — the entire erf-Morse penalty is its
bisection for `u`, which has to run once per bond and not once per energy call. MG needs no cache
at all (its `a` is a closed form). That is the practical content of the review's "MG brings closed
forms where erf-Morse needs a bisection".

**Recommendation**: if a new well form is adopted, take **MG**. erf-Morse buys nothing the data
can see, and it is the form that has to be cached to be affordable.

### 4.5 Method note — one ulp, 77 rebuilds

The first version of this change moved the 20-cell grid from 1186 to **1263** rebuilds with the
default `gauss` form, while every single-point energy stayed bit-identical to 12 digits. Cause:
the restructuring had written `(-2 alpha dr) * energy` as `dwell_dx * w * cshare` and the share's
`K` as `(K w) * (1/w)` — algebraically equal, one ulp apart, and a react MD amplifies one ulp
(Known Issue #33). Both are back to the delivered association, with a comment saying why.
**An energy-level identity check does not catch this; the trajectory fingerprint does.**

### 4.6 Known limitations of this step

- The table is keyed on the element **pair**, so it cannot distinguish C-C from C=C from C#C. The
  per-system fits differ substantially (s = 1.180 / 1.008 / 0.906; beta 0.564 / 0.488 / 0.934), and
  the median is used for all three. That is the stage-3b element factorisation and is why the
  class-A median rms lands at 19.5 rather than at the 2.07 the per-system fit reaches.
- The HB alpha modulation (`egbond_hb`, a softened alpha for a hydrogen-bond donor's X-H bond) is
  **not applied** in the new forms: it would enter through `a` resp. `u` and its chain rule has no
  counterpart there. Using `alpha_orig` keeps the well and its gradient exactly consistent and
  continuous; the cost is that such an X-H bond is not softened. Not measured against a reference —
  the class-A set contains no hydrogen bond.
- The inner side is capped at `y = 2` (C1, exact below y = 1.6, so the fitted region and the
  minimum are untouched), which bounds the well in [-D, 0] exactly as the Gaussian is in [k_b, 0].
  The repulsive wall stays the repulsion term's job, as it is for the Gaussian.
- Element pairs with no class-A data keep the Gaussian and say so at verbosity 2.
- Step (2) of the review's plan — freeing the curvature with an r0 re-solve — is NOT done.

`ctest`: 95/98 and 19/22, the same three known failures.

---

## Package 5 — calibration, tests, documentation

Done on the defaults as they stand: **`rev_budget_fix_h true`, `rev_share_form delivered`,
`rev_well_form gauss`**. If the operator flips the share or the well form, `cli_simplemd_16/18/20`
and `cli_gfnff_03` need another pass — their numbers are calibrations of the *current* default and
each test says so in its own header.

### 5.1 A harness finding that invalidates every "ctest unchanged" line above

`test_cases/cli/test_utils.sh` picks the binary from `release, debug, build, release_rocm,
release_cuda, release_vulkan, build_rocm` — **`build_rev` is not in that list**. So every
`cli_*` ctest run in packages 1-4 measured `release/curcuma` (Sep 14, md5 `e2f4a72a`), not the
package binary, unless `CURCUMA` is exported. The C++ tests (`gfnff_val_*`, `sqm_val_*`,
`react`, ...) are executables built in `build_rev` and were always testing the right thing; only
the bash `cli_*` tests were not.

**Run correctly** (`CURCUMA=build_rev/curcuma ctest -R "gfnff|sqm_val|react|cli_simplemd_"`):
**111 of 113 pass.** The two failures are `cli_simplemd_08/09` (acetic-acid dimer, `-method
gfnff`, "Total energy drift 0.446525 Eh >= 0.10"), which

- `HBUDGET_STATUS` section 5 already recorded on Sep 15 with the **same** 0.4465 Eh, before any
  of this work, and
- **pass** with `release/curcuma`, i.e. they are a property of this branch's build that predates
  these five packages. `-method gfnff` energies are bit-identical throughout (12 digits and the
  `dump_params` md5), so it is not an energy change.

The three tests that the package brief listed as "known failures, recalibrated in package 5"
(`cli_simplemd_16/18/19`) now **pass** against the build_rev binary.

### 5.2 Recalibrations

| test | what was wrong | what it is now |
|---|---|---|
| `cli_simplemd_16` | "T never exceeds 1.6x the target" measured the SAMPLING: at the default print frequency the run prints **17 rows**, and the max over 17 samples of a 6-DOF Maxwell-Boltzmann came out 2.17 / 1.23 / 1.79 / 1.74 / 1.91 over 3600-4400 K — the test passed or failed on which sample was printed | `-md.print_frequency 10` (1601 rows, identical trajectory and event counts), gate the **mean** at 1.5x (measured 1.14-1.17, 28 % margin) and keep a loose max at 8x (measured 5.5-6.0) as an instability guard. dE_jump thresholds untouched: measured median **0.00** and max **0.50-1.50** kJ/mol against 2.0 and 60.0 |
| `cli_simplemd_18` | the operating point moved, not the criterion: at 11500 K the 12-H2 bath now makes **568-592** rebuilds where the Sep-13 calibration measured 16-52, and the NVE slope scales with the event count (real topology changes, each with its own dE_jump) | temperature **11500 -> 8000 K**, where the count is back at 38-52, and the floor re-derived by the same rule — 4x the largest |slope_rev| (6.37e-04) and 7x the largest SE (3.6e-04): **1.6e-3 -> 2.5e-3 Eh/ps**. MAX_T_K keeps its 1.3x meaning, 15000 -> 10400 |
| `cli_simplemd_19` | its reference arm was `-method gfnff` at a 0.01 kcal/mol tolerance, but stage 3a (i)'s r0 fix gives revgfnff a constant equilibrium offset — measured **-0.5555 kcal/mol** on this very dimer | reference arm becomes `-method revgfnff -gfnff.topology_mode static` at the same tolerance. Measured: revgfnff **react == revgfnff static to 0.0000 kcal/mol** (-0.66200626 Eh both), i.e. the join is still exactly free, which is what the test exists to assert. The gfnff arm stays as a loose 2.0 kcal/mol bound on the offset |

Full calibration sweep for test 18 (10 ps, 0.1 ps buckets, seed 42, same bath):

| T / K | dt | slope_rev Eh/ps | SE | rebuilds | slope_gfnff |
|---:|---:|---:|---:|---:|---:|
| 8000 | 0.25 | **6.372e-04** | 3.6e-04 | 38 | 5.5e-05 |
| 8000 | 0.125 | **6.230e-04** | 3.6e-04 | 52 | 4.2e-06 |
| 9500 | 0.25 | 1.766e-03 | 7.6e-04 | 3274 | -1.1e-05 |
| 10500 | 0.25 | 2.521e-03 | 1.2e-03 | 188 | -7.6e-05 |
| 11500 | 0.25 | 4.311e-03 | 1.7e-03 | 568 | -9.5e-06 |
| 12500 | 0.125 | 4.375e-03 | 1.7e-03 | 422 | 8.8e-06 |

### 5.3 New tests

| test | asserts | measured |
|---|---|---|
| `cli_simplemd_20_gfnff_rev_h_budget` | the c2h6 runaway cell, 5 ps: max per-step abs(dEpot) outside rebuild intervals <= **150 kJ/mol** and min r(H-H) >= **1.0 a0**, AND a `-gfnff.rev_budget_fix_h false` control that must violate both | default **59.26** kJ / **1.770** a0; control **2593.60** kJ / **0.559** a0 |
| `cli_gfnff_03_rev_adduct_falsifier` | the class-C CH4 + H approach vs r2SCAN-3c: the delivered profile within 5 kcal/mol of its recorded min dev **-87.0**, and `conserving` above a **-10** kcal/mol floor | delivered -87.0, conserving **-1.5** |
| `cli_gfnff_04_rev_well_form` | gauss == default == the recorded 12-digit caffeine energy, gfnff untouched, and mg/erfmorse each differ from gauss and from each other by > 1e-6 Eh | mg - gauss +0.1266 Eh, erfmorse - gauss +0.1344, erfmorse - mg 7.81e-03 |

The second arm of tests 20 and 03 is the point: without it both would pass on a build where the
share is switched off entirely, and test 04's liveness clause is what stops a switchable well form
from dying silently. Test 04 reads the energy from the `-batch` JSONL because the printed
"Single Point Energy" carries 8 decimals and cannot express a 1e-11 identity statement.

### 5.4 Documentation

- **new** `docs/REV_GFNFF_STAGE3A.md` — the share, the hydrogen budget, the conserving mode and
  the two well forms, with every measured number, an explicit **what was NOT tested** list (no
  metals, no periodic systems, no hydrogen-bonded X-H under the new forms, no MD beyond 5 ps per
  cell, not every flag combination) and a **not implemented** list (free curvature + r0 re-solve,
  a bond-order-resolved well table, a donor rule for the conserving budget).
- `docs/REV_GFNFF_ROADMAP.md` — a stage-3a status table, the restated smoothness falsifier (three
  columns, with the reason the `dE_jump` statistic alone is not one) and the 20-cell correction.
- `CLAUDE.md` — links the new document from the rev-gfnff bullet.
- `AIChangelog.md` — one entry per delivered fact (H budget default, conserving share, well forms,
  the multi-gpu merge).

### 5.5 Open for the operator

1. **`rev_share_form conserving`** — flip or not. It is the only thing measured that removes the
   artificial radical adducts (package 3.3) and it costs 73-110 kcal/mol on dative/ylidic neutrals
   (package 3.6). My recommendation: not yet, and adopt it together with a donor rule.
2. **`rev_well_form mg`** — flip or not. MG and erf-Morse are indistinguishable on the data
   (package 4.4) and MG is the cheaper one; the guard does not open but equilibrium bond lengths
   move by up to 0.0067 A.
3. `cli_simplemd_08/09` — a pre-existing 0.4465 Eh NVE drift of plain `gfnff` on the acetic-acid
   dimer in this branch's build, recorded since Sep 15 and never chased.
4. Adding `build_rev` to `test_cases/cli/test_utils.sh`'s search list would not help (it comes
   after `release`, which exists); the reliable form is `CURCUMA=<path> ctest ...` and it is worth
   putting in the project's test instructions.


---

## Package 6 — the donor rule and the two default flips (2026-09-19)

Operator decision of 2026-09-19: build the donor rule, flip `rev_share_form` to `conserving` if it
closes the dative/ylide regression, and flip `rev_well_form` to `mg` unconditionally. **Both flips
were taken.** Four local commits: `118b289c` (donor rule), `a4712e83` (share default),
`1327aaee` (well default), `33199345` (four ctests re-pointed).

Binaries, all frozen before use: **`e53baabb`** (donor rule, defaults still delivered+gauss) ->
**`890457ed`** (share default flipped) -> **`916847ff`** (both flipped; every "final" number below
is from that copy). `-threads 1` everywhere, `*.topo.json` never reused, fresh directory per
structure.

### 6.1 The donor rule (`-gfnff.rev_share_donor_rule`, default true, `conserving` only)

An atom is granted `X_i >= 1` if, in the corner being evaluated, it has a partner that is either a
**group-13** element (empty p orbital — no bond count can reveal it) or an atom carrying **fewer
partners than its own nominal sigma valence** (a free coordination site). `max()` against the
charge rule, never a replacement, and only in the charge branch. Per-corner constant like `Q_i`,
so no new chain-rule term. Implemented in `FFWorkspace::prepareConservingShare`.

The package-3.6 table, re-measured arm by arm (dev against `-gfnff.rev_valence_share false`,
kcal/mol; binary `e53baabb`, i.e. the gauss well, so it is directly comparable with package 3):

| system | chg | `gfnff` (pinned) | rev share OFF | delivered | conserving, no rule | **conserving + rule** |
|---|---:|---:|---:|---:|---:|---:|
| H3N-BH3 | 0 | -0.897038980 | -0.765756400 | +0.02 | **+94.44** | **+0.00** |
| H3N-O | 0 | -0.667945910 | -0.600342710 | +0.00 | **+73.36** | **+0.00** |
| H3N-CH2 | 0 | -1.008339380 | -0.941139510 | +0.00 | **+109.50** | **+0.00** |
| H5O2+ | 1 | 0.708382040 | 0.528414250 | +137.68 | +120.52 | +120.52 |
| N2H7+ | 1 | 0.266793700 | 0.044595910 | +174.06 | +157.78 | +157.78 |

The "no rule" column reproduces WORK_STATUS 3.6 to the printed digit (+94.4 / +73.4 / +109.5 /
+120.5 / +157.8), which is the independent check on both. All three dative/ylide neutrals become
**bit-identical to the share-off energy**; the two proton-shared dimers are untouched by the rule
and stay 16-17 kcal/mol closer to the pinned `gfnff` value than `delivered` is.

**Nothing else moves.** Every falsifier is bit-identical to conserving without the rule:

| falsifier | conserving, no rule | conserving + rule |
|---|---|---|
| NH4+ / H3O+ / CH5+ / ClO4- / BF4- 1.394 / BF4- 1.143 | 0.827301258361 / 1.207219745650 / 0.781425328593 / 0.039597767376 / -1.470184050200 / 0.147985923798 | identical, 12 digits |
| class-C adducts, dev min / rms | -1.5/12.1, -1.4/5.1, -3.0/13.3, +0.0/7.0 | identical |
| rkt06, 11 points | rms 2.7614 | 2.7614 |
| 20-cell grid | 902 reb / 216.62 kJ / 72 / 0 of 449 / jump 21.3 / T 11139 | identical in every column |
| equilibrium toggle set (20 combinations x 5 molecules) | dE = 0.000000000 | identical |
| `gfnff` caffeine / benzene + dump md5 | -4.672737068614 / -2.362725526194, 4013d6fc / 6c3a87c8 | identical |
| FD gradient (rkt06 pt10 / c2h6 runaway / ch4_H f10 / adduct d=1.20) | 1.360e-08 / 1.218e-07 / 1.421e-08 / 3.815e-09 | identical |
| FD gradient at the rule's own geometry, H3N-BH3 | — | **1.713e-08** Eh/A |

**Scope checked while there**: DMSO and Me3P=O are inert in EVERY arm (+0.00 kcal/mol) — the
period >= 3 octet expansion already caps a sulfoxide's / phosphine oxide's donor, so the worry
recorded in package 3.8 was unfounded. Metals are still untouched by this mode (the d-block keeps
the delivered growth) and no rev-gfnff reference set contains one.

### 6.2 `rev_share_form` -> `conserving` (four sites: PARAM, GFNFF member, `setupRevSettings`
fallback, `RevSettings::share_conserving`)

Verified on `890457ed`: the DEFAULT arm reproduces the explicit `conserving` arm bit for bit on
the whole battery above, and `gfnff` is untouched.

### 6.3 `rev_well_form` -> `mg` (same four sites)

Verified on `916847ff`. **Scoping proven rather than assumed**: `-method gfnff` gives caffeine
-4.672737068614, benzene -2.362725526194 and the two `dump_params` md5s 4013d6fc / 6c3a87c8, i.e.
the whole rev path stays behind `rev_enabled`. The three forms are live and distinct on caffeine:
gauss -4.673521653477, mg (default) -4.546943898047, erfmorse -4.539136588545.

**MG alone** (i.e. against `-gfnff.rev_share_form delivered`, so it is comparable with package 4),
measured on `916847ff` — every row reproduces package 4.2:

| row | gauss+delivered (pkg 4) | mg+delivered (pkg 4) | mg+delivered (final binary) |
|---|---:|---:|---:|
| class-A median rms | 24.49 | 19.50 | **19.50** |
| class-A median dev D_e | -25.53 | -12.74 | **-12.74** |
| class-A median dev r90 | -0.330 | -0.058 | **-0.058** |
| guard, pooled MAD (167 reactions) | 1.0341 | 1.0439 | **1.0439** |
| class D dE_MAD / grad_RMS | 4.980 / 16.599 | 5.152 / 16.376 | **5.152 / 16.376** |
| rkt06 rms | 2.7140 | 2.67 | **2.67** |
| adducts dev min (CH4/NH3/H2O/N2H4 + H) | -87 / -107 / -55 / -90 | -89 / -110 / -89 / -94 | **-89.4 / -110.1 / -88.5 / -93.6** |
| 20-cell grid: reb / step max / n>=50 / jump max | 1186 / 391.4 / 487 / 48.6 | 1073 / 396.9 / 458 / 49.7 | **1073 / 396.88 / 458 / 49.7** |

### 6.4 BOTH flips together — one real interaction, and it is the grid

Packages 3 and 4 only ever measured the two pieces alone. On the 20-cell grid the combination is
**worse than either one**:

| arm | rebuilds | step max / kJ | n >= 50 | hard | jump max / kJ | n >= 50 | T_max / K |
|---|---:|---:|---:|---:|---:|---:|---:|
| gauss + delivered (pre-flip default) | 1186 | 391.44 | 487 | 0/591 | 48.6 | 0 | 8 306 |
| gauss + conserving | 902 | 216.62 | 72 | 0/449 | 21.3 | 0 | 11 139 |
| mg + delivered | 1073 | 396.88 | 458 | 0/535 | 49.7 | 0 | 8 178 |
| **mg + conserving (NEW DEFAULT)** | **1169** | **800.12** | **271** | **0/583** | **190.9** | **5** | **28 392** |

Four cells carry all of it — `c2h6/T2000_f0` 800.1 kJ (jump 190.9, T 28392), `c2h6/T2000_f8`
366.6 (101.8, 11505), `ch4_H/T2000_f0` 341.6 (183.6, 16919), `ch4_H/T1000_f10` 184.5 (132.1) —
while `ch4_H/T2000_f10`, the cell that motivated the share flip, is the BEST of the three arms in
the combined state (75.0 kJ, 5 events, jump 0.4, T 7445). Hard swaps stay 0 of 583, so this is
smooth-window overrun and not a discrete topology event, and both hot cells recover: mean T over
the last 0.5 ps is 2247 and 2347 K at a 2000 K setpoint. Net against the pre-flip default the
per-step event count improves 487 -> 271 and the tail gets worse. **Not root-caused.**
`-gfnff.rev_well_form gauss` recovers 216.62 / 72 / 21.3, `-gfnff.rev_share_form delivered`
recovers 396.88 / 458 / 49.7.

**Everything else is unaffected by the combination** (all on `916847ff`, DEFAULT arm):

| falsifier | value | vs mg alone |
|---|---|---|
| guard, pooled MAD / RMSD | **1.0439 / 2.0588** (ACONF 0.1966, ICONF 3.4077, MCONF 0.5868, PCONF21 1.6060, S66 0.8276) | identical |
| class D, 500 frames | dE_MAD **5.152**, grad_RMS **16.376** | identical |
| rkt06, 11 points | rms **2.72** (mg alone 2.67) | package 4.3's mg+conserving value |
| class-C adducts, dev min / rms | **-1.5/11.8, -1.4/5.1, -3.0/13.3, +0.0/2.9** | mg+delivered is -89.4/43.1, -110.1/48.7, -88.5/39.0, -93.6/42.8 |
| equilibrium toggle set | 20/20 at dE = **0.000000000** | — |
| `gfnff` identity | **-4.672737068614 / -2.362725526194**, md5 4013d6fc / 6c3a87c8 | — |
| FD gradient, 5 geometries | 1.539e-08 / 1.238e-07 / 1.404e-08 / 4.033e-09 / 1.725e-08 Eh/A | — |
| dative/ylide table on the final binary | +96.78 / +68.77 / +106.50 -> **+0.00 / +0.00 / +0.00** | the donor rule still closes it under mg |

**The strongest claim of the prior work is re-confirmed on the current code**: the conserving
share fully compensates the deeper, wider MG well on the adducts, which the delivered share does
not (the two adduct rows above).

**Class A is the one other mover**, and it is the share, not the well and not the donor rule:

| arm | median rms | median dev D_e | median dev r90 | mean rms |
|---|---:|---:|---:|---:|
| gauss + delivered | 24.49 | -25.53 | -0.330 | 28.39 |
| mg + delivered | 19.50 | -12.74 | -0.058 | 22.48 |
| **mg + conserving (DEFAULT)** | **20.52** | **-16.09** | **-0.120** | **22.51** |

**30 of the 32 bond types are bit-identical** between the last two rows. Only `ncl3_N-Cl` (rms
20.10 -> 23.57, dev D_e +17.98 -> -30.21, r90 +0.076 -> -0.319) and `hocl_O-Cl` (35.86 -> 33.30,
an improvement) move, so the median shift is rank re-ordering under one large mover, not a broad
degradation — the mean moves 22.48 -> 22.51. Attributed by ablation: identical with
`-gfnff.rev_share_donor_rule false`, and `delivered` gives 20.10 / +17.98 back. `ncl3_N-Cl` is
left open.

Also recorded, because it will surprise anyone comparing absolute numbers: **an MG well moves
`revgfnff`'s absolute energy away from `gfnff` by construction** — caffeine -4.5469 against
-4.6727, and -133.56 kcal/mol on the acetic-acid dimer of `cli_simplemd_19` where the Gaussian
gave -0.56. Relative energies do not move (the guard is 1.0439). And the compressed-BF4- probe
costs +33.4 kcal/mol under mg against +18.6 under gauss (that geometry is a perception question,
`FABLE_BOND_STATE.md`, not a share gate).

### 6.5 Tests — 113/113

`export CURCUMA=$PWD/build_rev/curcuma; cmake .; ctest -R "gfnff|sqm_val|react|cli_simplemd_|cli_gfnff_"`
gives **113 of 113**. (Package 5's baseline was 111/113; the two failures there,
`cli_simplemd_08/09`, were fixed at HEAD by `63ec8de3`.) Four tests needed re-pointing, none by
weakening a threshold — commit `33199345`:

| test | what changed |
|---|---|
| `cli_gfnff_03` | arms swapped: the default (no flag) must clear the -10 kcal/mol floor (measured -1.5), `-gfnff.rev_share_form delivered` carries the old-behaviour pin, which MOVED -87.0 -> **-89.4** because the delivered share is now evaluated with the MG well. Tolerance unchanged at 5.0 |
| `cli_gfnff_04` | identity is now default == explicit `mg`; the 12-digit PIN stays on `gauss` (-4.673521653477); liveness re-aimed at gauss/erfmorse vs the default |
| `cli_simplemd_20` | the negative control stopped firing. **EITHER flip removes this runaway alone**: under `conserving` `-gfnff.rev_budget_fix_h` is a NO-OP (H's cap is 0 by element — both arms 70.77 kJ / 1.536 a0), and mg+delivered gives the control 53.38 / 1.839. Both arms now pin `gauss` + `delivered` and reproduce 59.26 / 1.770 and 2593.60 / 0.559 exactly; a third arm asserts the shipped default (70.77 / 1.536) is inside the same bounds. Thresholds unchanged (150 kJ, 1.0 a0) |
| `cli_simplemd_19` | its three real assertions passed (react == static to 0.0000 kcal/mol, 0 formations under `order`, the `weight` control fires). Only the loose 2.0 kcal/mol rev-vs-gfnff offset bound broke, at -133.56 — the MG depth. A fourth sub-run with `-gfnff.rev_well_form gauss` now carries that bound at its original threshold |

### 6.6 Open after package 6

1. **The grid interaction of 6.4** — 800 kJ/mol per-step against 216 / 397 for the single flips,
   and 5 rebuild jumps >= 50 kJ where both single flips had 0. Not root-caused. This is the one
   thing a user could notice in hot react MD.
2. **`ncl3_N-Cl`** — the single class-A bond type the conserving share makes worse.
3. **The stage-3b bond-order-resolved well table** — the reason the class-A median is 19.5 and not
   the 2.07 a per-system fit reaches; the remaining half of the MG flip's benefit.
4. **H5O2+ / N2H7+** — the charge is split over two groups, so neither the charge rule nor the
   donor rule grants a full budget; `conserving` is only the smaller of two large errors there.
5. `-gfnff.rev_budget_fix_h` is now dead code in the default configuration. Worth either removing
   or documenting as delivered-only; it is documented here and in the PARAM help, not removed.

### 6.7 Method note

The zsh trap of the `cij` session recurred and was caught by a sanity check rather than by
discipline: a probe loop that bundled `-gfnff.rev_well_form gauss` into ONE shell variable made
all four well forms report the same energy, because the shell does not word-split it and the
parser drops a single argv element it cannot parse. Every harness call in this package passes flag
and value as separate literal tokens (or through a Python list). **If a switch appears to have no
effect, check the argv before the code.**

---

## Package 7 — the grid interaction of 6.4, root-caused (2026-09-19)

Task: reproduce the `mg + conserving` grid tail, localize it, attribute it per term and per bond,
isolate what the COMBINATION does that neither flip does alone, and fix only a minimal genuine
defect. **Result: it is not an interaction of the two flips.** The tail belongs to the
`conserving` share alone; the 20-cell grid was too small to show that. The energy jump itself is a
stage-1b blend-window resolution failure at `dt = 0.25 fs`. One commit, diagnostic only.

Binaries, both frozen before use: **`24a57b1c`** (HEAD `589443c3`, no source change — reproduces
the shipped table exactly, see 7.1) and **`fa09bcac`** (same plus the new `CURCUMA_BLENDDUMP`
log line). `-threads 1`, fresh directory per run, `*.topo.json` never reused, `-md.print_frequency 1`.

### 7.1 Reproduction, and one correction to the shipped table

The 20 grid cells, all four arms, binary `24a57b1c`:

| arm | rebuilds | step max / kJ | n >= 50 | T_max / K | package 6.4 said |
|---|---:|---:|---:|---:|---|
| gauss + delivered | 1186 | **222.53** | 485 | 8 305.9 | 1186 / 391.44 / 487 / 8 306 |
| gauss + conserving | 902 | **216.62** | 72 | 11 139.3 | 902 / 216.62 / 72 / 11 139 |
| mg + delivered | 1073 | **190.23** | 456 | 8 177.8 | 1073 / 396.88 / 458 / 8 178 |
| mg + conserving (DEFAULT) | 1169 | **800.12** | 271 | 28 391.5 | 1169 / 800.12 / 271 / 28 392 |

Rebuild counts and T_max reproduce exactly in all four arms; the two conserving arms reproduce in
every column. The two DELIVERED arms' `step max` does not: 222.53 / 190.23 here against the
recorded 391.44 / 396.88, with the event counts matching to 2 (485/487, 456/458). Checked by
re-running `ch4_H/T2000_f10` (the cell that carries both) with the orchestrator's own
`stepmax.py`: **222.53**, and the max is the same with and without the rebuild-interval exclusion.
The 391/397 figures come from an older binary (package 3.5 records 391.44 for the same cell), so
they are not comparable with the post-donor-rule build. **Every number in this package comes from
one binary.**

### 7.2 The event, localized

`c2h6/T2000_f0`, the worst cell: one event at t = 1.9595-1.9610 ps carries the whole cell
(+393.4, +707.1, +214.5, -800.1 kJ/mol on four consecutive steps, T_max 28 392 K at 1.96125).
Per-term (`terms2.py`): the first step is **bond** +320.7 kJ and angle +65; the +707 and the -800
are **bonded repulsion** (0.0576 -> 0.2657 -> 0.3251 -> 0.0791 Eh) plus the over-coordination term
(0 -> 0.15 Eh), i.e. the aftermath of a collision, not its cause. Every rebuild in the window
reports `dE_jump = 0.000000 Eh`; hard swaps are 0 in the whole run. The geometry (frames extracted
from the MD by truncated re-runs, exact to 1.8e-05 a0 against the status rows): H3 and H4 both sit
on C1 and have formed a transient H2, r(H3-H4) = 2.09 a0; two steps later it is **0.86 a0** —
the same artificial-H2 collapse as `RUNAWAY_STATUS.md`, but reached without any budget jump
(under `conserving` a hydrogen's cap is 0, so `Val_H = 1` always).

### 7.3 Per bond: what the two share rules do to that H2

Same geometry, fresh single point, all four arms (`CURCUMA_SHAREDUMP`), bond H3-H4 at r = 1.98 a0:

| arm | D (well) / Eh | c | E_pair / Eh |
|---|---:|---:|---:|
| gauss + delivered | 0.19724 | **0.0000** | **-0.00000** |
| mg + delivered | 0.20431 | **0.0000** | **-0.00000** |
| gauss + conserving | 0.19724 | 0.2531 | **-0.04993** |
| mg + conserving | 0.20431 | 0.2531 | **-0.05172** |

`f_H3 = f_H4 = Val/S = 1/1.9875 = 0.5031`, `c = f_H3 f_H4 = 0.2531`. The delivered left-over rule
gives the same pair `c = 0.000000` exactly. Both rules give the two C-H bonds the same
`c ~ 0.50`, so **the entire difference between the share rules on this motif is that one pair**.
Consequence, measured as the bond-term difference between the corner WITH the H3-H4 bond and the
corner without it, at the same geometry:

| arm | corner gap (bridged - unbridged) / kJ/mol |
|---|---:|
| gauss + delivered | **+124.4** |
| mg + delivered | **+128.6** |
| gauss + conserving | **-10.6** |
| mg + conserving | **-11.2** |

Under `delivered` the artificial geminal-H2 state costs 124 kJ/mol and the dynamics is pushed out
of it; under `conserving` it is isoenergetic. The well form moves this by 0.6 kJ/mol.

### 7.4 The jump itself: a blend window two steps wide

New `CURCUMA_BLENDDUMP=1` prints the stage-1b transition table per energy call (pair, forming,
tight, window `[w_a, w_b]`, r, coordinate c, corner weight s). Across the +393 kJ step:

| t / ps | pair | w_a -> w_b | r / a0 | c | s |
|---|---|---|---:|---:|---:|
| 1.959500 | 3-4 break | 0.32670 -> 0.02000 | 2.09317 | 0.326696 | **0.000000** |
| 1.959750 | 3-4 break | 0.32670 -> 0.02000 | 2.18296 | 0.170122 | **0.515777** |
| 1.960000 | 3-4 break | 0.32670 -> 0.02000 | 1.44075 | 0.999353 | 0.000000 (revert) |

**Half the break's blend window is traversed in one 0.25 fs step for a distance change of
0.09 a0.** The window is a fixed interval in the bo3 ORDER (0.327 -> 0.02), and the bo3 switch is
steepest for the smallest covalent sum, so for an H-H pair it is only ~0.175 a0 = 0.093 A wide in
DISTANCE - about two steps for a hydrogen at 2000 K. The energy it has to carry over that window
is the corner gap, measured per step over the whole trajectory (`spread.py`, max-min bond term
over the corners evaluated in that step): median **153.2**, p90 209.2, max **389.7** kJ/mol for
this arm. Reconstructing the corner weights from the corner sums and the blended bond term gives
s(C2-H3 formation) 0.0159 -> 0.9682 across the same step, i.e. **0.95 x 353 kJ/mol = 337 kJ of the
measured +321 kJ bond change**. The rest of the event is the collision that energy causes.

**The window, not the potential, is what fails.** Halving the time step on the same cell:

| dt / fs | rebuilds | step max / kJ | n >= 50 | T_max / K |
|---:|---:|---:|---:|---:|
| 0.25 | 76 | **800.1** | 62 | **28 392** |
| 0.125 | 102 | **38.7** | **0** | 5 579 |
| 0.0625 | 44 | **13.5** | **0** | 5 441 |

A discontinuity of the potential would survive a smaller step; this does not. The static potential
along the same path is smooth there: the fresh-perception single point changes by **+5.9 kJ/mol**
between the two frames where the MD jumps by +393.

### 7.5 The 20-cell grid was too small — 130 cells change the conclusion

Same three systems, ALL frames (c2h6 25, ch3nh2 25, ch4_H 15) at both temperatures = **130 cells
per arm**, one binary, same protocol:

| arm | sum reb | H-H formations | median step | p90 | max / kJ | cells > 100 kJ | T_max / K |
|---|---:|---:|---:|---:|---:|---:|---:|
| gauss + delivered | 9473 | 402 | 45.5 | 88.4 | 222.5 | 9 | 8 306 |
| gauss + conserving | 9022 | 138 | 42.9 | 116.2 | **1071.3** | 16 | **46 601** |
| mg + delivered | 8636 | 385 | 46.8 | 80.6 | 362.6 | 8 | 11 035 |
| mg + conserving | 9708 | 204 | 49.4 | 213.5 | 800.1 | **27** | 28 392 |

**`gauss + conserving` reaches 1071.3 kJ/mol and 46 601 K — worse in the maximum than the shipped
default.** Its worst cell is `c2h6/T2000_f7` (reproduced: 1071.3 kJ, 309 events, 172 rebuilds,
T_max 26 265); the 20-cell grid does not contain it. So the heavy tail is a property of the
**conserving share**, not of the combination. What MG adds is FREQUENCY: 27 cells above 100 kJ
against 16, p90 213.5 against 116.2, and on c2h6 at 2000 K (25 cells each) 20/25 cells above
100 kJ against 9/25 with 125 H-H formations against 70. Its mechanism contribution is ~0: the
corner gap moves 0.6 kJ/mol (7.3) and the H-H well itself is within 2.0-3.5 % of the Gaussian over
the whole perceived range (isolated-H2 scan, `mg/gauss` 1.020-1.035 at every r) because the fitted
`beta` for H-H is 1.135 /A^2 - the MG tail is wide for C-H, not for H-H.

The two share rules have **disjoint failure modes**, and each arm's worst cells say which:
every worst cell of the two delivered arms is `ch4_H` (the artificial radical adduct - 282 -> 20
H-H formations and 570 -> 76 rebuilds when the share is flipped, which is exactly what `conserving`
was built for), every worst cell of the two conserving arms is `c2h6` (the intramolecular geminal
H2). At `dt = 0.125 fs` over the same 130 cells the four arms are nearly indistinguishable in the
median (21.0 / 20.7 / 22.3 / 20.9) and the tail shrinks by 3-9x (max 134.9 / 312.5 / 124.0 / 591.6,
cells > 100 kJ 4 / 7 / 4 / 3).

### 7.6 The two structural hypotheses of the task brief — both falsified by measurement

- **(a) "the conserving share's smooth `min` is steeper than the delivered clip"**. True in
  principle (`shareMinOne` joins over `rev_share_min_width` = 0.1 in `Val/S`, `shareClip` over the
  whole unit interval), **but it never fires here**: the bridging hydrogen has `Val/S = 0.5031`,
  deep inside `shareMinOne`'s exactly-linear branch; the cubic join lives only on
  0.9 < Val/S < 1. Measured over the same 130 cells: `-gfnff.rev_share_min_width` 0.02 / 0.10
  (default) / 0.30 gives median step **49.4 / 49.4 / 48.9** and cells > 100 kJ **27 / 27 / 32**.
  No effect.
- **(b) "the new well forms lost the `w` damping, so a transition is smoothed once instead of
  twice"**. `gauss + conserving` keeps `w` on the well and produces the worst event of the whole
  ensemble, so `w` is not what is missing. **Latent, though**: for the new forms `E = well * c`
  and `calcBonds` skips the share when `w <= 1e-12`, so `c` would jump from `f_i f_j` to 1 on a
  pair whose own `w` dies while its atom stays over-claimed by others. **Not reached**: over
  126 290 `shareD` rows of the worst cell's trajectory, 6 956 have `w < 1e-3` and **0 of those have
  `c < 0.999`** - because the over-claim is carried by the pair's own `w`, so `f -> 1` as `w -> 0`.
  One cell, so the scope of that check is one trajectory.

### 7.7 Context: the share is worth a factor 15, either rule

Same 130 cells, `-gfnff.rev_valence_share false` (the only switch that removes the share):
median step **744.4** kJ, p90 2234.3, max **7738.5**, 87/130 cells above 100 kJ, T_max 9.96e+07.
So this is a choice between two second-order failure modes inside a large win, not between a good
and a bad option. Also measured and **negative**: `-gfnff.rev_share_onethree true`, which would
treat the geminal H-H as a 1,3 contact, is catastrophic (median 973.3 / 980.5 kJ, max 22 642 /
6 271, 84/130 cells above 100) - it removes the pair's claim on the hydrogen's valence AND leaves
its well at full depth.

### 7.8 Verdict and the option list (nothing implemented)

**Design tension, no minimal defect.** The share rules and the well forms do exactly what their
specifications say; the conserving rule's `c = f_i f_j` on a pair whose two partners are both
fully committed elsewhere IS the specification, and the delivered rule's exact zero there was an
accident of the left-over rule that happened to suppress the geminal-H2 artefact. The options, each
with what it would cost against the falsifiers that currently pass — **none of these was built**:

1. **Widen the transition window, or define it in distance instead of bond order.** This is what
   actually fails (7.4). The window `[c_seen, tr_begin]` is set at detection precisely because a
   hot X-H moves the coordinate by up to 0.15 per scan, and it is clamped to at least 0.10 in `c` —
   both deliberate. A distance-based window would be ~constant in steps across element pairs
   instead of narrowest for H. Cost: touches every rev-gfnff trajectory, so the whole class-D /
   NVE / react ctest calibration (tests 13-20) has to be re-measured; the equilibrium and
   single-point falsifiers are unaffected by construction.
2. **Require `dt <= 0.125 fs` when `conserving` is active**, or warn. Measured: removes the worst
   cell entirely and cuts the ensemble tail 3-9x (7.4, 7.5). Cost: 2x wall time, and it contradicts
   the documented "reactive MD needs dt <= 0.25 fs" — which is exactly the step at which the tail
   lives. No falsifier moves. Cheapest honest mitigation.
3. **Deny the share's product to a pair whose BOTH partners are over-claimed** (e.g. multiply `c`
   by the pair's settled weight `sig_p`; at the event `sig = 0.0000` while its well is credited
   with 25 % of full depth, so the share's own bookkeeping already disagrees with itself). Cost:
   `sig = 0` is also what a half-formed bond has at an exchange TS by construction (that is what
   the `2b - 1` rescaling was designed for, 3.1), so this is very likely to undo package 3's
   rkt06 path and the class-C adducts. Needs the full class-A/C + rkt06 set re-measured.
4. **Fix the perception rather than the share**: the geminal H-H is admitted as a bond at
   `r/rcov = 1.6` with a tight order of 0.147. Narrowing `rev_bo_center` for H-H would remove the
   artefact at its source. Cost: every bond-order-dependent quantity moves; largest re-validation.
5. **Do nothing.** The tail is bounded (thermostat recovers, hard swaps 0, no crash), it is 15x
   better than share-off, and it needs 2000 K react MD of a saturated hydrocarbon to appear.
   Cost: a user doing hot react MD sees 100-1000 kJ/mol per-step excursions.

### 7.9 Verification of the one commit (diagnostic only)

`CURCUMA_BLENDDUMP=1` in `FFWorkspace::updateTransitions()`, log-only, env-gated, inside a
rev-only function. Binary `fa09bcac` against `24a57b1c`: `gfnff` caffeine **-4.6727370686** /
benzene **-2.3627255262**, `revgfnff` caffeine **-4.5469438980** / benzene **-2.4079121860**,
all identical; the 20-cell grid x 4 arms is identical in **80 of 80** cells (rebuilds, step max,
n >= 50, T_max); `export CURCUMA=$PWD/build_rev/curcuma; ctest -R "gfnff|sqm_val|react|cli_simplemd_|cli_gfnff_"`
gives **113/113**.

### 7.10 Method note

The finding that the two flips "interact" came from a 20-cell grid with one seed per cell. Six and
a half times that many cells of the same three systems put `gauss + conserving` above the shipped
default in the maximum and moved the whole conclusion. Before attributing a tail to a specific
combination of options, check how many independent cells the tail actually rests on — here it was
one cell per arm. Second: a per-step energy jump that disappears when the time step is halved is
not a discontinuity of the potential; that one command separates "the model has a step in it" from
"the integrator cannot resolve a smooth switch" and should be the first thing run on any future
react-MD smoothness complaint.

---

## Package 8 — the time-step recommendation, warned and documented (2026-09-20)

Operator decision on package 7's option 2: **adopt `dt <= 0.125 fs` as an operating recommendation
for react-mode MD with the conserving share — warn and document, do NOT change the default
`-md.time_step`.** The structural fix (redefining the transition window in distance rather than in
bond order, option 1) stays deferred; the window code was not touched. Binary `4c323b80`
(HEAD `c72ee621` + this change).

### 8.1 The effective default time step — it is 0.25 fs, not 1.0

`-md.time_step` has PARAM default **1.0 fs** (`simplemd.h`), but `LoadControlJson` clamps it for
`-method revgfnff|gfnff-rev` to `-md.rev_dt_cap`, PARAM default **0.25 fs** (`simplemd.cpp:309-315`,
stage 1), printing its own warning. So a react run with no time-step flag runs at **0.25 fs** —
exactly the step at which package 7's tail lives, and the warning below therefore fires for the
common no-flag case rather than for an exotic one. It is worded accordingly.

### 8.2 The warning

One-time, at MD initialisation (`SimpleMD::Initialise`, right after the thread verbosity is
re-asserted — earlier is unsafe, Known Issue #31), `CurcumaLogger::warn`, i.e. verbosity >= 1,
plain ASCII. Gate, all four required:

`m_method in {revgfnff, gfnff-rev}` AND `topology_mode == react` AND `rev_share_form == conserving`
AND `m_dT > 0.125` (the measured safe point, a named constant next to the test).

`topology_mode` and `rev_share_form` are read from `ec_config` exactly the way the existing RATTLE
refusal reads `topology_mode` — the `gfnff` sub-scope first, the flat key as fallback — so both
`-gfnff.rev_share_form X` and the auto-routed flat `-rev_share_form X` are seen (verified, 8.3 g/h).
**The well form is deliberately NOT in the gate**: package 7.5 showed `mg` changes how often the
tail is visited, not the mechanism, and `gauss + conserving` reaches the larger maximum.

Text (one line in the log):

> rev-gfnff react MD at time_step 0.250 fs with rev_share_form=conserving (the default): a topology
> transition is blended over a fixed window in BOND ORDER, which for an H-H pair is only ~0.09 A
> wide in distance - about 2 MD steps for a hot hydrogen. A single step can then carry the whole
> corner gap (median 153, max 390 kJ/mol) as a large but SMOOTH energy spike; the thermostat
> recovers and no rebuild is discontinuous. Measured on c2h6 at 2000 K, dt 0.25 -> 0.125 fs takes
> the per-step maximum from 800.1 to 38.7 kJ/mol and T_max from 28392 to 5579 K. For quantitative
> react-mode runs use -md.time_step 0.125 (2x wall time). Detail: docs/REV_GFNFF_STAGE3A.md
> section 2.1

### 8.3 Verification — 10 cases, `c2h6` frame 16, 2000 K, `-maxtime 10`, binary `4c323b80`

| case | invocation (beyond `-md input.xyz`) | warning |
|---|---|---:|
| a | `-method revgfnff -gfnff.topology_mode react -md.time_step 0.25` | **1** |
| a2 | same, NO time-step flag (clamped 1.0 -> 0.25) | **1** |
| b | `... react -md.time_step 0.1` | 0 |
| c1 | `... react -gfnff.rev_share_form delivered -gfnff.rev_well_form mg -md.time_step 0.25` | 0 |
| c2 | `... react -gfnff.rev_share_form delivered -gfnff.rev_well_form gauss -md.time_step 0.25` | 0 |
| d | `-method gfnff -md.time_step 0.25` (no react, not rev) | 0 |
| e | `... react -gfnff.rev_well_form gauss -md.time_step 0.25` (conserving default) | **1** |
| f | `-method revgfnff -md.time_step 0.25`, topology_mode auto | 0 |
| g | `... react -rev_share_form delivered -md.time_step 0.25` (flat flag) | 0 |
| h | `... react -rev_share_form conserving -md.time_step 0.25` (flat flag) | **1** |

Exit code 0 in all ten; the count is the number of matching lines in the whole run log, i.e. the
warning is printed once per run, not per step.

### 8.4 Inertness and regression

Logging and documentation only, no numerical path touched: `gfnff` caffeine **-4.6727370686** Eh
(`-sp -verbosity 3`), unchanged. `export CURCUMA=$PWD/build_rev/curcuma;
ctest -R "gfnff|sqm_val|react|cli_simplemd_|cli_gfnff_"` gives **113/113** — no test greps the MD
startup block strictly enough to notice a new warning line (tests 16/18/20 run exactly the
warning's trigger combination and pass unchanged).

### 8.5 Files

`src/capabilities/simplemd.cpp` (the warning), `docs/REV_GFNFF_STAGE3A.md` (new section 2.2
"Recommended MD settings" + the open-item bullet of section 3 now points at it), `AIChangelog.md`
(one line), this file.

---

## Package 9a — stage 3a(iii) step 2: free curvature + r0 re-solve (`rev_well_form mg2`), opt-in

Operator decision 2026-09-20: build the two pieces the roadmap lists as "not started", each as an
**opt-in** path that does not change what `gauss`/`mg`/`erfmorse` compute and flips no default.

Binary for every number below: **`89a587bb0eb7bb60ba86608f8aed8129`** (`cur_v6`), frozen before
use. The fit that produced the table ran on `04fd1882`/`dc1bbd90` (same gauss/mg behaviour).

**Provenance caveat, measured rather than assumed**: the binary's md5 is NOT reproducible across
commits. `build_rev/curcuma` embeds `git describe` (`curcuma -version`), so re-linking the SAME
source after a commit changes the md5 — here `89a587bb` (`...-117-ge688254b`) became `7f2e2504`
(`...-125-gfa161a06`) with no source change at all. Checked, not inferred: all five well forms
agree to 12 digits on the four reference single points, and 9 of 9 react-MD cells (c2h6 T2000
f0/f10/f23 x {mg, mg2, mg3}) reproduce rebuild count, per-step max, n >= 50 and T_max exactly.
So "same md5" is the wrong identity test for this tree; the trajectory fingerprint is the right
one (memory `revgfnff-run-provenance-fingerprint`).

### 9a.0 The guardrail, verified at the level that can actually see it

`gauss`, `mg` and `erfmorse` are **bit-identical to the pre-session binary `4c323b80`**:

- 12 digits on the four reference single points (`gfnff` caffeine -4.672737068614 / benzene
  -2.362725526194; `revgfnff` caffeine -4.546943898047 / benzene -2.407912185968 under `mg`,
  -4.673521653477 / -2.363224128930 under `gauss`, -4.539136588545 / -2.404586308482 under
  `erfmorse`), and
- **8 react-MD cells** (`c2h6` T2000 frames 0/7/15/23 x {gauss, mg}, 5 ps, dt 0.25, seed 42,
  `-threads 1`) reproduce the OLD binary's rebuild count, per-step maximum, n >= 50 kJ and T_max
  **exactly, all 8 of 8**. This is the check package 4.5 says an energy identity cannot replace:
  a one-ulp change survives 12 printed digits but not a react trajectory. The new `Bond` member
  and the reordered form dispatch pass it.

### 9a.1 What the form is

Same MG well `E = -D (2y - y^2)`, `y = exp(-(a x + beta x^2))`, but `x = r - (r0_model + dr0)` and

| parameter | 'mg' (pinned) | 'mg2' (this step) |
|---|---|---|
| depth `s = D/\|k_b\|` | fitted | fitted |
| tail `beta` | fitted | fitted |
| curvature `ca` (`K_well = ca^2 * 2 alpha \|k_b\|`) | **1 by construction** | **solved from the reference** |
| minimum offset `dr0` | **0 by construction** | **solved from the reference** |

`ca = 1, dr0 = 0` reproduces `mg` bit-for-bit, which is how the code path was verified before any
fit was run (`mg2` == `mg` to 12 digits with a seeded table).

### 9a.2 The curvature is SOLVED, not fitted — and why that had to be measured first

Freeing `ca` as a fit parameter on the break-side objective **does not work**: the well's second
derivative at its own minimum is one point of a 15-point branch, so the simplex uses it as a spare
tail knob and drives it to zero. 8 of 32 bonds came out with a flat, QUARTIC bottom
(`ca -> 0` means `E''(r_min) = 0`) — an excellent break-side rms (median 1.10 vs MG's 2.36) and a
destroyed vibrational force constant. Fitting on the whole grid instead constrains it but drags in
the five COMPRESSED points, where the repulsion term carries the error and the well cannot absorb
it: `f2_F-F` went from a 1.42 break-side rms to 13.10.

The well is the only term whose curvature the form changes and that curvature is exactly
`ca^2 K_gauss` (beta enters only at `x^3`), so

    ca^2 = 1 + (k_ref - k_model) / (2 alpha |k_b|)

with `k_ref`/`k_model` the measured 3-point curvatures of the reference and of the delivered model
at their own sampled minima — the same quantity the class-A harness prints as "k ref/model". Zero
knobs, one closed form, both inputs measured. Fitted range `ca` = 0.30-1.96, one clamp
(`ncl3_N-Cl` at the 0.30 floor).

### 9a.3 The r0 re-solve: the rigid scan cannot do it, the relaxed geometry can

**First attempt, rejected by measurement.** Solving `dr0` so the candidate curve's sub-grid minimum
sits on the reference's (bisection inside the fit, exact to the digit *within the fit*) gives
`dr0(H-O) = -0.052 ... -0.166 A` and moves the optimised water O-H **0.9727 -> 0.9179 A** against a
reference 0.9644 — further from the reference than the pinned form and **5x** the equilibrium
shift of the package-4 precedent. Cause: a class-A scan is RIGID (the rest of the molecule is held
at the r2SCAN-3c geometry) and its minimum responds to `dr0` with `d r*/d dr0 = K_well/k_total`,
which is ~0.1 there, while a relaxed optimisation responds with ~0.33-1.0. The rigid-scan minimum
and the equilibrium bond length are simply not the same quantity.

**Shipped solve.** The class-A r_eq frame IS the r2SCAN-3c minimum, so the distance of the
stretched pair in that frame is the reference equilibrium bond length. Per class-A bond type,

    dr0 = (b_ref - b_model) / response,   response = d b_model / d dr0 measured with a +0.05 probe

(`optbond.py`: optimise the molecule, measure the same pair; `mkdr0.py`: one Newton step, median
per key). Measured response 0.30-1.0 for 31 of 32 bonds (`h2o2_O-O` 3.7), fitted `dr0` range
-0.174 ... +0.235 A, none at the +-0.40 A bound. `(D, beta)` are then re-fitted at that fixed
offset.

**Result, over all 32 class-A systems, |b_model - b_r2SCAN-3c| in Angstrom:**

| arm | median | mean | max |
|---|---:|---:|---:|
| gauss | 0.0293 | 0.0390 | 0.124 |
| mg | 0.0246 | 0.0395 | 0.207 |
| **mg2** | **0.0060** | **0.0170** | 0.179 |
| **mg3** | **0.0044** | **0.0087** | 0.043 |

The equilibrium bond lengths move MORE than `mg`'s and they move **towards** the reference: the
median error against r2SCAN-3c drops by a factor 5-7. That is the honest way to read the guard row
below.

### 9a.4 / 9b.4 Acceptance, every row, one binary

Class-A harness (`--mode kept` + `-gfnff.topology_mode react`, `-method revgfnff`, 32 bond types;
the rms is over the WHOLE common grid, which is why the inner branch matters):

| arm | rms median | rms max | dev D_e median | dev r90 median | dev k median |
|---|---:|---:|---:|---:|---:|
| gauss | 24.68 | 58.80 | -25.79 | -0.341 | -169 |
| mg (shipped default) | 22.15 | 63.56 | -14.32 | -0.066 | -158 |
| erfmorse | 22.16 | 63.47 | -14.86 | -0.079 | -163 |
| **mg2** | **15.83** | 64.62 | -16.17 | -0.101 | **-44** |
| **mg3** | **13.22** | 59.07 | **-6.86** | -0.131 | -82 |

Note the mechanism: the offline BREAK-SIDE fit rms is unchanged by the free curvature (MG 2.158,
MG2 2.157 median over 32), while the harness rms drops 22.15 -> 15.83. The free curvature buys the
INNER branch (dev k -158 -> -44), which the break-side objective never saw and the harness rms
does. The bond-order split (9b) then buys the depth (dev D_e -14.32 -> -6.86).

Everything else, same binary:

| row | gauss | mg | erfmorse | mg2 | mg3 |
|---|---:|---:|---:|---:|---:|
| guard, pooled MAD (167 reactions, ACONF+ICONF+MCONF+PCONF21+S66) | **1.0341** | 1.0439 | 1.0427 | 1.0543 | 1.0547 |
| max equilibrium bond-length shift vs gauss (4 molecules, opt) | - | 0.0063 | 0.0067 | **0.0212** | **0.0212** |
| median \|b - b_r2SCAN-3c\| over 32 class-A systems [A] | 0.0293 | 0.0246 | - | **0.0060** | **0.0044** |
| class D dE_MAD (20 systems, mean) | 4.980 | 5.152 | 5.184 | **4.590** | **4.555** |
| class D grad_RMS (mean) | 14.575 | 14.428 | 14.635 | **11.397** | **11.313** |
| rkt06 path rms, ref = image 0 (conserving) | 2.3840 | 2.3788 | 2.3932 | **2.2665** | **2.2665** |
| rkt06 path rms (delivered) | 2.3505 | 2.3447 | 2.3589 | **2.2599** | **2.2599** |
| FD gradient, worst of the 4 standard points [Eh/A] | - | 1.02e-07 | - | **1.14e-07** | **1.14e-07** |

Class-C adducts, dev min / rms per scan, **under the shipped `conserving` share**:

| scan | mg | mg2 | mg3 |
|---|---|---|---|
| CH4 + H | -1.53 / 11.79 | -1.53 / 11.81 | -1.53 / 11.81 |
| NH3 + H | -1.36 / 6.29 | -1.33 / 6.49 | -1.33 / 6.49 |
| H2O + H | -3.04 / 14.19 | -3.05 / **13.47** | -3.05 / **13.47** |
| N2H4 + H | +0.00 / **2.91** | +0.00 / 3.43 | +0.00 / 3.78 |

Under `delivered` all three arms are the known -90 / -109 / -94 / -91 kcal/mol: the well form does
not fix that, the share does (package 3.3), and mg2/mg3 do not make it worse
(mg -89.4/-110.1/-88.5/-93.6 vs mg2 -90.2/-108.8/-94.2/-91.8).

`ctest -R "gfnff|sqm_val|react|cli_simplemd_|cli_gfnff_"` with `CURCUMA=build_rev/curcuma`:
**113/113**, unchanged.

**The guard opens by 2x the package-4 precedent and this is the one row that is worse.**
gauss -> mg was 1.0341 -> 1.0439 (+0.0098); gauss -> mg2/mg3 is 1.0341 -> 1.0543/1.0547
(+0.0202/+0.0206), i.e. +2.0 % on a 1.03 kcal/mol MAD. No subset collapses; the worst reaction is
unchanged at -20.2 kcal/mol (ICONF) in every arm. The r0 re-solve is what pays for it: the
curvature-only variant (same table with `dr0 = 0`) gives 1.0508/1.0490 but leaves the relaxed bond
lengths 0.0135-0.024 A from gauss and 0.020 A from the reference, so the cheaper guard buys a worse
equilibrium. Both variants are measured; the operator picks.

### 9a.5 A defect in `scripts/revgfnff_wellfit.py`, found and fixed here

The scan that feeds the fit did not pass `-gfnff.rev_well_form gauss`. `prepare()` recovers
`(r0, alpha, k_b)` by fitting a quadratic to `ln(D/w)`, which is exact **only if D is the
Gaussian**. When the script was written `gauss` was the default and the omission was invisible; the
Sep 19 flip to `mg` made it silently wrong. Symptom: the reported "delivered" break rms came out
**9.4** (that is MG's) instead of 19.2, and the first free-curvature table built on it made the
class-A harness WORSE (17.88/14.84 -> the corrected fit gives 15.83/13.22). Fixed; the script now
forces `gauss` in the scan and says why.

**Cross-check of the whole pipeline**: with the fix, regenerating `rev_well_table.h` from scratch
reproduces the COMMITTED `mg` table **exactly** (21 element pairs, max |ds| = 0.0000,
max |dbeta| = 0.0000), and the per-system medians reproduce package 4.6's quoted C-C values
(s 1.180 / 1.008 / 0.906, beta 0.564 / 0.488 / 0.934 for single / double / triple).

### 9a.6 Known limitations of this step

- `dr0` is one Newton step from a linearised response; the residual median |b - b_ref| is 0.0060 A
  (mg2), not 0. A second iteration was not run.
- The `ca` solve uses the class-A curvature of a RIGID scan; for a bond whose molecule relaxes
  strongly that is not exactly the vibrational force constant. No frequency was measured.
- `ncl3_N-Cl` sits on the `ca = 0.30` floor: its reference curvature is far below the model's and
  the closed form asks for a smaller one than the floor allows.
- The pair table takes an INDEPENDENT median per parameter, so a pair with n > 1 can get a
  combination no single system had. Inherited from package 4; 9b removes it for the pairs that
  have more than one order.
- The HB alpha modulation is still not applied (package 4.6), the y-cap at 2 is unchanged, and an
  element pair with no class-A data still keeps the Gaussian.

---

## Package 9b — stage 3b: the bond-order-resolved well table (`rev_well_form mg3`), opt-in

### 9b.1 The key is a CONTINUOUS order, and it is the one the model already computes

The design risk named in the task brief is a discrete single/double/triple classification that can
flip during a trajectory — the failure mode of `RUNAWAY_STATUS.md` and of package 7. It is avoided
by construction:

    order = 1 + pibo * (hyb_i == 1 && hyb_j == 1 ? 2 : 1),   clamped to [1, 3]

`Bond::rev_order`, `GFNFF::continuousBondOrder()`. This is `refreshReactBondOrders()`'s expression
with its ROUNDING and its `pi > 0.5` threshold removed: the second, degenerate pi system of a
linear sp-sp bond enters as a smooth FACTOR on `pibo`, so the order goes to exactly 1 as the pi
order vanishes instead of stepping by 1. The table is then interpolated **linearly between the
fitted orders of that pair and clamped at the outermost**, so a benzene C-C (order **1.666**,
measured) gets a well between the single and the double fit rather than being forced onto one.

`pibo` is a topology quantity (per rebuild, per stage-1b corner), exactly like `fc` and `alpha`, so
the selected parameters are constants inside one energy call: no geometry derivative, no new term
in the gradient (confirmed by the FD row above).

**The table is keyed on the RUNTIME order, not on a chemical label.** Measured at every class-A
r_eq frame: c2h6 1.000, c2h4 2.000, c2h2 3.000, ch2nh 1.993, hcn 2.992, co 2.989, n2 2.998,
h2co 1.994, n2h2 1.997 — and **o2 3.000, not 2**, because GFN-FF gives both oxygens hyb = 1 and the
sp-sp doubling applies. Keying on the label would have put O2's fit at an order the force field
never asks for.

### 9b.2 What the class-A set can and cannot resolve

| pair | orders with data | members |
|---|---|---|
| C-C, C-N, C-O, N-N | 1, 2, 3 | c2h6/c2h4/c2h2, ch3nh2/ch2nh/hcn, ch3oh/h2co/co, n2h4/n2h2/n2 |
| O-O | 1 and 3 | h2o2, o2 (order 2 is interpolated between them) |
| C-H | 1 only, n = 2 | ch4 + hcn (sp3 vs sp C-H; a HYBRIDISATION difference, not an order one) |
| H-O | 1 only, n = 2 | h2o + ch3oh |
| the other 15 pairs | 1 only, n = 1 | exact per-system fit |

So 30 of the 32 class-A bonds now have their OWN key (n = 1) and only 4 bonds still share two
keys. The two that remain (C-H, H-O) are averaged over a hybridisation difference, which a bond
ORDER dimension cannot separate — that is why the class-A median lands at 13.22 and not at the
per-system floor.

### 9b.3 The smoothness re-verification (130 cells, not 20)

Package 7's method: the same 3 systems x 2 temperatures x EVERY frame = 130 cells per arm, 5 ps,
dt 0.25 fs, CSVR 10, seed 42, `-threads 1`, `-md.print_frequency 1`.

| arm | cells | sum rebuilds | median step | p90 | max / kJ | cells > 100 kJ | T_max / K | n >= 50 | max dE_jump | n(jump) >= 50 |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| gauss | 130 | 9713 | 44.1 | 154.3 | 417.2 | 19 | 11 964 | 2355 | 259.4 | 35 |
| mg (default) | 130 | 9881 | 50.3 | 215.7 | 446.8 | 27 | 16 919 | 2567 | 272.3 | 14 |
| mg2 | 130 | 8684 | 51.2 | 214.9 | 432.6 | **21** | 22 002 | 2412 | 359.8 | 26 |
| mg3 | 130 | 8656 | 47.4 | 254.2 | 441.3 | 24 | 22 002 | 2718 | **456.3** | **89** |

Hard swaps (`CURCUMA_VERB=3`, `begin_form`/`begin_break` with `s >= 0.99`), 20-cell grid:
**0 of 623 (mg), 0 of 497 (mg2), 0 of 469 (mg3)**.

Reading: the per-step statistic is comparable to `mg` (median and p90 within noise, MAX slightly
lower, 3-6 fewer cells above 100 kJ than `mg`, 12 % fewer rebuilds), T_max is ~30 % higher, and the
**rebuild `dE_jump` tail is clearly worse for mg3** — 89 events >= 50 kJ/mol against 14 for `mg`.

### 9b.4 Is that tail the order variable? Measured: no.

Two independent checks.

**(a) The response to the order is exactly linear.** New DIAGNOSTIC PARAM
`-gfnff.rev_well_order_override` (Double, default -1 = off, read only by `mg3`; verified inert -
all five forms bit-identical with and without the PARAM in the build) forces every bond's order, so
`dE/d(order)` is measurable. Scanning c2h6 frame 0 from order 1.00 to 3.00 in steps of 0.05:
E rises monotonically by 26.33 kcal/mol, `|dE/d(order)|` is 10.5-15.9 kcal/mol per unit order, and
the **largest single step is 0.7930 kcal/mol against a linear prediction of 0.7929** — i.e. the
well's response to a drifting pi order is as smooth as the drift itself, with no step anywhere in
the interval.

**(b) The large rebuild jumps carry NO order change.** 5 independent cells (c2h6 f4/f10,
ch3nh2 f7/f19, ch4_H f10 at 2000 K, 1.5 ps each) re-run with `CURCUMA_BONDDUMP=1`, correlating each
rebuild's `dE_jump` with the maximum |change of `rev_order`| over the bonds that survive that
rebuild: **259 rebuilds, 10 with `dE_jump` >= 50 kJ/mol, 0 of them with any order change at all**
(max |d order| = 0.0000). Conversely the 4 rebuilds that DO move an order carry a maximum
`dE_jump` of **0.1 kJ/mol**.

So the mg3 tail is not the new dimension; it is the ordinary "a deeper, wider well makes a topology
change cost more" mechanism, with mg2 (same form, no order dimension) already at 359.8/26 between
`mg`'s 272.3/14 and mg3's 456.3/89.

### 9b.5 Known limitations of this step

- C-H and H-O still average two systems each; that difference is hybridisation, not bond order.
- O-O has no order-2 datum; order 2 is a straight interpolation between h2o2 and o2.
- Every pair with a single class-A member has an n = 1 fit and no cross-validation whatsoever.
- The interpolation is C0 in the order (piecewise-linear, knots at the fitted orders). The order is
  constant within an energy call, so this never enters a force; it would matter only if a future
  design made the order a continuous function of the geometry.
- No metal, no charged system, no hydrogen bond, no periodic system was tested; the class-A set is
  H/C/N/O/F/Cl only.
- `rkt06`'s absolute number here (2.27-2.38) does not reproduce package 4's 2.71 — that table used
  `-batch_reuse_topology true`, which keeps the reactant image's topology for the whole path and
  gives 8.06 on this binary. All arms above use the same (fresh-perception) protocol, so the
  comparison between them is valid; the absolute value is not comparable with package 4's.

---

# Package 10 — the MD time-step unit fix, and everything rev-gfnff measured re-derived in TRUE fs

AI-generated, machine-tested. Branch `reactff2-llm`, start HEAD `f1ea50d5` (clean). Binaries frozen
per phase: `curcuma_pre` md5 `7f2e2504` (HEAD, pre-fix), `curcuma_post` md5 `46dbfbef` (core fix),
`curcuma_final` (fix + warning re-threshold + PARAM help). Harness in
`.../scratchpad/p10/` (`rundt.sh`, `sweep.sh`, `agg.py`, `jump.py`, `probe_clock.py`, `nve.py`,
`t20fix.sh`). Every number below was taken with a frozen copy; the md5 stands in each run
directory's `wall.txt`.

## 10.1 The bug and the constant

`SimpleMD` integrates in Angstrom / amu / Hartree. Velocities are `sqrt(kb_Eh*T/m)`, i.e.
`sqrt(Eh/amu)`, so the step multiplying a velocity to give an Angstrom carries the unit
`A*sqrt(amu/Eh) = sqrt(amu*A^2/Eh)`. Derived three independent ways from curcuma's own CODATA-2018
constants, all agreeing to 10 digits:

| route | value |
|---|---|
| `sqrt(ATOMIC_MASS_UNIT * 1e-20 / (HARTREE_TO_KJMOL*1000/AVOGADRO))` | 1.9516144204 fs |
| `sqrt(1.66053906660e-27 * 1e-20 / 4.3597447222071e-18)` (CODATA Eh) | 1.9516144204 fs |
| `sqrt(AMU_TO_AU) * ANGSTROM_TO_BOHR` = 80.68242 atomic time units | 1.95161442 fs |

The user's `-md.time_step` was handed to the integrator unconverted, so **every curcuma MD ran at
1.9516144x the requested step and `-MaxTime` was stretched by the same factor**, for every method
and every caller of `SimpleMD`. Found on `origin/feature/multi-gpu` (`ef462fcf`); re-implemented
here directly rather than cherry-picked, because that commit's diff also touches an adaptive
step-rejecting integrator this branch does not have.

(Noticed while there, **not** changed: `CurcumaUnit::Constants::ATOMIC_TIME_TO_FS` and
`CurcumaUnit::Time::ATOMIC_TIME_TO_FS` are both `24.188843265857`, which is aut -> ATTOseconds, not
femtoseconds — off by 1000. Neither has a single use site anywhere in `src/`, so nothing is wrong
today; it is a loaded gun for the next person who reaches for it.)

## 10.2 The evidence (A.2 / A.3 / A.4)

**Frequency cross-check** — water relaxed with the method under test (`-opt.optimizer lbfgs`),
Hessian frequencies (a path validated against xtb 6.7.1 to 0.13 %, Known Issue #28), one O-H
stretched by 0.006 A, the local-mode period read off an NVE trajectory at ~0 K:

| method | local mode / cm-1 | period expected / fs | MD reports, PRE | ratio | MD reports, POST | dev |
|---|---:|---:|---:|---:|---:|---:|
| gfnff | 3635.9 | 9.1742 | 4.6961 | 1.9536 | 9.1663 | -0.09 % |
| gfn2 | 3646.6 | 9.1474 | 4.6861 | 1.9520 | 9.1466 | -0.01 % |
| gfn1 | 3730.8 | 8.9409 | 4.5825 | 1.9511 | 8.9428 | +0.02 % |

Method-independent, as it must be: the integrator only ever sees a gradient in Eh/Angstrom.

**The scheme is untouched, proven not argued**: the pre-fix binary at `-dt 0.25` and the fixed one
at `-dt 0.4879036051` (= 0.25 x 1.9516144204) give **bit-identical** trajectories — CH4 / gfnff /
NVE / 202 frames, max |dx| = **0.000e+00 Angstrom**.

**NVE conservation still scales as dt^2** (rms fluctuation of Etot over 1 ps, `-thermostat none`):
CH4 gfnff 4.14e-06 / 1.03e-06 at dt 0.5 / 0.25 (ratio 4.01), CH4 gfn2 7.01e-06 / 1.68e-06 (4.17),
H2O gfn2 1.16e-05 / 4.06e-06 (2.87). Below dt 0.125 the series flattens against the 1e-6 Eh print
precision. H2O + gfnff shows a dt-independent ~1.3e-5 Eh floor — **pre-existing**, reproduced on
the pre-fix binary at matched native steps (1.65e-05 / 1.40e-05 / 9.3e-06 / 1.33e-05), not related
to this fix.

**Every static path is byte-identical**: `-dump_gradient` files (12-digit energy, 15-digit
gradients) for gfnff / revgfnff / gfn2 on caffeine, C6H6 and H2O — 9 of 9, `diff -r` clean.

**New ctest `md_time_axis`** (`test_cases/check_md_time_axis.py`) locks the clock to the Hessian
frequency for gfnff and gfn2; passes on the fixed binary (0.23 s), fails on the pre-fix one with
the factors above.

## 10.3 ctest blast radius

Full suite, same machine, same `build_rev`, `CURCUMA` pointed at the frozen binary:

| | failures | which |
|---|---:|---|
| PRE (baseline) | 12 / 304 | 11 known pre-existing + `md_time_axis` (fails BY DESIGN — that is the test working) |
| POST | 12 / 304 | the same 11 pre-existing + **`cli_simplemd_20_gfnff_rev_h_budget`** |

The 11 pre-existing: `confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`,
`cli_curcumaopt_07_opt_multixyz`, `cli_confscan_01..07`. `cli_simplemd_16/18/19` — the three the
brief flagged as candidates — **pass unchanged**; their criteria (slope, mean, event counts) turned
out to be time-scale robust.

**`cli_simplemd_20` is FLAGGED, not recalibrated** (see the new header block in its `run_test.sh`).
Its negative control (`-gfnff.rev_budget_fix_h false`) must violate both bounds and no longer does:

| arm | PRE (real 0.4879 fs) | POST (real 0.25 fs) |
|---|---|---|
| gauss+delivered, fix_h on | 59.26 kJ / 1.770 a0 | 25.40 / 2.080 |
| gauss+delivered, fix_h false (control) | **2593.60 kJ / 0.559 a0** | 32.68 / 2.120 |
| shipped default | 70.77 / 1.536 | 31.48 / 1.994 |

The time-step fix is not implicated: restoring the OLD physical regime on the FIXED binary
(`-md.rev_dt_cap 0 -md.time_step 0.4879036051 -maxtime 9758.072102 -md.coupling 19.516144204
-md.remove_com_motion 195.16144204`, i.e. every dt-derived setting x 1.9516144204) reproduces all
three arms **exactly**: 59.26 / 2593.60 / 70.77 kJ and 1.770 / 0.559 / 1.536 a0.

It is **not** a mechanical recalibration, because there is no dt to recalibrate TO. Holding the
discrete dynamics fixed (CSVR per-step ratio 0.025, COM removal every 400 steps, 20000 steps) and
varying only the true step length, the control arm gives max per-step dEpot / min r(H-H):

| 0.125 | 0.20 | 0.25 | 0.30 | 0.35 | 0.40 | 0.45 | 0.4879 | 0.50 | 0.55 | 0.60 |
|---|---|---|---|---|---|---|---|---|---|---|
| 17.19/1.757 | 19.81/2.082 | 32.68/2.120 | 34.36/2.144 | 45.21/1.986 | 52.39/2.098 | 49.31/1.969 | **2593.60/0.559** | 65.99/1.930 | 91.20/1.941 | **641.62/0.458** |

0.45 and 0.50 bracket 0.4879 and both stay inside, so the violation is a rare event of ONE chaotic
trajectory that the step length reshuffles, not a threshold in dt. (`-md.seed` does not help: it
does not change the initial velocities on this path — 8 seeds, all bit-identical.) **Operator
question**: the falsifier for `-gfnff.rev_budget_fix_h` now exists only ABOVE the shipped 0.25 fs
cap, so the default needs either a new cell/temperature that exposes the budget at 0.25 fs, or a
re-scoping of what test 20 asserts.

## 10.4 The react-MD time step in TRUE femtoseconds (C.1 / C.2 / C.3)

Protocol anchored first: the pre-fix binary at nominal dt 0.25 over package 7's 130 cells
reproduces package 7's shipped-default row **exactly** — sum rebuilds 9708, median 49.4, p90 214.2,
max 800.1 kJ, 27 cells above 100 kJ, T_max 28 392 K, worst cell `c2h6_T2000_f0`.

Re-measured on the fixed binary, 130 cells (c2h6 / ch3nh2 / ch4_H, every frame, 1000 and 2000 K),
**9.758 ps each = the same physical exposure as package 7**, shipped defaults, `-threads 1`,
`-md.print_frequency 1`, fresh directory per cell:

| true dt / fs | sum reb | median / kJ | p90 | max | cells>100 | cells>200 | steps>=50 per 1e4 | T_max / K |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| (old nominal 0.25 = real 0.488) | 9708 | 49.4 | 214.2 | 800.1 | 27 | 15 | 10.2 | 28 392 |
| 1.0 | 9209 | 97.8 | 347.1 | 1133.1 | 64 | 31 | 159.6 | 5.5e7 |
| 0.5 | 10475 | 51.8 | 218.4 | 453.1 | 25 | 16 | 13.1 | 18 928 |
| 0.25 | 9677 | 27.5 | 64.8 | 483.0 | 11 | 9 | 1.2 | 11 948 |
| 0.125 | 13597 | 13.1 | 61.5 | 204.2 | 10 | 1 | 1.9 | 10 405 |
| **0.0625** | 20922 | **6.3** | **16.1** | **78.5** | **0** | **0** | **0.2** | 9 321 |

Median, p90, cells>100, cells>200 and T_max are all monotone in dt; the single-cell maximum and the
raw event count are not (chaotic scatter of one cell). **0.0625 fs is the first true step at which
no cell of the 130 exceeds 100 kJ/mol in a single step**, which is the direct analogue of what
package 8's 0.125 claimed on its ONE cell. It also equals the old threshold rescaled
(0.125 / 1.9516 = 0.0640) — the arithmetic prior and the measurement agree, but the recommendation
rests on the measurement.

Worst cell of package 7, same discrete dynamics, only the step length varying:
0.4879 fs -> 76 reb / **800.1** kJ / 62 events / **28 392** K; 0.25 -> 102 / 47.8 / 0 / 5 206;
0.125 -> 170 / 20.1 / 0 / 5 520; 0.0625 -> 82 / 13.2 / 0 / 5 602.

**C.3 — what `rev_dt_cap` = 0.25 means now.** It clamps to a genuine 0.25 fs where it used to
integrate 0.488 fs, so the shipped default improved by roughly a factor two in every robust
statistic for free (median 51.8 -> 27.5, cells>100 25 -> 11, T_max 18 928 -> 11 948 K). It is still
**outside** the band where the whole sample is bounded (11 of 130 cells above 100 kJ, max 483), so
the warning still fires for a plain react run with no flags. **The PARAM default was NOT changed**;
only its help text now says what the number means. Operator decision.

**Shipped change**: `rev_react_dt_advice` 0.125 -> **0.0625** in `SimpleMD::Initialise`, and the
warning text now carries the 130-cell table plus an explicit note that every time in it is a real
femtosecond. Gate re-verified: fires at 0.25 and 0.125, silent at 0.0625 and 0.05, silent for
`-gfnff.rev_share_form delivered` and for plain `gfnff`.

## 10.5 mg / mg2 / mg3 at the corrected clock (C.4) — package 9's ordering does NOT survive

Same 130 cells, same 9.758 ps, one frozen binary per phase. `n>=50` is the count of rebuild
`dE_jump` events above 50 kJ/mol; the rate normalises it by the rebuild count, which differs
between arms.

| clock / arm | rebuilds | max dE_jump | n>=50 | rate / 1000 reb | step median | step p90 | step max | cells>100 | T_max |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| **old (real 0.488 fs)** mg | 14551 | 304.4 | 24 | 1.65 | 49.4 | 214.2 | 800.1 | 27 | 28 392 |
| old mg2 | 12882 | 352.8 | 37 | 2.87 | 48.3 | 183.3 | 460.9 | 27 | 23 040 |
| old mg3 | 12459 | 440.5 | 72 | **5.78** | 47.3 | 182.7 | 700.1 | 22 | 56 270 |
| **true 0.25 fs** mg | 14506 | 260.5 | 9 | **0.62** | 27.5 | 64.8 | 483.0 | 11 | 11 948 |
| true 0.25 mg2 | 11801 | 300.6 | 21 | **1.78** | 25.2 | 120.8 | 447.8 | 14 | 15 480 |
| true 0.25 mg3 | 11780 | 283.5 | 16 | 1.36 | 25.7 | **60.3** | **316.5** | **9** | 15 738 |
| **true 0.0625 fs** mg | 31371 | 317.2 | 11 | 0.35 | 6.3 | 16.1 | 78.5 | 0 | 9 321 |
| true 0.0625 mg2 | 23695 | 254.4 | 2 | **0.08** | 6.7 | 17.8 | 81.6 | 0 | 9 380 |
| true 0.0625 mg3 | 26712 | 354.0 | 6 | 0.22 | 6.9 | 18.6 | 79.5 | 0 | 9 380 |

The old-clock rows reproduce package 9's **ordering** (mg < mg2 < mg3, mg3 ~3.5x mg's rate) on this
binary; the absolute values differ from package 9's table (272.3/14, 359.8/26, 456.3/89) because
that table was taken on a different binary and does not reproduce package 7's own shipped-default
row either, while the measurement here does, exactly.

**At the corrected clock the ordering breaks.** At true 0.25 fs: mg (0.62) < mg3 (1.36) < mg2
(1.78) by jump rate, and on the per-step statistics **mg3 is the best of the three** (max 316.5
against 483.0 / 447.8, 9 cells above 100 kJ against 11 / 14). At true 0.0625 fs: mg2 (0.08) < mg3
(0.22) < mg (0.35), and **all three have zero cells above 100 kJ/mol per step**. mg3's penalty
relative to mg falls from 3.5x to 2.2x and it is never the worst arm again.

**Consequence for the pending adoption decision**: the "mg3 has a clearly worse smoothness tail"
argument was an artefact of the too-coarse clock and should not weigh against `mg3`. The other two
costs recorded in package 9 — the conformer/S66 guard opening 1.0341 -> 1.0543/1.0547 and the
0.0212 A equilibrium bond-length shift — are time-step independent and stand unchanged.

## 10.6 Files touched

- `src/core/units.h` — `MD_TIME_UNIT_FS` / `FS_TO_MD_TIME` with the derivation in the comment.
- `src/capabilities/simplemd.cpp` — conversion in `Verlet()`, `Rattle()` (incl. the `1/dt` of the
  constraint velocity correction) and `NoseHover()`; `m_dt2` removed; `rev_react_dt_advice`
  0.125 -> 0.0625 and the warning text rewritten.
- `src/capabilities/simplemd.h` — `m_dt2` member removed; `rev_dt_cap` help text.
- `test_cases/check_md_time_axis.py` + `test_cases/CMakeLists.txt` — new ctest `md_time_axis`.
- `test_cases/cli/simplemd/20_gfnff_rev_h_budget/run_test.sh` — header note only, **no threshold or
  flag changed**, left failing on purpose.
- `docs/REV_GFNFF_STAGE3A.md` 2.1 / 2.2 / 2.3, `docs/REV_GFNFF_STAGE1.md`,
  `docs/REV_GFNFF_ROADMAP.md`, `AIChangelog.md`.

## 10.7 What was NOT done

- No default value changed (`rev_dt_cap` 0.25, `rev_well_form mg`, `rev_share_form conserving`,
  `time_step` 1.0 all untouched).
- Nothing imported from `origin/feature/multi-gpu` beyond the idea: no adaptive integrator, no GPU
  eigensolver, no `docs/MD_LARGE_SYSTEMS.md`, no electrostatics-cutoff or optimizer change.
- `cli_simplemd_20` not recalibrated (10.3).
- The mg/mg2/mg3 comparison re-ran only the smoothness axis; the class-A / guard / class-D
  falsifiers are time-step independent by construction and were not re-run.
- Everything packages 1-9 wrote in "fs" outside the files listed in 10.6 still carries the old
  scale. The two `docs/` pages most read (STAGE1 2.x, STAGE3A 2.1/2.3) now carry an explicit
  conversion box; `test_cases/revgfnff/_log/*.md` do not.

---

# Package 11 — the mg / mg2 / mg3 react-MD tail, re-measured with paired replicates (2026-09-21)

Measurement only. No source change, no default change, no adoption recommendation.

## 11.0 PRE-REGISTRATION (written and committed BEFORE the sweep was run)

**Why a pre-registration.** The react-MD tail is a rare-event statistic of a CHAOTIC trajectory
(a 1-ulp change is amplified — Known Issue #33 — and `-md.seed` does not vary the initial
velocities on this path), so ONE run per (cell, arm) is ONE sample. Package 9 and package 10
each drew an ordering of `mg`/`mg2`/`mg3` from a single 130-cell aggregate, and package 10's
summary was then contradicted by the orchestrator's direct spot-check on the single
most-scrutinised cell. With several plausible summary statistics available after the fact, the
statistic that decides is fixed here, in advance, in writing.

**Anchor (the identity check for the harness, reproduced before anything else was run).**
`c2h6/T2000_f0`, true `-md.time_step 0.25`, unperturbed start, `-maxtime 5000` (5 ps true),
`-threads 1`, shipped defaults otherwise, binary `build_rev` at HEAD `d60fae2b`
(md5 `cc6aace6dcbc4800aaf4b6c5cbcc7a13`): **mg 28 rebuilds / 28.3 kJ / 0 steps >= 50 / 5206 K;
mg2 26 / 32.2 / 0 / 5069; mg3 42 / 256.8 / 12 / 7030** — all three reproduce exactly through the
harness `scripts/revgfnff_tail_sweep.py` (md5 `d4b3b0fe96d6f37e81d5f61b88be558f` at the time of
this pre-registration; `revgfnff_tail_sweep.py anchor` re-checks it).

**Design.**
- Cells: the 130-cell set of packages 7/10 — `{c2h6, ch3nh2, ch4_H}` x `{1000, 2000 K}` x every
  frame (25 / 25 / 15).
- Replicates: **6 per cell** — replicate 0 is the UNPERTURBED start (so the historical numbers
  are rows of the CSV), replicates 1-5 displace every atom by a vector of exactly **1e-5 A** in a
  uniformly random direction. The seed depends only on (cell, replicate), so the SAME six start
  geometries are used for every arm and every dt: every comparison is **paired** on
  (cell, replicate). Verified in advance that the perturbation does change the trajectory
  (`c2h6/T2000_f0`, mg3, dt 0.25: rebuilds 42 / 98 / 34 / 64 / 52 / 80 over the six replicates).
- Arms: `mg` (reference), `mg2`, `mg3`, `gauss`, and `gauss + -gfnff.rev_share_form delivered`
  (the pre-2026-09-19 default) as historical context. Everything else stays at the shipped
  default (`rev_share_form conserving`, `rev_share_donor_rule true`, `rev_budget_fix_h true`),
  confirmed via `CURCUMA_REVDUMP=1`.
- dt (true fs): 0.25 (the shipped `rev_dt_cap`), 0.125, 0.0625, each at a FIXED true duration of
  5 ps, so the step count scales with 1/dt and exposure is not confounded with dt.
- 130 x 6 x 5 x 3 = **11 700 trajectories**, `-threads 1`, 8 concurrent processes.

**Per-trajectory quantity.** `max_step_kj` = the maximum over steps of |Epot(i+1) - Epot(i)| in
kJ/mol, with every interval whose end row is preceded by a `REACT rebuild` line excluded —
byte-for-byte the quantity packages 7/9/10 reported as "per-step max |dEpot|".

**PRIMARY decision statistic.** A trajectory is **bad** iff `max_step_kj >= 100 kJ/mol` (the same
threshold package 8/10 used for "cells above 100 kJ/mol", chosen because it is the criterion the
existing dt recommendation rests on). Per (arm, dt): `p = #bad / #trajectories` over all
130 x 6 = 780 trajectories.

**PRIMARY decision rule.** Arm X is declared **worse than `mg`** at a given dt iff the paired
difference `dp = p_X - p_mg` is positive AND its 95 % **cluster-bootstrap** confidence interval
excludes 0. The bootstrap resamples the **130 cells** with replacement (10 000 draws, fixed
RNG seed 20260921), keeping all replicates of a drawn cell and keeping the pairing — cells are
the independent unit, replicates inside one cell are not. Anything else is reported as **"not
distinguishable at 6 replicates"**. The exact McNemar binomial p-value on the discordant
(cell, replicate) pairs is reported alongside as descriptive only (it ignores the clustering).

**SECONDARY, also fixed here**: the same bad-fraction at a 50 kJ/mol threshold; the distribution
of `max_step_kj` per arm (median / p90 / p99 / max); steps >= 50 kJ/mol per ps of trajectory and
the fraction of trajectories with at least one; the rebuild `dE_jump` tail (max, n >= 50);
`T_max`; and the **overlap of the failing-cell sets** between arms (a cell fails for an arm if at
least one of its six replicates is bad) with its Jaccard index.

**Raw data.** One CSV row per (cell, replicate, arm, dt) in
`test_cases/revgfnff/_log/tail_remeasure.csv`. No raw MD log is stored (a verbosity-1 log of a
5 ps / 0.0625 fs run is ~19 MB): stdout is parsed while it streams and the run directory is
deleted immediately. Logs are kept only for the handful of runs used in the 11.C attribution.

**What this measurement cannot answer.** It is the smoothness axis only. The two other costs
recorded in package 9 — the conformer/S66 guard opening 1.0439 -> 1.0543/1.0547 and the 0.0212 A
equilibrium bond-length shift — are static, time-step-independent quantities and are not re-run
here. The adoption decision is the operator's.

## 11.1 What was run

11 700 trajectories, all `rc = 0`, 19 750 s of CPU in 2 471 s wall at 8 concurrent single-thread
processes. Binary `cur11` = a frozen copy of `build_rev/curcuma` rebuilt at HEAD `d60fae2b`
(`MAKE_EXIT=0`, md5 `cc6aace6dcbc4800aaf4b6c5cbcc7a13`); harness
`scripts/revgfnff_tail_sweep.py`; raw data `tail_remeasure.csv` (11 700 rows, 1.2 MB). **Nothing
was cut** — the full design of 11.0 was executed (130 cells x 6 replicates x 5 arms x 3 dt, 5 ps
true each). No raw MD log was written: stdout was parsed as it streamed and each run directory
deleted. The `rev dump (flags)` line confirms the shipped defaults were active
(`share_form conserving`, `share_donor_rule true`, `budget_fix_h true`).

**Anchor: reproduced exactly**, through the harness itself (`revgfnff_tail_sweep.py anchor`):
mg 28 rebuilds / 28.3 kJ / 0 steps >= 50 / 5206 K, mg2 26 / 32.2 / 0 / 5069, mg3 42 / 256.8 / 12 /
7030 — the orchestrator's three rows to the digit. It reproduces again after the parser was
optimised, and the three rows are also IN the CSV as `c2h6_T2000_f0, replicate 0, dt 0.25`.

**Validity check on the perturbation.** It does change the trajectory (`c2h6/T2000_f0`, mg3,
dt 0.25: rebuilds 42 / 98 / 34 / 64 / 52 / 80 over the six replicates) and it does not bias the
statistic: at dt 0.25 the unperturbed replicate 0 alone gives p = 0.0385 / 0.0615 / 0.0385 /
0.0385 / 0.0231 for mg / mg2 / mg3 / gauss / gauss_delivered against 0.0415 / 0.0492 / 0.0431 /
0.0462 / 0.0231 for the five perturbed replicates. So replicate 0 is a sample from the same
distribution as the others — which is precisely why a single one of them cannot decide anything.

## 11.2 PRIMARY result — no well form is distinguishable from `mg`

Fraction of the 780 trajectories per arm whose max per-step |dEpot| (rebuild steps excluded)
reaches 100 kJ/mol; `dp` is the paired difference against `mg` with its 95 % cluster-bootstrap
interval over the 130 cells (10 000 draws); `*` marks an interval that excludes 0.

| dt / fs | arm | bad / 780 | p | dp vs mg [95 % CI] | McNemar | median | p90 | p99 | max |
|---:|---|---:|---:|---|---:|---:|---:|---:|---:|
| 0.25 | **mg** | 32 | 0.0410 | (reference) | | 22.4 | 51.0 | 300.5 | 434.0 |
| 0.25 | mg2 | 40 | 0.0513 | +0.0103 [-0.0167, +0.0346] | 0.396 | 22.1 | 52.8 | 307.3 | 493.3 |
| 0.25 | mg3 | 33 | 0.0423 | +0.0013 [-0.0256, +0.0256] | 1.0 | 21.9 | 51.3 | 303.0 | 390.8 |
| 0.25 | gauss | 35 | 0.0449 | +0.0038 [-0.0231, +0.0333] | 0.788 | 21.6 | 49.3 | 294.6 | 530.0 |
| 0.25 | gauss+delivered | 18 | 0.0231 | -0.0179 [-0.0513, +0.0167] | 0.065 | 23.3 | 41.4 | 127.8 | 129.7 |
| 0.125 | **mg** | 35 | 0.0449 | (reference) | | 9.8 | 29.3 | 152.3 | 204.2 |
| 0.125 | mg2 | 28 | 0.0359 | -0.0090 [-0.0359, +0.0167] | 0.419 | 9.9 | 26.9 | 162.8 | 228.5 |
| 0.125 | mg3 | 40 | 0.0513 | +0.0064 [-0.0192, +0.0295] | 0.615 | 9.9 | 30.5 | 163.3 | 220.4 |
| 0.125 | gauss | 19 | 0.0244 | -0.0205 [-0.0474, +0.0064] | 0.037 | 9.3 | 23.8 | 143.5 | 232.9 |
| 0.125 | gauss+delivered | 12 | 0.0154 | **-0.0295 [-0.0564, -0.0026]\*** | 2.8e-06 | 10.0 | 19.9 | 136.6 | 136.8 |
| 0.0625 | all five | 0 | 0.0000 | +0.0000 [0, 0] | 1.0 | 5.4-5.8 | 9.7-13.0 | 51-73 | 72-100 |

**By the pre-registered rule: neither `mg2` nor `mg3` (nor `gauss`) is distinguishable from `mg`
at any of the three time steps.** Every interval contains 0; the smallest McNemar p among them is
0.0365 (`gauss`, dt 0.125), and that one is exactly the case the clustering matters for — the
unclustered test would call it significant at 0.05 while the cluster-bootstrap CI
([-0.0474, +0.0064]) does not, because `gauss`'s 12 failing cells are not 780 independent draws.
The
one arm that IS distinguishable is `gauss + delivered` at dt 0.125, and it is **better** than `mg`
(dp = -0.0295, CI [-0.0564, -0.0026]) — the pre-2026-09-19 default, i.e. the share rule, not the
well form, is the only lever this measurement resolves. That is the same conclusion package 7
reached by a different route ("the tail belongs to `conserving`, the well form only changes how
often it is visited") and it survives replication.

At the secondary 50 kJ/mol threshold the picture is identical (mg 0.1064, mg2 0.1141, mg3 0.1077
at dt 0.25; all CIs contain 0; only `gauss + delivered` separates). A post-hoc `mg3` vs `mg2`
comparison (NOT pre-registered) is also null at every dt and threshold: the largest effect is
dp = +0.0154 [-0.0013, +0.0321] at dt 0.125 / 100 kJ.

**At true dt 0.0625 fs no arm produces a single trajectory above 100 kJ/mol** in 780 tries — this
confirms package 10's 0.0625 recommendation on 30x the sample it was measured on, and it holds for
every well form and both share rules.

## 11.3 Why the aggregate and the single cell disagreed — the failures are sporadic

Number of the 130 cells that fail (>= 1 of their 6 replicates above 100 kJ/mol), and how many
replicates of a failing cell fail:

| dt | arm | failing cells | k=1 | k=2 | k=3 | k=4 | k=5 | k=6 | shared with mg (Jaccard) |
|---:|---|---:|---:|---:|---:|---:|---:|---:|---|
| 0.25 | mg | 21 | 15 | 4 | 1 | 0 | 0 | 1 | — |
| 0.25 | mg2 | 31 | 24 | 5 | 2 | 0 | 0 | 0 | 9 / 43 (0.21) |
| 0.25 | mg3 | 26 | 19 | 7 | 0 | 0 | 0 | 0 | 8 / 39 (0.21) |
| 0.25 | gauss | 22 | 17 | 3 | 0 | 0 | 0 | 2 | 9 / 34 (0.26) |
| 0.25 | gauss+delivered | 3 | 0 | 0 | 0 | 0 | 0 | **3** | 0 / 24 (0.00) |
| 0.125 | mg | 19 | 11 | 4 | 1 | 2 | 1 | 0 | — |
| 0.125 | mg2 | 18 | 12 | 3 | 2 | 1 | 0 | 0 | 3 / 34 (0.09) |
| 0.125 | mg3 | 31 | 24 | 6 | 0 | 1 | 0 | 0 | 11 / 39 (0.28) |
| 0.125 | gauss+delivered | 2 | 0 | 0 | 0 | 0 | 0 | **2** | 1 / 20 (0.05) |

Two things follow.

**(a) The failing cells barely overlap between arms** (Jaccard 0.09-0.28). "Which cell is the
worst cell" is mostly a property of the trajectory, not of the arm — which is exactly how package
10's aggregate and the orchestrator's single-cell check could both be right and still contradict
each other.

**(b) Under `conserving`, ~75 % of failing cells fail in exactly ONE of six replicates.** No cell
fails in 5 or 6 replicates for `mg2` or `mg3`, and only one does for `mg` (`c2h6_T2000_f6` at
0.25, `c2h6_T2000_f24` at 0.125). Under `delivered` it is the opposite: 3 cells fail, and each
fails in **all six** replicates (`ch4_H_T1000_f10`, `ch4_H_T2000_f10`, `ch4_H_T2000_f14` — the
artificial radical-adduct cells package 7 named). So the two share rules fail in qualitatively
different ways: `delivered` has a small, reproducible, identifiable set of bad configurations;
`conserving` has a diffuse, chaotic tail spread over a fifth of all cells. A single run can find
the first reliably and the second only by luck.

**The anchor cell is the textbook case.** `c2h6/T2000_f0` at dt 0.25, max per-step |dEpot| over
the six replicates:

| arm | rep 0 (unperturbed) | rep 1 | rep 2 | rep 3 | rep 4 | rep 5 |
|---|---:|---:|---:|---:|---:|---:|
| mg | 28.4 | 37.1 | 29.6 | 33.6 | 42.8 | 55.0 |
| mg2 | 32.2 | 35.9 | 29.1 | 47.4 | 33.1 | 37.2 |
| mg3 | **256.8** | 43.4 | 35.8 | 49.9 | 44.7 | 38.5 |
| gauss | 36.1 | 28.6 | 50.8 | 50.5 | 49.4 | 57.2 |

The 256.8 kJ that the orchestrator's spot-check found is a **singleton**: a 1e-5 A displacement of
the start removes it, five times out of five, and every other arm stays in the same 28-57 kJ band
on the same cell. It is a real event of a real trajectory — it is not evidence about `mg3`.

## 11.4 The `mg3` anchor event, attributed (n = 1)

`c2h6/T2000_f0`, replicate 0, mg3, dt 0.25, the +256.8 kJ step t = 2.98825 -> 2.98850 ps. All
twelve of that trajectory's >= 50 kJ steps lie in one 6 fs window (2.9845-2.9905 ps); this is one
event, not twelve.

**Which term.** Per-step decomposition (`terms2.py` on a verbosity-3 replay, which reproduces the
event to the digit): **bond +199.8 kJ/mol**, angle +54.1, non-bonded repulsion +19.4, bonded
repulsion -16.5, everything else < 0.1 — sum +256.8. Over-coordination is exactly 0 across the
step; every rebuild in the window reports `dE_jump = 0.000000 Eh`, so there is no hard swap.

**Which pair, and what the blend does.** `CURCUMA_BLENDDUMP=1`: the transition in flight is the
**break of a transient geminal H3-H4** (both hydrogens on C1; C2-H4 had broken 59 fs earlier),
window `[w_a, w_b] = [0.32344, 0.02000]` in bond ORDER. Across the jump step the corner weight `s`
goes **0.000000 -> 0.756416** — 76 % of the window in one 0.25 fs step — for a distance change of
only 2.09478 -> 2.22454 a0.

**What the corner gap is made of.** `CURCUMA_SHAREDUMP=1` prints both corners of the same energy
call. Bond-term difference (unbridged minus bridged) at the pre-jump geometry, per pair:

| pair | r / a0 | E bridged | E unbridged | gap / kJ |
|---|---:|---:|---:|---:|
| C1-H3 | 2.19426 | -0.13570535 | -0.08268342 | **+139.2** |
| C1-H4 | 2.53478 | -0.11871910 | -0.07242078 | **+121.6** |
| C1-H5 | 2.10546 | -0.17413552 | -0.16590890 | +21.6 |
| C1-C2 (+0.19) + the three C2-H (+0.72) | — | — | — | +0.91 total |
| **H3-H4 (the transitioning pair)** | 2.09478 | -0.05107869 | -0.05107869 | **0.00** |
| | | | **total** | **+283.3** |

`0.756416 x 283.3 = 214 kJ` against the measured bond change of +200 (the rest is the geometry
moving within the step). **The gap is the re-derivation of the force constant `fc` of the three
C-H bonds on the carbon that hosts the geminal pair** — the same quantity `BREAK_TAIL_STATUS.md`
isolated — and **not** the H-H pair's own well, which is bit-identical in both corners by
construction (a transitioning pair's well belongs to every corner, `ff_workspace.cpp`).

**How much of this is `mg3`?** Single points at the *same* geometry with the three well forms
(`CURCUMA_SHAREDUMP=1`): the six C-H wells are **bit-identical between `mg2` and `mg3`** and
differ from `mg` by **-0.19 % to +0.37 %** (C1-H3 0.16150247 mg / 0.16188876 mg2 = mg3; C1-H4
0.14192413 / 0.14164877). The only pair `mg3` changes appreciably is **C-C** (0.13395966 mg /
0.13470578 mg2 / 0.15957860 mg3, +19.1 % over mg) — and C-C contributes **0.19 kJ/mol** to the
283 kJ corner gap. The `mg3`-specific order dimension therefore touches nothing that carries this
event: on ethane the order table differs from the element-pair table only for C-C (H-H and C-H
have a single order-1 row, verified in `rev_well_table_v2.h` and against a live dump).

**Verdict: mechanism (a), without the "amplified by a deeper well" part.** It is package 7's
`conserving` mechanism verbatim — a transient geminal H-H whose fixed-width bond-order transition
window is about two MD steps wide in distance at 0.25 fs and 2000 K, so one step carries three
quarters of a 283 kJ corner gap. At the event geometry the well form is worth **under 1 kJ/mol of
that 283**, so it cannot be the amplifier; what `mg3` changed is which trajectory arrives at that
configuration. **n = 1** — this is one event of one trajectory, and 11.3 shows five sibling
trajectories of the same cell and arm that never reach it.

One caveat on the classic "halve dt and the event disappears" check: on this cell it does
(256.8 -> 12.0 kJ at 0.125), but a different dt is a *different* trajectory after the first
divergence — at 0.0625 the same cell/arm gives 1096 rebuilds and 128 steps >= 50 kJ with a max of
79.5. The dt statement is only safe as the 780-trajectory statement of 11.2, not per cell.

## 11.5 What this does and does not settle

- **Does**: on the react-MD smoothness axis, at 780 paired trajectories per arm and dt, `mg2` and
  `mg3` are **not distinguishable from `mg`** — neither better nor worse. Package 9's "`mg3` has a
  clearly worse tail" and package 10's "`mg3` is the best of the three" are BOTH unsupported; each
  read an ordering out of a sample of one trajectory per cell. The orchestrator's spot-check is a
  correct measurement of a singleton.
- **Does**: `gauss + delivered` is measurably smoother than any `conserving` arm at dt 0.125
  (the only result whose CI excludes 0), and its failures are reproducible where `conserving`'s
  are sporadic.
- **Does**: true 0.0625 fs bounds every arm (0 of 780 above 100 kJ/mol each).
- **Does not**: touch the two other costs of `mg2`/`mg3` — the conformer/S66 guard (1.0439 ->
  1.0543/1.0547) and the 0.0212 A equilibrium shift — or their class-A/class-D benefits. Those are
  static and stand as package 9 measured them. **The adoption decision is the operator's**; this
  package only removes the smoothness tail from the list of arguments in either direction.
- **Method note.** The decisive number was not any aggregate but the per-cell replicate table:
  ~75 % of failing cells fail in 1 of 6 replicates, and the failing-cell sets of two arms overlap
  by a Jaccard of 0.2. Any future react-MD tail claim should report that ratio before reporting a
  maximum — a maximum over one trajectory per cell is a draw from a distribution whose spread this
  package measures for the first time.
