# WORK_STATUS — rev-gfnff work packages 1-7 (2026-09-18 / 19)
Packages done: 7/7

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
