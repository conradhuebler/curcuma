# WORK_STATUS — rev-gfnff work packages 1-5 (2026-09-18)
Packages done: 4/5

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
