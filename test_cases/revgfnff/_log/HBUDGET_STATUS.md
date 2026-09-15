# rev-gfnff 3a(ii): hydrogen keeps its nominal valence in the share budget (2026-09-15)

AI-generated, machine-tested. Worktree `curcuma-head`, branch `fix/h-valence-budget` off `30c186a4`,
uncommitted. Binary `build/curcuma` md5 **a4b9de6eb2bf89fbf2866f310daea678** (MAKE_EXIT=0, no new
warning in the touched files); every run dir carries that md5 in its `wall.txt`. Harness:
`.../scratchpad/hb/` (`run_one.sh`, `filt.py`, `perstep.py`, `analyse.py`, `fals.sh`, `pathB.py`,
`fdchk3.py`).

## 1. The change - `PARAM(rev_budget_fix_h, Bool, false, ...)`, DEFAULT OFF

    Pass 4 of prepareValenceShare, per atom:
      Z != 1 : Val_i = Val_Z + softplus_50(settled_i - Val_Z), dval_i = sigmoid_50(...)  (unchanged)
      Z == 1 : Val_i = Val_Z = 1 exactly,                      dval_i = 0                (flag on)

`dval` is read in exactly one place - the `dVal_i/d(settled)` factor of `applyValenceShareGradient`
- so a constant budget must contribute exactly 0 there; that is the whole chain rule, confirmed by
the FD in 4. 4 files, **+38 / -2 lines**: `gfnff.h` (PARAM + member), `gfnff_method.cpp`
(`setupRevSettings()` read, `rv.budget_fix_h`, `rev dump (flags)`), `ff_workspace.h`
(`RevSettings::budget_fix_h`), `ff_workspace_gfnff.cpp` (Pass 4).

**Off-arm provenance**: caffeine `-sp` revgfnff **-4.67352165**, gfnff **-4.67273707**; the 22-cell
grid off-arm reproduces `QP_STATUS` A.2 exactly - **960 rebuilds / max 471.1 / 3 events >= 50 kJ /
5 of 478 `begin_*` at s >= 0.99**. **Mechanism** (`CURCUMA_SHAREDUMP=1` on the delivered binary's
own off-arm frames of c2h6/T2000_f16 at t = 3206.5 / 3207.0 / 3207.5 / 3208.0 fs): Val(H4) off =
1.011 / 1.058 / **1.993** / 1.469, f(C1-H4) off = 1.000 / 0.0097 / **0.9999** / 0.454 - budget and
share swing over their whole range in 1 fs with no topology event; on, Val(H4) = **1.0000** at all
four.

## 2. (b) Per-step smoothness, the two runaway cells (`-md.print_frequency 1 -md.dump_frequency 1`)

max |Epot(t+dt) - Epot(t)| over adjacent printed steps, intervals containing a `REACT rebuild` line
excluded. **This is the real metric**: the rebuild `dE_jump` books 0.0 while the pair is crushed.

| cell | arm | per-step max / kJ | >= 50 kJ | T_max / K | min r(H-H) ever / a0 |
|---|---|---:|---:|---:|---:|
| c2h6/T2000_f16 | off | **2593.60** @3.20825 ps | 22 | 62 098 | **0.559** |
| c2h6/T2000_f16 | on | **59.26** @4.02625 ps | 17 | 5 086 | 1.770 |
| ch4_H/T2000_f10 | off | **1088.96** @1.42425 ps | 180 | **1.594e8** | **0.124** |
| ch4_H/T2000_f10 | on | 391.44 @step 1 (same in both arms), then **222.53** | 398 | 8 306 | 1.594 |

The residual 222.53 kJ is a 6500 K H2/CH4 vibration, not a discontinuity: Epot runs -0.6205/-0.6336/
-0.6425/-0.5896/-0.5678/-0.6526/-0.6647 Eh with T anti-correlated (2631 -> 6476 K).

## 3. (c) 22-cell grid, both arms, same binary

| arm | rebuilds | median | p99 | max dE_jump | >= 50 kJ | < 1 kJ | T_max / K | hard swaps |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| off (= delivered) | 960 | 0.00 | 1.8 | **471.1** | **3** | 98.1 % | **1.594e8** | **5 of 478** |
| on | 1186 | 0.00 | 1.7 | **48.6** | **0** | 97.7 % | **8 306** | **0 of 591** |

Per-step metric over the grid (~399 000 intervals per arm, t = 0.00025 ps init step excluded from
both): off **max 2593.60** kJ, 322 events >= 50 (143 outside ch4_H/T2000_f10); on **max 222.53**,
485 (**88** outside that cell). The all-22 count goes the wrong way for one measured reason: the off
arm's ch4_H/T2000_f10 **destroys itself at 1.42575 ps** (T > 1e5 K, Epot flat at -0.30..-0.37 Eh for
3.6 ps) and stops doing chemistry - 129 of its 180 events fall before the blow-up, 91 per ps against
the on arm's 80 per ps over the full 5 ps, and the on arm runs 126 rebuilds there against 27.
10 of 22 cells differ between the arms at all; the other **12 are identical**.

## 4. (d) Falsifiers, flag on vs off

| falsifier | verdict |
|---|---|
| NH4+ / H3O+ / CH5+ / ClO4- (12 digits) | **dE = 0.00000000 kcal, dG_max diff 0.00e+00** |
| BF4- compressed 1.143 A (q -1) | +0.118269033480 both, **bit-identical** (F untouched) |
| BF4- realistic 1.394 A (q -1) | -1.470184046880 both, **bit-identical** |
| rkt06_h_h2, 11 points | rms **2.714** both, max dE **0.000e+00 Eh**; barrier +3.41@4, pt5/pt10 -3.06/-3.79 |
| equilibrium 2x2 (share x fix_h) x 5 systems | **20/20 dE = +0.000000000 kcal, dG_max 0.00e+00** |
| `gfnff` caffeine / benzene | -4.672737068614 / -2.362725526194, `dump_params` md5 d297bc3b… / 77c134bb… = the records |

The softplus floor saturates as expected: an H with one partner has Val = 1.0139 off, 1.0000 on, and
`f = Val/w > 1` clips to exactly 1 either way. revgfnff caffeine -4.673521653477, benzene
-2.363224128930 = the CIJ/VALFIX records. **FD gradient** (central, dx 1e-4 A, fresh dir per
displacement, `-gfnff.cache_topology false`, `-batch_reuse_topology false`), worst |analytic - FD|
in Eh/A: rkt06 pt10 **1.194e-08** on (1.194e-08 off = the committed value); c2h6 off-arm frame
t = 3207.0 fs **1.210e-07** on (7.014e-06 off). Both far below 1e-6. The off arm's 7.0e-6 is FD
truncation across the beta = 50 softplus: dx 3e-4 / 1e-4 / 5e-5 give 6.313e-05 / 7.014e-06 /
1.753e-06 - exact dx^2 scaling.

## 5. Tests

`ctest -R gfnff` **60/65**, `ctest -R cli_simplemd_` **17/22**, the same 5 failures in both:
`cli_simplemd_16/18/19` (the known pre-existing rev ones) plus `cli_simplemd_08/09` (acetic-acid
dimer, NVE drift 0.447 Eh >= 0.10). 08/09 run **`-method gfnff`**, bit-identical to 12 digits and to
the `dump_params` md5 per section 4, so they are pre-existing in this worktree's build.

## 6. Open, and two harness caveats

- Default OFF; making it the default is an operator decision. On this evidence the on arm is
  strictly better on every smoothness measure and costs nothing on any falsifier.
- Not fixed: the shared first-step 391 kJ jump of ch4_H/T2000_f10, and the 88 remaining >= 50 kJ
  steps of the other 21 cells (within the thermal amplitude at 5000-7000 K and dt = 0.25 fs).
- **`bt/run_one.sh`'s verbosity-3 filter silently drops almost every status row**: it tests `$1` of
  the RAW line, which can carry a leading ANSI reset - 20 000 rows collapsed to 2 and the per-step
  metric came out empty. `filt.py` strips ANSI first.
- **A parallel launch overwrote a run directory** (an aborted launch's processes were still writing),
  and the "on" cell then reproduced the "off" trajectory exactly. Caught by fingerprinting the
  rebuild count against the grid arm (135/204 for c2h6/T2000_f16, 27/126 for ch4_H/T2000_f10).
