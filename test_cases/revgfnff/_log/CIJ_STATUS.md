# rev-gfnff stage 3a(ii): valence share of the bond well - status (2026-09-13)

AI-generated measurement job, machine-tested. src changed, uncommitted, branch `reactff2-llm`, binary
`build_rev/curcuma` (MAKE_EXIT=0, no new warnings in the touched files). Harness + raw JSON in the
session scratchpad `.../scratchpad/cij/` (`pathB.py`, `decomp.py`, `sweepB.py`, `toggle2.sh`,
`fdcheck2.py`, `jump/`). `gfnff` bit-identical throughout.

## 1. What changed (5 files, 0 new element data, 2 new rev PARAMs)

- `ff_workspace_gfnff.cpp` `prepareValenceShare()` (main thread, per topology corner, before the
  partitions): `sum_i = Σ_k w_ik` over the corner's own bond list, `w = revWeight` = the term-weight
  switch the well is already multiplied with (bo_center 2.0 / k -7.5, w_join-shifted).
- `calcBonds()`: `c_ij = (f_i + f_j)/2`, `f_i = shareClip((Val_i - sum_i + w)/w)`. `shareClip` is the
  C1 smoothstep `x^2(3-2x)` on [0,1], **exactly 0 and exactly 1 outside -> no new parameter** (the
  clip width is the unit interval of f; that is the "one soft-clip width" of the review). Multiply the
  well by `c`; accumulate `acc.dEdshare(i) += -K f_i'/2` (`K = k_b e^{-a dr^2}`); the own-`w`
  derivative enters through `wfac = c + w dc/dw`.
- `applyValenceShareGradient()` (main thread after the reduction): `(Λ_i + Λ_j) dw/dr dr/dx` over the
  same bond list, `Λ_i = Σ_pairs dE_pair/d sum_i` - the `calcOverCoordination` pattern.
- `Val_i = revValence(Z)` (the existing sigma-valence table of the over-coordination term).
- `gfnff_method.cpp`: **"an H is never sp"** - the bridging-atom rule (`sp H or group 7 -> bstrength
  = bstren[1]*0.30`, `gfnff_ini.f90:1170-1201`) no longer fires for hydrogen in rev mode.
- PARAMs: `rev_valence_share` (Bool, true), `rev_h_not_sp` (Bool, true); both gated on `rev_enabled`,
  so `-method gfnff` is untouched.

## 2. Falsifier - class-B `rkt06_h_h2`, 11 r2SCAN-3c points, kcal/mol vs the reactant

| arm | path rms | barrier | TS points 5 / 10 (ref +2.53 / +2.57) |
|---|---:|---|---|
| neither (pre-change) | 15.79 | +39.80 @pt5 | +39.80 / +39.22 |
| c only | 29.69 | +72.40 @pt5 | +72.40 / +71.90 |
| H-never-sp only | 48.70 | +3.41 @pt4 | **-111.1 / -112.2** |
| **both (final)** | **2.71** | **+3.41 @pt4** (ref +2.57 @pt10) | -3.06 / -3.79 |

**rms 15.79 -> 2.71, target <= 10 met.** The two rules are complementary, neither works alone: without
the share the H-not-sp wells are two full wells at the TS (-111 each vs the reactant's -102 -> 120 kcal
too low); without H-not-sp the wells are 3.3x too shallow (fc -0.0538 vs -0.1788 Eh for the same bond
in H2) and the share halves them again (+72).

## 3. Class-B set, rms before -> after (kcal/mol)

| system | rms before | rms after | bar_ref | bar before | bar after |
|---|---:|---:|---:|---:|---:|
| rkt06_h_h2 | 15.79 | **2.71** | 2.57 | 39.22 | -3.79 |
| rkt01_h_hcl_h2_cl | 28.65 | 25.48 | 6.83 | 6.77 | 6.77 |
| hclhts_h_hcl | 14.77 | 17.16 | 13.88 | 48.34 | 53.95 |
| hfhts_h_hf | 11.93 | 22.51 | 32.46 | 59.86 | 84.98 |
| n2h_h_n2h2 | 41.67 | 39.24 | 7.96 | 46.03 | -39.03 |
| rkt14_h_oh_h2_o | 35.16 | 35.16 | 3.37 | 52.61 | 52.61 |
| rkt10_f_h2_hf_h | 39.53 | 39.53 | 22.50 | 6.93 | 6.93 |
| n2_h_n2h | 47.27 | 47.27 | 6.52 | 6.77 | 6.77 |
| hf2ts_h_f2_hf_f | 69.95 | 69.95 | 0.01 | -0.23 | -0.23 |
| px13_h2o_2 / hf_2 / nh3_2 | 604/685/235 | 627/709/416 | 1885/1862/1694 | - | - |

(`n2_h2_n2h2`, `n2h2_h2_n2h4`, `n2h4_h2_2nh3` have no `points.xyz` in the ref dir -> not measured.)
The other 11 systems are unchanged or a few kcal worse - they are the positive-TS-error regime where a
factor <= 1 can only raise the transition state (see 6).

## 4. Equilibrium + bit-identity (final binary)

- c on/off x H-not-sp on/off, **2x2**: caffeine, benzene, 2h2, n2_3h2, ch4_H frame 0 -
  **dE = +0.00000000 kcal/mol and dG_max = 0.00e+00 Eh/Angstrom for every combination** (12 digits).
  Every bond has c = 1 exactly: `sum_{k != j} w < Val_i - w` for a saturated atom because w <= 1.
- Absolute values equal the pre-change build: caffeine -4.673521653477, benzene -2.363224128930 Eh
  (pre-change doc: -4.67352165 / -2.36322413).
- `gfnff`: caffeine -4.672737068614, benzene -2.362725526194 Eh and every gradient component
  identical to the pre-change binary (12 digits, `curcuma_base`); `-gfnff.dump_params` md5 identical
  (caffeine d297bc3b91ad75224f224b0cb4c3d189, benzene 77c134bbfd37284059c15d89c6125efc).

## 5. Rebuild smoothness (22 cells: c2h6/ch3nh2/ch4_H x 1000/2000 K x 3 start frames, 5 ps, dt 0.25)

| arm (same binary) | cells | rebuilds | median | max | p99 | <1 kJ | <5 kJ | >=50 kJ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| before (both toggles off) | 22 | 862 | 0.00 | **641.3** | 430.6 | 0.883 | 0.917 | 47 |
| final | 22 | 1076 | 0.00 | **10.70** | 1.20 | 0.986 | 0.998 | **0** |

Documented baseline (`POLY_JUMP_BASELINE.md`, older binary, 45 cells incl. 3000 K): pooled median 0.00,
max 423.5. **No regression - 5x below the baseline max, and the 1000/2000 K max falls 60x against the
same-binary pre-change arm.** Worst cell: ch3nh2 2000 K f16, 10.7 kJ.

## 6. Gradient

FD of the full energy (central, dx = 0.002 A, fresh dir per displacement, `-batch_reuse_topology
false`), final binary: worst |analytic - FD| **9.8e-06** (rkt06 TS, c = 0.5) and **2.05e-05** Eh/A
(a 3-atom frame with f in the middle of the clip band). With `applyValenceShareGradient()` disabled the
same geometry degrades to **6.45e-05** - the sum chain rule is present and correct.

## 7. Tests

- `ctest -R gfnff` in `build_rev/`: **65/65 passed, 91.12 s, exit 0** (pre-change: 65/65, 91.6 s).
  Includes `gfnff_rev_fd`, the 8 react tests and `cli_simplemd_19_gfnff_rev_form_refuses_hbond`
  (**passed**, 0.72 s: water dimer forms zero bonds in react mode and stays within 0.01 kcal/mol of
  static). No test file or `test_cases/cli/CMakeLists.txt` touched.

## 8. Open / limits (measured, not fixed)

- **`Val_i` is the nominal `revValence(Z)`.** The review's `max(Val_Z, settled topology bonds)` is
  **vacuous**: with `Val_i = max(Val_Z, n_i)` and `sum_{k != j} S_ik <= n_i - 1` one gets
  `f_i >= 1/S >= 1` for every atom, i.e. `c == 1` always - which also contradicts the review's own
  worked example (`c_CH = 0.54` needs `Val_H = 1` for an H with two partners). Chosen reading: nominal.
  **Consequence, measured**: a hypervalent equilibrium loses ~half of every bond to that atom -
  NH4+ +0.345 Eh (+216 kcal/mol), H3O+ +0.206 Eh (+129) vs the same run with the share off (>6-coordinate
  metals, ClO4-, ClF3 likewise). A "settled" count would fix those but would have to switch inside the
  1.20-1.30 r_eq band where exchange transition states live, so it was not adopted.
- `c` is effectively a step of the topology: a pair enters the bond list with `w = 0.96`, so `f` lands
  at ~0.001 and never crosses the middle of the clip. All smoothness therefore rests on the stage-1b
  corner blend - measured above, it holds (0 events >= 50 kJ).
- Untouched: the well form (3a iii) - at the rkt06 TS the Gaussian is still on its left wall; the
  remaining class-B systems (hfhts, rkt14, n2_h_n2h, rkt10) are unchanged by this rule and are the
  3a(iii)/3b targets. Stage 2 charges, rev on GPU, `docs/` - not touched.
- Reproduce: `python3 .../scratchpad/cij/pathB.py build_rev/curcuma rkt06_h_h2 out.json`;
  `sweepB.py BIN out.json`; `jump/run_one.sh` (see `POLY_JUMP_BASELINE.md` section 6 for the grid).
