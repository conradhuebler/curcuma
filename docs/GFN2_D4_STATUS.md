# GFN2-D4 status (alignment to tblite)

Status of the native GFN2-xTB D4 dispersion path (`-method gfn2`, canonical
backend in `release/`) vs the tblite reference, as of 2026-05-28.

## What was fixed

| Commit  | Fix |
|---------|-----|
| `07d0420` | **D4 CN alignment**: GFN2 now uses dftd4's EN-weighted covalent CN (`dispersion/d4_ncoord.{h,cpp}`, port of cpp-d4 `ncoord_d4`, no log-cap). Switched on via `D4ParameterGenerator::setD4CovalentCN(true)`. GFN-FF keeps its log-capped erf-CN. |
| `0a56abf` | **Step-3a diagnostic** (`test_cases/sqm_reference/diag_curcuma_d4_potential.cpp`): evaluates `dE_D4/dq` at tblite's converged Mulliken charges, compares to `vat_tblite − gf2.vat_extra` (implied tblite D4 SCF potential). Localises whether the residual is in the D4 model. Now part of `ctest -L d4_diag`. |
| `19ce842` | **`set_refalpha_gfn2` zeta scaling**: curcuma's `computeC6Reference` was missing the inner-shell `zeta(ga, hardness(refsys)·gc, zeff(refsys), refh(ir,num)+zeff(refsys))` factor that scales `secaiw` in dftd4's `set_refalpha_gfn2_num` (reference.f90:343-377). New table `d4_refh_charges` (118×7, verbatim from dftd4 reference.inc) + gating by `m_use_d4_covalent_cn` so GFN-FF stays on its unscaled C6 table. Cache invalidation on flag toggle. |

## Where we are now

**Single-point energy vs tblite** (`release_tblite/dumps/*_gfn2.json`):

| Molecule | Original | After CN | After CN + α-zeta | Notes |
|----------|----------|----------|--------------------|-------|
| H₂       | -3.3e-9  | -3.3e-9  | -3.3e-9            | exact (refh=0 for H ⇒ zeta=1) |
| H₂O      | +4.68e-5 | +3.41e-5 | **-0.67 µEh**      | sub-µEh (70× ↓) |
| NH₃      | +8.64e-5 | +7.02e-5 | **-1.42 µEh**      | sub-µEh (60× ↓) |
| CH₄      | -3.68e-5 | +1.43e-5 | -1.79e-5           | **C-path overcorrection** |
| triose   | n/a      | n/a      | -3.21 mEh          | 66 × CH₄-style residual |
| complex  | n/a      | n/a      | DIVERGES           | DIIS charge sloshing, not D4 |

`d4_diag` ctest (regression for `max|diff|` per atom in the D4 SCF potential at
tblite's charges; tolerances tighten as we improve):

```
H2     : tol 1.5e-6
H2O    : tol 1.0e-6
NH3    : tol 1.5e-6
CH4    : tol 1.5e-5   ← to be tightened once C-path is fixed
triose : tol 1.5e-4   ← to be tightened once C-path is fixed
```

Run: `cd release && ctest -L d4_diag --output-on-failure`. Tighten the
tolerances in `test_cases/sqm_reference/CMakeLists.txt` as the residual closes.

## Open issues

### 1. C-path overcorrection (CH₄ +18 µEh, triose +3.2 mEh)

After the α-zeta fix, curcuma's `dE_D4/dq` is still ~2–5 % too high on C-H
references (CH₄ C ratio 1.04, H ratio 1.02; triose ratios up to 1.08).
The α-zeta correction is strong (~0.50) on O/N references (large refh) and
weak (~0.92) on C references (small refh, e.g. CH₄ refh ≈ 0.03); on the C side
the fix slightly overshoots.

The per-atom error sums for triose:
- max|dq| (Mulliken vs tblite) = 1.55e-4 on a C atom
- mean|dq| ≈ 3.1e-5
- D4-vat max|diff| = 1.20e-4 per atom
- 66 atoms × ~50 µEh/atom = ~3 mEh total — exactly the observed ΔE.

Likely candidates for the remaining ~5 % on the C side:
- a missing layer in the C reference data (per-ref `ascale` / `sscale` subtleties
  for the CH/CH₂/CH₃/CH₄ refs);
- the `max(α_corrected, 0)` clamp interacting differently for small refh;
- an alphaiw extraction discrepancy specific to C reference states.

Probe: per-pair `C6(C,H)` and `C6(C,C)` at tblite's charges — extend
`diag_curcuma_d4_potential` to also print per-pair C6 and diff against tblite's
`d4%c6` (would need a small dump-tblite patch to expose `model%c6`, or compute
the per-pair C6 from `gwvec ⊗ c6ref` independently).

### 2. complex (231 atoms) — SCF divergence, **not** a D4 issue

```
iter 0:  -328.10 Eh   (close to tblite -329.53)
iter 1:  -326.60 Eh
iter 2:  -307.31 Eh   (+20 Eh)
iter 3:   -94.35 Eh   (+213 Eh)
iter 5: +2419.24 Eh   (sign flip after DIIS engages)
iter 149: +2764.9 Eh, max|dq|~2–6 (oscillating)
```

Classic charge-sloshing pattern: the bare-`H0` guess + Coulomb iterations
diverge when DIIS engages at `diis_start = 5`. The same mechanism the
`diis_start=5` warmup fixed for HCN reappears at this size with many polar
groups. Independent of the D4 model — orthogonal next investigation: longer
damping warmup, larger damping, level shifting, or a better initial guess
(e.g. EEQ-charge-derived diagonal H0 shift).

## Architecture pointers

- D4 entry points for GFN2: `xtb_native.cpp::calcDispersionEnergy` and
  `addDispersionPotential` — both set `setD4CovalentCN(true)`,
  `setUseD4SingleShotEEQ(false)`, `setTopologyCharges(m_wfn.q_at)`.
- D4 evaluator (per-reference path): `dispersion/d4_evaluator.cpp:188-228` ⇒
  uses `D4ParameterGenerator::weightedC6Gfn2`.
- C6 reference matrix: `d4param_generator.cpp::computeC6Reference` (zeta-scaled
  when `m_use_d4_covalent_cn` is on).
- New reference table: `d4_refh_charges` (in `d4_reference_data_fixed.cpp`).
- Diagnostic tool: `test_cases/sqm_reference/diag_curcuma_d4_potential.cpp`
  (`ctest -L d4_diag`).
- Reference Fortran: `release_tblite/_deps/dftd4-src/src/dftd4/{model,reference,
  reference.inc}.f90` — `set_refalpha_gfn2_num`, `set_refq_gfn2_num`,
  `weight_references`, dispmat construction in `tblite/disp/d4.f90`.
