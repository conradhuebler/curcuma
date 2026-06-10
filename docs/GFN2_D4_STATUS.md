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

disp+elec residual at fixed tblite density (`diag_curcuma_energy_components`):

| Molecule | + α-zeta | + refcovcn (2-body) | + ATM (3-body) | Notes |
|----------|----------|---------------------|----------------|-------|
| H₂       | -3.3e-9  | -3.3e-9   | -3.3e-9   | exact |
| H₂O      | -0.67 µEh| 3.0e-10   | **2.7e-10** | essentially exact |
| NH₃      | -1.42 µEh| 1.0e-9    | **1.4e-10** | essentially exact |
| CH₄      | -1.79e-5 | 5.9e-9    | **2.4e-11** | C-path resolved |
| triose   | -3.21 mEh| -0.78 mEh | **5.9e-9**  | refcovcn + ATM → sub-nEh |
| complex  | DIVERGES | DIVERGES  | DIVERGES  | DIIS charge sloshing, not D4 |

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

**What we ruled out (cross-checked line-by-line against `dftd4-src/src/dftd4/`):**
- `weight_references`, `zeta`, `dzeta` formulas — identical.
- Constants `ga=3`, `gc=2`, `wf=6` — identical.
- `effective_nuclear_charge`, `chemical_hardness` (curcuma's `zeta_zeff`/`zeta_c`)
  — bit-identical for H/C/N/O (and likely for all 118).
- `refq` (GFN2 reference charges) — bit-identical (verified for H/C/N/O).
- `refh` (newly-added `d4_refh_charges`) — verbatim from dftd4 reference.inc.
- `set_refgw` multi-gaussian `ngw` — formula identical.
- `refsys` (curcuma `d4_refsys_data`) — bit-identical for C/N/O.
- Casimir-Polder integrator (trapezoidal over 23-point freq grid) — identical
  closed-form weights.
- The zeta-correction is gated by `m_use_d4_covalent_cn` with cache invalidation
  on toggle (verified GFN-FF unaffected).

**Reference-C6 matrix probe — DONE (2026-05-28).** Built a standalone Fortran
dumper (`test_cases/sqm_reference/dump_dftd4_c6.f90`) that constructs the exact
GFN2 model (`new_d4_model(..., ref=d4_ref%gfn2)`) and writes `model%c6(iref,
jref,isp,jsp)`, plus a C++ differ (`diag_curcuma_d4_c6.cpp`) against curcuma's
`m_c6_flat_cache`. `ctest -L d4_c6` (10 tests). Findings:

1. **Reference C6 WAS wrong** — bit-identical on hcount=0 references, but up to
   **57 % off (O) / 18 % off (C)** on the `hcount≠0` (α-zeta-corrected)
   intermediate references. Root cause: `precomputeC6ReferenceMatrix()` runs in
   the ctor with `m_use_d4_covalent_cn=false` (unscaled). Native GFN2 sets the
   flag true afterwards, which *invalidates* the cache, but `GenerateParameters`
   never rebuilt it (only `GenerateDispersionPairsNative`, the GFN-FF path, did).
   So `set_refalpha_gfn2` was effectively dead for the GFN2 energy path.
2. **Fix**: added `if (!m_c6_reference_cached) precomputeC6ReferenceMatrix();`
   to `GenerateParameters` (`d4param_generator.cpp`). Reference C6 now
   bit-identical to dftd4 (~1e-10). GFN-FF unaffected (flag stays false → cache
   never invalidated → no-op; GFN-FF total energy unchanged, verified).
3. **But this does NOT resolve the residual.** The dispersion energy is
   *byte-identical* before and after the fix (CH4 disp = -0.0006801449,
   triose = -0.0713917953, both). The corrected (hcount≠0) references carry
   **~zero CN/charge weight** for H/C/N/O at typical geometries, so the
   reference-C6 error never reached the energy. It was a genuine latent bug
   (would bite high-CN systems), now fixed and regression-locked, but it is
   **ruled out as the CH4/triose residual cause.**

**So the residual is NOT in the reference C6 matrix.** The next probe found it.

**Weighted per-atom C6 probe — ROOT CAUSE FOUND (2026-05-28).** Built
`test_cases/sqm_reference/dump_dftd4_atomic_c6.f90` (reproduces tblite's
`weight_references` + `get_atomic_c6` at given charges) + C++ differ
`diag_curcuma_atomic_c6.cpp` (curcuma's `weightedC6Gfn2`). For triose: the D4
covalent **CN is bit-identical** (max diff 6.75e-12), but the **weighted
per-atom C6 diverges** — O-O ~0.5%, O-C ~3-6%, **C-C up to 14%, curcuma
systematically higher**, scaling with reference count (carbon = 7 refs, worst).

Drilling into the carbon `gwvec` (CN-Gaussian × charge-zeta): refcn, refq, the
zeta formula, eta, zeff, ga/gc/wf, ngw, and the reference C6 ALL match — yet
curcuma's `gwvec` was far less peaked than tblite's (CH4 carbon ref4: curcuma
0.74 vs tblite 0.99). **The bug: dftd4 uses TWO distinct reference-CN tables —
`refcn` for the `ngw` bucketing (`set_refgw`, `nint(refcn)`) and `refcovcn`
(covalent CN) for the actual Gaussian weighting (`set_refcn` → `model%cn` →
`weight_references`).** Carbon `refcovcn` = [0, 0.919, 1.908, 2.831, **3.749**,
2.918, 0.856] vs `refcn` = [0, 0.987, 1.998, 2.999, **3.984**, 3.142, 1.000].
Curcuma's `m_refcn` is loaded from `refcn` and used for BOTH the ngw bucketing
AND the Gaussian — but the Gaussian must use `refcovcn`. At CN≈3.70 (sp3 C),
`refcovcn₄=3.749` (d=−0.05, sharply peaked) vs `refcn₄=3.984` (d=−0.29, broad),
which is exactly the 14% over-broadening. **Proof:** recomputing the carbon
`gwvec` with `refcovcn` for the Gaussian (keeping `refcn` for ngw) reproduces
tblite's `gwvec` to ~1e-5 across all 7 references.

This is the C-path / triose residual: the affected references DO carry CN/charge
weight (unlike the reference-C6 error), so the wrong `refcovcn` directly shifts
the per-atom C6 → the dispersion energy → the 3.2 mEh.

**Fix LANDED (2026-05-29).** Ported the `refcovcn` (118×7) table from dftd4
`reference.inc` (auto-extracted, 372 entries) into `d4_reference_data_fixed.cpp`
as `d4_refcovcn` → `m_refcovcn`, and switched the CN-Gaussian in
`weightedC6Gfn2::buildRefW` to `weight_cn(cn − refcovcn)` while keeping `m_refcn`
only for the `ngw` rounding (exactly mirroring dftd4 set_refgw/set_refcn). The
GFN-FF path (`getChargeWeightedC6`/`precomputeGaussianWeights`, wf=4) still uses
`m_refcn` — untouched, GFN-FF total energy byte-unchanged (cli_gfnff_01
−37.24064873, identical before/after).

Results (combined `disp+elec` vs tblite, `diag_curcuma_energy_components`):

| Molecule | before fix | after fix |
|----------|-----------|-----------|
| H₂O      | 0.67 µEh  | **3.0e-10** |
| NH₃      | 1.42 µEh  | **1.0e-9**  |
| CH₄      | 1.79e-5   | **5.9e-9** (C-path overcorrection RESOLVED) |
| triose   | 3.21 mEh  | **0.78 mEh** (4× ↓) |

Weighted per-atom C6 (`diag_curcuma_atomic_c6`) now bit-identical to tblite
`get_atomic_c6` (triose max|rel| 3.6e-11, was 14%). Regression: `ctest -L
"d4_c6;d4_diag;gfn2_align"` all green; `gfn2_align` tolerances tightened
(H₂O/NH₃ 1e-8, CH₄ 1e-7, triose 1e-3).

**Remaining triose 0.78 mEh — AUDITED & EXPLAINED (2026-05-29): missing ATM
three-body term.** curcuma's `D4Evaluator` omits the ATM 3-body term entirely
(`D4Params::s9` is documented "not used by this evaluator yet"), while tblite's
GFN2 uses **s9=5.0** (`xtb/gfn2.f90:54`, disp3 cutoff 25.0). Built
`dump_dftd4_atm.f90` (dftd4 `get_dispersion3` at q=0 C6 — the ATM is
charge-independent, `disp.f90:110`). tblite's pre-SCF `dispersion` container
equals the ATM term **exactly** (10 digits), and the ATM exactly equals the
remaining `disp+elec` residual:

| Molecule | tblite ATM = pre-SCF disp container | = remaining residual |
|----------|-------------------------------------|----------------------|
| H₂O      | +3.05e-11 | yes |
| NH₃      | +8.78e-10 | yes |
| CH₄      | +5.87e-9  | yes |
| C₆H₆     | +7.07e-6  | — |
| triose   | **+7.7945e-4** | **= the 0.78 mEh** |

So tblite splits GFN2 D4 as: ATM 3-body → pre-SCF `dispersion` container;
SCC 2-body → `eelec`. curcuma puts the (now-correct) 2-body in `m_E_dispersion`
and had **no ATM at all**. Tool: `dump_dftd4_atm` (built under USE_TBLITE).

**ATM LANDED (2026-05-29).** `D4Evaluator::computeATM` ports dftd4
`get_atm_dispersion` (energy + radial gradient + CN chain rule), evaluated at the
q=0 reference C6 (charge-independent → no q-response, matches tblite
`get_dispersion_nonsc`). Wired into `xtb_native::calcDispersionEnergy` (GFN2 only;
GFN-FF uses the per-pair path and is structurally untouched — total byte-unchanged).
Results: ATM energy matches tblite to 10 digits; **triose disp+elec 0.78 mEh →
5.9e-9** (H₂O 2.7e-10, NH₃ 1.4e-10, CH₄ 2.4e-11). ATM gradient validated against
finite differences to machine precision (CH₄ 7e-20, C₆H₆ 5e-16, triose 6e-13
Eh/Bohr). `gfn2_align` tolerances tightened (triose 1e-6). **From the original
3.2 mEh triose residual, the native GFN2 D4 is now sub-nEh vs tblite at fixed
density** (refcovcn 2-body fix + ATM 3-body).

NB: `test_xtb_gradient` FAILs (~6-8e-4) for C₆H₆/CH₃OCH₃/butane on **both ngfn1
and ngfn2** — a PRE-EXISTING full-gradient inaccuracy unrelated to D4/ATM (GFN1
has no D4); the automated set (H₂O/CH₄/NH₃) passes. Separate investigation.

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
  when `m_use_d4_covalent_cn` is on); cache (`m_c6_flat_cache`) rebuilt in
  `GenerateParameters` when invalidated (2026-05-28 fix).
- New reference table: `d4_refh_charges` (in `d4_reference_data_fixed.cpp`).
- Diagnostic tools: `test_cases/sqm_reference/diag_curcuma_d4_potential.cpp`
  (`ctest -L d4_diag`) and `diag_curcuma_d4_c6.cpp` + `dump_dftd4_c6.f90`
  (reference-C6 vs dftd4 `model%c6`, `ctest -L d4_c6`).
- Reference Fortran: `release_tblite/_deps/dftd4-src/src/dftd4/{model,reference,
  reference.inc}.f90` — `set_refalpha_gfn2_num`, `set_refq_gfn2_num`,
  `weight_references`, dispmat construction in `tblite/disp/d4.f90`.
