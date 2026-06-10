# WP2 — GFN1 to 1e-8 Eh vs tblite (localize & close residuals)

> Part of [SQM_VALIDATION_ROADMAP.md](SQM_VALIDATION_ROADMAP.md).
> Status: **ADD** (specified, not started). Depends on **WP3** (per-container diff).
> 🤖 AI-generated plan.

## Target

Bring native GFN1 (`curcuma::xtb::XTB(GFN1)`, `xtb_native.cpp` + `xtb_h0.cpp` /
`xtb_coulomb.cpp` / `xtb_scf.cpp`) to **|E_native − E_tblite| ≤ 1e-8 Eh** on the
full validation set, removing its 12 xfails from the suite (WP1). Also re-examine
the single GFN2 `complex` 6.95e-5 residual.

## Measured residuals to close (2026-05-30)

| molecule | \|dE\| (Eh) | sign (native−tblite) |
|---|---:|:--:|
| He2 | 1.5e-8 | − |
| H2O | 2.1e-6 | + |
| H2 | 1.3e-6 | − |
| NH3 | 1.6e-5 | − |
| LiH | 4.3e-5 | + |
| HCN | 4.4e-5 | + |
| CH4 | 4.7e-5 | − |
| C6H6 | 4.8e-4 | − |
| acetic_acid_dimer | 1.4e-3 | + |
| caffeine | 1.8e-3 | − |
| triose | 8.4e-3 | − |
| complex | 3.5e-2 | − |

Two regimes:
- **Small covalent (H2…HCN), ~1µEh–50µEh:** likely H0 self-energy / repulsion /
  CN-counting precision differences vs tblite (parameter or formula detail), not
  dispersion (too small, mixed sign).
- **Larger (C6H6, acetic, caffeine, triose, complex), 0.5–35 mEh, growing with
  size and mostly over-binding:** a size-extensive term is off. Candidates:
  isotropic/3rd-order Coulomb, GFN1 D3 (BJ) scaling/cutoff, or the halogen-bond /
  repulsion container. **Cause is open** — do not assume; localize first (WP3).

## Method (localize before fixing)

1. **Per-container diff (uses WP3).** Run `diag_curcuma_energy_components` (GFN1)
   on triose + complex (the largest residuals) with tblite's injected density, to
   attribute the mEh to repulsion / dispersion / electronic / halogen. The
   container carrying the residual names the bug. Extend the `gfn2_align`-style
   audit to GFN1 (WP3 delivers this).
2. **Cross-check the existing Phase-0 diagnostics for GFN1** already present in
   `test_cases/sqm_reference/` (`test_xtb_overlap`, `test_xtb_h0`,
   `test_xtb_coulomb`, `test_xtb_scf_snapshot`) — currently only run on the small
   Phase-0 dump set (H2/He2/LiH/H2O/CH4/NH3/C6H6). Confirm overlap/H0/Coulomb are
   bit-clean for GFN1 there; if a small-molecule residual (e.g. CH4 47µEh) shows
   up in H0 or Coulomb, that isolates the parameter/formula gap cheaply.
3. **Dispersion check.** GFN1 uses D3(BJ) (`D3ParameterGenerator::createForGFN1`,
   `calcDispersionEnergy` GFN1 branch). Diff native GFN1 D3 energy against tblite's
   D3 container on triose/complex. The over-binding sign means if D3 is the cause
   it is *too attractive* (wrong s8/a1/a2 or cutoff), not missing.
4. **Fix the localized term**, re-measure, tighten the corresponding `sqm_val_*`
   xfail to a normal 1e-8 test (WP1 mechanism flips it automatically).

## Suggested order (cheap → expensive)

1. Small-molecule H0/Coulomb diff (steps 2) — fast, isolates parameter gaps.
2. triose container audit (step 1) — one molecule, names the dominant term.
3. D3 container diff (step 3) if dispersion implicated.
4. complex last (largest, slowest) to confirm size-extensivity of the fix.

## Done when

- All 12 GFN1 `sqm_val_*` tests pass at 1e-8 (xfails removed in WP1's CMake).
- The GFN1 row of `docs/SQM_VALIDATION.md` shows ≤1e-8 across the set.
- Each closed residual has a one-line root-cause note in `AIChangelog.md`.

## Caveats / not in scope

- No `✅ TESTED`/`APPROVED` by AI; machine-tested only.
- 1e-8 is the agreement-with-tblite target, not a statement of physical
  correctness — both implement the same model; this measures *implementation*
  fidelity. Genuinely different but defensible choices (if any are found) must be
  documented, not silently tolerance-loosened.
- Elements/regimes outside the set (open-shell, charged, heavy elements,
  solvation) are not covered here.

## Results (2026-05-30) — 🤖 AI-generated, machine-tested only

**Localized (WP3 tooling + new D3 probes):** the entire GFN1 dispersion residual
was the **two-body D3 C6 coefficient**. Every other D3 term was verified
bit-matching s-dftd3 (tblite's GFN1 D3 engine): the molecular CN, the C8/C6 ratio
(`10.72·r4r2_i·r4r2_j` ≡ s-dftd3 `3·sqrt_z·sqrt_z`), the BJ critical radius
`a1·sqrt(c8/c6)+a2`, and the pair energy `-s6·c6·t6 - s8·c8·t8`; He2 (CN=0) was
exact. New probes: `dump_dftd3_atomic_c6.f90` (s-dftd3 `get_atomic_c6`) and
`diag_d3_c6.cpp` (curcuma interpolated C6).

**Root cause:** `src/core/.../ff_methods/d3_reference_cn.cpp`
(`reference_cn_data_complete`, 721 values) was **scrambled** — carbon read
`[-1,-1,0.96,-1,-1,9.81,0]` instead of s-dftd3's `[0,0.987,1.999,2.999,3.984,-1,-1]`,
corrupting the Gaussian reference weights → carbon C6 +61% (C-C) / +26% (C-H),
H −3.6%, O ≈exact. The **C6 reference table** (262444 values) and
`m_number_of_references` were correct.

**Fix:** regenerated the CN table from s-dftd3 `reference.f90 reference_cn(7,103)`
(column-major = curcuma's `[elem*7+ref]` layout). Also pinned `d3_s9=0.0` in
`createForGFN1` (tblite `gfn1.f90:53`; already effectively 0). Outcome:

| metric | before | after |
|---|---:|---:|
| interpolated C6 vs s-dftd3 (worst pair) | +61% | 0.00000% |
| injected-density \|disp+elec\| triose | 9.28e-3 | 1.44e-5 |
| GFN1 total dE triose | +8.4e-3 | +8.96e-4 |
| GFN1 total dE complex | +3.5e-2 | +3.36e-3 |

Shared with UFF-D3 / GFN-FF-D3 → **regression-checked, none** (GFN-FF 21/22 as
baseline, gfn2/gradient/cpscf/d4 green). `gfn1_align` tolerances tightened.

**Electronic residual — FIXED (2026-05-30).** The size-extensive electronic
difference (triose ~8.8e-4, complex ~3.36e-3) was a **double-counted GFN1
third-order potential**. `XTB::addThirdOrderPotential` (`xtb_thirdorder.cpp`)
added the atom-resolved `v_at(i)=q_i²·Γ_i` into **both** `pot.v_at` *and*
(broadcast) `pot.v_sh`; `expand_potential` forms `v_ao = v_sh + v_at`, so the
Fock saw the third-order at 2×. Over-penalised charge → less polarised SCF fixed
point (triose Coulomb-ES2 0.660 vs tblite 0.676), electronic +8.8e-4. The same
double-broadcast was mirrored in the CPSCF response (`xtb_response.cpp`). Fix:
keep the GFN1 third-order in `v_at` only (its intended home) in both places.

How it was localized: Phase-0 diffs proved **S and H0 bit-perfect** for GFN1
(max|Δ| ~1e-13) and γ structurally exact; the **n0 reference occupations exact**
(`q_at` from curcuma-n0 − tblite-n_sh matched tblite to ~1e-16); the **ref-P
one-shot Fock** (now extended to GFN1 in `test_xtb_scf_snapshot.cpp`) reproduced
tblite orbital energies to ~5e-6. All three SCF modes (plain/diis/broyden)
converged to the *same* wrong fixed point → not the mixer, the Fock/potential.

Result (curcuma − tblite total, after fix): triose 8.96e-4 → **1.45e-5**, complex
3.36e-3 → **5.09e-5**, caffeine 1.8e-3 → 5.3e-6, acetic_acid_dimer 1.4e-3 →
2.1e-6, C6H6 4.8e-4 → 2.2e-6, small molecules ~6e-8..2e-7. triose electronic now
matches tblite to 8 decimals (−131.33260976). Gradient FD-validated
(`test_xtb_gradient` ngfn1) and CPSCF green — the fix also corrects the GFN1
gradient/response, which double-counted too.

**D3 C8 tail — FIXED (2026-05-31).** The dispersion residual was a systematic
**+0.06% error in the D3 `C8/C6` ratio**. `D3ParameterGenerator::getR6`
(`d3param_generator.cpp`) used an empirical `C8/C6 = 10.72·r4r2_i·r4r2_j` with a
non-standard `r4r2` table; the authoritative s-dftd3 form is
`C8/C6 = 3·r4r2_i·r4r2_j` with `r4r2(z)=√(½·⟨r⁴⟩/⟨r²⟩·√Z)` (s-dftd3
`data/r4r2.f90:69`, `damping/rational.f90:173`). The empirical fit was +0.055…
+0.069% high on every pair → a size-extensive dispersion bias. Fix: exact form +
the verbatim s-dftd3 raw `⟨r⁴⟩/⟨r²⟩` table (118 elements; identical to
`D4ParameterGenerator::m_r4_over_r2`). This also corrects the BJ radius
`R0=a1·√(C8/C6)+a2`.

Result: triose dispersion −0.04199628 → **−0.04201074** (= tblite, exact); GFN1
total now bit-matches tblite. **10/12 molecules at 1e-8** (H2, LiH, H2O, CH4,
NH3, C6H6, HCN, acetic_acid_dimer, caffeine, triose). xfail list reduced to
`He2` (~1.5e-8, dispersion-only numerical floor) and `complex` (~2.6e-7, 231
atoms — residual accumulates; cause not yet localized).

**Scope of validation (important).** This is validated **only** against tblite
(genuine s-dftd3) via the GFN1 path — that pins the shared D3 kernel (C6 tables,
CN, `C8/C6`, BJ `R0`, damping sum) as correct at GFN1's params
(s6=1, s8=2.4, a1=0.63, a2=5.0, s9=0). **UFF-D3, `gfnff-d3`, and standalone-D3
remain UNVALIDATED** — they use different damping params / code paths and have no
authoritative s-dftd3 reference. The `getR6` fix improves their `C8/C6` to exact,
but their full energies are not checked against anything authoritative. The older
"D3 <1% on 10/11" table does not establish correctness (its references were never
tied to s-dftd3).
