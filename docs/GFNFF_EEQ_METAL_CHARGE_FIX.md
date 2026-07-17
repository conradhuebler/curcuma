# GFN-FF Phase-1 EEQ Metal Chi-Shift Fix (RhCl(CO)2 / ED39)

Status: 🤖 AI-generated, ⚙️ machine-tested. Human production testing pending.
Date: 2026-07-17

## Summary

Curcuma's native GFN-FF Phase-1 EEQ **topology charges** for transition-metal (TM)
complexes were wrong: the TM was given too little positive charge and its ligands too
much. This propagated into the charge-dependent bond factor (`fqq`) and the Coulomb term.
For ED39 = RhCl(CO)2 the native single-point was **+12.9 mEh above** xtb 6.6.1 `--gfnff`.

Root cause: the Fortran reference applies a **metal chi-shift** (`mchishift = -0.09`,
i.e. `chieeq += 0.09`) to the Phase-1 EEQ electronegativity of every TM (`imetal==2`),
making the metal *less* electronegative so it carries *more* positive topology charge
(gfnff_ini.f90:413-418). Curcuma's actually-used Phase-1 solver
(`EEQSolver::calculateTopologyCharges`) never applied this shift.

After the fix the Phase-1 charges match the Fortran reference **to 6 decimals** for every
TM complex tested, and the ED39 bond + Coulomb terms match xtb exactly.

## Root cause (file:line)

`src/core/energy_calculators/ff_methods/eeq_solver.cpp`, function
`EEQSolver::calculateTopologyCharges` (the Phase-1 topology-charge path called from
`gfnff_method.cpp:9146/9176/9181/9428/11085`).

Before the fix, line ~2735 built the Phase-1 RHS as:

```cpp
chi(i) = -params_i.chi + dxi(i) + cnf_term;   // no metal chi-shift
```

The separate `EEQSolver::solveEEQ` function *did* contain a metal-shift branch (line 2599),
but that function is only reached from line 1158 (a different, non-Phase-1 code path), so the
shift was effectively dead for the topology charges that feed `fqq` and Coulomb.

Fortran reference:
- `external/gfnff/src/gfnff_ini.f90:417` — `topo%chieeq(i) = topo%chieeq(i) - gen%mchishift` inside `if (imetal(i).eq.2)`
- `external/gfnff/src/gfnff_param.f90:756` — `gen%mchishift = -0.09_wp`
- `imetal(i) = param%metal(at(i))` (gfnff_ini.f90:273); the line-274 caveat (`nb<=4 & group>3 → imetal=0`) can only clear main-group (`imetal==1`) flags — every `imetal==2` element has periodic group ≤2 or negative — so testing curcuma's `metal_type[z-1]==2` is faithful.

## The fix

Added the Phase-1-only metal chi-shift in `calculateTopologyCharges`, gated on
`metal_type[z_i-1] == 2` (curcuma's `metal_type` array in `gfnff_par.h:370` is a verified
mirror of Fortran `param%metal(86)`):

```cpp
chi(i) = -params_i.chi + dxi(i) + cnf_term;
if (z_i >= 1 && z_i <= 86 && metal_type[z_i - 1] == 2) {
    constexpr double MCHISHIFT = -0.09;  // gfnff_param.f90:756
    chi(i) -= MCHISHIFT;                 // effectively += 0.09
}
```

Because the gate is `metal_type==2` (zero for H/C/N/O and all non-TM elements), the change is
a mathematical **no-op for every neutral-organic system**. Phase 2 deliberately overwrites
`chieeq` without the shift (gfnff_ini.f90:715), so the shift stays Phase-1-only as in Fortran.

## Verification — ED39 RhCl(CO)2

Phase-1 topology charges (verbosity 3), vs Fortran analyzer
(`gfnff-gfnff_analyze-test .../ED39/mol.xyz - 4`):

| Atom | qa before | qa after | Fortran ref |
|------|-----------|----------|-------------|
| Rh (Z=45) | 0.261693 | **0.355794** | 0.355794 |
| Cl (Z=17) | 0.045102 | **0.003080** | 0.003080 |

The Phase-1 RHS `chieeq(Rh)` moved -1.070317 → -0.980317 (exactly +0.09), matching Fortran.
The A-matrix was already bit-identical; only the Rh RHS entry differed.

ED39 energy components (native after fix vs xtb 6.6.1 `--gfnff`):

| Term    | before      | after (native) | xtb ref     |
|---------|-------------|----------------|-------------|
| Bond    | -0.98363    | **-0.9965630** | -0.99656    |
| Coulomb | -0.13213    | **-0.1357345** | -0.13573    |
| Angle   | +0.00063    | +0.0006301     | +0.00410    |
| **Total** | **-1.04046782** | **-1.05685502** | **-1.05338689** |

Bond (via `fqq`) and Coulomb now match the reference to 5+ digits. The residual
**-3.47 mEh** is entirely the **Angle** term (native +0.00063 vs xtb +0.00410), which is
**unchanged by the charge fix** — it is the known unimplemented GFN-FF metal-coordination
angle correction (`feta` = 1.0 for all, see ff_methods/CLAUDE.md "Metal coordination").
|residual| improved 12.9 → 3.5 mEh.

## Spot-check — other MOR41 metal complexes

Phase-1 charges now match the Fortran reference **exactly (6 decimals)** for all:
Fe (ED02) 0.331561, Ni (ED03) 0.300502, Cr (PR01) 0.345743, Rh (ED39) 0.355794.

Totals vs Fortran-analyzer reference (`- 0`):

| System | before | after | ref | \|before-ref\| | \|after-ref\| |
|--------|--------|-------|-----|-----------|-----------|
| ED39 RhCl(CO)2 | -1.04046782 | -1.05685502 | -1.05338689 | 12.9 mEh | **3.5 mEh** |
| ED03 Ni(CO)3   | -1.41115994 | -1.41148205 | -1.41148176 | 32.2 mEh | **0.0003 mEh** |
| PR01 Cr(CO)6   | -2.70835095 | -2.70825132 | -2.71178260 | 3.43 mEh | 3.53 mEh |
| ED02 Fe(CO)4   | -1.88289483 | -1.88330502 | -1.87614184 | 6.75 mEh | 7.16 mEh |

Interpretation (honest):
- **ED39, ED03 — large improvement.** ED03 becomes essentially exact (0.3 µEh).
- **PR01, ED02 — total essentially unchanged (±0.1–0.4 mEh).** The charges are now correct
  (they match the reference exactly), but ED02 over-binds the reference by ~7 mEh **both
  before and after**, so the ~7 mEh gap is a **non-charge** issue (metal bond parameters /
  angle `feta` / repulsion — separate documented TODOs). Before the fix, the *wrong* TM
  charges were partially masking these other errors; correcting the charges removes that
  accidental cancellation, which is why the total moves by a fraction of a mEh. This is a
  correctness improvement in the charges even where the total does not shrink.

The fix is correct-by-construction: it makes the Phase-1 topology charges equal the Fortran
reference. Remaining metal-complex residuals are in the angle/bond/repulsion metal
corrections, not in the charges.

## Neutral-organic regression check

`cd release && ctest -R gfnff` → **75/76 pass**. The single failure `gfnff_gpu_vs_cpu_h2`
is pre-existing and unrelated: `std::runtime_error: File not found:
../../external/gfnff/test/h2.xyz` (a missing-fixture / working-dir issue, GPU test, H2 has no
metals). Passing: 21/21 gfnff CPU, 37/37 validation (caffeine/triose/monosaccharide/…),
28/28 gfnff solvation. The `metal_type==2` gate guarantees byte-identical results for all
neutral-organic systems.

## Not tested / caveats
- Transition-metal GFN-FF remains AI-implemented and not human production tested.
- Only Fe/Ni/Cr/Rh carbonyl-type complexes checked; the fix is general (any `imetal==2`
  element) but broader TM validation (charged complexes, open shells, 5d metals) is pending.
- The residual metal angle/bond/repulsion corrections (`feta`, TM bond parameters) are
  separate unimplemented items and still leave a few mEh on some complexes.
