# WP Phase 3b-4 — GFN2-D4 Mulliken charge-response: OUTCOME

**Status:** Resolved for the D4 target; raw ∂q/∂x deep-dive → [PHASE3B5](PHASE3B5_MULTIPOLE_RESPONSE_WP.md)
**Updated:** 2026-05-25
**Branch:** `feature/sqm-implementation`
**Build:** always in `release/` (`cd release && make -j4`)

---

## What this WP originally assumed (and why it was wrong)

The original plan held that the remaining error in the `d4_charge_source="mulliken"`
gradient came from a **missing multipole-interaction linearization** (+ an explicit
overlap term), to be added to `computeMullikenChargeResponse`. Investigation
disproved this:

1. The multipole-interaction linearization was implemented (two-pass bilinear
   δ-substitution of the energy-gradient SD/DD/SQ terms, `mpPairGrad` in
   `xtb_response.cpp`). It contributes only ~1e-6 — two orders below the
   discrepancy. **Not the fix.**
2. The explicit-overlap term (`EXPL_FAC·lam_pair·P`) **worsens** agreement at any
   nonzero value, both signs. **Dead end — removed.**

## The actual bug: inverted response sign

A decoupled diagnostic (Test D, `test_xtb_cpscf.cpp`) calls
`computeMullikenChargeResponse` with `dEdq = e_A`, yielding the raw `∂q_A/∂x`,
FD-validated against Mulliken charges with a tightened SCF (1e-9). It showed the
analytical response was **≈ −(FD truth)** for H₂O, HCN **and** CH₄.

Root cause: `RHS_SIGN = −1` was "validated" using a **sign-blind max-abs ratio**
(`|Δanal|/|Δnum|`), so a clean sign flip read as "ratio ≈ 1.0, correct." It
inverted the entire response.

**Fix: `RHS_SIGN = +1`** (`xtb_response.cpp`). Effect:

| metric | RHS=−1 (old) | RHS=+1 (fixed) |
|--------|--------------|----------------|
| raw ∂q/∂x err, H₂O / HCN / CH₄ | 1.7e-1 / 1.6e-1 / 8.0e-2 | 1.1e-2 / 2.1e-2 / 6.9e-3 |
| D4 isolation err (target <5e-5), H₂O | 4.19e-5 | **2.46e-5** |
| D4 isolation err, HCN | **5.66e-5 (FAIL)** | **7.75e-6** |

The D4-mulliken gradient now meets `< 5e-5` for H₂O **and** HCN. All 6 regression
tests pass (`ctest -R "xtb_gradient|d4_dedq|energy_methods|xtb_cpscf"`).

## What shipped

- `RHS_SIGN = +1`, locked `constexpr` with a comment explaining the FD validation.
- Removed the dead explicit-overlap term and all temporary env calibration knobs.
- Multipole charge-response terms (3b-4a integral-Pulay + 7b interaction) **kept
  but gated off** (`MP_RESPONSE_ENABLED = false`): correct for low-multipole
  systems (CH₄ → ~1% with them on) but mis-contracted for polar molecules. Off by
  default; immaterial to the D4 target.
- New gated tests in `test_xtb_cpscf.cpp`:
  - **Test D** — raw `∂q_A/∂x` vs FD (tight SCF), tol 3e-2, guards `RHS_SIGN`
    (a sign regression gives ~1.5e-1).
  - **Test C** — D4 isolation now **gated** at `< 5e-5` for H₂O/HCN.
- Added `XTB::setScfThreshold` (test support for low-noise FD).

## What remains (separate WPs)

- **Raw ∂q/∂x to < 1% for polar molecules** — the multipole property-gradient
  contraction needs re-derivation. See
  [PHASE3B5_MULTIPOLE_RESPONSE_WP.md](PHASE3B5_MULTIPOLE_RESPONSE_WP.md).
- **Pre-existing ~4e-5 Eh/Å GFN2 energy-gradient baseline error** (separate from
  the response) — still open; see the note in the AP6b WP.

---

## References

- Handy & Schaefer, J. Chem. Phys. 81 (1984) 5031 — Z-vector.
- tblite `coulomb/multipole.f90::get_multipole_gradient_0d`.
- [docs/D4_Q_RESPONSE.md](D4_Q_RESPONSE.md), [AIChangelog.md](../AIChangelog.md).
