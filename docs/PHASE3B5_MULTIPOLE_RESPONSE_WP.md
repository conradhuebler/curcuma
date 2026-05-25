# WP Phase 3b-5 — Multipole charge-response: raw ∂q/∂x to <1% for polar molecules

**Status:** Open (deep-dive). Prereq done: [PHASE3B4](PHASE3B4_MULLIKEN_RESPONSE_WP.md) (RHS_SIGN=+1).
**Created:** 2026-05-25
**Branch:** `feature/sqm-implementation`
**Build:** always in `release/` (`cd release && make -j4 test_xtb_cpscf`)

---

## Goal

Make the **raw** Mulliken charge response `∂q_A/∂x` agree with FD to **< 1%** for
strongly multipolar molecules (H₂O, HCN), then re-enable the multipole response
terms by default. This is a correctness goal beyond the D4-gradient target (which
is already met without these terms, see PHASE3B4).

## Where things stand (RHS_SIGN=+1, reliable FD: SCF 1e-9, h=2e-4 Å)

Raw `∂q_0/∂x` err (Test D), with the multipole terms toggled (`MP_RESPONSE_ENABLED`
and per-family `MP_RESP_MASK` in `xtb_response.cpp`):

| config | H₂O | HCN | CH₄ |
|--------|-----|-----|-----|
| multipole OFF (shipped default) | 1.09e-2 | 2.12e-2 | 6.85e-3 |
| 3b-4a integral-Pulay only | 2.08e-2 | 2.46e-2 | 8.27e-3 |
| 7b interaction only (SD+DD+SQ) | 1.95e-2 | 2.45e-2 | 1.89e-3 |
| both ON | 1.97e-2 | 2.92e-2 | **4.72e-4** |

`|num|`: H₂O 8.95e-2, HCN 8.60e-2, CH₄ 3.96e-2.

**Key facts:**
- With both multipole terms ON, **CH₄ reaches ~1.2%** — the machinery is
  structurally sound for low-multipole systems.
- For **polar molecules the multipole terms make it worse**. 7b family
  decomposition (3b-4a off): SD *helps* H₂O (1.09e-2→8.8e-3) but DD/SQ *hurt* it,
  while **CH₄ needs SQ** (mask 3→7: 8.8e-3→1.9e-3). So the same SQ/DD term is
  correct for CH₄ and wrong for H₂O — the error scales with multipole magnitude.
- Multipole OFF, the residual is **directional**: H₂O's symmetry-axis (y)
  component is at ~90% of FD but the perpendicular (x, antisymmetric) at ~68%.

## Hypothesis (most likely)

The implemented multipole response is **not the correct property-gradient
skeleton-Fock contraction**. The CPSCF property gradient is
`dL/dx = Tr(D_z·Fˣ_skel) − Tr(W_z·Sˣ) + explicit`, where `Fˣ_skel` includes
`∂V_mp/∂x` at **frozen** SCF moments/potential (the density-induced moment change
is already in `D_z` via the Z-vector kernel, which includes multipoles — Test B is
machine-precise).

- **3b-4a** correctly contracts `∂(multipole integrals)/∂x` with `D_z` and the
  converged `v_dp/v_qp`. This is one half of `∂V_mp/∂x`.
- **7b** currently linearizes the *interaction energy* via δ-moments
  `G(δM_i,M_j)+G(M_i,δM_j)`. This is the **moment-variation**, which is likely
  **already represented in `D_z`** (double counting) rather than the frozen-moment
  `∂A_mp(x)/∂x` skeleton term. For CH₄ (tiny moments) the mistake is negligible
  and the net is right; for polar molecules it dominates and is wrong.

**Conjecture:** replace 7b with the frozen-moment distance/damping derivative of
the multipole interaction (Section 5 `∂A_mp/∂x` at fixed SCF moments) contracted
appropriately with `D_z`-derived potentials, and verify there is no double count
against the Z-vector kernel. Derive from the tblite property-gradient form, not by
linearizing the energy gradient.

## Diagnostic harness (already in place)

- **Test D** (`test_xtb_cpscf.cpp`): set `dEdq = e_A`, get `∂q_A/∂x`, FD vs Mulliken
  charges, tight SCF (`XTB::setScfThreshold(1e-9)`), `h=2e-4`. Set `CPSCF_DUMP=1`
  to print per-component `anal | num | residual`.
- **Levers** in `computeMullikenChargeResponse` (`xtb_response.cpp`): flip
  `MP_RESPONSE_ENABLED=true`; `MP_RESP_MASK` bits 1/2/4 = SD/DD/SQ. For finer
  decomposition, re-introduce temporary env gates (the prior session used
  `CPSCF_*` getenv knobs; removed from the shipped code).
- **FD reliability:** SCF must be ≥1e-9 tight or the FD explodes (Mulliken-charge
  SCF noise / h). Confirmed stable for h ∈ [5e-5, 5e-4] at SCF 1e-9.

## Pitfalls (confirmed)

- The max-abs ratio `|Δanal|/|Δnum|` is **sign-blind** — never validate the
  response with it (this masked the original sign bug). Use signed per-component
  residuals (Test D).
- The explicit-overlap term `−EXPL·lam_pair·P` is a dead end (worsens both signs).
- The (mull−eeq) isolation (Test C) is a D4-level check; it does **not** isolate
  raw `∂q/∂x` (dEdq weighting + D4 energy-difference confound). Use Test D.

## Verification / acceptance

```bash
cd release && make -j4 test_xtb_cpscf && ./test_cases/test_xtb_cpscf
```
- Test D err < 1e-2 (stretch <1% of |num|) for H₂O/HCN/CH₄ with multipole ON.
- Re-enable `MP_RESPONSE_ENABLED`, tighten Test D tol, keep Test C <5e-5.
- Regression `ctest -R "xtb_gradient|d4_dedq|energy_methods|xtb_cpscf"` stays 6/6.

## References

- Handy & Schaefer, J. Chem. Phys. 81 (1984) 5031 — Z-vector property gradient.
- tblite `coulomb/multipole.f90::get_multipole_gradient_0d`.
- [PHASE3B4_MULLIKEN_RESPONSE_WP.md](PHASE3B4_MULLIKEN_RESPONSE_WP.md).
