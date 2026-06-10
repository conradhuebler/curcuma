# WP3 — Component-level validation (per-container dE vs tblite)

> Part of [SQM_VALIDATION_ROADMAP.md](SQM_VALIDATION_ROADMAP.md).
> Status: **ADD** (specified, partially-existing infrastructure). Feeds **WP2**.
> 🤖 AI-generated plan.

## Why

The total-energy gate (WP1) says *whether* native matches tblite; it does not say
*which term* is off when it does not. WP2 needs a per-container breakdown
(repulsion / dispersion / electronic / interactions / halogen) to localize the
GFN1 residual (and the GFN2 `complex` 7e-5). The committed `*.ref.json` were
intended to carry this but currently **do not** (see gap below).

## Current state

- **Tool exists:** `diag_curcuma_energy_components` (label `gfn2_align`) injects
  tblite's converged density into native `XTB` and prints per-container
  |E_native − E_tblite|. It reads the `e_components_tblite` block from the full
  forensic dumps (`dump_tblite_multipole`, `release_tblite/dumps/*_gfn2.json`).
  Today it runs only on H2/H2O/NH3/CH4/triose, **GFN2 only**.
- **Patch exists:** `patches/tblite/curcuma-tblite.patch` exposes the five
  `tblite_get_result_energy_component_*` getters (repulsion / dispersion /
  electronic / interactions / halogen).
- **GAP:** the compact `dump_tblite_reference` (WP-A) calls those getters via
  `pullComponent()`, but the committed `*.ref.json` came out **without** an
  `energy_components` block (the getters returned an error / empty at the time —
  the `ok` guard dropped the block). So WP3's data is currently missing from the
  committed references.

## Tasks

### 1. Fix `energy_components` capture in `dump_tblite_reference.cpp`

Diagnose why `pullComponent()` returns false for the committed references:
- Confirm the `release_tblite` build actually has the diagnostic patch applied
  (the getters are no-ops / error without it). `dump_tblite_multipole` writes
  `e_components_tblite`, so on a patched build the getters work — compare the two
  dumpers' call sequence (error handle reuse, call ordering, container presence
  for H/C/N/O where `halogen` is legitimately absent).
- Make the per-getter failure non-fatal *individually* (record the components
  that ARE present rather than dropping the whole block if one — e.g. halogen —
  is absent).
- Regenerate the committed `*.ref.json` so each carries `energy_components`
  (molecular totals: repulsion, dispersion, electronic, interactions, halogen).

### 2. Component gate in `validate_sqm.py`

Add an optional component comparison when `energy_components` is present in the
reference:
- Compare per-container totals native-vs-reference.
- **Gating rule:** use **combined `dispersion + electronic`** as the gated number,
  not per-container max. tblite distributes D4/D3 across the dispersion and
  electronic containers differently than curcuma lumps it (curcuma puts native D4
  in `m_E_dispersion`); the per-container split is not 1:1 but the **sum** is
  comparable. This matches the rationale already documented in
  `diag_curcuma_energy_components`'s header and memory
  [[gfn2-component-audit-phase-abc]].
- Report repulsion separately (should be bit-identical — a good sanity bit).
- Component gate is **reported**, and gated at a WP-appropriate tolerance
  (looser than 1e-8 on a self-converged density; tighten as WP2 closes residuals).

### 3. Extend the injected-density audit to GFN1 + full set

- Generate GFN1 forensic dumps (`dump_tblite_multipole gfn1 ...`) for the full
  set so `diag_curcuma_energy_components` can run for GFN1.
- Register `ecomp_<mol>_gfn1` ctests (label `gfn1_align`) alongside the GFN2 ones,
  extended to the full molecule set (currently 5 molecules, GFN2 only).

## Done when

- Every committed `*.ref.json` contains an `energy_components` block.
- `validate_sqm.py` reports (and gates, at a stated tolerance) the combined
  disp+electronic and repulsion components when present.
- The injected-density audit (`diag_curcuma_energy_components`) runs for **both**
  methods on the full set; WP2 can read off which container carries the GFN1
  triose/complex residual.

## Caveats

- Per-container numbers on a **self-converged** native density differ more from
  tblite than on the **injected** density — keep the two paths distinct: the
  injected-density audit (`*_align`) is the precise localizer; the `.ref.json`
  component gate is the coarse end-to-end check.
- `halogen` is legitimately absent for H/C/N/O molecules — treat "absent" as
  "0 / not applicable", not a failure.
