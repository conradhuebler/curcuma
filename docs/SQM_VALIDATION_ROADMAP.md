# SQM Validation & Consolidation Roadmap

> ⚠️ **AI-generated, machine-tested only — not human production tested.**
> Author: Claude. Baseline measured 2026-05-30 on `release/` (USE_TBLITE=OFF)
> against committed tblite references in `test_cases/sqm_reference/reference_data/`.

## Goal

**Honest 1e-8 Eh agreement** between the native GFN1/GFN2-xTB implementation
(`curcuma::xtb::XTB`, `src/core/energy_calculators/qm_methods/xtb_native.cpp`)
and tblite, across the validation molecule set — and a test suite that **reflects
that target truthfully**: the gate is 1e-8 everywhere, and molecules that do not
yet reach it are marked as tracked expected-failures (not hidden by loosened
tolerances).

This supersedes the loosened, partly-fabricated tolerances written into
`docs/SQM_VALIDATION.md`, `CLAUDE.md`, `README.md`, `AIChangelog.md` and
`test_cases/sqm_reference/CMakeLists.txt` during a tooling outage — correcting
those is **WP1**.

## Measured baseline (2026-05-30, dE = E_native − E_tblite)

CLI prints energies to 8 decimals, so the validator's effective resolution is
1e-8 Eh. "≤1e-8" means native matches tblite to all printed digits.

### GFN2 — 11/12 already meet 1e-8

| molecule | E_native (Eh) | E_tblite (Eh) | \|dE\| (Eh) | 1e-8 |
|---|---:|---:|---:|:--:|
| H2 | -0.98202304 | -0.98202304 | ≤1e-8 | ✓ |
| He2 | -3.48627274 | -3.48627274 | ≤1e-8 | ✓ |
| LiH | -0.77927818 | -0.77927818 | ≤1e-8 | ✓ |
| H2O | -5.07036982 | -5.07036982 | 1.9e-9 | ✓ |
| CH4 | -4.17507458 | -4.17507458 | ≤1e-8 | ✓ |
| NH3 | -4.42618246 | -4.42618246 | ≤1e-8 | ✓ |
| C6H6 | -15.87853118 | -15.87853118 | ≤1e-8 | ✓ |
| HCN | -5.50406617 | -5.50406617 | ≤1e-8 | ✓ |
| acetic_acid_dimer | -28.89721517 | -28.89721517 | ≤1e-8 | ✓ |
| caffeine | -42.14723025 | -42.14723025 | ≤1e-8 | ✓ |
| triose | -119.90699419 | -119.90699419 | ≤1e-8 | ✓ |
| **complex** (231 at) | **-329.52707823** | **-329.52714777** | **6.95e-5** | ✗ |

### GFN1 — 0/12 meet 1e-8

| molecule | E_native (Eh) | E_tblite (Eh) | \|dE\| (Eh) | 1e-8 |
|---|---:|---:|---:|:--:|
| H2 | -1.03616163 | -1.03616294 | 1.3e-6 | ✗ |
| He2 | -3.25174863 | -3.25174861 | 1.5e-8 | ✗ |
| LiH | -0.88140542 | -0.88144840 | 4.3e-5 | ✗ |
| H2O | -5.76863902 | -5.76864110 | 2.1e-6 | ✗ |
| CH4 | -4.27428599 | -4.27423856 | 4.7e-5 | ✗ |
| NH3 | -4.83026580 | -4.83024948 | 1.6e-5 | ✗ |
| C6H6 | -15.89466767 | -15.89419025 | 4.8e-4 | ✗ |
| HCN | -5.78258461 | -5.78262866 | 4.4e-5 | ✗ |
| acetic_acid_dimer | -31.57700486 | -31.57843892 | 1.4e-3 | ✗ |
| caffeine | -44.51168034 | -44.50985543 | 1.8e-3 | ✗ |
| triose | -130.56700873 | -130.55861510 | 8.4e-3 | ✗ |
| complex | -343.21472744 | -343.17980350 | 3.5e-2 | ✗ |

**Sign note (important, corrects an earlier wrong claim):** on the large systems
native GFN1 is *more* negative than tblite (triose −8.4, caffeine −1.8, complex
−35 mEh → native **over-binds**), while acetic_acid_dimer is *less* negative
(+1.4 mEh). A missing attractive dispersion term would make native uniformly
*less* negative, so the residual is **not** simply "missing D3" (native GFN1 D3
exists since commit `063a1d2`). Root cause is **open** and must be localized
per-container — this is exactly why WP3 (component validation) gates WP2.

## Test-design principle (how 1e-8 is reflected)

`test_cases/sqm_reference/CMakeLists.txt` registers one ctest per (molecule ×
method), gating total energy via `validate_sqm.py` at **`--tol-energy 1e-8`**
for **all** of them. Molecules that do not yet meet 1e-8 are marked
`WILL_FAIL TRUE` (ctest expected-failure) with the measured residual in an
adjacent comment. Consequences:

- The honest gate is 1e-8 everywhere — no tolerance is loosened to mask a gap.
- The suite stays green today (xfails are expected), so it is usable in CI.
- When a method reaches 1e-8 for a molecule, ctest reports the xfail as a
  **failure** ("passed but was expected to fail"), forcing removal of the xfail —
  the gap cannot be silently re-opened, and progress is self-documenting.

xfail inventory at baseline: GFN2 `complex`; GFN1 all 12.

## Work packages

| WP | Title | Status | File |
|----|-------|--------|------|
| **WP1** | Honest 1e-8 tolerances + correct fabricated numbers | ADD | [SQM_WP1_honest_tolerances.md](SQM_WP1_honest_tolerances.md) |
| **WP2** | GFN1 → 1e-8 (localize & close residuals) | ADD | [SQM_WP2_gfn1_accuracy.md](SQM_WP2_gfn1_accuracy.md) |
| **WP3** | Component-level validation (per-container dE) | ADD | [SQM_WP3_component_validation.md](SQM_WP3_component_validation.md) |
| **WP4** | Performance benchmark native vs Fortran | DEPOSITED (below) | — |
| **WP5** | SQM consolidation (shared GFN1/GFN2 paths) | DEPOSITED (below) | — |

Dependency: **WP1** first (turns the suite honest). **WP3** delivers the
per-container diff tool that **WP2** needs to localize the GFN1 residual. WP2 also
re-checks the single GFN2 `complex` 7e-5 residual. WP4/WP5 follow once the methods
are honestly green at 1e-8.

---

## WP4 (deposited) — Performance benchmark: native vs Fortran

**Goal:** quantify native GFN1/GFN2 single-point + gradient wall-time against the
Fortran reference (tblite) on the same molecules, to track the "as fast as / faster
than Fortran" objective and find hotspots.

**Scope:**
- New `test_cases/sqm_reference/bench_sqm.cpp` (or a ctest-driven timing harness):
  best-of-N wall time for `ngfn1`/`ngfn2` SP and SP+grad over the molecule set.
- When `USE_TBLITE`, time tblite `get_singlepoint` on the same geometries in the
  same binary for an apples-to-apples comparison; emit a `native | tblite | ratio`
  table.
- Record a baseline table in this roadmap. Note hotspots (Eigen solve, D4
  reference build, STO integrals) for later optimization WPs.
- Reference point already measured: GFN2 `complex` (231 at) SP+grad ~4.4 s after
  the D4-regeneration fix (was ~7.1 s) — see [[gfn2-d4-scf-regeneration-fix]].

**Not in scope:** GPU; threading rework (separate initiative).

## WP5 (deposited) — SQM consolidation (shared GFN1/GFN2 paths)

**Goal:** reduce duplication and per-call overhead now that the legacy standalone
`gfn1.cpp`/`gfn2.cpp` are gone and everything routes through `xtb_native.cpp`.

**Candidate items (to be scoped after WP1–WP3):**
- Unify the GFN1/GFN2 dispersion entry (`calcDispersionEnergy`) — D3 vs D4 select
  behind one prepared-once guard, mirroring the WP-C `m_d4_prepared` pattern, so
  GFN1 D3 is not regenerated per SCF iteration either.
- Shared CN / EEQ / reference-data paths between the GFN1 and GFN2 code so a fix
  lands once (the per-container diffs from WP3 will show where they must differ).
- Audit remaining per-SCF recompute in `xtb_native.cpp` (gamma matrix, multipole
  setup) for geometry-fixed quantities that can move out of the SCF loop.
- Fold the analytical-D3-gradient plan (`plans/analytical_d3_gradient_plan.md`)
  into the consolidated dispersion path.

**Not in scope:** new methods; changing the public CLI surface.
