# SQM Validation Suite — Native GFN1/GFN2 vs tblite

> ⚠️ **AI-generated, machine-tested only — not human production tested.**
> Baseline measured 2026-05-30 on `release/` (USE_TBLITE=OFF). Re-measure after
> any change to the native xTB code or the tblite pin.
>
> **Roadmap & open work:** [SQM_VALIDATION_ROADMAP.md](SQM_VALIDATION_ROADMAP.md)
> (WP1 honest tolerances, WP2 GFN1→1e-8, WP3 component validation, WP4 perf, WP5
> consolidation).

## What this is

An end-to-end validation suite for the **native** GFN1/GFN2-xTB implementation
(`curcuma::xtb::XTB` in `src/core/energy_calculators/qm_methods/xtb_native.cpp`),
modeled on the GFN-FF suite. It runs the real `curcuma` binary via `gfn1`/`gfn2`
(the canonical native methods) and compares the single-point energy against
**committed** tblite reference JSONs,
so a native-only build (no tblite) can validate.

## The target and how the tests reflect it

The target is **honest 1e-8 Eh** agreement with tblite. Every `sqm_val_*` ctest
gates total energy at exactly `1e-8` (`validate_sqm.py --tol-energy 1e-8`). The CLI
prints energies to 8 decimals, so 1e-8 Eh is the validator's resolution; "≤1e-8"
means native matches tblite to every printed digit.

Molecules that do **not** yet reach 1e-8 are marked `WILL_FAIL` (ctest
expected-failure) with their measured residual recorded in
`test_cases/sqm_reference/CMakeLists.txt`. The gate is **never loosened** to hide a
gap. When a method reaches 1e-8 for an xfail molecule, ctest reports it as a
**failure** ("passed but expected to fail"), which forces removal of the xfail —
progress is self-documenting and a closed gap cannot silently re-open.

## Components

| File | Role |
|---|---|
| `test_cases/sqm_reference/molecules/*.xyz` | 12 test molecules (H2, He2, LiH, H2O, CH4, NH3, C6H6, HCN, acetic_acid_dimer, caffeine, triose, complex) |
| `test_cases/sqm_reference/reference_data/<mol>_<method>.ref.json` | committed tblite references (24 files) |
| `test_cases/sqm_reference/dump_tblite_reference.cpp` | compact reference generator (USE_TBLITE build only) |
| `test_cases/sqm_reference/validate_sqm.py` | runs `curcuma`, gates total energy vs reference at 1e-8 |
| `test_cases/sqm_reference/CMakeLists.txt` | registers `sqm_val_<mol>_<method>` ctests (labels `gfn1_validation`/`gfn2_validation`), 1e-8 gate + xfail list |

The `*.ref.json` schema: `total_energy`, `charges`, `dipole`, `homo_energy`,
`lumo_energy`, `gradient` (norm + per-atom values), and — once WP3 lands —
`energy_components` (repulsion/dispersion/electronic/interactions/halogen; requires
the tblite diagnostic patch). **Note:** the currently committed references do *not*
carry the `energy_components` block (see WP3).

## Running the validation (native, no tblite)

```bash
cd release && make -j4
ctest -L gfn2_validation --output-on-failure   # 11 pass + 1 xfail (complex) -> green
ctest -L gfn1_validation --output-on-failure   # 12 xfail -> green
```

## Regenerating the references

References are committed and only need regeneration when the tblite pin changes or
the molecule set grows. Requires a tblite build (`-DUSE_TBLITE=ON`, e.g.
`release_tblite/`); the `energy_components` block additionally needs the diagnostic
patch (`patches/tblite/curcuma-tblite.patch`).

```bash
cd release_tblite && make dump_tblite_reference -j4
cd ..
DUMP=release_tblite/test_cases/sqm_reference/dump_tblite_reference
MOL=test_cases/sqm_reference/molecules
OUT=test_cases/sqm_reference/reference_data
for m in H2 He2 LiH H2O CH4 NH3 C6H6 HCN acetic_acid_dimer caffeine triose complex; do
  for meth in gfn1 gfn2; do
    $DUMP $meth $MOL/$m.xyz $OUT/${m}_${meth}.ref.json
  done
done
```

## Measured residuals (2026-05-30 baseline)

`dE = E_native − E_tblite`, in Hartree. "≤1e-8" = matches to all 8 printed digits.

### GFN2 — meets 1e-8 on 11/12

| molecule | \|dE\| (Eh) | 1e-8 | | molecule | \|dE\| (Eh) | 1e-8 |
|---|---:|:--:|---|---|---:|:--:|
| H2 | 3.3e-9 | ✓ | | HCN | 4.7e-9 | ✓ |
| He2 | 9.7e-10 | ✓ | | acetic_acid_dimer | 1.8e-9 | ✓ |
| LiH | 4.0e-9 | ✓ | | caffeine | 3.7e-10 | ✓ |
| H2O | 1.9e-9 | ✓ | | triose | 9.7e-10 | ✓ |
| CH4 | 4.3e-9 | ✓ | | **complex** (231 at) | **6.95e-5** | ✗ |
| NH3 | 1.4e-9 | ✓ | | C6H6 | 4.2e-10 | ✓ |

Native GFN2 reproduces tblite to numerical noise on every small/medium molecule,
including the 66-atom `triose`. The single open case is the 231-atom `complex`
(6.95e-5 Eh ≈ 0.07 mEh) — re-examined in WP2.

### GFN1 — does not yet meet 1e-8 on any

| molecule | \|dE\| (Eh) | 1e-8 | | molecule | \|dE\| (Eh) | 1e-8 |
|---|---:|:--:|---|---|---:|:--:|
| H2 | 1.3e-6 | ✗ | | HCN | 4.4e-5 | ✗ |
| He2 | 1.5e-8 | ✗ | | acetic_acid_dimer | 1.4e-3 | ✗ |
| LiH | 4.3e-5 | ✗ | | caffeine | 1.8e-3 | ✗ |
| H2O | 2.1e-6 | ✗ | | triose | 8.4e-3 | ✗ |
| CH4 | 4.7e-5 | ✗ | | **complex** (231 at) | 3.5e-2 | ✗ |
| NH3 | 1.6e-5 | ✗ | | C6H6 | 4.8e-4 | ✗ |

He2 is borderline (1.5e-8); the rest range from ~1µEh to 35 mEh, growing with
size. **Sign matters:** the large systems (triose, caffeine, complex) are *more*
negative than tblite — native GFN1 **over-binds** — while acetic_acid_dimer is
*less* negative. A missing attractive dispersion term would push native uniformly
*less* negative, so the residual is **not** "missing D3" (native GFN1 D3 exists,
`xtb_native.cpp::calcDispersionEnergy` GFN1 branch). **Root cause is open**; it is
localized per energy container in WP3 and closed in WP2.

## What is tested / not tested

- **Tested**: total single-point energy of native GFN1/GFN2 vs tblite on 12
  molecules (1–231 atoms), self-converged native SCF (the real user dispatch via
  MethodFactory).
- **Reported, not yet the automatic gate**: Mulliken charges, gradient (native
  gradients are FD-validated separately by `test_xtb_gradient`); per-container
  energy decomposition (WP3; precise injected-density audit is `ctest -L gfn2_align`).
- **Not tested**: open-shell, charged species, elements beyond the set
  (H/C/N/O/Li/He), solvation, excited states, numerical stability across all
  geometries.

## Caveats

- The 1e-8 gate is the **agreement-with-tblite** target (implementation fidelity:
  both implement the same model), not a statement of physical correctness.
- xfail entries carry the *measured* residual; they are tracked defects, not
  accepted tolerances. Tightening (removing xfails) as WP2/WP3 close residuals is
  the expected trajectory.
