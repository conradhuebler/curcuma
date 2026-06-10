# tblite diagnostic patches

`external/tblite/` is a FetchContent `SOURCE_DIR` (pinned to upstream commit
`9ca469a77cb9433dc6b0469ba5f6effc48a1bda8`). It is **not** version-controlled by
curcuma, so the diagnostic C-API extensions curcuma relies on live only in that
fetched tree and are dropped by a clean re-fetch. This directory keeps them as a
re-applicable patch.

## `curcuma-tblite.patch`

Adds **diagnostic-only** C-API exports to tblite (no change to the physics or the
public production API). Two groups:

1. **AP6b multipole diagnostics** — converged multipole AO integrals
   (`tblite_get_result_{dipole,quadrupole}_integral`) and atomic moments +
   potentials (`..._atomic_dipole/quadrupole`, `..._potential_vat/vdp/vqp`).

2. **GFN2 component audit** (2026-05) — per-container atom-resolved energies
   (`tblite_get_result_energy_component_{halogen,repulsion,dispersion,
   interactions,electronic}`); the `electronic` lump is Tr(P·H0) + Coulomb-shell
   + 3rd-order + multipole, the ATM 3-body lands in `dispersion` (pre-SCF), the
   SCC 2-body in `electronic`.

Touched files: `src/tblite/results.f90`, `src/tblite/xtb/singlepoint.f90`,
`src/tblite/api/result.f90`, `include/tblite/result.h`.

### Who needs it
- `test_cases/sqm_reference/dump_tblite_multipole` (writes `*_tblite` multipole
  fields + the `e_components_tblite` block + `shell_charges`).
- transitively, `diag_curcuma_energy_components` (`ctest -L gfn2_align`), which
  reads those dump fields.

The standalone dftd4 tools (`dump_dftd4_c6`, `dump_dftd4_atomic_c6`,
`dump_dftd4_atm`) call dftd4's **public** API directly and need **no** patch —
`dftd4-src` stays pristine.

## Apply (after a fresh configure that re-fetched tblite)

```bash
git -C external/tblite apply patches/tblite/curcuma-tblite.patch
# or, without git:
patch -p1 -d external/tblite < patches/tblite/curcuma-tblite.patch
# then rebuild the USE_TBLITE targets (release_tblite/)
```

Check first with `git -C external/tblite apply --check patches/tblite/curcuma-tblite.patch`.

## Regenerate (after extending the patch)

```bash
git -C external/tblite diff > patches/tblite/curcuma-tblite.patch
```

Base commit: `9ca469a` (see top-level `CMakeLists.txt` FetchContent_Declare(tblite)).
These are AI-generated diagnostic patches; they only add read-only getters used by
the curcuma↔tblite alignment test tooling.
