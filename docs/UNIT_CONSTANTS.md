# Unit constants: one system (CODATA 2018), old values behind a build switch

⚠️ AI-generated, machine-tested only — not human production tested.

## Summary

Until Sep 25, 2026 curcuma converted between Bohr and Ångström (and between Hartree and
kcal/mol / eV) with **several different constants**, each hard-coded where it was used. The
differences are small (up to 4e-7 relative), but they are systematic, and they surfaced as a
measurable CPU/GPU disagreement: native GFN-FF on the CPU built its coordination-number radii
with one Bohr radius and the GPU with another, which put **+1…2.5e-10 Eh per bond** between the
two (+5.0e-7 Eh on `polymer_2x`, 7320 atoms), always with the same sign.

Now every conversion goes through `src/core/units.h`:

| build | every site uses |
|---|---|
| default | **CODATA 2018** (`CurcumaUnit::Length::BOHR_TO_ANGSTROM = 0.529177210903`, `Energy::HARTREE_TO_KCALMOL = 627.5094740631`, `Energy::HARTREE_TO_EV = 27.211386245988`) |
| `-DUSE_LEGACY_UNIT_CONSTANTS=ON` | exactly the value that site used before, bit for bit |

Each site keeps its old value visible as the argument of
`CurcumaUnit::Length::bohr_radius_or_legacy(...)`, `angstrom_to_bohr_or_legacy(...)`,
`Energy::hartree_to_kcalmol_or_legacy(...)` or `hartree_to_ev_or_legacy(...)`.
The legacy build exists for reproducing earlier results and for exact comparison with reference
codes that carry their own constants.

## The values that were in use

Bohr radius in Å (for Å→Bohr factors the reciprocal is listed); "rel." is relative to CODATA 2018.

| value | origin | rel. | where (before Sep 2026) |
|---|---|---:|---|
| 0.529177210903 | CODATA 2018 | 0 | `units.h`; GFN-FF input geometry; native xTB `AA_TO_AU`; D4 generator (1.8897261246257702); react radii |
| 0.5291772105638411 | ≈ CODATA 2014 | −6.4e-10 | `gfnff.h` Bohr→Å of the GFN-FF numerical gradient (the forward conversion used 2018, so the round trip was off) |
| 0.52917721092 | CODATA 2010 | +3.2e-11 | `global.h` `au` (used for gradient unit conversion by many methods), `STOIntegrals.hpp`, `xtb_fragment_scf.cpp` |
| 1.88972612546 | 2018, truncated | −4.4e-10 | D3 parameter generator (two sites) |
| 1.889726125 | 2018, truncated | −2.0e-10 | UFF/QMDFF workspace `m_au` |
| 1.8897259886 = 1/0.529177249 | **CODATA 1986** | +7.2e-8 | GFN-FF CPU CN radii (`cn_calculator.cpp`, `gfnff_method.cpp`), GBSA, D4 charge model, dead GPU/ROCm code |
| 0.52917726 | GFN-FF Fortran literal (`gfnff_param.f90:380`, not a CODATA value) | +9.3e-8 | GFN-FF reference radius table `covalent_rad_d3` (also used by the GPU CN), Phase-1 EEQ topological distance, ALPB |
| 0.529177 | 6 digits | −4.0e-7 | ANCopt step control, EHT basis coordinates, one GFN-FF printout |

Energy:

| value | origin | where |
|---|---|---|
| 627.5094740631 | CODATA 2018 | most of the code |
| 627.503 | older | UFF torsion/inversion/vdW force defaults (`forcefieldgenerator.cpp`), rel. −1.0e-5 |
| 27.211386245988 | CODATA 2018 | `units.h` |
| 27.2113957 | older (xtb) | GFN-FF Hückel orbital-energy scaling (`huckel_solver.h`), rel. +3.4e-7 |
| 27.211 | truncated | four HOMO-LUMO printouts |

For comparison, the reference programs themselves disagree: the pprcht/gfnff **library** uses
0.52917726 for its parameter tables, its **app** reads coordinates with 0.529177249 (CODATA 1986),
and curcuma reads coordinates with CODATA 2018. On `polymer_2x` that input convention alone
accounted for most of the remaining per-term deviation against pprcht (bond +0.0241 → +0.0036,
repulsion −0.0296 → +0.0011 kcal/mol when pprcht was temporarily rebuilt with curcuma's constant).

## Effect of the switch

Measured Sep 25, 2026 with `scripts/refset_regression.py` (every structure of the set).

**Legacy build vs the binary before the change** — must be, and is, bit-identical:

| method | MOR41 (95) | GMTKN55 (2462) |
|---|---|---|
| gfnff | 0 differences | 0 differences |
| gfn2 | 0 | 0 |
| gfn1 | 0 | 0 |

**Default (CODATA 2018) vs the binary before the change**:

| method | MOR41 | GMTKN55 |
|---|---|---|
| gfnff | 54 changed, max 5.0e-5 kcal/mol | 1027 changed, max 1.5e-4 kcal/mol |
| gfn2 | 0 changed | 0 changed (geometry was already CODATA 2018; `au` only scales gradients) |
| gfn1 | 0 changed | 2 changed, max 6.3e-6 kcal/mol (D3 constant) |

GFN-FF against the port reference pprcht (per-structure MAD over MOR41 + GMTKN55, 2555
structures): **0.00006 kcal/mol in both builds** — the CODATA shift is far below the known
port residuals (max 0.0445, `MB16-43/04`). CPU and GPU agree to all ten printed decimals in both
builds (triose, polymer 1410 atoms, polymer_2x 7320 atoms); before the CN radii were unified
they differed by up to 5.0e-7 Eh.

## When to use the legacy build

- Reproducing a number computed before Sep 25, 2026 exactly.
- Bit-level comparison with a reference code that uses one of the old constants, e.g. the
  pprcht/xtb GFN-FF tables (0.52917726). The default is within 1.5e-4 kcal/mol of it on the
  reference sets, so this only matters for sub-1e-6 Eh work.

Configure with `cmake -DUSE_LEGACY_UNIT_CONSTANTS=ON ..`. It is a compile-time switch on purpose:
the constants are `constexpr` and several feed static parameter tables.
