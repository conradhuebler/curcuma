# Solvation Models in Curcuma

Implicit solvation support for GFN and MNDO methods via TBLite and Ulysses interfaces.

## Overview

Curcuma supports implicit solvation models for quantum chemistry calculations:

- **TBLite Interface**: CPCM, GB (Generalized Born), ALPB for GFN methods. ⚠️ **Pending**:
  this path needs `USE_TBLITE=ON`; the reference dev build (`release/`) has TBLite off, so
  the TBLite Interface solvation models are currently unverified in that build. Native ALPB/GBSA/COSMO/CPCM below are
  unaffected.
- **Ulysses Interface**: GBSA (Generalized Born) for GFN and MNDO methods
- **Native GFN1/GFN2**: ⚠️ self-consistent **ALPB** and **GBSA** (June 2026) —
  AI-implemented, machine-tested vs tblite (≤1e-8 Eh on the validation set), human
  production testing pending. See below.
- **Native GFN1/GFN2 COSMO**: 🤖 AI-generated, ⚙️ machine-tested only, human production testing pending:
  self-consistent **ddCOSMO** (Nov 2026), a port of xtb 6.7.1 `--cosmo` (section "Native GFN1/GFN2 COSMO").
- **Native GFN1/GFN2 CPCM**: 🤖 AI-generated, ⚙️ machine-tested only, human production testing pending: tblite 0.4.0 `--cpcm`
  on the same ddCOSMO solver (Nov 2026, section "Native GFN1/GFN2 CPCM").
- **Native GFN1/GFN2 CPCM-X**: 🤖 AI-generated, ⚙️ machine-tested only, human production testing pending: xtb 6.7.1 `--cpcmx`
  (conductor COSMO SCF plus the post-SCF COSMO-RS/SMD step of CPCM-X), single points only (Nov 2026, section "Native GFN1/GFN2 CPCM-X").
- **Native GFN-FF**: ⚠️ self-consistent **ALPB** (June 2026, WP5) — the Born reaction
  field couples into the EEQ charge solve (`A_eeq += B`), so the EEQ charges polarize in
  the solvent. Matches **xtb 6.7.1** (`--gfnff --alpb`) to ≤1e-8 Eh (7 mol × 4 solvents,
  `ctest -L gfnff_solvation`); analytic gradient FD-validated. GFN-FF has no separate GBSA
  model (the reference always uses ALPB, `gfnff_alpb.F90:23-24`), so `-gfnff.solvent_model
  gbsa` maps to ALPB. `-method gfnff -gfnff.solvent water -gfnff.solvent_model alpb`.
  Machine-tested only.

## Native GFN1/GFN2 ALPB / GBSA (no external dependency)

The native `gfn1`/`gfn2` backends include self-consistent ALPB and GBSA solvation
models matching the tblite parameterization (Born electrostatics + CDS surface
tension/H-bond + state shift, see "Reference states"). GFN1 uses CM5 charges, GFN2 Mulliken charges, as in
tblite. The solvent reaction field enters the SCF (the wavefunction polarizes), the
energy is added to the total, and the nuclear gradient is analytic. ALPB uses the
P16 Born kernel; GBSA is ALPB with the analytical-PB shape term off and the
classical Still kernel (exactly tblite's ALPB=11/12 vs GBSA=21/22 distinction).

```bash
# Dotted form required (the flat -solvent is ambiguous across providers):
curcuma -sp mol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model alpb
curcuma -sp mol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model gbsa
curcuma -opt mol.xyz -method gfn1 -xtb.solvent dmso                          # defaults to alpb
```

- `-xtb.solvent <name>`: water, dmso, acetone, chloroform, methanol, … (tblite set)
- `-xtb.solvent_model`: `alpb` (default when a solvent is given), `gbsa`, `cosmo` (ddCOSMO, below), `cpcm` (tblite's CPCM, below), `cpcmx` (CPCM-X, below, `-sp` only), `none`.
  Legacy numeric codes (3=alpb, 2=gbsa, 1=cpcm, 0=none) are still accepted; `cosmo` has no numeric code.
- ROCm (`-gpu rocm`): the reaction field v_at += B q_at enters the device SCF potential and the resident loop (`-xtb.rocm_device_solvation`, default true; energy and gradient on the host); see [GPU_TUNING.md](GPU_TUNING.md).

**Validation (machine-tested):** native total energy matches tblite to ≤1e-8 Eh for
7 molecules × {water, dmso, acetone, chloroform} × {gfn1, gfn2} — for **both** models
(ALPB: `ctest -L _solvation`, 56 tests; GBSA: `ctest -L _gbsa`, 56 tests); analytic
gradients FD-validated for both (`ctest -R xtb_solvation_numgrad`). **Not tested:**
other solvents, ions, metals, large systems, long MD/opt stability. See
[SQM_SOLVATION_WP.md](SQM_SOLVATION_WP.md).

## Supported Solvents

Native GFN1/GFN2 (ALPB and GBSA) and native GFN-FF (ALPB) share one name handling (`src/core/solvation/solvent_names.h`):
names are case-insensitive (`Water`, `CHCl3`), spaces and underscores are ignored, and the aliases of xtb 6.7.1 are
accepted. A name that has no parameter set for the chosen method and model is an **error** (exit code 1), never a silent
gas-phase run. `none`, `vacuum` and `gas` mean no solvent.

| Canonical name | Accepted spellings (besides the canonical one) |
|---|---|
| `water` | `h2o` |
| `ch2cl2` | `dichloromethane`, `dichlormethane`, `methylenechloride`, `dcm` |
| `chcl3` | `chloroform`, `trichloromethane`, `tcm` |
| `acetonitrile` | `methylcyanide` |
| `cs2` | `carbondisulfide`, `carbon disulfide` |
| `ether` | `diethylether`, `diethyl ether` |
| `dmf`, `dmso` | `dimethylformamide`, `dimethylsulfoxide` |
| `furane` | `furan` |
| `hexane` | `n-hexane`, `nhexane`, `n-hexan`, `nhexan` |
| `thf` | `tetrahydrofuran` |
| `ethylacetate` | `ethyl acetate` |
| `woctanol` | `octanol wet`, `wet octanol`, `wetoctanol` |

Also available under their own names: `acetone`, `aniline`, `benzaldehyde`, `benzene`, `dioxane`, `ethanol`, `hexadecane`, `methanol`,
`nitromethane`, `octanol`, `phenol`, `toluene` (25 ALPB solvents; the GBSA sets are smaller: GFN1 12, GFN2 14, as in tblite and
xtb, e.g. `dmf` and `hexane` have no GFN1 GBSA parameters). `dichloromethane`, `diethyl ether`, `carbon disulfide`, `ethyl acetate`
and `octanol wet` are curcuma additions. The alias list is that of the xtb source tree in `external/xtb` (`normalizeSolventName`,
newer than the 6.7.1 binary): the 6.7.1 binary accepts only `furane`, `h2o`, `n-hexane` and the canonical names (measured: it rejects `dcm`,
`dichloromethane`, `methylenechloride`, `furan`, `diethylether`, `tetrahydrofuran`, `trichloromethane`, `methylcyanide`), curcuma accepts the longer list.

### Reference states (`-solvent_state`)

One flat argument for GFN1, GFN2 and GFN-FF (`-solvent_state gsolv|bar1M|reference`, case-insensitive; the dotted
`-xtb.solvent_state` and `-gfnff.solvent_state` also work). The state shifts the solvation free energy as `xtb --alpb water [state]`:

| State | Meaning | Shift (298.15 K) |
|---|---|---|
| `gsolv` (default) | gas and solution at the same concentration | 0 |
| `bar1M` | 1 bar ideal gas to 1 mol/L solution | +3.0312e-3 Eh (+1.902 kcal/mol), all solvents |
| `reference` | 1 bar ideal gas to the pure solvent, as in xtb | solvent dependent, e.g. water +1.8138e-3 Eh |

One implementation serves all methods (`solutionStateShift`): ideal-gas molar volume 24.79 L/mol as in xtb and pprcht/gfnff,
CODATA-2018 kB. Against xtb 6.7.1 the shifts agree to 3e-9 Eh (`ctest -R cli_sqm_18`, 3 solvents x 3 methods). tblite derives
24.78960 L/mol from R instead of 24.79 and therefore differs by 1.4e-8 Eh in `bar1M`; it was not ported separately. The temperature is fixed at 298.15 K.

### Not implemented natively

The ionic strength of xtb's GBSA, a choice of the Born kernel (ALPB always P16, GBSA always Still), the Lebedev grid size of ALPB/GBSA (230)
and, for ALPB/GBSA, a solvent given only by its dielectric constant (COSMO and CPCM accept one). GFN-FF with `solvent_model gbsa` runs ALPB with a
warning, as pprcht/gfnff does; xtb 6.7.1 runs a Still-kernel GBSA with the GFN2 parameters there. GFN-FF has no COSMO and no CPCM.

## Native GFN1/GFN2 COSMO (ddCOSMO)

🤖 AI-generated (Nov 2026), ⚙️ machine-tested only; human production testing pending. Code: `src/core/solvation/ddcosmo.{h,cpp}`
(class `DdCosmoSolvation`), a line-by-line port of xtb 6.7.1 `src/solv/cosmo.f90` and its ddCOSMO backend `src/solv/ddcosmo/{core,solver}.f90`
(Lipparini, Stamm, Cances, Maday, Mennucci; Cances et al., J. Chem. Phys. 139, 054111 (2013); Lipparini et al., J. Chem. Theory Comput. 9,
3637 (2013)). Lebedev grids: `src/core/solvation/lebedev_xtb_grids.h`, generated from xtb `lebedev.f90` by `scripts/extract_lebedev_grids.py`.

```bash
curcuma -sp mol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cosmo       # = xtb mol.xyz --gfn 2 --cosmo water
curcuma -sp mol.xyz -method gfn1 -xtb.solvent 37.5 -xtb.solvent_model cosmo        # = xtb --cosmo 37.5 (dielectric constant)
curcuma -opt mol.xyz -method gfn2 -xtb.solvent_model cosmo -xtb.solvent_epsilon 20 # epsilon alone also switches COSMO on
```

- **Model (as xtb)**: cavity of spheres with the COSMO radii of xtb (`vdwradcosmo.f90`, radius scale 1.0); surface charge per sphere in real
  spherical harmonics up to `lmax` 6; Jacobi/DIIS solver to 1e-8; switching width eta 0.2; cavity grid = the Lebedev grid nearest to half of
  `-xtb.cosmo_grid` (230 -> 110 points), surface-area grid = 230 points. Energy E = 1/2 (1 - 1/eps) <X, Psi> + sum_i gamma_i A_i (solvent
  accessible surface, probe 1.3 A). Atomic charges: Mulliken for GFN2 and for GFN1 (xtb uses CM5 only with GBSA/ALPB, not with COSMO).
  No hydrogen-bond term, no empirical shift, no reference-state shift: `-solvent_state` other than `gsolv` is an error with cosmo.
- **Parameters**: a solvent name takes epsilon and the surface tensions of its ALPB set of the method (GFN1 or GFN2), exactly as xtb does
  (xtb `model.f90` with cosmo); all 25 ALPB solvents. A number, or `-xtb.solvent_epsilon` without a name, takes the surface tensions of xtb's
  default COSMO set (`param_cosmo.fh`, identical to the GFN1/ALPB/water set, checked). `-xtb.solvent_epsilon` with a name keeps the name's
  surface tensions and replaces epsilon. xtb's `--cosmo inf` (conductor, own parameter set `param_cosmo_inf.fh`) exists only as the
  COSMO step of CPCM-X (`DdCosmoSolvation::Model::XtbCosmoInf`), not as a `cosmo` solvent.
- **SCF coupling**: per SCF iteration one direct (L X = g) and one adjoint (L* S = Psi) solve, restarted from the previous solution; the
  potential is v = keps sqrt(4 pi) X_00 - J^T zeta. Gradient: surface term plus the ddCOSMO force kernels (fdoka/fdokb/fdoga) plus the two
  field terms; the adjoint is re-solved to 1e-3 of the threshold first, as in xtb.
- **PARAMs** (module `xtb`): `cosmo_grid` (230), `cosmo_lmax` (6), `cosmo_conv` (1e-8), `cosmo_eta` (0.2), `solvent_epsilon`.
- **GPU**: no device path. The model has no Born matrix, so `-gpu cuda|rocm` keep the host-driven SCF loop for the solvation potential.
  Run once with `-gpu cuda` (methanol/water, gfn2, RTX 5080): energy equal to the CPU run within 5e-12 Eh (n = 1); ROCm and Vulkan not run.

**Validation (Nov 22, 2026, xtb 6.7.1 `/opt/bin/xtb --gfn N --cosmo <solvent> --sp --acc 0.01`, curcuma `-scf_threshold 1e-9`):**
20 single points, {water, methanol, caffeine, NH4+, acetate} x {gfn1, gfn2} x {water, acetone}: max |E_curcuma - E_xtb| = 2.8e-8 Eh
(caffeine; its gas-phase difference without solvent is already 3.1e-8 Eh, the other molecules <= 5.2e-9 Eh); `--cosmo 37.5` and
`-xtb.solvent_epsilon 37.5` on methanol 3.1e-9 Eh. Analytic gradient against central differences of curcuma's energy, max over all coordinates:
h = 5e-4 A, 8 cases ({water, methanol, NH4+, acetate} x {gfn1, gfn2}, water or acetone): <= 4.7e-7 Eh/A (step truncation; acetate gfn2
falls from 4.5e-7 to 7.2e-8 at h = 2e-4 A); h = 2e-4 A: water gfn2 3.8e-8, acetate gfn2 7.2e-8, caffeine gfn2 1.6e-7, NH4+ gfn1/acetone
3.3e-7 Eh/A (5.3e-8 at h = 5e-4 A: the ddCOSMO threshold 1e-8 leaves an energy noise of about 1e-10 Eh that limits small steps).
Against xtb's analytic gradient: <= 1.7e-8 Eh/A for water, methanol, NH4+, caffeine; acetate differs by 5e-4 Eh/A, but already in the gas
phase (8.6e-4 Eh/A), where xtb's own finite difference (-3.62442e-3) agrees with curcuma's analytic value (-3.62446e-3) and not with xtb's
analytic one (-2.76e-3), so this is an xtb gradient issue outside COSMO (not investigated further). Test: `ctest -R cli_sqm_19`.
**Not tested**: other solvents than water/acetone and epsilon 37.5, metals and d elements (radius stub 2.223 A as xtb), large systems,
molecules with grid points on a sphere pole in the force kernel, geometry optimisation and MD stability, GPU, periodic systems (not
supported). Not implemented: `--cosmo inf` as a `cosmo` solvent (only inside CPCM-X), ionic strength, the `.cosmo` file output of xtb
(`--tmcosmo`).

## Native GFN1/GFN2 CPCM (tblite model)

🤖 AI-generated (Nov 2026), ⚙️ machine-tested only; human production testing pending. tblite 0.4.0 (commit 9ca469a)
`src/tblite/solvation/cpcm.f90` on the ddCOSMO class above (`DdCosmoSolvation::Model::TbliteCpcm`); solvent table
`src/core/solvation/tblite_solvent_data.h` (epsilon column of tblite `data.f90`, `get_solvent_data`).

```bash
curcuma -sp mol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcm      # = tblite run --method gfn2 --cpcm water
curcuma -sp mol.xyz -method gfn1 -solvent 37.5 -solvent_model cpcm                # = tblite --cpcm 37.5 (flat flags also work)
curcuma -md mol.xyz -method gfn2 -xtb.solvent_model cpcm -xtb.solvent_epsilon 20  # epsilon alone also switches CPCM on
```

tblite's `cpcm_dd.f90` is the ddCOSMO core of xtb (`core.f90` + `solver.f90`), read side by side: same switching function, neighbour
list, Y_lm basis, Jacobi/DIIS (25 vectors, drop 10), H^-1/2 norm, 200 iterations, adjoint re-solved to conv * 1e-3 for the gradient;
it returns gradients instead of forces (fdoka/fdokb/fdoga with the opposite sign), the same total. So no second core was written;
CPCM differs from COSMO only in the parameters:

| | xtb COSMO (`cosmo`) | tblite CPCM (`cpcm`) |
|---|---|---|
| keps | 1/2 (1 - 1/eps) | 1/2 (1 - 1/eps) / (1 + 0.571412) |
| eta (switching width) | 0.2 | 2.0 (`cpcm_input` default; the CLI does not change it) |
| cavity grid | nearest to 230/2 = 110 points | nearest to 74 points |
| epsilon of a name | ALPB parameter set (water 80.2) | tblite solvent table (water 78.36), 112 names incl. `inf` |
| non-polar term | SASA, surface tensions of the ALPB set | none |
| radii, lmax, conv, charges | COSMO radii, 6, 1e-8, Mulliken (GFN1 too) | same |
| reference state | gsolv only (error otherwise) | gsolv only (tblite: "Solution state shift not supported for CPCM") |

`cosmo_grid`, `cosmo_lmax`, `cosmo_conv`, `cosmo_eta` do not apply to `cpcm`, which keeps tblite's fixed defaults; no new PARAMs.

**Validation (Nov 22, 2026, `tblite run --method gfnN --cpcm <solvent> --acc 0.0001`, curcuma `-scf_threshold 1e-10`):** 24 single points
{water, methanol, caffeine, NH4+, OH-, H2S} x {gfn1, gfn2} x {water, hexane}: max |E_curcuma - E_tblite| = 7.8e-10 Eh (caffeine gfn2), equal to
the gas-phase difference (caffeine gfn2 7.5e-10, methanol gfn2 5.0e-10 Eh). `--cpcm 37.5` (name and `-xtb.solvent_epsilon`) 4.9e-10,
pyridine (no ALPB set) 8e-11, `inf` 3e-14 Eh. Analytic gradient against `tblite --grad`: max 1.9e-9 Eh/A over the 24 cases. Against central
differences of curcuma's energy (h = 2e-4 A): methanol gfn2 3.7e-8, NH4+ gfn1/hexane 4.8e-9, caffeine gfn2 1.2e-7, OH- gfn2 3.1e-8,
H2S gfn1/hexane 1.4e-7, water gfn1 3.5e-8 Eh/A. Test: `ctest -R cli_sqm_20`.

**Properties of tblite's model that curcuma reproduces (not curcuma bugs):** for small neutral molecules the CPCM energy lies *above*
the gas phase (water gfn2: +2.0e-3 Eh in tblite and curcuma; methanol +8.9e-4 Eh), ions are stabilised (NH4+ -0.10 Eh). With eta = 2 the
switching function is 0.5 at the centre of a neighbour sphere, so grid points inside neighbour spheres stay partly in the cavity. The surface
is stiff: in a 500 fs NVE run of methanol (gfn2, water, n = 1) the total energy spans 4.2e-4 Eh (gas, ALPB, COSMO: 1.5-1.6e-4), with single-step
spikes; at the spike frame curcuma's energy equals tblite within 3e-11 Eh and its gradient tblite's within 1.6e-10 Eh/A, and the FD error
there falls as h^2 (3.9e-5 at h = 1e-3 A, 9.8e-8 at 5e-5 A), i.e. large curvature, not a wrong gradient.
**Not tested**: other solvents than water/hexane/pyridine/37.5/inf, metals and d elements (radius stub 2.223 A), large systems, long MD,
ROCm and Vulkan, periodic systems. CUDA keeps the host-driven SCF loop as for COSMO; one run (caffeine gfn2/water, RTX 5080, n = 1)
equals the CPU energy within 3e-14 Eh and the gradient within 3.5e-11 Eh/A. Not implemented: tblite's `cpcm_input` options `nang`, `rscale`, `eta`, `lmax` (tblite's CLI exposes
none of them either).

## Native GFN1/GFN2 CPCM-X

🤖 AI-generated (Nov 2026), ⚙️ machine-tested only; human production testing pending. Port of CPCM-X 1.1.0 (github.com/grimme-lab/CPCM-X,
LGPL-3.0; Stahn, Ehlert, Grimme, J. Phys. Chem. A 127, 7036 (2023)) and of the SMD part of numsa 0.2.0, as xtb 6.7.1 calls them
(`xtb/src/solv/cpx.F90`, `prog/main.F90`). Code: `src/core/solvation/cpcmx.h`, `cpcmx_data.cpp` (parameters, tables, database, solute input),
`cpcmx_crs.cpp` (charge averaging, rings, restoring free energy, Broyden solver, state correction), `cpcmx_smd.cpp` (SMD CDS, numsa surface);
driver `NativeXtbMethod::calculateCpcmx()`.

```bash
curcuma -sp mol.xyz -method gfn2 -xtb.solvent water -xtb.solvent_model cpcmx -xtb.cpcmx_db ~/CPCM-X/DB/xtb   # = xtb --gfn 2 --cpcmx water --sp
CURCUMA_CPCMX_DB=~/CPCM-X/DB/xtb curcuma -sp mol.xyz -method gfn2 -solvent hexane -solvent_model cpcmx
```

- **Steps (as xtb)**: (1) SCF in xtb's conductor set `gfn_cosmo_inf` (epsilon = inf, element radius scales, probe 2.949 A, its surface
  tensions), E_COSMO and the segments (position, charge, area, potential); (2) gas-phase SCF, E_gas (also all reported properties);
  (3) CPCM-X: sigma averaging (r_av, 2 r_av, orthogonalised), ring atoms, ideal-state energy E_COSMO - E_gas - eta RT, averaging correction,
  self-consistent solvent segment potentials (Broyden, 500 cycles, 1e-4 kcal/mol), solute chemical potential, state correction
  -RT ln(rho V_m / M), SMD CDS (probe 0.4); E = E_gas + dG_solv at 298.15 K. Verbosity 2 prints the contributions as xtb -v does;
  `getEnergyDecomposition()` has them as `cpcmx_dG_*`.
- **Solvent database**: read at run time from `-xtb.cpcmx_db <dir>`, else `$CURCUMA_CPCMX_DB`, else `<exe dir>/cpcmx_db`,
  `<exe dir>/../share/curcuma/cpcmx_db`; the files are `DB/xtb/<solvent>.cosmo` of CPCM-X (88 solvents, 8.9 MB, not shipped; the ctest
  carries water, hexane, acetonitrile). Their numbers are read in single precision, as xtb compiles them (default-real literals of the
  `.fh` include files; checked: `.fh` = `.cosmo` for all 88). Names: CPCM-X `solvent_name` (case-insensitive, e.g. `dmso`, `thf`, `dcm`);
  `methanol`, `xylene`, `benzonitrile`, `heptanol` have a name but no database file; an unknown name or a missing file is an error.
  xtb 6.7.1 rejects `dmso` (its CPCM-X is older); `dimethylsulfoxide` works in both.
- **Reproduced Fortran details**: single-precision constants (BtoA 0.529177, jtokcal, 1.15, 0.3, 0.43, 1.8); both covalent radii of a
  bond taken from the first atom (`det_bonds`); the Gauss-Seidel order of `iterate_solvent`; the SMD coordinates of the first n atoms with
  n = number of atoms owning a segment; xtb's file round trip of the solute (`-xtb.cpcmx_reference_rounding true`: positions, charges,
  areas, potentials to 9 decimals, energy to 10; false changes dG by <= 2e-10 Eh on four molecules).
- **SMD oxygen radius** (`-xtb.cpcmx_smd_oxygen_radius`, solvents with Abraham acidity < 0.43, i.e. not water): the CPCM-X in xtb 6.7.1
  enlarges `rad(8)`, the radius of the **8th atom whatever its element** (`atom8`, default, nothing below 8 atoms, so the result depends on
  the atom order); CPCM-X 1.1.0 enlarges every oxygen (`oxygen`); `none` disables it. Both add 1.8 (0.43 - alpha) (an Angstrom amount in
  the SMD paper) to a radius in Bohr. Found by comparison: the Fortran CPCM-X 1.1.0 with this one line changed to `rad(8)` reproduces
  xtb 6.7.1's CDS term to 6e-11 Eh (caffeine, glucoside, trisaccharide in hexane and acetonitrile), unchanged it does not (up to 2e-2 Eh); water in hexane differs between the rules by 4.1e-3 Eh, the trisaccharide in hexane by 2.1e-2 Eh.
- **No gradient** (xtb has none either): `hasGradient()` is false, `-sp -dump_gradient` is an error, every other command (`-opt`, `-md`,
  `-hessian`, ...) is rejected at start. CPU only (`-gpu` is an error), no `-large_system_mode`, `-solvent_state` gsolv only (CPCM-X applies
  its own state correction). GFN1 warns as xtb does (CPCM-X is parametrised for GFN2). Temperature, probe and solver settings are fixed
  at xtb's values (298.15 K, 0.4, 500, 1e-4).

**Validation (Nov 22, 2026, xtb 6.7.1 `/opt/bin/xtb mol.xyz --gfn N --cpcmx <solvent> --sp --acc 0.01`, curcuma `-scf_threshold 1e-10`):**
56 single points, {water, methanol, caffeine, OH-, NH4+, methyl alpha-D-glucopyranoside, trisaccharide C18H32O16 (66 atoms)} x {gfn1, gfn2}
x {water, hexane, acetonitrile, dmso}. Total energy: max |E_curcuma - E_xtb| = 4.8e-8 Eh except trisaccharide/gfn2 (4.8e-6 Eh, its gas-phase
D4 ATM difference; with `-xtb.d4_atm_cutoff 40.0` 6.1e-8). The differences equal the gas-phase differences of the native methods (max
6.1e-8). Contributions against CPCM-X 1.1.0 run on xtb's own `xtb.cosmo` (full precision, linked as a Fortran harness; it reproduces xtb's
total to <= 7e-10 Eh for water solvent): dG_IS 1.5e-8 (it contains E_COSMO - E_gas), dG_av 2.3e-9, dG_res 4.8e-9, dG_ss 1.3e-10, dG_CDS
4.6e-9 (`atom8` vs xtb) / 3.2e-10 (`oxygen` vs CPCM-X 1.1.0). Segments of water/gfn2 against xtb.cosmo: same 224 segments and atoms,
charges 1.0e-9, areas exact, positions 3.2e-7 Bohr and potentials 2e-8 (curcuma's CODATA-2018 Bohr vs xtb's 0.52917726). Run time
0.3-0.8 s for 24-66 atoms. Test: `ctest -R cli_sqm_21`.
**Not tested**: solvents other than water, hexane, acetonitrile, dmso; elements other than H, C, N, O; open-shell molecules; molecules
with buried atoms (no segment); temperatures other than 298.15 K. **Not implemented**: CPCM-X's own programs (SAC models, sigma
profiles, isodensity cavities, `.prop` and `smd_*` files in the working directory, which numsa reads when present), other parameter
sets than `xtb` (gtb, isodens, orca, tm), a gradient.

## Usage

### Basic Syntax

```bash
# Single point calculation in water
curcuma -sp molecule.xyz -method gfn2:tblite -solvent water

# Geometry optimization in DMSO
curcuma -opt molecule.xyz -method gfn2:tblite -solvent dmso

# Semi-empirical PM6 in acetone (requires Ulysses)
curcuma -sp molecule.xyz -method pm6:ulysses -solvent acetone
```

### Advanced Options

#### TBLite Solvation Models

```bash
# CPCM (Conductor-like Polarizable Continuum Model)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent water -solvent_model 1

# GB (Generalized Born)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent water -solvent_model 2

# ALPB (Analytical Linearized Poisson-Boltzmann)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent water -solvent_model 3
```

#### Custom Dielectric Constant

```bash
# Use custom epsilon (overrides solvent name)
curcuma -sp molecule.xyz -method gfn2:tblite \\
        -solvent_epsilon 78.4  # Water's dielectric constant
```

## Method Availability

| Method | Provider | Solvation | Models Available |
|--------|----------|-----------|------------------|
| **GFN2** | TBLite | ✅ Yes | CPCM, GB, ALPB |
| **GFN2** | Ulysses | ✅ Yes | GBSA |
| **GFN2** | Native | Yes (AI-generated, machine-tested) | ALPB, GBSA, COSMO (ddCOSMO), CPCM (tblite), CPCM-X (`-sp`) |
| **GFN1** | TBLite | ✅ Yes | CPCM, GB, ALPB |
| **GFN1** | Ulysses | ✅ Yes | GBSA |
| **GFN1** | Native | Yes (AI-generated, machine-tested) | ALPB, GBSA, COSMO (ddCOSMO), CPCM (tblite), CPCM-X (`-sp`) |
| **PM6** | Ulysses | ✅ Yes | GBSA |
| **PM6** | Native | ❌ Not yet | Planned |
| **AM1** | Ulysses | ✅ Yes | GBSA |
| **PM3** | Ulysses | ✅ Yes | GBSA |
| **MNDO** | Ulysses | ✅ Yes | GBSA |

## Compilation Requirements

### Enable TBLite Support

```bash
cd build
cmake -DUSE_TBLITE=ON ..
make -j4
```

### Enable Ulysses Support

```bash
cd build
cmake -DUSE_ULYSSES=ON ..
make -j4
```

### Check Compilation Status

```bash
curcuma --version  # Shows enabled features
```

## Output Examples

### TBLite with Water Solvation (Verbosity 1)

```bash
$ curcuma -sp h2o.xyz -method gfn2:tblite -solvent water -verbosity 1

[INFO] Initializing TBLite quantum chemistry method
[PARAM] accuracy: 1
[PARAM] SCF_maxiter: 100
[PARAM] temperature: 300.0 K
[PARAM] solvation_model: GB (Generalized Born)
[PARAM] solvent: water
[INFO] TBLite implicit solvation enabled

[ENERGY] TBLite Energy: -5.123456 Eh
```

### Ulysses with DMSO Solvation (Verbosity 1)

```bash
$ curcuma -sp molecule.xyz -method pm6:ulysses -solvent dmso -verbosity 1

[INFO] Initializing Ulysses quantum chemistry method
[PARAM] method: PM6
[PARAM] solvation_model: GBSA (Generalized Born)
[PARAM] solvent: dmso
[INFO] Ulysses implicit solvation enabled

[ENERGY] Ulysses Energy: -42.789012 Eh
```

## Solvation Models Explained

### CPCM (Conductor-like PCM)
- Treats solvent as conductor
- Fast and robust
- Good for polar solvents
- Reference: Cossi et al., J. Comput. Chem. 24, 669 (2003)

### GB (Generalized Born)
- **Default model for TBLite** (`-solvent_model 2`)
- Pairwise atom-based approximation
- Fast analytical gradients
- Excellent for optimization
- Reference: Still et al., J. Am. Chem. Soc. 112, 6127 (1990)

### ALPB (Analytical Linearized Poisson-Boltzmann)
- Advanced GB variant with parameter fitting
- High accuracy for small molecules
- Reference: Ehlert et al., J. Chem. Theory Comput. 17, 4250 (2021)

### GBSA (Generalized Born + Surface Area)
- **Default model for Ulysses**
- GB model + non-polar cavity term
- Includes atomic surface areas
- Reference: Hawkins et al., J. Phys. Chem. 100, 19824 (1996)

## Energy Decomposition

With `-verbosity 2` or higher:

```bash
[PARAM] electronic_energy: -42.123 Eh
[PARAM] core_repulsion: 15.456 Eh
[PARAM] solvation_energy: -0.234 Eh  # Implicit solvation contribution
[ENERGY] Total Energy: -26.901 Eh
```

## Performance Considerations

- **GB/GBSA**: Fastest (analytical gradients)
- **CPCM**: Medium speed (iterative solver)
- **ALPB**: Similar to GB, highly accurate

For geometry optimization, **GB is recommended** for speed.

## Troubleshooting

### Solvation Not Applied

**Problem**: Solvent parameter ignored
**Solution**: Specify method provider explicitly:
```bash
# Correct (explicit provider)
curcuma -sp mol.xyz -method gfn2:tblite -solvent water

# May fall back to native (no solvation)
curcuma -sp mol.xyz -method gfn2 -solvent water
```

### Compilation Errors

**Problem**: TBLite/Ulysses not found
**Solution**: Check CMake configuration:
```bash
grep "USE_TBLITE\|USE_ULYSSES" build/CMakeCache.txt
```

## References

1. **TBLite Library**: Grimme et al., J. Chem. Phys. 143, 054107 (2015)
2. **Ulysses Package**: Menezes et al., ChemRxiv (2023)
3. **GBSA Model**: Grimme, J. Comput. Chem. 32, 1456 (2011)
4. **ALPB Model**: Ehlert et al., J. Chem. Theory Comput. 17, 4250 (2021)

## See Also

- [Parameter Registry Guide](archive/PARAMETER_MIGRATION_GUIDE.md)
- [TBLite Interface Documentation](../src/core/energy_calculators/qm_methods/tbliteinterface.h)
- [Ulysses Interface Documentation](../src/core/energy_calculators/qm_methods/ulyssesinterface.h)
