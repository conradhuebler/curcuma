# SQM Solvation — Work Packages (native GFN-FF / GFN1 / GFN2)

> Status: 🤖 AI-generated, in progress. Started 2026-06-07. Operator-approved scope:
> **all three models (ALPB, GBSA, CPCM)**, **self-consistent in the SCF**, **full
> device-resident on GPU**. Plan: `~/.claude/plans/was-wir-noch-brauchen-linked-porcupine.md`.
>
> This file is the executable hand-off. Each WP lists concrete files, functions,
> line anchors (approximate — re-grep), and acceptance criteria.

## ✅ DONE (2026-06-07): WP0 + WP1 + WP2 — self-consistent ALPB & GBSA in native GFN1/GFN2 (CPU)
> ⚙️ Machine-tested, **not** human production tested. Operator decided: full total ΔG
> (Born + CDS + shift), validated per-WP against `release_tblite/`.
> - **WP0**: `scripts/extract_tblite_solvation_params.py` → `src/core/solvation/tblite_solvation_params.h`
>   (76 entries: alpb {epsv,c1,soset,sx} + cds {rprobe,gamscale,sqrtGhbond} + shift {smass,rhos,gshift}).
>   `scripts/extract_tblite_cm5.py` → `src/core/solvation/cm5_data.h` (CM5 a0 + atomic radii).
>   `dump_tblite_reference --solvent` extended; 56 ALPB refs in `reference_data/`.
> - **WP1**: `ALPBSolvation` (`ff_methods/alpb_solvation.{h,cpp}`) method-gated tblite path
>   (D3 radii kept; tblite params; CDS tension `gamscale·4π·1e-5` + H-bond folded into the
>   Born diagonal; CM5 for gfn1). Wired into `xtb_native.{h,cpp}` (member + InitialiseMolecule
>   create + per-geometry `update` + in-SCF `addSolvationPotential` + `energySolvation` into
>   `m_E_total` + gradient before `/=au`). PARAMs in `xtbinterface.h`, `applyXtbScfConfig`.
>   Usage: `-method gfn2 -xtb.solvent water -xtb.solvent_model 3`.
> - **KEY FIX**: curcuma `m_sasa` lacks the 4π that tblite's surface carries → CDS tension
>   needs `·4π`, CDS H-bond `-kcaltoau·g²` (the 4π cancels). Born/keps/radii were already exact.
> - **Validation**: 56/56 refs ≤1e-8 Eh vs tblite (worst 6.4e-9); FD gradient PASS (gfn1 ~1e-8,
>   gfn2 ~1e-4 = base gfn2 D4-response floor). `ctest -L _solvation` (56 = 28 gfn1 + 28 gfn2;
>   the per-method labels are `gfn1_solvation`/`gfn2_solvation`, but ctest ANDs multiple `-L`
>   so use the shared `_solvation` substring) + `ctest -R xtb_solvation_numgrad`; gas-phase
>   byte-identical (gfn{1,2}_validation 24/24).
> - **WP2 (GBSA, 2026-06-07)**: `-xtb.solvent_model gbsa` (descriptive names
>   `none`/`cpcm`/`gbsa`/`alpb`, case-insensitive; legacy numeric 0..3 still accepted —
>   `solventModelCode()` in xtb_native.h; the PARAM is now a String). GBSA = ALPB with `alpbet=0`
>   (so `keps=1/eps-1`, shape term drops) + the classical **Still** Born kernel
>   (tblite api/solvation.f90:191 sets `born_kernel%still` for GBSA versions 21/22,
>   vs P16 for ALPB 11/12). `ALPBSolvation::setUseAlpb(false)` branches the param
>   load (`getTbliteSolvationParam(...,alpb=false)` → gbsa entries, already in the
>   header), `buildBornMatrix`, and a new `addGradientStill`; CDS/SASA/HB/CM5/shift
>   shared. **Bug fixed en route**: `applyXtbScfConfig` read `solvent_model` only via
>   `is_number_integer()`, but the CLI routes it as a float (`2.0`) — so the value was
>   silently dropped and ALPB-by-default ran regardless (ALPB tests passed only because
>   the default is 3). Now accepts any number/numeric-string. Dumper `--model gbsa`
>   (version 21/22), `validate_sqm.py` maps `solvent_model` → `-xtb.solvent_model 2/3`.
>   **Validation**: 56 gbsa refs ≤1e-8 Eh (`ctest -L _gbsa`); FD gradient 16/16
>   (alpb+gbsa, `ctest -R xtb_solvation_numgrad`); ALPB 56/56 + gas-phase 24/24 unchanged.
> - **WP4a (GPU correctness, 2026-06-07)**: GFN1+solvent on GPU was already correct
>   (host-driven path folds `addSolvationPotential` into the uploaded `v_ao`), but
>   **GFN2+solvent+GPU was wrong by up to 3.8 mEh** (ALPB) / 1.8 mEh (GBSA): the
>   device-potential / resident-loop path builds the potential on-device and omitted
>   the reaction field, so the SCF converged to the gas-phase fixpoint and the
>   post-SCF `energySolvation()` was evaluated at unpolarized charges. Fix: gate the
>   device-potential path off when `m_solvation` is active (`xtb_native.cpp:~722`),
>   so GFN2+solvent falls back to the host-driven device loop (in-SCF, GPU eigensolve).
>   Now bit-identical to CPU (0.0 µEh, H2O→caffeine). New `ctest -L gpu_*_solvation /
>   gpu_*_gbsa` (28, water × 7 mol × {gfn1,gfn2} × {alpb,gbsa}); gpu component
>   suite 121/121, gas-phase gpu 24/24 unchanged.
> - **WP4b (full GPU residency for GFN2+solvent, 2026-06-07, DONE)**: the device
>   potential build now adds `v_at += B·q_at` in-SCF (cuBLAS dgemv on the resident
>   `dQat`), so GFN2+solvent uses the fully device-resident loop instead of the WP4a
>   host fallback. New seam `GpuScfBackend::{supportsDeviceSolvation,beginSolvation}`
>   (`XtbGpuContext::beginSolvation` uploads the nat×nat Born matrix once per geometry;
>   reset in `beginPotential`); `ImplicitSolvationModel::deviceBornMatrix()` returns
>   `m_born_mat.data()` for GFN2 (Mulliken) and nullptr for the GFN1 CM5 path (which
>   keeps the host loop). No device **energy** change — the host recomputes the total
>   from the downloaded charges; only the in-SCF charges needed the reaction field.
>   The WP4a gate now blocks the device path only when the device can't take solvation
>   (CM5 / no backend support). Verbosity-2 logs "…+solvation) built on the device".
>   **Validation**: GFN2+solvent+GPU bit-matches CPU (0.0 µEh, H2O→caffeine, both
>   models); `ctest -L gpu_*_solvation|gpu_*_gbsa` 28, gpu component 121/121, gas-phase
>   gpu 24/24, CPU suite 112+1+24 unchanged. NOTE: residency, not a measured `-sp`
>   speedup (eigensolve-bound).
> - **WP5a (GFN-FF routing + model select, 2026-06-07, DONE)**: `-gfnff.solvent water`
>   now works. Root cause: the native `gfnff` PARAM module declared no `solvent`, so the
>   CLI value never reached `GFNFF::InitialiseMolecule` (`m_parameters["solvent"]`) and
>   solvation was silently gas-phase. Registered `solvent` + `solvent_model` (`alpb`
>   default / `gbsa`; numeric 2/3 accepted) in `gfnff.h`; wired `solvent_model` into the
>   ALPB creation (`setUseAlpb`, Still kernel for gbsa); fixed the hard-coded "ALPB" log
>   label. Energy now changes sensibly (H2O gfnff: gas −0.32726 → water ALPB −0.33982,
>   GBSA −0.33920). `test_gfnff_numgrad --solvent NAME [--solvent-model gbsa]` added.
>   **SUPERSEDED by the "✅ DONE WP5" section at the top of this file:** the WP5a "no
>   external reference / energy unvalidated / WP5b-deferred" status no longer holds —
>   the reference is xtb 6.7.1 (`--gfnff --alpb`), the self-consistent EEQ coupling
>   (A_eeq += B, what WP5a called "WP5b") is implemented, and the ALPB energy matches
>   xtb to ≤1e-8 Eh. GBSA remains approximate (params not extracted).
> - **Remaining**: CPCM (WP3); gfnff GBSA param extraction; GPU-solvation test coverage.
>   Wrapper `tblite-gfn2 -tblite.solvent_model 3` segfaults (pre-existing TBLiteInterface
>   bug, unrelated; native path unaffected).

## 0. Context & key facts (read first)

- **Reference lives in-tree:** `external/tblite/src/tblite/solvation/` contains the
  full reference: `alpb.f90`, `born.f90`, `cpcm.f90`, `cpcm_dd.f90`, `surface.f90`,
  `cds.f90`, `shift.f90`, type defs under `data/` and the **parameter data** under
  `data/alpb/param_{alpb,gbsa}_<solvent>.fh`. The `.fh` files contain **both**
  `gfn1_alpb_<solvent>` **and** `gfn2_alpb_<solvent>` (and gfnff) variants.
- **tblite `alpb_parameter` = `{epsv, c1, soset, sx(94)}`** (`data/alpb.f90`):
  epsv = dielectric, c1 = Born radius scaling, soset = Born offset, sx(94) =
  per-element dielectric descreening. **CDS (surface tension) and shift are
  separate, charge-independent containers** (`cds.f90`, `shift.f90`) — they enter
  the total ΔG_solv but NOT the wavefunction/potential.
- **tblite ALPB kernel == curcuma's existing `alpb_solvation.cpp`** (verified):
  - P16 kernel default (`born_kernel%p16`), Still optional.
  - `keps = (1/ε − 1)/(1 + αβ)`, `αβ = 0.571412/ε` (`alpb.f90:189-190`).
  - ALPB shape correction `+ keps·αβ/adet` to every matrix element (`alpb.f90:273-277`).
  - OBC Born radii `obc = [1.0, 0.8, 4.85]`, descreening default 0.8, born_scale
    default 1.0, born_offset default 0.0 (`born.f90:52-55`); `rho = vdwr·descreening`,
    `svdw = vdwr − born_offset`.
  - ⇒ The curcuma Born machinery is **reusable for GFN1/GFN2** by feeding tblite's
    `(epsv, c1, soset, sx)` for the chosen method. Only parameter VALUES differ.
- **In-SCF coupling math (exact):** charge-dependent solvation energy is
  `E = ½ qᵀ B q` with symmetric `m_born_mat` B (already scaled by keps and including
  Born self-energy, HB diagonal, ALPB shape). Therefore the Fock potential is
  `v_i = dE/dq_i = (B·q)_i`. SASA/gshift are charge-independent → no potential.
- **Double-count avoidance:** the SCF energy is `Tr(P·H0) + energyCoulombShell() +
  energyThirdOrder() + energyMultipole()` (`xtb_native.cpp:1147`), H0 = bare. Add a
  separate `energySolvation()` (½ form) to `m_E_total` — do NOT fold it into the band
  term. This mirrors the Coulomb ES2 handling exactly.
- **CLI routing gap (must fix):** `-solvent` is registered in 4 modules
  (gfnff_external/tblite/ulysses) but NOT in native `gfnff`/`xtb` → the flat flag is
  ambiguous → it stays in the command module and never routes. Confirmed: native
  `gfnff -solvent water` returns the gas-phase energy (silently ignored). Native
  GFN1/GFN2 must register the PARAMs in the `xtb` module and be used via the dotted
  form `-xtb.solvent water` (flat `-solvent` stays ambiguous by design).
- **Validation build:** the default `release/` build has **no TBLite** (`gfn2:tblite`
  errors). Use `release_tblite/` for live tblite comparison, or the committed
  reference suite (`dump_tblite_reference`, `scripts/validate_sqm.py`, ctest
  `gfn{1,2}_validation`).

### Already done (2026-06-07, release/ green, default path byte-unchanged)
- `src/core/solvation/implicit_solvation.h` — abstract `ImplicitSolvationModel`
  `{init(Z,solvent,method), update(Z,xyz), addPotential(q,v_at), energy(q), addGradient}`.
- `ALPBSolvation` (`ff_methods/alpb_solvation.{h,cpp}`) now implements it:
  - `addPotential(q_at, v_at)` → `v_at += m_born_mat * q_at` (after `getEnergyParts`).
  - `energy()` wrapper → `getEnergy()`.
  - `init(..., const std::string& method = "gfnff")` + `m_host_method` member.
  - GFN-FF call site unchanged (default argument).

---

## ✅ DONE (2026-06-07): WP5 — GFN-FF self-consistent ALPB ("the riddle")
> ⚙️ Machine-tested. The GFN-FF ALPB was a frozen-charge post-hoc add-on (energy at
> gas EEQ charges) and unvalidated. **Fix:** couple the Born reaction field into the
> *linear* EEQ solve — `A_eeq(:n,:n) += B` (the symmetric Born matrix) before solving,
> exactly the reference `gfnff_engrad.F90:1346-1350` (`A(:n,:n)+=gbsa%bornMat`). EEQ is
> linear, so one solve gives the self-consistent (solvated) charges; the GFN-FF Coulomb
> term then uses them and the gsolv parts (gborn+ghb+gsasa+gshift) are added separately
> — no double count. The reference EEQ gradient is itself frozen-charge ("without dq/dr",
> :308), so curcuma matches it by using the solvated charges in the existing addGradient.
> - Reference source: `/home/conrad/src/curcuma/external/gfnff` (read in place; network
>   clone failed). Reference numbers: `/opt/bin/xtb` 6.7.1 (`--gfnff --alpb SOLVENT`).
> - `EEQSolver::setReactionField(const Eigen::MatrixXd*)` + injection in
>   `calculateFinalCharges` before `dispatchSolve`; `ALPBSolvation::bornMatrix()`;
>   `gfnff_method.cpp` sets the field in `prepareCNAndEEQ` (update→setReactionField) and
>   forces the Phase-2 re-solve when solvated (`eeq_charges_current && !m_solvation`).
> - **Validated:** native gfnff+ALPB total = xtb to ≤1e-8 Eh, 7 mol × {water,dmso,acetone,
>   chloroform} (`ctest -L gfnff_solvation`, refs `reference_data/gfnff_alpb_xtb.ref.json`);
>   gradient FD `gfnff_numgrad_alpb_water`; gfnff gas byte-identical; native gfn1/gfn2 untouched.
> - **GBSA approximate** (~1-3 mEh): reuses the ALPB param set + Still kernel; the dedicated
>   gfnff GBSA params (`gfn1_/gfn2_<solvent>` in the reference) are not yet extracted (warns).
> - Remaining: extract gfnff GBSA params; self-consistent gradient adjoint (optional — the
>   reference is frozen-charge); GPU residency for gfnff solvation.

## WP0 — Parameters + reference harness

**WP0.1 Verify tblite field units before extracting.** Read
`external/tblite/src/tblite/solvation/data/alpb.f90` (`get_alpb_param`) and
`alpb.f90` (`create_alpb_input`, `new_alpb`). Confirm how `epsv/c1/soset/sx` are
consumed (units!). curcuma's `loadSolventParameters` does `born_offset = soset *
0.1 * aatoau` and `born_scale = c1` — **the `0.1*aatoau` is gfnff-flavoured; tblite
likely uses soset directly (Bohr or Å).** Nail this or every non-water solvent will
be subtly wrong.

**WP0.2 Extraction script.** Write `scripts/extract_tblite_solvation_params.py` to
parse `data/alpb/param_{alpb,gbsa}_*.fh` → C++ header
`src/core/solvation/tblite_solvation_params.h`. Layout: a lookup keyed by
`(method ∈ {gfn1,gfn2,gfnff}, model ∈ {alpb,gbsa}, solvent)` returning
`{epsv, c1, soset, sx[94]}`. ~24 ALPB + ~14 GBSA solvents × 3 methods. Mechanical,
verifiable (spot-check `gfn2_alpb_water` epsv=80.2, c1=1.47438678, sx[0]=0.18678116).

**WP0.3 CDS + shift parameters.** Extract surface-tension (CDS, `data/cds.f90` /
`param_*` cds blocks) and the per-(method,solvent) `gshift` (`data/shift.f90`).
Needed for total-ΔG matching (not for the wavefunction). Defer if only the
electrostatic ΔG is being validated first.

**WP0.4 Reference dump.** Extend `dump_tblite_reference` + `scripts/validate_sqm.py`
to emit tblite **ALPB** and **CPCM** energies + gradients for the 12-molecule suite
× {gfn1,gfn2} × {water,dmso,acetone,chloroform}. Skip GB (C-API broken). Build
`release_tblite/` to generate.

**Acceptance:** generated header compiles; spot-checked values match `.fh`; reference
JSON written for the suite.

---

## WP1 — Self-consistent ALPB in native GFN1/GFN2 (CPU)

**WP1.1 Method-aware parameters.** In `alpb_solvation.cpp::loadSolventParameters`
(~:263), branch on `m_host_method`: for `gfn1`/`gfn2` use the WP0 tblite header
`(epsv,c1,soset,sx)` and set **`m_lhb = false`** (tblite ALPB has no HB term in the
Born matrix; HB-like effects live in CDS) and gfnff SASA gamscale → 0; for `gfnff`
keep the current path. Verify `m_born_offset`/`m_born_scale` unit handling per WP0.1.

**WP1.2 Wire into the native SCF** (`xtb_native.{h,cpp}`):
- Add member `std::unique_ptr<ImplicitSolvationModel> m_solvation;` and
  `std::string m_solvent = "none"; int m_solvent_model = 0;`.
- Read config in `InitialiseMolecule`/ctor from `m_parameters["xtb"]["solvent"]`
  (+ top-level fallback, mirror `native_xtb_method.cpp` d4_charge_source pattern,
  ~:91-96). If solvent != none → create model (model 3→ALPB now; 2→GBSA WP2; 1→CPCM
  WP3), `init(Z, solvent, method)`.
- **Geometry update:** call `m_solvation->update(Z, m_geometry_bohr)` once per
  geometry (right after CN / γ build, near `buildGammaMatrix`).
- **In-SCF hook:** in the `if (!use_device_potential)` block (`xtb_native.cpp:803-819`),
  after `addThirdOrderPotential(m_pot)` (:810), add a private
  `addSolvationPotential(m_pot)` that does `if (m_solvation) m_solvation->addPotential(m_wfn.q_at, m_pot.v_at);`.
  (v_at is folded into v_ao by `expand_potential` → enters Fock AND the host-GPU
  upload automatically.) Mirror in the device-post rebuild (:1138-1142) and the
  gradient rebuild (:1188-1191).
- **Energy:** add `energySolvation()` returning `m_solvation ? m_solvation->energy(m_wfn.q_at) : 0.0`;
  add `m_E_solvation` to `m_E_total` (`xtb_native.cpp:1154`) and to the verbose
  decomposition (:1158-1169). Do NOT add into `Tr(P·H0)`.
- **Gradient:** in the `if (gradient)` block (~:1185), after the potential rebuild,
  `if (m_solvation) m_solvation->addGradient(Z, m_geometry_bohr, m_wfn.q_at, grad)`.
  Note units: ALPB gradient is Eh/Bohr; `xtb_native` converts at the end — add the
  solvation gradient in the SAME space as the other terms (before the `/= au`/unit
  step) — verify against the existing D4 gradient injection point.

**WP1.3 PARAM registration.** In `xtbinterface.h` `BEGIN_PARAMETER_DEFINITION(xtb)`
(:33) add (names/semantics mirror tblite for consistency):
`PARAM(solvent, String, "none", ...)`, `PARAM(solvent_model, Int, 0, "0=none,1=CPCM,2=GBSA,3=ALPB")`,
`PARAM(solvent_epsilon, Double, -1.0, ...)`. Run `make GenerateParams`, check no
validation warnings. Usage: `-xtb.solvent water -xtb.solvent_model 3` (dotted; flat
`-solvent` stays ambiguous). Optionally also register `solvent` in the native gfnff
module to fix the GFN-FF routing gap (WP5).

**WP1.4 Validation.** `release_tblite/`: compare native `gfn1`/`gfn2` + ALPB ΔG_solv
(electrostatic Born part) vs tblite ALPB on the suite. Target ≤1e-6 Eh for the Born
term once units (WP0.1) and `m_lhb=false` are correct. Full total needs CDS+shift
(WP0.3). FD gradient test (new ctest, analog `test_gfnff_gradients`) on
H₂O/CH₄/methanol, tol 1e-5 Eh/Bohr.

**Acceptance:** `solvent=none` byte-identical to current; `-xtb.solvent water`
changes the energy and polarizes charges; Born ΔG matches tblite ≤1e-6; FD gradient
passes; regression ctest green.

---

## WP2 — Self-consistent GBSA  ✅ DONE (2026-06-07, see top of file)
> Implemented by extending `ALPBSolvation` (setUseAlpb + Still kernel), not by
> resurrecting `gbsa.{h,cpp}`. The notes below are the original spec.

- Consolidate the second, unused GBSA implementation `src/core/solvation/gbsa.{h,cpp}`
  (GBOBC-II, energy+gradient) under `ImplicitSolvationModel`: add `addPotential`
  (= Born-matrix·q, same math as ALPB), method-aware params from WP0 (`param_gbsa_*`),
  `update`/`energy`/`addGradient` to the interface signatures.
- Selected via `solvent_model = 2`. GBSA = ALPB with `αβ = 0` (no shape term).
- **Validation:** vs Ulysses GBSA and, where available, xtb `-gbsa`.

**Acceptance:** model selectable; ΔG matches reference; FD gradient passes.

---

## WP3 — Native CPCM (new)

- New `src/core/solvation/cpcm.{h,cpp}` implementing `ImplicitSolvationModel`, ported
  from `external/tblite/src/tblite/solvation/cpcm.f90` + `cpcm_dd.f90` + `surface.f90`.
- Molecular surface from the existing Lebedev/SASA machinery (`lebedev_grid.h`,
  ALPB SASA). Build the surface interaction matrix S; solve `S·q_surf = −f(ε)·V_solute(grid)`;
  the reaction-field potential at the atoms is the in-SCF `v_at` (self-consistent).
- Selected via `solvent_model = 1`; only ε needed (`solvent_epsilon`, or auto from name).
- **Validation:** vs tblite CPCM (WP0.4). CPCM tol may be looser if surface
  discretization differs.

**Acceptance:** converges; ΔG within agreed tol of tblite CPCM; FD gradient passes.

---

## WP4 — GPU device-resident solvation

- **Free win:** WP1 already works on the default host-GPU path — the solvation term
  is folded into host `v_at`→`v_ao`, which is uploaded per iteration anyway
  (`xtb_native.cpp:~831-866`). Confirm `gpu_gfn{1,2}_validation` stays bit-stable with
  a solvent.
- **Fully-resident Stage-6 loop** (the `-scf_gpu_potential` path, where the potential
  is built entirely on device): add a CUDA kernel for the Born potential `B·q`
  (O(N²)) + energy at the `GpuScfBackend` seam (`xtb_native.h:~314+`,
  `cuda/xtb_gpu_context.{h,cu}`): `beginSolvation()` (upload Born matrix / radii once
  per geometry), `residentAddSolvationPotential()` (inside the device potential
  build), solvation energy scalar in the device energy reduction, gradient in Stage 4
  or host-post.
- CPCM device surface solve is heavier → do last.

**Acceptance:** device vs CPU energy bit-stable; gradient ~1e-8 Eh/Å at
`-scf_threshold 1e-8`; `ctest -L gpu_scf|gpu_gfn1_validation|gpu_gfn2_validation` green.

---

## WP5 — GFN-FF self-consistent + validation

- Fix the GFN-FF CLI routing gap: register `solvent` in the native gfnff PARAM module
  so `-solvent`/`-gfnff.solvent` reaches `gfnff_method.cpp:802`.
- Optional upgrade: couple ALPB self-consistently into the GFN-FF EEQ solve
  (`eeq_solver` + reaction field) instead of the current post-hoc add-on, so the EEQ
  charges feel the solvent. At minimum, validate the existing GFN-FF ALPB and document
  status in `docs/GFNFF_STATUS.md`.

**Acceptance:** GFN-FF `-solvent water` actually changes the energy; documented.

---

## Cross-cutting: docs & changelog
- Update `docs/SOLVATION.md` (add native GFN1/GFN2/GFN-FF columns), keep this WP file
  current, add the feature to the top-level CLAUDE.md capability list + README, one
  line in `AIChangelog.md`.
- Honesty per CLAUDE.md: only the operator sets ✅ TESTED/APPROVED. Document what was
  tested (which molecules/solvents) and what was not (system classes, CDS/shift if
  deferred, metals).

## Open questions for the operator
1. Validate each WP against `release_tblite/` as we go, or finish the native wiring
   first then validate in a batch?
2. Match tblite's **total** ΔG_solv (Born + CDS + shift) for "100%", or is the
   self-consistent **electrostatic** ΔG the immediate target (CDS/shift as a later WP)?
