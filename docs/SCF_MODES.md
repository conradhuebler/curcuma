# Native GFN SCF Convergence Modes

Selectable SCF strategies and initial guesses for the native GFN1/GFN2-xTB
implementation (`curcuma::xtb::XTB`, the canonical `gfn2` / `gfn1` backends).

> Status: AI-implemented, machine-tested. The default mode is now `broyden`
> (charge-vector quasi-Newton mixing, like tblite); it converges the 231-atom
> `complex` (previously divergent) from the bare guess and reaches the same
> energies as the historic DIIS on the established set (all native-GFN ctests
> pass). Human production testing pending.

## Why

The native SCF used a single fixed strategy — linear density mixing + Pulay DIIS
on the Fock matrix, with a hardcoded `diis_start = 5`. On large, polar systems
this extrapolates the density out of the convergence basin (classic charge
sloshing): the 231-atom `complex` diverged from a near-correct iter-0 energy to a
sign-flipped blow-up. tblite/xtb avoid this by mixing the low-dimensional SCC
**charge vector** with modified Broyden (quasi-Newton). Curcuma now does the same
by default (`broyden`); the historic DIIS path remains available as `-scf_mode
diis`. See [GFN2_SCF_STATUS.md](GFN2_SCF_STATUS.md) for the diagnosis.

Additionally the SCF parameters from the user config were never wired through
(`GFN2Method::updateGFN2Parameters()` was an empty stub) — the native object
always used its header defaults. That is now fixed: the settings below route
through the `xtb` config scope to the native object.

## SCF modes (`-scf_mode`)

| Mode          | What is mixed       | Behaviour                                                                 | Use when |
|---------------|---------------------|--------------------------------------------------------------------------|----------|
| `broyden`     | SCC charge vector   | Modified-Broyden (Johnson 1988) quasi-Newton mixing of the charge/multipole vector (`q_sh`, +`dp_at`/`qp_at` for GFN2). **The tblite/xtb mixer. Default.** Most robust. | Default; everything. |
| `diis`        | Fock matrix         | Damped warmup (`diis_start` iters) then Pulay DIIS. The historic path.    | Reproducing pre-2026-05-29 behaviour. |
| `plain`       | density matrix      | Linear density mixing every iteration, no DIIS. Slowest, simple.         | Debugging; baseline. |
| `level-shift` | density matrix       | Saunders-Hillier virtual-orbital shift on the Fock matrix + density mixing, shift faded near convergence. | Alternative for stiff systems. |

Aliases: `normal` -> `plain`, `levelshift`/`ls` -> `level-shift`, `pulay` -> `diis`.

**Why `broyden` is different.** `diis`/`plain`/`level-shift` mix the density or
Fock matrix; `broyden` mixes the low-dimensional self-consistent *charge* vector
(shell charges, plus atomic dipoles/quadrupoles for GFN2) with an approximate
inverse Jacobian built from history — exactly what tblite/xtb do by default
(`tblite/scf/mixer/broyden.f90`). This is the architectural reason tblite
converges stiff systems like `complex` where Fock-DIIS diverges.

## Initial guess (`-scf_guess`)

| Guess  | Behaviour                                                                  |
|--------|---------------------------------------------------------------------------|
| `h0`   | Bare Hamiltonian (zero shell charges) at iter 0. **Default.**             |
| `eeq`  | Single-shot dftd4 EEQ atomic charges (`curcuma::dispersion::D4ChargeModel`), partitioned across shells by reference occupations. Starts the SCF in the correct basin. |

## Controlling parameters

All route through the `xtb` config scope (flat CLI flags auto-route there);
defaults equal the historic native defaults, so unset = no change.

| Flag             | Default | Meaning |
|------------------|---------|---------|
| `-scf_mode`      | `broyden` | Strategy (see above). |
| `-scf_guess`     | `h0`    | Initial charge guess (see above). |
| `-scf_damping`   | `0.4`   | Density mixing factor: `P = damp*P_new + (1-damp)*P_old`. Lower = stronger damping. |
| `-scf_threshold` | `1e-6`  | Convergence threshold on max\|dq_shell\| (and dE). |
| `-diis_start`    | `5`     | Damped warmup iterations before DIIS (diis/level-shift). |
| `-diis_subspace` | `6`     | DIIS history depth (Fock matrices kept). |
| `-level_shift`   | `0.2`   | Virtual-orbital shift magnitude (Eh) for `level-shift` mode. |

## The `complex` (231 atoms) — converged

The default (`broyden`) converges it from the bare guess with no flags; all of the
following reach the tblite reference (E approx -329.527 Eh). The historic `diis`
path diverges (eigen solve fails):

```
# DEFAULT — Broyden from the bare guess, no flags (the way tblite/xtb does it):
curcuma -sp complex.xyz -method gfn2                                            # 34 iters

# alternatives (no longer needed, but available):
curcuma -sp complex.xyz -method gfn2 -scf_mode level-shift -scf_damping 0.2     # 56 iters
curcuma -sp complex.xyz -method gfn2 -scf_guess eeq -scf_mode diis -scf_damping 0.1   # 43 iters

# historic path — diverges on complex:
curcuma -sp complex.xyz -method gfn2 -scf_mode diis                             # eigen solve fails
```

`broyden` is the principled fix: it is the same charge-vector quasi-Newton mixer
tblite uses, so it converges stiff systems from the bare-H0 guess at the default
damping — which is why it is now the default. The density/Fock-space modes
(`diis`, `plain`, `level-shift`) need stronger damping than the 0.4 default on hard
systems (use 0.1-0.2) and/or the EEQ guess; `broyden` does not.

## What was tested

- `complex` (231 atoms, C83H122N12O14): the default (`broyden`, no flags)
  converges to -329.52707823 Eh in 34 iters (tblite reference approx -329.53 Eh);
  `level-shift` (56 iters) and the EEQ guess (43 iters) also converge; the historic
  `-scf_mode diis` diverges (eigen solve fails), as documented.
- Consistency: `broyden` reaches the **identical** energy as `-scf_mode diis`
  on H2, H2O, CH4, NH3, C6H6, triose (GFN2) and H2O/NH3 (GFN1) — same fixed point,
  to all printed digits.
- Regression with `broyden` as default: the full native-GFN ctest set passes —
  `xtb_gradient_{H2O,CH4,NH3}`, `d4_dedq`, `xtb_cpscf`, `cli_sqm_08_ngfn1_baseline`,
  `cli_sqm_09_ngfn2_baseline` (7/7). (`cli_sqm_11` fails only because TBLite/`ipea1`
  is not compiled in this build — unrelated.)

## What was NOT tested

- Charged / open-shell systems, transition-metal complexes.
- Interaction of the EEQ guess with `-opt` / MD trajectories over many steps.
- Level-shift / EEQ guess for GFN1 (the code path is shared, but only GFN2 was
  validated on the `complex`).

## Implementation

- SCF loop + mode dispatch + charge pack/unpack: `src/core/energy_calculators/qm_methods/xtb_native.cpp`
- Broyden mixer (`BroydenMixer`, modified Broyden / Johnson 1988): `src/core/energy_calculators/qm_methods/broyden_mixer.h`
- Level-shift helper (`applyLevelShift`): `src/core/energy_calculators/qm_methods/xtb_scf.cpp`
- EEQ guess (`seedEEQGuess`): `xtb_native.cpp`, reusing `dispersion/d4_charge_model.{h,cpp}`
- Config wiring (`applyXtbScfConfig`): `xtb_native.h`, called from `gfn2_method.cpp` / `gfn1_method.cpp`
- Parameter registration: `xtbinterface.h` (`xtb` PARAM scope)

`broyden` mixes the SCC charge vector (`q_sh`, plus `dp_at`/`qp_at` for GFN2,
`q_at` derived from `q_sh`); the other modes mix the density / Fock matrix.
