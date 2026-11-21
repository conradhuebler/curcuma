# Spin polarisation of native GFN1/GFN2 (`-xtb.spin_polarization`)

Status: 🤖 AI-implemented, ⚙️ machine-tested (Nov 21, 2026). Not ✅ TESTED: human production testing is pending until the operator removes this line.

## What it is

`-spin N` alone (see [SCF_MODES.md](SCF_MODES.md), "Open shell") fills one orbital set with two Fermi fillings (n_alpha, n_beta). That is the
tblite model *without* spin polarisation: the triplet O2 lies 1.65 kcal/mol above the singlet. With `-xtb.spin_polarization true` (flat
`-spin_polarization true` is routed to the module `xtb`, as is `-xtb.spin_polarization_scale`) the SCF runs the tblite `spin-polarization`
container with two Fock matrices and the spin constants W (tblite `src/tblite/spin.f90`, nspin = 2):

- Shell magnetisation `m_l = sum_{mu in l, nu} (Pa - Pb)_{mu nu} S_{nu mu}` (alpha minus beta Mulliken population of shell l).
- `E_spin = 1/2 sum_A sum_{l,l' on A} W^A_{ll'} m_l m_l'` (all shell pairs of one atom; W < 0), added to the electronic energy.
- `u_l = sum_{l'} W_{ll'} m_{l'}`; `F^sigma = F + sigma * 1/2 S_{mu nu} (u_mu + u_nu)`, sigma = +1 (alpha), -1 (beta); `F` is the charge Fock matrix.
- Both spins have their own orbitals and Fermi filling (`fermiFillingChannel`); `P = Pa + Pb` feeds every other term unchanged. `m_sh` is
  appended to the Broyden vector (`[q_sh | dp | qp | m_sh]`) and to the warm-start state.
- Gradient: W has no geometry dependence. Two changes only (`xtb_gradient.cpp`): the energy-weighted density sums both orbital sets, and the
  overlap-gradient coefficient gets `+ Pm_{mu nu} (u_a + u_b)` with `Pm = Pa - Pb` (opposite in sign to the charge term `- P (v_a + v_b)`).
- Acts only with `-spin N > 0`; with `-spin 0` nothing changes. Without the option all results are bit-identical to before.
- `spin_polarization_scale` multiplies all W (tblite `wscale`); calibration only.

Code: `XTB::solveEigenUhf`, `spinPotential`, `energySpin`, `buildSpinFocks` (`xtb_scf.cpp`); constants `parameters/spin_constants.hpp`.

## Constants: origin and open points

86 elements x 6 values (ss sp pp sd pd dd), copied value for value from tblite **v0.8.0** (`src/tblite/data/spin.f90`), indexed by angular momentum, not
by shell position (transition metals are ordered d,s,p). 85 of the 86 rows differ from the tblite commit pinned in `CMakeLists.txt` (9ca469a, Dec 2024).
- **TODO operator: licence.** tblite is LGPL-3.0-or-later; whether the table may be copied into this GPL-3.0 tree is not decided.
- **TODO operator: literature reference.** The tblite source gives none (only the PR title "Updated spGFN2 spin constants (#155)"); not looked up, not guessed.
- The tblite pin in `CMakeLists.txt` still gives the old constants for `USE_TBLITE` runs with spin polarisation; raise it before comparing with `tblite-gfn2`.

## Not supported (hard error, no silent fallback)

`-gpu cuda|rocm|vulkan`, `-large_system_mode fragments|dc|sparse`, `-scf_mode` other than `broyden`, `-scf_precond frag`, `-eigensolver` other than `mkl`,
`-scf_extrapolation_apply xlbomd`. The analytic response (`xtb_response.cpp`) already rejects open shells. `-spin` stays global (no spin per fragment).
Not implemented: UHF without the spin constants, spin contamination analysis, orbital output for the beta set (`C`/`eps` hold the alpha set).

## What was tested (Nov 21, 2026, release build, CPU, `-threads 1`)

Reference: tblite 0.8.0 Python wheel, `add("spin-polarization", 1.0)`, same electronic temperature.

- **Energies** against tblite, GFN2 and GFN1, 300 K and 2000 K: H, O, OH (uhf 1, 3), H2 (0, 2), O2 (0, 2, 4), Pt7 (0, 2, 4, 6, 8): all differences <= 2.8e-8 Eh
  (most < 1e-8; the 2000 K Pt7 values 1e-8 to 2.8e-8 at the default `-scf_threshold`). O2 triplet GFN2: -7.931969836872 vs -7.931969836906 Eh.
- **Gradients** against tblite (Eh/Bohr): O2, O2 tilted, OH, H2, Pt7 (uhf 4, 6, GFN1 and GFN2), Pt7 + 2 H2 + O2 (uhf 2, 4, 2000 K): max difference 4e-8, mostly < 1e-8.
- **FD gate** (`-dump_gradient` against central differences of `-sp`, Angstrom steps, `-scf_threshold 1e-10`): O2 (uhf 2, 4), O2 GFN1, OH, H2 (uhf 2), Pt7 uhf 4, 6, GFN1 uhf 4: error scales
  with h^2 (O2: 2.1e-6 at h = 1e-3, 5.3e-7 at 5e-4); smallest 5e-8 Eh/A (Pt7, h = 5e-4). Pt7 uhf 2 at 300 K is bistable: displaced geometries land on another SCF solution (3.5 mEh apart); only
  the 12 of 21 components that stay on one branch are comparable (6e-8 Eh/A); at 2000 K all 21 are (2.4e-7).
- **Regression**: without the option 11 cases (water, O2 spin 0 and 2, OH, H2, Pt7, Pt7 + 2 H2 + O2; GFN1 and GFN2) give bit-identical energies and `-dump_gradient` files against the binary of `267dc325`;
  `ctest -R "cli_sqm_|xtb_gradient|gradient_unit_contract|sqm_"`: 343 of 343; `ctest -L structinv`: 168 of 172, the 4 failures (`reactions_doebereiner_control_gfn2`) fail identically with the old binary (SCF not converged).
- **ctest** `cli_sqm_16_spin_polarization`: triplet O2 energy and gradient (GFN1, GFN2) against tblite, option off, `-spin 0`, three excluded modes.
- **Pt7 + 2 H2 + O2** (structure `reactions-doebereiner-pt7`, GFN2, 2000 K): see the spin ladder below.

## Spin ladder of the Doebereiner Pt7 + 2 H2 + O2 start geometry (GFN2, `-sp`, spin_polarization true)

| uhf | native 2000 K | tblite 2000 K | native 300 K | tblite 300 K |
|---|---|---|---|---|
| 0 | -41.94941114 | -41.94941116 | not converged (500 it) | not converged |
| 2 | -41.95057830 | -41.95057830 | -41.88990742 | not converged |
| 4 | -41.95253772 | -41.95253773 | -41.90091402 | not converged |
| 6 | -41.96038415 | -41.96038416 | -41.90283926 | -41.89334208 (other solution) |
| 8 | -41.96143396 | -41.96143398 | -41.90507124 | -41.90507123 |
| 10 | -41.95432148 | -41.95432150 | -41.90141416 | -41.90141416 |

The minimum is uhf 8 at both temperatures (not 4 or 6), 0.66 kcal/mol below uhf 6 at 2000 K. Short MD (cold Pt7 input, spin 8, 1000 steps, 1 thread, machine under foreign load 25): 38.3 s with the option
against 30.1 s without (26 against 33 steps/s, not half), 20.4 against 20.0 SCF iterations per step (200 steps), no SCF error.

## Known limits

- The SCF has several solutions (different magnetisation patterns). The result depends on guess and damping: GFN1 full system uhf 4, 2000 K: `-scf_guess eeq` (default) -42.50317, `-scf_guess h0` -42.51162 (tblite -42.51162).
  Multi-solution cases were not resolved; compare energies only where both programs converge to the same solution. At 300 K the 13-atom system often needs more than 500 Broyden cycles (also in tblite).
- Not tested: MD beyond 1000 steps, `-opt` of open-shell metal clusters (only a smoke run on O2), solvation, `-hessian`, periodic or charged systems, elements other than H, O, Pt, GFN1 d-shell systems other than Pt.
- GPU: not implemented, not tested.
