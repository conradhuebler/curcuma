# What runs on the ROCm device and what stays on the host (Oct 17, 2026)

AI-generated from a code reading by four read-only agents; three claims were re-checked by the main session (multipole
potential on the host, no `supports*` override on the ROCm side, no solvation code in the GFN-FF wrapper), the rest is
taken from the agent reports and not checked line by line. Nothing was run for this list; all timings come from
`ROCM_GFN_POTENTIAL_2026-10-14.md` (peo201-chain, 1410 atoms, 1 host thread, n = 1) and `GPU_STEP_PAUSES_2026-10-14.md`.
Paths: `qm_methods/` = `src/core/energy_calculators/qm_methods/`, H = `qm_methods/rocm/xtb_hip_context.hip`.

## GFN1 / GFN2 (`-gpu rocm`)

The ROCm backend does not set `supportsResidentLoop`, `supportsDevicePotential`, `supportsOnTheFlyMultipole` or
`supportsDeviceSolvation` (only the CUDA class in `xtb_gpu_method.cpp:400,500,511,523` does). The SCF therefore runs as a
host loop (`xtb_native.cpp:1493-1600`) with a device eigensolver. Older documents that call the ROCm GFN path
"fully device-resident" refer to the matrices, not to the loop.

| Part | Where | Note |
|---|---|---|
| Basis tables, S, H0, gamma (`k_gamma`), dipole/quadrupole integrals, L = chol(S) | device, then downloaded once per geometry | 3 nao^2 + nsh^2 (+ 9 nao^2 for GFN2) cross the bus (`xtb_native.cpp:578-660`) |
| Hamiltonian CN and self-energy | host and device (computed twice) | host copy feeds D4 and the gradient |
| Fock, eigensolve, density, Mulliken partial sums, multipole moments | device | per iteration only length-nao/nat vectors cross the bus |
| Coulomb ES2 potential (gemv on the downloaded gamma), third order, Fermi occupations, population to shell/atom charges, Broyden, convergence test | host | O(nsh^2) to O(nao), 1 thread |
| GFN2 multipole potential and energy (18 nat^2 interaction matrices held on the host) | host | 29 of 44 ms host time per iteration in the doc |
| D4 potential in the SCF (GFN2) | host weights (`buildRefWFlat`), device pair contraction `k_d4_dedq` | 2*7*nat doubles up per iteration |
| GFN1 dispersion (D3) and halogen bond, energy and gradient | host | no device kernel |
| After the SCF: P and C download (2 nao^2), host Fock `m_F`, band energy | host | also when the gradient runs on the device |
| Gradient: repulsion, H0/Pulay with multipole integrals, CN on-site, ES2 | device | `H:279-461` |
| Gradient: D4 pair terms (`k_d4_grad`), ATM (`k_d4_atm`, FP64) | device, weights and C6 tables built on the host | ATM 1.2 s of the 15 s run (doc) |
| Gradient: GFN2 multipole interaction, Hamiltonian CN chain rule, D4 covalent-CN fold, solvation | host | O(nat^2), threaded |
| Bond orders | not computed (`getBondOrders()` returns an empty vector) | |

CUDA has device versions of most host rows (`cuda/xtb_gpu_context.cu`): `k_multipole_potential` (:2594, `_otf` :2664),
`k_d4_build_refw` (:2440), `k_occupations` (:2955), `runBroydenUpdate` (:5899), the whole step `residentScfStep`
(:6083-6185).

## GFN-FF (`-gpu rocm`)

All energy terms run as kernels (bond, angle, torsion, inversion, repulsion, D4 dispersion with ATM, HB/XB, Coulomb, EEQ
with projected PCG). The wrapper never calls `GFNFF::Calculation()`. Host work: topology and parameters (setup, again after a
displacement over 0.5 Bohr), list rebuilds (repulsion list about 10.5 of 15.1 ms, D4 CN and weights 16 of 47 ms, XB
detection and the bond-HB cross reference), the MD integrator and thermostat (3N doubles per step each way), the host
Schur step of the dense Cholesky EEQ path.

- **GPU fallbacks to the CPU method:** net charge on several fragments with `frag_charge_model ensemble`; **implicit solvent
  (new, Oct 17, 2026)**.
- **Defect found while writing this list:** `-gfnff.solvent water -gpu rocm` returned the gas-phase energy without a message
  (caffeine, one run each: ROCm -4.67273706 Eh, CPU with solvent -4.68976245 Eh), because the ALPB term exists only in
  `GFNFF::Calculation()`. The wrapper now requests the CPU method (counted fallback, `-gpu_strict true` exits with 3). The CUDA
  wrapper shares `gfnff_gpu_method_impl.h`; its behaviour was not run here.
- **PBC:** the periodic kernels exist but do not change the energy (`TODO.md`).
- **Corrected statements:** `TODO.md` and `ROCM_EVALUATION_2026-10-13.md` listed device ppcg, localinv, K1, far field and
  staged upload as missing; the code has them (`eeq_solver.cpp:939-941`, `gfnff_hip_method.h:56-74`, `H:783`).

## Candidates to move (agent estimate, not measured)

1. GFN2 multipole potential and energy (medium): largest single host item; CUDA kernels are portable one to one, the 18 nat^2
   matrices need an upload (about 286 MB at nat = 1410) or the on-the-fly variant.
2. D4 reference weights on the device (medium); needs device charges.
3. Coulomb gemv, third order, charge scatter, energies (small to medium).
4. Fermi occupations (small), Broyden and convergence (medium); worth it only together with 1 to 3.
5. Complete `residentScfStep` for ROCm (large); leaves only dq and a few energies on the bus.

GFN1 gains only from 3 to 5 (no multipole, no D4). Each step needs a correctness check against the CPU.
