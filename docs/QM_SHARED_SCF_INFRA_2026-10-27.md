# Native RHF and master's SCF infrastructure (Oct 27, 2026)

> 🤖 AI-generated, ⚙️ machine-tested. One machine (Core Ultra 9 285K, 24 threads, RTX 5090), n = 1-3 per number,
> not reviewed by the operator. Human production testing is pending.

After the merge of master into `feature/qm` (`fde64d16`) the SCF infrastructure that master built for the
semi-empirical methods was checked for what the native RHF (`QMEngine`, `qm_scf.cpp`) can use instead of
duplicating it. CLI and parameter routing of the shared `scf` module: [SCF_MODES.md](SCF_MODES.md),
[PARAMETER_SYSTEM.md](PARAMETER_SYSTEM.md).

## Candidates

| # | master piece | reusable for HF | decision, measurement |
|---|---|---|---|
| 1 | FP32 eigensolver / pseudo-diagonalisation (`-scf_pseudo_diag`, `eigensolveStep`, xtb device loop) | technically yes, not worth it now | device RHF loop, `-gpu cuda`, def2-SVP: eigensolve (X^T F X + syevd) 1032 ms of a 15.4 s loop (6.7 %, 25.3 s wall) for C80H162 (1930 functions, 86 ms per iteration), 2546 ms of 36.4 s (7.0 %, 67.8 s wall) for C120H242 (2890). J/K builds are 89 %. A pseudo step could save at most ~4 % of the loop (~2.5 % of the wall time); port of the kernels into `cuda/qm_gpu_scf.cu` not done, TODO entry |
| 2 | `gpu_device_pool.cpp`, `-gpu_devices`, `-gpu_device` | already used | `EnergyCalculator` puts the leased (or first `-gpu_devices`) device into the method config as `gpu_device`; `QMEngine` passes `m_gpu_device` to every QM device context (J/K, direct, SCF, gradient, RI, MP2). Checked: `-gpu_devices 0` runs on the card; `-gpu_device 1` on the one-GPU machine falls back with one counted `qm_direct` fallback; with `-gpu_strict true` it exits 3. Multi-GPU not testable here (one GPU) |
| 3 | SCF extrapolation (`-scf_extrapolation aspc\|gauss`, `src/core/charge_extrapolation.h`) | yes, implemented | `-qm.scf_extrapolation none\|aspc\|gauss` (same names, same weights as xtb; flat `-scf_extrapolation` routes to `qm` for native QM) extrapolates the converged AO densities of the earlier geometries and projects the prediction onto an idempotent closed-shell density in the new metric. Default `none` (as xtb). Measured below |
| 4 | fragment Broyden preconditioner (`-scf_precond frag`), `-scf_guess fragments` | no | the preconditioner acts on the SCC charge vector of a Broyden mixer; the RHF mixes Fock matrices with DIIS. A fragment guess (block density of separately converged molecules) could help water/urea clusters; not implemented, TODO entry |
| 5 | GPU fallback counting, `-gpu_strict` | already used, one gap closed | every QM device failure path calls `curcuma::reportGpuFallback` (26 call sites in `qm_*.cpp`, 27 with the new one); `lmp2` with `-gpu cuda` ran its CPU-only part silently, now counted as `qm_lmp2` like `qm_ccsd` (`-gpu_strict true` exits 3, checked on (H2O)4) |
| 6 | D3 device tables (`D3ParameterGenerator::buildDeviceTables`, ROCm `d3Dispersion`) for hf-3c | no | the hf-3c D3(BJ)/gCP/SRB terms are O(N^2) host loops of milliseconds next to the O(N^4) HF; the merge kept `createForHF3C` next to master's `buildDeviceTables`, no parameter changed; `ctest -L qm` (qm_hf3c terms vs PySCF + simple-dftd3) passes |

## Density extrapolation (candidate 3)

`QMEngine::buildExtrapolatedGuess` (qm_scf.cpp): P_pred = sum_j w_j P_j (j = 0 the previous geometry), weights from
`curcuma::extrapolation::weights`; natural orbitals of (S X)^T P_pred (S X) with X = S^-1/2, the n_occ largest kept,
P = 2 C C^T. Falls back to the orbital warm start when fewer than two densities exist or the predicted occupations
have no gap of at least 1 (exact: 2 and 0). The history is kept in `finishConvergedSCF` (one entry per geometry,
`m_geom_serial`), cleared by `resetForNewMolecule`.

SCF iterations, HF/def2-SVP methanol, 8 threads (all SCFs tight, since a gradient follows):

| run | SCFs | iterations | mean |
|---|---|---|---|
| MD 0.5 fs, 500 K, 30 fs, `none` | 61 | 615 | 10.08 |
| MD, `aspc` (order 3) | 61 | 562 | 9.21 |
| MD, `gauss` (order 3) | 61 | 568 | 9.31 |
| `-opt` from `methanol.alt1`, `none` | 49 | 491 | 10.02 |
| `-opt`, `gauss` | 49 | 543 | 11.08 |
| `-opt`, `aspc` | 49 | 572 | 11.67 |

ASPC lowers the orbital gradient of the start density by about a factor 10 (4e-2 -> 3e-3..6e-3 per MD step) but saves
only one iteration of ten: from there the DIIS contraction (0.25-0.4 per iteration) to the tight 3e-7 decides. In the
optimisation the steps are too irregular and both predictors cost iterations, so the default stays `none`; use
`aspc` for fixed-step MD. The converged results do not depend on the guess (first MD step identical, later steps
drift apart at the 1e-7 Eh level after 30 steps because the trajectory is chaotic). Not tested: MD energy
conservation over long runs, the device loop with extrapolation (it starts from the same host density), larger
molecules.

## SCF iteration count on the merged base (Task 3)

### Verified stop (`-qm.scf_verified_stop`, default true, energy-only SCFs)

Before, the loop stopped at iteration n when |dE_n|, the orbital gradient g_n of P_n and dP_n were small, then
diagonalised once more and returned P_(n+1) after one more full J/K build. Now the test runs right after a full build
of P_n: |dE_n|, g_n and dP_(n-1) small -> P_n, F(P_n) and the orbitals C_(n-1) it was built from are the result, and
the extra build is skipped (`ScfVerifiedStop`, host and device loop; the device loop downloads J/K with the new
`QMGpuSCF::downloadJK`). In integral-direct runs the build on which convergence is predicted (g_(n-1)^2 / g_(n-2) and
dP_(n-1) below the thresholds) is made a full one. The reported iteration count is unchanged; one J/K build less.

Tight SCFs (analytic gradients, mp2/ri-mp2/lmp2/ccsd) keep the old rule: the old rule returns a density one DIIS step
(g about 4x smaller) better than the tested one. With the verified stop the methanol HF gradient deviated 5.3e-8
Eh/A from the tight reference instead of 4.7e-9, so it was restricted to energy-only runs, where the energy is second
order in g. Checked after the restriction: methanol, anthracene (HF) and water (HF-3c) gradients identical to the
previous binary (4.7e-9, 6.1e-9, 8e-14 Eh/A from the tight reference `-qm.scf_convergence density -qm.scf_threshold 1e-10`).

Wall time, `-sp -method hf`, def2-SVP, default flags; before = the binary of the merge `fde64d16` (the routing change before this one left the SCF unchanged), min of n runs:

| system | nbf | threads | n | before (s) | after (s) | iterations | E_elec after - before (Eh) |
|---|---|---|---|---|---|---|---|
| anthracene | 246 | 24 | 3 | 2.03 | 1.94 | 12 / 12 | +1e-10 |
| C20H42 | 490 | 24 | 3 | 7.76 | 7.30 | 11 / 11 | 0 |
| C40H82 | 970 | 24 | 2 | 36.14 | 34.01 | 11 / 11 | 0 |
| anthracene | 246 | 1 | 1 | 33.19 | 31.10 | 12 / 12 | +1e-10 |
| C20H42 | 490 | 1 | 1 | 107.44 | 100.67 | 11 / 11 | 0 |
| C80H162, device loop (`-gpu cuda`) | 1930 | 24 | 1-2 | 25.3 | 24.4 | 12 / 12 | +2e-10 |

All energies are within 1e-10 Eh of the tight references of [QM_SCF_CONVERGENCE_2026-10.md](QM_SCF_CONVERGENCE_2026-10.md)
(anthracene -1307.2299507207, C20H42 -2163.0355541443, C40H82 -5029.1769178829, C80H162 -11508.2026479164 Eh).
The GPU "before" is the run of the same day before the change (the copied reference binary did not find the CUDA
plugin, its `-gpu cuda` runs fell back to the CPU and are not used). 5-6 % on the CPU, 3-4 % on the GPU.

### Correlated and tight runs (unchanged)

| case | threads | iterations | wall before / after (s) |
|---|---|---|---|
| ri-mp2 C20H42 | 24 | 13 / 13 | 11.77 / 11.86 |
| ccsd (H2O)4 | 24 | 14 / 14 | 8.41 / 8.34 |
| ccsd (H2O)6 | 24 | 14 / 14 | 56.64 / 56.64 |
| lmp2 C40H82 | 24 | 14 / 14 | 53.66 / 53.79 |

No change was found that shortens the tight phase without loosening it: the 2-3 extra iterations take the orbital
gradient from 1e-5 to 3e-7 at the DIIS contraction of 0.25-0.4 per iteration. Returning the tested density instead
of the next one (verified stop) makes the gradient 10x worse (above). Not tried: a second-order (Newton/CPHF) final
step, which needs an orbital-Hessian solve costing several J/K builds.

### Other observations

- Device loop C80H162: FP32 phase ends at dP 3.6e-4; the FP32 SCF energy is 1.1e-4 Eh off the FP64 one, so the first
  FP64 iteration jumps to dP 6e-3 and the FP64 phase needs 4 builds (9.9 s of 18 s incl. the final build).
  `-qm.gpu_fp32_switch 1e-4` / `3e-4`: 14 / 13 iterations, 27.7 / 26.3 s (worse); default unchanged.
- Setup outside the loop at 1930 functions with `-gpu cuda`: about 2.9 s to the end of the 1e integrals, 2.5 s
  orthonormaliser and guess, 1.5 s direct-J/K setup (host tables 1.3 s); not optimised here.
- Warm starts (MD, -opt) run tight SCFs, so they do not get the verified stop; `-scf_extrapolation aspc` saves about
  one iteration per MD step (above).
