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
