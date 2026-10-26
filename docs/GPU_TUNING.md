# GPU and Large-System Tuning Options

> 🤖 AI-generated options and measurements, machine-tested only. Human production testing pending.

Reference for manual tweaking of multi-GPU runs and large systems (thousands of atoms). Every option below lists its default, what it changes, when to touch it, and what was measured. Measurements: 4x RTX A4500 (20 GB, Ampere sm_86, PCIe), 36 cores, CUDA 13.3, Sep 2026, unless noted. Background and history: [MULTI_GPU.md](MULTI_GPU.md).

CLI flags are global keys (`-gpu_device 2`) unless written with a scope (`-gfnff.gpu_coulomb_implicit false`). Environment variables are read at start-up.

## 1. Device selection and batch distribution

| Option | Default | Effect |
|---|---|---|
| `-gpu cuda\|rocm\|vulkan\|cpu` | `none` (CPU) | Backend. `auto` was removed Oct 26, 2026 (warning, counted fallback, CPU). |
| `-gpu_device N` | -1 (runtime default, device 0) | Pins a single calculation to device N. Indices are after `CUDA_VISIBLE_DEVICES` / `ROCR_VISIBLE_DEVICES` (SLURM sets these). Invalid index: error + CPU fallback. |
| `-gpu_devices 0,2` | all visible | Device pool for batch workers; since Sep 28, 2026 also the device set of a single calculation (it runs on the first listed device and splits only over the listed ones). |
| `-gpu_workers_per_device N` | 1 | Batch workers sharing one device. |
| `-gpu_strict true\|false` | false | Every GPU degradation (CPU fallback of a worker, plugin missing or declined, device init failure, slower eigensolver library or split only partly used, GFN-FF EEQ device-solve fallback, split flags ignored in a batch worker) ends the run with **exit code 3**. Without it they are counted and summarised after the run as warnings, whatever the verbosity. **When to set it**: timing runs and production jobs on unfamiliar hardware, where a silent CPU fallback costs hours before anyone notices. Costs nothing; a clean run is unaffected (checked Sep 28, 2026). |

- **Batch workers** (`-sp`/`-opt` on multi-XYZ, ConfSearch, Hessian) lease a device slot each; surplus workers wait. The pool is active when more than one device is visible, or when `-gpu_devices`/`-gpu_workers_per_device` is given.
- **Host threads per GPU worker** = cores / GPU slots (CPU batch workers stay single-threaded).
- **When to change `-gpu_workers_per_device`**: small/medium molecules where one context does not saturate the card. Measured, 8x polymer (1410 atoms) GFN2: 4 GPUs 28.5 s, 8 slots (2 per GPU) 25 s. Watch device memory: each worker holds its own context.
- **When GPU batching does not pay**: small molecules. 40 conformers of 90 atoms: 4 GPUs 2.0 s vs CPU with 36 workers 0.58 s.
- `curcuma -methods` lists the devices a plugin sees.

## 2. Native GFN1/GFN2 on the GPU (CUDA)

| Option | Default | Effect / when to change |
|---|---|---|
| `-gpu_memory_check true\|false` | true | Refuses a basis whose estimated device memory exceeds free memory (warning + CPU run). Set false only if the estimate is too pessimistic for a card that just fits. |
| `-gpu_sparse_integrals auto\|on\|off` | auto | Screened pair storage for S, H0 and the GFN2 multipole integrals. `auto` = when fewer than 50 % of AO pairs survive the distance screen (integrals < 1e-20 dropped). `on` forces it (used for validation: identical energies on small molecules), `off` forces dense nao^2 storage. polymer_2x: 6 % of pairs, device memory ~50 -> 17 GB. |
| `-scf_fp32_threshold X` | 1e-5 on the GPU, 1e-3 on the CPU | Mixed precision: FP32 eigensolves until max\|dq\| < X, then FP64. Larger X = more FP64 iterations (slower on FP64-weak consumer cards, closer CPU gradients). polymer (1410 atoms): 1e-3 -> 4 FP64 solves, SCF 6.7 s, gradient vs CPU 2.6e-6; 1e-5 -> 1 FP64 solve, 4.5 s, 7.0e-6 Eh/A. Below ~3e-6 the FP32 noise floor adds iterations. |
| `-scf_mixed_precision true\|false` | **false on a device with full-rate FP64** (A100, H100/H200/GH200, B200, V100 - compute capability 6.0/7.0/8.0/9.x/10.x), true elsewhere | Mixed precision is a consumer-GPU trick: it pays where FP64 runs at 1/32 to 1/64 of FP32 (A4500: FP64 iteration 38 s vs FP32 9.4 s) and hurts where FP64 is half of FP32. Measured on an H200 NVL, polymer_2x: **FP32 iteration 5.29 s, FP64 3.36 s** - and the FP32 phase converged to a fixed point **1.1 kcal/mol off** (it reported max\|dq\| 7.1e-6 while the true residual was 9.4e-3), which the SCF only discovered through its FP64 check, at the cost of 42 iterations instead of 15. |
| `-scf_pseudo_diag true\|false` | **true since Oct 15, 2026** (operator decision; only acts in FP32 SCF steps, i.e. on GPUs with mixed precision on) | FP32 SCF steps rotate the previous eigenvectors once between the occupied and virtual blocks and re-orthonormalise them instead of diagonalising (only GEMMs; docs/GFN2_GPU_COST_PLAN.md stage 2). Used after the previous max\|dq\| < `-scf_pseudo_diag_max_dq` (0.05), at most `-scf_pseudo_diag_max_steps` (8) in a row, closed shell, HOMO-LUMO estimate > 30 kT, all rotation angles < 0.1 - a full solve otherwise. Measured (4x A4500): polymer_2x 113 -> **89 s**, 12 iterations either way, energy identical to 8 decimals; polymer on 1 GPU SCF 11.8 -> 10.6 s, gradient 2.8e-9 Eh/A from the full-solve run; MOR41 on the GPU 95/95 and GMTKN55 2459/2460 identical to 8 decimals (one 1e-8 Eh). Only FP32 steps use it, so it does nothing where `-scf_mixed_precision` is off (the default on full-rate-FP64 GPUs such as the H200). The step runs on the calculation's own GPU even when the eigensolve is split. |
| `-scf_pseudo_diag_fp64 true\|false` | false (EXPERIMENTAL) | With `-scf_pseudo_diag`: also in FP64 steps, in FP64. The gradient's energy-weighted density is then built as C_o (2 G_oo) C_o^T from a rebuilt Fock matrix, because the rotated vectors do not diagonalise the occupied block (without that: 8.6e-4 Eh/A error on polymer; with it 3.6e-9). **Slower on FP64-weak GPUs**: polymer_2x on 4x A4500 117 -> 167 s (the FP64 step on one card 25 s for F*C alone, the W correction 29 s). Meant for full-rate-FP64 GPUs (H200: mixed precision off, so without this flag no step uses the pseudo path) - not measured there. |
| `-scf_threshold 1e-8` | 1e-5 | Tight SCF. Needed for GPU gradients that must agree with the CPU beyond ~1e-4 Eh/A (the loose default lands the device Broyden at a different point inside the tolerance). |
| `-xtb.rocm_device_potential true\|false` (ROCm, GFN2) | **true** since Oct 24, 2026 | The host SCF loop sends the mixed charges, moments and D4 reference weights; the device builds the SCF potential (gamma q_sh, third order, multipole with device-built interaction matrices, D4 dE/dq), the Fock matrix and the eigensolve. A failed setup or step moves the run to the host potential build with one warning (`CURCUMA_ROCM_FAIL_POT=begin\|solve` forces it). gfx1200, 14 small structures and peo201-chain: dE vs CPU <= 8.2e-12 Eh, dG <= 2.7e-6 Eh/A, on vs off <= 5.7e-14 Eh; peo201-chain `-sp` 14.33 -> 14.18 s (n = 3), potential 10.9 -> <1 ms per iteration. With it on, `-xtb.rocm_device_multipole` has no effect. |
| `-xtb.rocm_resident_loop true\|false` (ROCm, GFN1 and GFN2, closed and open shell) | **true** since Oct 24, 2026 | GFN2 needs `rocm_device_potential`. The whole SCF iteration runs on the device (D4 weights, potential, Fock, eigensolve, Fermi occupations, density, charges and moments, SCC energies, max\|dq\|, Broyden); the host receives one integer and five scalars per step. GFN1 (since Oct 24, 2026, also off with `-scf_gpu_gfn1_resident false`): SCC vector = shell charges, v_sh = gamma q_sh, v_at = Gamma_A q_A^2; gfx1200, caffeine, fi-pyridine, 12 MOR41, peo201-chain: dE vs CPU <= 1.1e-10 Eh (peo201-chain; others <= 3.4e-12), dG <= 2.6e-6 Eh/A, on vs off <= 1.1e-10 Eh, SCF iterations equal to off; peo201-chain `-sp` 18.64 -> 18.51 s, PR24 (120 atoms) `-md` 100 fs 37.16 -> 37.66 s (median of 3; per SCF iteration 35.0 -> 37.4 ms at 120 atoms). Open shell follows `-xtb.rocm_resident_open_shell`, `-scf_precond frag` follows `-xtb.rocm_resident_precond`, an implicit solvent follows `-xtb.rocm_device_solvation`; `CURCUMA_SCF_PRECOND_FILE` and `-keep_diis true` keep the host-driven loop. A failed setup takes the host-driven loop, a failed step restarts that SCF on it; one warning (`CURCUMA_ROCM_FAIL_RES=begin\|step` forces it). gfx1200, caffeine, 12 MOR41, peo201-chain: dE vs CPU <= 5.5e-12 Eh, dG <= 2.7e-6 Eh/A, on vs off <= 1.0e-11 Eh, SCF iterations equal to off; peo201-chain `-sp` 14.18 -> 14.06 s, PR24 `-md` 100 fs 54.50 -> 54.01 s (median of 3). |
| `-xtb.rocm_resident_open_shell true\|false` (ROCm, GFN1 and GFN2, open shell) | **true** since Oct 24, 2026 | Needs `rocm_resident_loop`. An open shell (n_alpha != n_beta, `-spin`) stays on the device-resident loop: the occupation step fills the same eigenvectors with n_alpha and n_beta electrons (aufbau, then a Fermi level per spin channel) and sums them, as the host loop. gfx1200, caffeine +1, 16 open-shell GMTKN55 structures (doublets to quartets, charged and neutral, Pb, Te), ED01 triplet (Cr), PR07 +1 (W), PR24 +1: dE vs CPU <= 8.7e-11 Eh and dG <= 8.9e-6 Eh/A (ED01 triplet, GFN1, T = 0; the same with the loop off), on vs off <= 2.0e-13 Eh, SCF iterations equal to off; also at `-xtb.electronic_temperature` 0, 300, 1000, 5000. false keeps open shells on the host-driven loop. |
| `-xtb.rocm_resident_precond true\|false` (ROCm, GFN1 and GFN2, `-scf_precond frag`) | **true** since Oct 24, 2026 | Needs `rocm_resident_loop`. The device Broyden of the resident loop applies the host-built M (M*F for alpha*F, M*dF for alpha*dF) and the `-scf_precond_guard` rule on \|\|F\|\|_2, as the CPU loop ([SCF_PRECOND.md](SCF_PRECOND.md)). gfx1200, 8 systems (urea/water clusters, segmented PEO, caffeine, AcCl, Bz, PR07) plus peo201-chain, GFN1 and GFN2: dE vs CPU with the same switch <= 1.0e-10 Eh (peo201-chain GFN1; others <= 1.1e-11), dG <= 5.2e-6 Eh/A, SCF iterations equal to the CPU except PR07 GFN2 (8 against 10), guard hits at the same iteration; peo201-chain `-sp` (segmented, `-scf_precond_segment_atoms 3 -scf_precond_max_fragment_atoms 20`) GFN2 14.50 -> 11.21 s, GFN1 19.16 -> 15.26 s against frag off (median of 3). false keeps `-scf_precond` other than none on the host-driven loop, where it is ignored with a warning. |
| `-xtb.rocm_device_solvation true\|false` (ROCm, GFN1 and GFN2, native ALPB/GBSA) | **true** since Oct 24, 2026 | The Born matrix B (nat^2 doubles, 16 MB at 1410 atoms, limit 2 GB) goes up once per geometry and the device SCF potential adds v_at += B q_at (GFN2: `rocm_device_potential` step and resident loop; GFN1: resident loop, plus the CM5 shift B cm5), so solvated runs keep the resident loop. Solvation energy and nuclear gradient stay on the host (post-SCF, converged charges). A failed upload or step returns the solvation to the host-driven loop with one warning (`CURCUMA_ROCM_FAIL_SOLV=begin\|step` forces it, `CURCUMA_ROCM_SOLV_MAX_MB` lowers the limit). gfx1200, caffeine, caffeine +1, acetate -1, 6 MOR41 x {water ALPB/GBSA, methanol ALPB, toluene GBSA, hexane ALPB} x GFN1/GFN2 (90 pairs) plus peo201-chain in water: dE vs CPU <= 7.3e-11 Eh (peo201-chain GFN1; small systems <= 2.7e-11), dG <= 4.8e-6 Eh/A, on vs off <= 6.7e-11 Eh, SCF iterations equal to off; peo201-chain water `-sp` GFN2 17.39 -> 16.96 s, GFN1 16.55 -> 16.42 s (median of 3). false keeps solvated runs on the host-driven loop. |
| `-xtb.rocm_d4_postscf_device`, `-xtb.rocm_d3_device`, `-xtb.rocm_cn_device` `true\|false` (ROCm; GFN2 D4, GFN1 D3(BJ), both CN) | **true** since Oct 24, 2026 | Post-SCF dispersion on the device in FP64: GFN2 D4 weights at the converged charges, 2-body, ATM with its q = 0 C6 matrices (2 nat^2 doubles) and neighbour list, dE/dCN fold; GFN1 D3 CN, weights, pair terms and CN chain rule; the Hamiltonian CN is taken from the device integral build and the D4 covalent CN computed on the device. The D4 path needs the reference tables of the resident loop (`rocm_resident_loop`); after a host-driven SCF the host code runs. A device failure returns that path to the host with one warning (`CURCUMA_ROCM_FAIL_POSTSCF=d4\|d3\|cn` forces it). gfx1200, 64 small structure/method pairs (caffeine, caffeine +1, fi-pyridine, 13 MOR41, 16 open-shell GMTKN55) plus peo201-chain, gas phase and water ALPB: on vs off <= 1.4e-12 Eh, dG <= 1.6e-14 Eh/A (D4 or D3 alone: dE 0, dG <= 6.2e-16), SCF iterations equal; peo201-chain `-sp` GFN2 14.33 -> 13.41 s (ATM kernel 1.22 -> 0.36 s), GFN1 19.33 -> 18.64 s; PR24 `-md` 100 fs GFN2 54.15 -> 51.29 s, GFN1 37.62 -> 37.67 s (median of 3). AI-generated, machine-tested only ([ROCM_RESIDENT_LOOP_WP.md](ROCM_RESIDENT_LOOP_WP.md) WP8a). |
| `-xtb.rocm_grad_resident`, `-xtb.rocm_grad_tail_device`, `-xtb.rocm_defer_host_integrals` `true\|false` (ROCm, GFN1 and GFN2) | **true** since Oct 25, 2026 | After the resident loop: the gradient potentials are built on the device at the converged charges (gas phase, GFN2 with `d4_charge_source mulliken`); the gradient tail (GFN2 multipole interaction, CN chain rule) and the repulsion energy run on the device in FP64; the host keeps lazy copies (no multipole-integral download, no post-SCF host potential rebuild, S/H0 downloaded on demand). A device failure returns that path to the host code with one warning (`CURCUMA_ROCM_FAIL_GRAD=resident\|tail\|defer` forces it). gfx1200, 64 small structure/method pairs plus peo201-chain and peo201-water500: on vs off <= 3.6e-12 Eh, dG <= 1.2e-15 Eh/A, SCF iterations equal; peak host RSS peo201-chain GFN2 2172 -> 914 MB, peo201-water500 GFN2 6700 -> 681 MB; peo201-chain `-sp` with gradient GFN2 14.04 -> 12.86 s, GFN1 19.25 -> 17.87 s, PR24 `-md` 100 fs GFN2 51.27 -> 49.65 s (median of 3). AI-generated, machine-tested only ([ROCM_RESIDENT_LOOP_WP.md](ROCM_RESIDENT_LOOP_WP.md) WP8b). |
| `-scf_pseudo_diag true\|false` on ROCm (GFN1 and GFN2, resident loop) | **false on ROCm** (Oct 24, 2026; the PARAM default true is the CUDA one) | Port of the CUDA pseudo step (same gates, `-scf_pseudo_diag_max_dq`, `_max_steps`, `_fp64`), closed shell only. peo201-chain `-sp` GFN2 14.11 -> 12.35 s, GFN1 18.60 -> 16.90 s (median of 3), 5 of 12 / 3 of 10 iterations pseudo; energy within 2.7e-12 Eh, gradient 1.6e-7 Eh/A of the full solve. Off by default because with `-scf_precond frag` caffeine needs 6 instead of 4 iterations and PR24 MD gains nothing; details in [ROCM_RESIDENT_LOOP_WP.md](ROCM_RESIDENT_LOOP_WP.md) WP7 (AI-generated, machine-tested). |

### Multi-GPU eigensolve (one large molecule on several GPUs)

The per-iteration full-spectrum eigensolve of the device-resident SCF can be spread over several GPUs; integrals, Fock build, mixing and most of the gradient stay on the calculation's own device (`-gpu_device`, default 0; with `-gpu_devices` the first listed device). **On by default** above `gpu_eigensolver_min_nao` (4000) basis functions when more than one device is usable - the text said "off by default" until Sep 28, 2026, contradicting the table and the code (`xtb_gpu_method.cpp`). **Since Sep 2026 the gradient's energy-weighted density W also goes through the distributed density path** (`-gpu_density_devices`), because it is the same SDDMM as the density itself: 3309 -> 855 ms on 4 GPUs at nao 15444.

| Option | Default | Effect / when to change |
|---|---|---|
| `-gpu_eigensolver_devices all\|0,1,2,3\|none` | **all usable devices** when more than one is usable, the calculation is not pinned with an explicit `-gpu_device`, and it does not run inside a batch worker; else off | Devices for the eigensolve (at least two); may include the calculation's own device. `all` means the `-gpu_devices` list when one is given, else every visible device (since Sep 28, 2026; before, `-gpu_devices` was ignored here). `none` switches it off. Batch workers never split, **also not with an explicit flag** (since Sep 28, 2026: ignored with a warning and counted as a GPU fallback; before, `all` made every worker use every GPU). |
| `-gpu_eigensolver_backend auto\|mp\|mg` | auto | `mp` = cuSOLVERMp over NCCL, `mg` = cusolverMg (CUDA toolkit, deprecated in CUDA 13). `auto` tries mp, then mg. |
| `-gpu_eigensolver_min_nao N` | 4000 | Below N basis functions the single-GPU solver is used: the column scatter/gather and NCCL collectives cost more than they save for small matrices. This gate is what keeps the default harmless for everyday molecules. |
| `-gpu_eigensolver_block N` | 128 | Column block size of the 1 x ndev distribution. Standalone n = 15444, 4 GPUs: nb 64/128/256 -> FP64 21.9/20.7/20.9 s, FP32 5.50/5.24/5.36 s. |
| `-gpu_eigensolver_fp32 true\|false` | true | Also distribute the FP32 iterations of the mixed-precision SCF (`mp` only; `mg` never gets FP32 solves, see below). Set false to keep FP32 on the calculation's own GPU. |

- On a failure the SCF continues with the single-GPU solver; at verbosity >= 1 a warning names the reason, verbosity 2 reports backend, GPUs and number of solves after each SCF.
- The solvers live in `libcurcuma_cuda_mgpu.so` next to `libcurcuma_cuda.so` and are loaded on first use, so a host without cuSOLVERMp/NCCL still runs all single-GPU calculations.
- What is distributed: FP64 iterations run the whole generalized solve on the GPUs (reduction `sygst`, `syevd`, back-transform; the Cholesky factor is distributed once per geometry). FP32 iterations distribute `syevd` and the back-transform: the cuSOLVERMp FP32 reduction is not faster than one GPU (0.56 vs 0.57 s), so it stays on the calculation's own device.
- **Back-transform C = L^-T Q as a column split** (since Sep 29, 2026, both precisions): every eigenvector column needs only its own column of Q and all of L, so each GPU keeps a full copy of L (1.9 GB FP64 / 0.95 GB FP32 at nao 15444, one peer copy per geometry; the calculation's own device uses its L in place) and solves its own columns with a plain `trsm`, without any communication. The distributed `cublasMpTrsm` it replaces did not scale over PCIe (10.6 / 20.8 / 10.9 s on 1 / 2 / 4 A4500; nsys: three devices idle while the fourth ran). polymer_2x FP64 back-transform 20.8 -> 5.6 s (2 GPUs), 10.9 -> 2.7 s (4 GPUs); FP32 0.35 -> 0.12 s per iteration (4 GPUs). The cuBLASMp trsm workspace is no longer allocated. `CURCUMA_GPU_EIG_MP_TRSM=1` restores `cublasMpTrsm` (FP64) and the device-0 FP32 back-transform, `CURCUMA_GPU_EIG_LOCAL_BACKXF=1` only the latter. A device that has no room for the copy of L falls back to `cublasMpTrsm` for that solve.
- Memory: the eigensolver workspaces and the per-GPU solver buffers are released as soon as the SCF has converged, before the post-SCF phase and the gradient. The FP64 verification copy of F is kept on a helper device (see `CURCUMA_GPU_EIG_LOCAL_BACKUP`): polymer_2x device-0 peak 1 / 2 / 4 GPUs 15.3 / 15.8 / 13.3 GB (Sep 29, 2026; before 15.3 / 18.5 / 15.8).
- **Every distributed solve is verified** (`-gpu_eigensolver_verify`, default true). A random vector is pushed through the matrix twice - once directly and once through the returned eigenpairs (for the generalized path with the metric S = L L^T) - and the eigenvalues are checked for a genuinely permuted spectrum. A backend that fails is dropped for that precision, the input is restored and the solve returns to the calculation's own device, with a warning naming the residual. Measured cost: polymer_2x 183.3 s with verification vs 184.9 s without, i.e. inside the noise; the gradients agree to 1.2e-14.
- Why it is not a one-off check: **cusolverMg passes the first solve and degrades afterwards** (polymer, nao 3222: solve 1 residual 1.4e-6, solve 2 8.7e-2 with buffer reuse, 1.2e-3 with everything rebuilt per call). A single gate at the start let a run converge to an energy **1.65 kcal/mol wrong**; with per-solve verification the same run is rejected and lands on the correct energy. cusolverMg is therefore useful only where it verifies: at nao 558 it passes every solve, at 3222 its FP32 is rejected and only its FP64 is used.
- **What FP32 costs on a datacenter GPU** (operator, H200 NVL, Sep 18, 2026): besides being slower
  per iteration there, the FP32 phase can converge to the wrong fixed point at this size. curcuma
  never accepts convergence on an FP32 step, so the result stayed correct
  (-11784.87804452 Eh), but the run needed 42 iterations. Two things now prevent that: mixed
  precision defaults to OFF on such a device (see above), and if an FP64 iteration finds a residual
  more than 10x above what FP32 last claimed, the SCF says so and stays in FP64 for the rest.
- **Second machine, same verdict** (operator, Sep 18, 2026, 2x RTX PRO 5000 Blackwell, a build
  where only cusolverMg was found): `cusolverMg: its FP32 eigenpairs failed verification here
  (relative residual 0.003513, 0 eigenvalues out of ascending order); FP32 iterations stay on
  device 0 and only FP64 is distributed`. The run then converged in 12 iterations to
  -11784.87804452 Eh, i.e. exactly the reference energy, in 91 s (FP32 iteration ~4.0 s on one
  Blackwell, the single FP64 one 25.3 s). So Mg's FP32 defect is not specific to the A4500 box,
  and the fallback does what it is supposed to do on hardware we cannot test here.
- `CURCUMA_GPU_EIG_VERIFY_ALWAYS=1` prints the residual of every verification, `CURCUMA_GPU_EIG_CORRUPT=1` deliberately damages the eigenvectors once to prove the check reacts (it reports residual 0.43 and falls back, and the run still gives the correct energy).
- Results: energies identical to the single-GPU run at the printed precision; with `-scf_threshold 1e-9` energies and gradients agree to <= 1.2e-9 (complex, 231 atoms, gfn1 and gfn2, mp and mg). At the loose default threshold the gradients differ by up to 1e-5 Eh/A, because the FP32 iterations take a slightly different path inside the tolerance (same as GPU vs CPU).

Measured, polymer_2x GFN2 single point + gradient (7320 atoms, nao 15444), 4x RTX A4500 (PCIe, P2P), SCF on device 0, `CURCUMA_GPU_PROFILE=1`, energy -11784.87804452 Eh in every row:

| Configuration | Wall | SCF iterations (FP32+FP64) | per FP32 iteration | per FP64 iteration | peak device 0 | peak other GPUs |
|---|---|---|---|---|---|---|
| single GPU | 419 s | 14+2 | 12.2 s | 85 s | 17.9 GB | - |
| `-gpu_eigensolver_devices 0,1` (syevd only, before the FP64 generalized path) | 327 s | 15+1 | 11.7 s | ~67 s | 17.7 GB | 5.4 GB |
| `-gpu_eigensolver_devices all` (syevd only, same build) | 260 s | 13+1 | ~10.6 s | ~49 s | 15.4 GB | 3.2 GB |
| `-gpu_eigensolver_devices all` (current: FP64 generalized, early release) | 242 s | 13+1 | 9.4 s | 38.2 s | 13.5 GB | 5.0 GB |
| `-gpu_eigensolver_devices 1,2,3` (device 0 not in the solver) | 244 s | 11+1 | 10.3 s | 48.5 s | 12.3 GB | 7.0 GB |

- Iteration counts differ between rows because the FP32 path lands at slightly different points; compare the per-iteration columns.
- Both features are ON by default from 4000 basis functions up when several GPUs are visible (operator decision, Sep 2026); the rows above with fewer devices need the explicit option.
- Measure on idle GPUs: the `mem` column of the profiler is device-wide, so a concurrent job (e.g. `ctest -L gpu`) inflates it and slows the run - one comparison in this file was repeated for that reason.
- FP64 iteration split (all 4 GPUs): sygst 3.8 s (single GPU: two trsm 21 s), syevd 19.5 s (49 s), back-transform 2.7 s (10.5 s; 10.9 s with the former `cublasMpTrsm`). FP32: syevd 4.6 s (7.5 s).
- Including device 0 is ~10 % faster per iteration; excluding it lowers device 0 by ~1.2 GB and raises the others by ~2 GB. Use `1,2,3` when device 0 is the one that does not fit.
- The pattern density is distributed separately, see the next section.
- Gradients at the default `scf_threshold` differ from the single-GPU run by up to 4.6e-4 Eh/A (different FP32 path inside the loose tolerance, same as GPU vs CPU); with `-scf_threshold 1e-9` they agree to <= 1.2e-9 (complex).
### Checklist: making ONE calculation use several GPUs

What is distributed inside a single point / optimisation step:

| part | device | share of an A4500 run (polymer_2x, 4 GPUs) |
|---|---|---|
| eigensolve (FP32 and FP64 iterations) | all listed GPUs | 86 of 194 s |
| screened-pattern density | all listed GPUs | 10 s |
| integrals, Fock build, potential, charges, post-SCF, and the gradient except its W build | the calculation's own device only | ~98 s |

So the ceiling for one calculation is set by the part that stays on one device - on the A4500 box the whole single point went 419 -> 194 s (2.16x on 4 GPUs), not 4x, and that is the honest expectation.

**Prerequisites, in the order they bite:**

1. **cuSOLVERMp, cuBLASMp and NCCL at build time.** Without them only `cusolverMg` (CUDA toolkit) is available, and curcuma gives Mg **FP64 solves only** - its FP32 eigenvectors are wrong at large n (measured: polymer_2x diverged from the first iteration). Since mixed precision is the default, that means roughly one of twenty iterations is distributed and the run looks single-GPU. Verify at configure time:
   `cmake .. -DCUSOLVERMP_ROOT=... -DCUBLASMP_ROOT=... -DNCCL_ROOT=...` must print
   `-- Multi-GPU eigensolver: cuSOLVERMp <path>, NCCL <path>`.
2. **The devices must be listed as INDICES.** `-gpu_devices 2` means "device number 2", not "two GPUs" - use `-gpu_devices all` or `0,1`. Same for `-gpu_eigensolver_devices` / `-gpu_density_devices` (`all`, `solver` or a list).
3. **Check what actually ran**, at `-verbosity 2`, after the SCF:
   `GPU multi-GPU eigensolver: cuSOLVERMp on 2 GPUs, 20 solves` - the backend name and the solve count are the test. `cusolverMg on 2 GPUs, 1 solves` means prerequisite 1 is missing.
   `GPU distributed density: pattern density on devices [0,1], 21 steps`.
4. **The libraries must be found at run time too**: their directories are baked into the RPATH of `libcurcuma_cuda_mgpu.so`, so a compute node needs the same paths (shared file system) or `LD_LIBRARY_PATH`.

Note that the status lines for the density and for a successful eigensolve verification are
`info` (verbosity >= 2); only a failed verification is a warning. A run at the default verbosity
that prints nothing about the density did not necessarily skip it - use `-verbosity 2`.

**Measuring it on a new machine** (nothing here is measured on H200 yet - PCIe A4500 numbers do not transfer to NVLink):

```
CURCUMA_GPU_PROFILE=1 curcuma -sp big.xyz -method gfn2 -gpu cuda -verbosity 2 \
    -gpu_eigensolver_devices none -gpu_density_devices none      # baseline, one GPU
CURCUMA_GPU_PROFILE=1 curcuma -sp big.xyz -method gfn2 -gpu cuda -verbosity 2   # defaults
```

The profile's `eig FP32 / eig FP64 / density P` rows say how much of the run is distributable at all; if they are not the majority, more GPUs cannot help much and the answer is a bigger share on the device instead (integrals and the Fock build are single-device; of the gradient only the W build is distributed, the H0-Pulay pair kernel is the next candidate at 454 ms and a 176 KB reduction).

### Distributed pattern density

`P(r,c) = sum_k Cw(r,k) C(c,k)` is a sum over the occupied columns, so each GPU can evaluate the whole screened pattern over its own slice of columns and the partials are added on the calculation's device. This is exact, not an approximation. Only a column slice of C (n x kn) and one partial pattern array travel per SCF step; the pattern indices are uploaded once per geometry.

| Option | Default | Effect / when to change |
|---|---|---|
| `-gpu_density_devices all\|0,1,2\|solver\|none` | **all visible devices** under the same condition as the eigensolve above, else off | Helper devices for the density (the calculation's own device is always one of the workers and is skipped in the list). `solver` reuses `-gpu_eigensolver_devices`, `none` switches it off. Needs the screened storage (`-gpu_sparse_integrals`, automatic for large systems). |
| `-gpu_density_min_nao N` | 4000 | Below N basis functions the single-device kernel is used (polymer, nao 2975: the split is 10 % slower). |
| `-gpu_multipole_otf auto\|on\|off` | auto | GFN2 multipole interaction matrices: stored, or rebuilt per iteration in the kernel. `auto` rebuilds above ~1 GB of matrices (~2700 atoms), which is what makes 7320 atoms fit a 20 GB card (18 nat^2 doubles = 7.7 GB, plus the same again as host upload copies). Costs ~50 ms/iteration at 1410 atoms, where storing is the better choice (measured: 9.13 s stored vs 9.87 s on the fly). Energies identical either way. Replaces `CURCUMA_GPU_MP_OTF`, which still works and wins over the flag. |

- Measured, polymer_2x GFN2 on 4x RTX A4500 (identical with the defaults and with explicit `-gpu_eigensolver_devices all -gpu_density_devices solver`, both 194 s): density 3.35 -> **0.85 s per SCF step** (46.9 -> 10.2 s total), FP32 iteration 9.4 -> 6.9 s, wall **242 -> 194 s**, energy identical. Device 0 peak 13.5 -> 13.9 GB, helpers 5.0 -> 5.8 GB.
- Too small to pay: polymer (1410 atoms, nao 2975) goes 50 -> 55 ms per step, i.e. the transfers cost more than the split saves. The split needs at least 16 columns per worker and falls back to the single device otherwise.
- On any failure it falls back to the single-device kernel for the rest of the run and warns (at verbosity >= 1).
- Rejected alternatives for the same kernel: a row-major variant (contiguous rows, occupations applied in the kernel) was **8x slower** (4.69 vs 0.60 s over 12 steps on polymer); a banded GEMM does not help because the atom order of polymer_2x is not spatially local - a column block spans 93-99 % of the rows although only 7.9 % of the atom pairs are inside the cutoff.

- Not measured: `mg` end to end on polymer_2x in FP64 (standalone FP64 26.4 s vs `mp` 20.7 s on 4 GPUs), GFN1 on polymer_2x, NVLink systems (H200).

Environment variables:

| Variable | Effect |
|---|---|
| `CURCUMA_GPU_PROFILE=1` | Stream-synchronised per-phase device timings (integrals, potential, Fock, reduce/syevd/back-transform per precision, density, charges, Broyden) plus the post-SCF host phases. Adds synchronisation cost; for diagnosis only. |
| `CURCUMA_GPU_MP_OTF=0\|1` | The old form of `-gpu_multipole_otf`. Still honoured and still wins over the flag, for reproducing older runs. |
| `CURCUMA_GPU_FULL_DENSE_S=1` | Screened integral storage only: compute the dense S/H0 for every atom pair again instead of skipping the pairs outside the screening cutoff (Sep 29, 2026). For comparison only - the skipped integrals are below `sparse_eps` (1e-20). polymer_2x: `k_overlap_h0` 11.5 -> 1.3 s. |
| `CURCUMA_GPU_THREAD_PER_ATOM=1` | Use the thread-per-atom originals of the D4 ATM, on-the-fly multipole potential and multipole energy kernels instead of the warp-per-atom versions (Sep 29, 2026). For comparison only; the originals launch only nat threads and were latency-bound (polymer_2x: D4 ATM 3.8 -> 1.6 s, multipole potential 0.33 -> 0.11 s per iteration). Results differ in the last bits (summation order). |
| `CURCUMA_GPU_EIG_LOCAL_BACKUP=1` | Keep the FP64 verification copy of F on the calculation's device instead of the helper device with the most free memory (Sep 29, 2026). The copy is only read when a distributed solve is rejected; moving it lowered the device-0 peak of polymer_2x by 1.8 GB (2 GPUs 17.6 -> 15.8 GB, 4 GPUs 15.1 -> 13.3 GB). |
| `CURCUMA_GPU_DIST_POTRF=0` | Keep the Cholesky factorisation of S on the calculation's device instead of the multi-GPU eigensolver's devices (cuSOLVERMp `potrf`, since Sep 29, 2026; same gates as the FP64 eigensolve, result checked with `S v = L (L^T v)`). |

Automatic large-system behaviour (no option; triggered above ~5000 basis functions on the CUDA device-resident path):
- the host does not build the dense dipole/quadrupole integrals or the post-SCF host Fock matrix (host fallbacks build them on demand; disabled for `-d4_charge_source cpscf`);
- a single point without gradient keeps P and C on the device and takes the SCC energies from the last device step (`MolecularOrbitals()` is empty in that case).

Reference numbers, polymer_2x GFN2 (7320 atoms, nao 15444), one A4500: 375 s, 16.1 GB device, 15.6 GB host, E identical to the CPU reference (4468 s). One H200 (operator, before these changes): 110 s.

## 3. GFN-FF on the GPU (CUDA)

**One molecule on several GPUs (Sep 29, 2026).** The implicit Coulomb term (tiles of atom pairs) and
the projected-PCG EEQ (column blocks of the EEQ matrix + every PCG matrix-vector product) are spread
over the devices; everything else stays on the calculation's device. Results equal the single-device
run up to summation order (ctest `gfnff_gpu_split_equals_single`: charges 3e-15 e, energies 3e-14 Eh).
MD step polymer_2x (7320 atoms, 16 host threads, A4500): **1 / 2 / 4 GPUs 0.139 / 0.099 / 0.076 s**
(the same step was 0.377 s before the Sep 29 work: frozen topology 0.38 -> 0.20 s, each-pair-once
Coulomb kernel 0.20 -> 0.14 s, split 0.14 -> 0.076 s). Energy call 120 -> 60 ms on 4 GPUs (EEQ 77 -> 35,
Coulomb 32 -> 14 ms). polymer (1410 atoms, one fragment): 16.0 -> 12.0 ms per step, Coulomb only.

**Coulomb far field without erf (Oct 10, 2026, no switch).** For gamma_ij r >= 6, tested as
r^2 >= 36 (alp_i + alp_j), erf(gamma_ij r) is exactly 1.0 in FP64 (CUDA erf saturates from x = 5.9007,
checked on 2e9 sample points) and the Gaussian of the force is below 2.3e-16 relative. The kernels then
evaluate q_i q_j / r (energy) and 1/r (EEQ matrix) instead of erf and exp. This is not a cutoff.
Affected: `k_coulomb_tiles`, `k_coulomb_implicit` (rsqrt when there is no `coulomb_r_cut`, otherwise sqrt
and division after the unchanged cutoff test), `k_coulomb` (sqrt and division), and the matrix builds
`k_eeq_build_matrix`, `k_eeq_build_augmented`, `k_eeq_build_fragment_matrices` (1.0/r). Not changed:
`k_eeq_build_cols` (multi-GPU split, keeps erf; same values since erf is 1.0 there). The ROCm kernels got the same far-field branch on Oct 13, 2026 (see [ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md](ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md)).
- EEQ matrix: bit-identical, 0 differing elements in 25 builds over 12 systems (4 to 7320 atoms, incl. MOR41
  PR16/PR24/PR25 and an 8-water cluster where 252 of 276 pairs take the shortcut).
- Single points (11 systems, every Coulomb kernel variant, with and without `coulomb_r_cut 60`
  and `eeq_distance_cutoff 40`): old/new |dE| <= 4.6e-13 Eh, |dg| <= 1.1e-15 Eh/Bohr; the old binary run
  twice differs by the same amount (atomics order).
- polymer_2x (7320 atoms, RTX 5080, 50 fs MD, n = 3 alternating): `gpu.coulomb` 15.4 -> 6.6 ms,
  `eeq_solve` 16.1 -> 12.3 ms, step 50.2 -> 37.6 ms; Epot after 50 fs unchanged (-902.161126 Eh), PCG
  iterations unchanged (16). Kernel times (CUPTI): `k_coulomb_tiles` 11.2 -> 2.4 ms (3.7 ms with sqrt and
  division instead of rsqrt), `k_eeq_build_matrix` 6.6 -> 2.8 ms. 99.9 % of the pairs are far-field there.
- Not checked: periodic cells (a `unitcell` comment did not change the -sp energy, so the PBC path was
  not exercised), multi-GPU, other GPUs than the RTX 5080 (FP64 at 1/64; on an FP64-strong card the gain is smaller).

| Option | Default | Effect / when to change |
|---|---|---|
| `-gfnff.gpu_split_devices auto\|0,1\|none` | auto | Devices for the split. `auto` = the `-gpu_devices` pool, else every visible device; off in a leased batch worker and when the run is pinned with an explicit `-gpu_device`. `none` = one device. |
| `-gfnff.gpu_split_min_atoms N` | 1000 | Coulomb split from N atoms on (measured gain at 1410 atoms, see above; below 1000 not measured). |
| `-gfnff.gpu_split_eeq_min_atoms N` | 4000 | EEQ split from N atoms on: at 1410 atoms the per-iteration peer latency eats the gain (EEQ 8.1 vs 8.2 ms), at 7320 atoms 77 -> 35 ms. Only the projected-PCG path is split. |
| `CURCUMA_GFNFF_COULOMB_ROWWISE=1` | unset | Diagnostic: the old row-wise Coulomb kernel (every pair twice). The tiled kernel is 3x (7320 atoms, 93.7 -> 31.6 ms) to 6x (1410 atoms, 17.4 -> 2.7 ms) faster; gradients equal to 9e-16. |
| `-gfnff.coulomb_implicit true\|false` | true | **CPU**: evaluate the N^2/2 Coulomb pairs on the fly from the per-atom EEQ charges and alpeeq instead of building and storing the pair list. The stored list is 128 bytes per pair - 3.4 GB and ~0.5 s of pure write bandwidth at 7320 atoms, and threading that loop changes nothing (measured). polymer_2x single point, 36 threads: wall 4.49 -> 3.93 s, host RSS 6.66 -> 3.72 GB, pair setup 496 -> 0.1 ms, and the Coulomb evaluation itself 1733 -> 1096 ms (summed over threads) because it no longer streams the list. MD (polymer, 100 steps) is unchanged at 7.7 s, small molecules unchanged. Energies identical, gradients within 4e-16. Set false for the stored list. Not used with `eeq_distance_cutoff > 0`. |
| `-gfnff.gpu_coulomb_implicit true\|false` | true | The device enumerates all Coulomb pairs itself instead of a host-built N^2/2 list (CUDA `k_coulomb_implicit`, ROCm `k_coulomb_tiles_hip`). polymer_2x: -1.5 s setup, -6 GB host memory, energy identical. Ignored with `eeq_distance_cutoff > 0`. Set false only to compare against the pair-list path. |
| `-gfnff.gpu_disp_pairs_on_device true\|false` | true | Builds the D4 dispersion pair list on the device. polymer_2x: dispersion pairs 630 -> 32 ms, wall 6.8 -> 5.6 s, energy identical, gradient 5e-16; device peak higher (list lives on the device).  Default true on CUDA since Oct 14, 2026 (27 structures from 2 to 7320 atoms, energy and gradient equal to the host build, MD step median equal, mean and longest step lower on bulk systems); the setup build ran before the geometry was on the device (kept all pairs) and the pair order was arbitrary, both fixed, the list is sorted like the host list (`CURCUMA_GPU_DISP_SORT=0` keeps arrival order). ROCm: default true as well since Oct 14, 2026 (ordered build with the gather CSR on the device, rebuild 427 -> 47 ms at 6200 atoms). See [GPU_STEP_PAUSES_2026-10-14.md](GPU_STEP_PAUSES_2026-10-14.md). |
| `-gfnff.gpu_hb_detect_on_device true\|false` | true (CUDA, since Oct 14, 2026) | CUDA only; ROCm ignores it. The H-bond (A, H, B) triples are enumerated on the device at every HB/XB refresh (one warp per AB pair, sorted list, SoA written on the device) instead of by the host detection plus upload. urea400-water1000: 36.2 -> 8.1 ms per refresh, same triple set as the host list on 21 structures (3 to 7320 atoms). `-gfnff.gpu_hb_detect_verify true` compares with the host detection at every refresh. See [GPU_STEP_PAUSES_2026-10-14.md](GPU_STEP_PAUSES_2026-10-14.md). |
| `-gfnff.gpu_disp_pairs_on_device true\|false` | **true on CUDA and ROCm** (trait `default_disp_pairs_on_device`) | Builds the D4 dispersion pair list on the device, also at every MD rebuild. polymer_2x: setup dispersion pairs 630 -> 32 ms, MD rebuild 209 -> 69 ms, energy identical, gradient 5e-16. Default ON since Oct 15, 2026 after 67 structures of `test_cases/structures` (energy difference 0, gradient <= 2.2e-16 against the host list). Since the same day the list is grouped by row and built after the geometry upload; before, the setup build ran with zero device coordinates and held all N(N-1)/2 pairs (energies right, `dc6dcn`/dispersion 1.9x/1.6x slower). |
| `-gfnff.gpu_rep_pairs_on_device true\|false` | **true on CUDA and ROCm** (ROCm since WP9, Oct 25, 2026) | Rebuilds the non-bonded repulsion pair list on the device when the Verlet skin rule fires (2 count/build kernels, per-pair alpha and repab as on the host). Only for distance-filtered lists (systems above `nb_cell_list_min_atoms`). polymer_2x GFN-FF MD, 4 A4500, 300 steps, one run each: rebuild 26.5 ms host build + 5.9 ms upload -> 7.4 ms total (37 rebuilds), mean step 23.3 -> 20.0 ms. Non-bonded repulsion identical to 10 printed digits on 121 output lines over 120 fs (15 device rebuilds), Coulomb differs by 1e-10 at most. Not covered on CUDA: other system classes (metals, charged fragments, many H-H 1,3/1,4 pairs were not singled out in the test). ROCm (RX 9060 XT, WP9): ordered row build with the gather CSR on the device (pairs in ascending i, j; CUDA: arrival order), host build from the same tables as counted fallback (`CURCUMA_ROCM_FAIL_REPDEV=1`); MD 200 steps, every rebuilt list equal to the host list (`CURCUMA_ROCM_REP_DUMP`), repulsion refresh per rebuild peo201-chain 3.66 -> 0.76 ms (15 rebuilds), urea400-water1000 21-31 -> 7.6 ms (27 rebuilds), median step unchanged (6.19 / 43.2 ms), energies per step on vs off <= 3.6e-12 Eh. |
| `CURCUMA_DC6DCN_TABLES=0` (env) | tables on (CUDA) | Disables the tabulated dc6dcn kernel (`k_d4_partner_tables` + `k_dc6dcn_per_pair_tab`): per atom and partner element the reference sums are tabulated once per step, the pair kernel then needs nref instead of nref x nref terms. polymer_2x GFN-FF MD, 4 A4500: kernel time per GPU and step 4.11 ms -> 1.86 ms (incl. table build 0.24 ms), mean MD step 20.2 -> 19.1 ms, median 14.3 -> 12.9 ms (296 steps, one run each). Against the old kernel: 95 MOR41 structures, gradient max difference 1.5e-16 Eh/A; 120 fs of MD, all printed energy terms equal to the print precision. Falls back to the old kernel above 48 distinct elements. |
| `-gfnff.rocm_dc6dcn_tables true\|false` | false (ROCm) | ROCm port of the tabulated dc6dcn (`k_d4_partner_tables` + `k_dc6dcn_per_pair_tab`, WP9); `CURCUMA_DC6DCN_TABLES=0\|1` overrides it on ROCm too. With tables the fused gather `k_dispersion_gather_uv` is off. MD 200 steps, median ms/step of 3 runs, fused uv (default) / uv per pair (`CURCUMA_ROCM_DISP_FUSED=0`) / tables: peo201-chain 6.11 / 6.14 / 6.18, urea400-water1000 43.0 / 47.0 / 47.7, PR24 within the scatter (0.91-1.0). Hence off. Energy against off <= 4.7e-12 Eh per MD step, MOR41 + 308 GMTKN55 `-sp` gradients <= 4e-17 Eh/A. |
| `-gfnff.rocm_device_schur true\|false` | **true** (ROCm) | Schur step of the dense EEQ (cholesky route; N > 128 or more than 8 fragments, below that the one-workgroup kernel) on the device (WP9): nfrag + nfrag^2 sums down instead of N (nfrag + 1). MD 200 steps, ms/step: urea-water-cluster197 with `eeq_rocm_cpu_fragment_threshold 0` 2.99 -> 2.90 (the default CPU fragment route takes 0.98), macrocycle-bgal 2.94 / 2.96, peo201-chain with `-eeq_solver.eeq_ppcg_min_nfrag 0` 19.5 -> 19.3. Fault injection `CURCUMA_ROCM_FAIL_SCHUR=1`. |
| `-gfnff.gpu_eeq_mixed_precision true\|false`, `-gfnff.gpu_eeq_mixed_rho` | **true**, 1e-5 (CUDA split EEQ only, i.e. several GPUs and >= `gpu_split_eeq_min_atoms`) | The PCG iterations multiply with an FP32 copy of the matrix; a round ends at rho x its starting residual, the true residual is formed with the FP64 matvec and the CG restarts from it, until the FP64 residual meets the unchanged tolerance. polymer_2x MD, 4 A4500, 3 runs each x 296 steps: step mean 18.8 -> 17.9 ms, median 13.4 -> 12.2 ms, EEQ median 9.5 -> 8.4 ms; per solve 3 rounds, 17.4 instead of 15.9 iterations. Energies equal to 1e-10 (printed) over 300 steps. rho: 1e-3 4 rounds, 1e-4 4, 1e-5 3 (best), 1e-6 3, 1e-7 3 rounds but more iterations. FP32 without restarts does not work: the true residual stays at 8e4 x the tolerance, the local-inverse guard falls back to Jacobi (67 iterations). Not covered: one GPU (dense symv stays FP64), other system classes, ROCm. `CURCUMA_EEQ_MIXED=0` switches it off as well. |
| `-scf_pseudo_diag_max_dq`, `-scf_pseudo_diag_max_steps` | **0.3**, **12** (were 0.05, 8; since Oct 16, 2026; also read by the opt-in CPU pseudo step) | Gate of the GPU pseudo-diagonalisation. Profile of polymer_2x GFN2 (4 A4500): iterations 0-2 full FP32 solves (5.9 s), 3-10 pseudo (2.9 s), 11 forced full refresh (5.9 s, directly before the final FP64 solve), 12 FP64 (27.4 s). max_steps 12 removes the refresh (-2.9 s), max_dq 0.3 turns iteration 2 into a pseudo step (-2.9 s more): wall 92.5 -> 86.5 s, still 13 iterations, energy equal to 8 decimals. Gradient at the default threshold 1e-5 changes by 6e-6 Eh/A (run-to-run 1e-14); against a `-scf_threshold 1e-8` run without pseudo steps the error goes 1.8e-5 -> 2.2e-5 Eh/A, i.e. within the threshold. MOR41 GPU (95): energies identical, gradient max 3.8e-6 Eh/A. CPU opt-in path, polymer GFN2: energy and 12 iterations unchanged, SCF 13.9 -> 13.75 s. |
| `-scf_gpu_gfn1_resident true\|false` (CUDA, GFN1) | **true** since Oct 16, 2026 (`CURCUMA_GFN1_RESIDENT=0` does the same) | GFN1 runs in the fully device-resident SCF loop of GFN2 (potential, Fock, eigensolve with pseudo-diagonalisation, occupations, density, charges, Broyden on the device; the host polls max\|dq\| and the energies): SCC vector = shell charges, v_sh = gamma q_sh, v_at = Gamma_A q_A^2 (atom-resolved third order), no D4 weights, no multipoles. peo201-water500 (2910 atoms, 8028 AOs, 4 A4500, 3 runs each): wall 44.7 -> 21.7 s, SCF 33.8 -> 13.1 s, 11 iterations either way; iterations 3-9 are pseudo steps (0.49 s instead of 1.4-2 s). Energy equal to 8 decimals; gradient differs by 4e-6 Eh/A (threshold 1e-5). MOR41 GFN1 on the GPU (95): energies identical, gradient max 7.9e-6 Eh/A. Not used with implicit solvation, without the device-built integrals, or for the device `-scf_precond` (GFN2 only). Memory: GFN1 has more AOs than GFN2; urea400-water1000 (17600 AOs) and polymer_2x (20056) do not fit on a 20 GB card (22.4 / 27.9 GB needed on device 0) and fall back to the CPU. |
| `CURCUMA_D4_POSTSCF_DEVICE`, `CURCUMA_D3_DEVICE`, `CURCUMA_MP_AMAT_DEVICE`, `CURCUMA_MP_INTS_EAGER`, `CURCUMA_REP_DEVICE`, `CURCUMA_GRAD_TAIL_DEVICE`, `CURCUMA_GRAD_DEVPOT`, `CURCUMA_D4_CN_DEVICE`, `CURCUMA_CN_HOST`, `CURCUMA_HOST_INTS_EAGER` (env, CUDA GFN1/GFN2) | all device paths **on**; `=0` restores the host code of that piece (`CURCUMA_MP_INTS_EAGER=1` the eager host build of the AO multipole integrals) | Host work removed from a GPU run on Oct 16-17, 2026 (reference: `CURCUMA_GPU_PROFILE=1`, `-threads 1` and `-threads 16` took equally long before and after, the host pieces were serial or memory bound): GFN2 post-SCF D4 (weights, 2-body, ATM with C6 matrices and neighbour list, CN fold), GFN1 D3(BJ) (CN, weights, pair loop, dE/dCN, CN chain rule), multipole interaction matrices built by a kernel (stored path) and the dense AO dipole/quadrupole integrals no longer built eagerly on the host, repulsion energy, the gradient tail (direct multipole gradient, dE/dr -> dE/dCN, CN chain rule) and the converged potentials of the gradient (gamma q_sh, third order, multipole, D4 dE/dq, v_ao) built on the device instead of by the host and uploaded. peo201-water500 (2910 atoms, 4 A4500): GFN2 dispersion phase 633 -> 454 ms (rest = ATM kernel on one card), gradient phase 713 -> 384 ms; GFN1 dispersion 581 -> 72 ms; polymer (1410 atoms) GFN2 total 4.9 -> 3.7 s, multipole setup 749 -> 8 ms. Against the all-host runs: MOR41 (95) energies identical, gradient max 1e-11 Eh/A (GFN2; one earlier run of the same comparison showed 4.9e-9 once, not reproduced in the rerun) and 7.9e-6 Eh/A (GFN1, the resident SCF path, see `-scf_gpu_gfn1_resident`). Second round (Oct 17): the D4 covalent CN (host loop in GenerateParameters, 46 ms at 2910 atoms) and the Hamiltonian CN (downloaded from the device build instead of a host loop) come from the device, the second nsh^2 download of gamma is gone (55 / 115 ms), and S, H0, L, gamma stay on the device and are fetched only when a host consumer asks (`ensureHostIntegrals`: host-driven loops, host gradient, responses, getters; if the device no longer holds the geometry the host builds them). peo201-water500 TOTAL: GFN2 11.54 -> 11.34 s, GFN1 18.88 -> 18.40 s (2 runs each); pre-SCF 223 -> 135 ms (GFN2), 214 -> 81 ms (GFN1). The host copies cost 5.7 GB of RAM and about 1 s at 7320 atoms before. Checked with solvation, `-d4_charge_source cpscf|eeq` (host fallbacks that read the lazy matrices): identical to the eager copies. Still on the host: the single-shot EEQ guess prepares O(N) tables only (the solve is on the device), the GFN1 halogen-bond term (loops over the precomputed X-B-A triples only), `getFock()` is empty on the resident path. PR18 (MOR41) shows a rare run-to-run jump of the GFN2 gradient (4e-7, 1 of 6 sweeps of host pieces against host pieces, i.e. without any device piece involved; consecutive runs with equal settings differ by <= 3e-12, n = 36 pairs): GPU non-reproducibility, not these changes. Not covered by this table: ROCm and Vulkan. |
| `CURCUMA_GPU_DENSITY_SYRK` (env, CUDA GFN1/GFN2, dense integral storage only) | off (`=1` on) | Density P = C diag(occ) C^T by `cublasDsyrk` on sqrt(occ)-scaled occupied columns plus a mirror of the lower triangle instead of the GEMM with Cw = occ C: half the flops, one extra n x ncol buffer. Applies only where the integrals are stored dense (the screened-pattern density of sparse storage, e.g. peo201-water500 and the GFN2 chain, is untouched). Measured on one RTX A4500, `-sp`, wall time median of 3 alternating runs: GFN1 peo201-chain (1410 atoms, dense) 7.74 -> 7.01 s (-9 %), energy equal to 8 digits, 11 SCF iterations both; GFN1/GFN2 urea-water (197 atoms) 1.41-1.46 s both (GFN1: 15 iterations against 14 with the GEMM, same energy); GFN2 peo201-chain 4.45 s both (sparse path); water500 GFN1 28.9 s, GFN2 14.8 s, no change. Not bitwise equal to the GEMM density (rounding), so the default stays off. Not run: MD/opt, open shell, multi-GPU, other cards. |
| `-gfnff.eeq_mixed_precision true` | false | FP32 EEQ factorisation + FP64 refinement (few-fragment paths). For FP64-weak cards; measure per card. |
| `-gfnff.eeq_rmsd_threshold X` (Bohr) | 0 | MD/opt: reuse the EEQ Cholesky factor while the geometry moved less than X per atom (0 = refactor every step, exact reference behaviour). |
| `-eeq_solver.solve_method NAME` | auto | EEQ solver by name, same names on CPU and GPU (`auto`, `cholesky`, `ppcg`, `pcg`, `batched`, `lu`, `ldlt`; `-gpu rocm`: `auto`, `cholesky` only; `auto` = `ppcg` from `eeq_ppcg_min_atoms` = 500 atoms), verbosity 1 prints which one runs: [EEQ_SOLVER_CONTROL.md](EEQ_SOLVER_CONTROL.md). |
| `-eeq_solver.eeq_ppcg_min_nfrag N`, `-eeq_solver.eeq_ppcg_min_atoms M` | 1, 500 | **CPU and GPU alike (GPU since Sep 27, 2026)**: a system with at least N EEQ fragments and M atoms solves the Phase-2 EEQ by projected PCG (GPU: WP7-E) instead of the exact dense Schur-Cholesky (GPU: WP7-A); before, only the CPU switched. polymer_2x (7320 atoms, 1502 fragments), MD step on one A4500: GPU "EEQ + finish" ~180 ms instead of ~1500 ms; CPU CN+EEQ ~975 instead of ~1415 ms (12 threads). Epot of a 4-step MD identical to all 6 printed decimals on CPU/GPU, exact/projected; single point -917.51616997 Eh for all. `N = 0` forces the exact solve on both sides (checked with `-eeq_solver.` and `-gfnff.eeq_ppcg_min_nfrag 0`). The GPU reads `eeq_ppcg_min_nfrag/min_atoms/max_iter` but not `eeq_ppcg_tol`: its WP7-E tolerance is `pcg_tolerance`, absolute on the residual, where the CPU's is relative. Below the thresholds (e.g. an 8-water cluster) the GPU keeps WP7-A. **Not alike for a single fragment**: the GPU takes projected PCG only for nfrag > 1; one fragment always gets the dense device Cholesky, while the CPU switches to projected PCG from 500 atoms ([MULTI_GPU_GAPS.md](MULTI_GPU_GAPS.md) F-4, noted Sep 28, 2026). `-verbosity 2` prints the chosen path (`EEQ GPU Phase 2: ... path=`). |
| `-gfnff.nonbonded_skin_bohr S` | 2 (since Sep 27, 2026; was 0) | Verlet skin for the non-bonded repulsion list (see [GFNFF_PAIR_LIST_REFRESH.md](GFNFF_PAIR_LIST_REFRESH.md)): built at 20 + S Bohr, rebuilt only when an atom moved more than S/2. With S = 0 the list is built on the host every step, which on the GPU path is most of the step. polymer_2x, 30 MD steps (dt 1 fs, CSVR 300 K, 24 host threads, one A4500), measured Sep 27, 2026: energy call 431 -> 182 ms with S = 1 or 2 (rebuilds 31 -> 6 / 3), whole run 24 -> 17 s; CPU only 1176 -> 1120 ms (-5 %). polymer (1410 atoms), 30 steps, 8 host threads: GPU 45.9 -> 26.6 ms (-42 %), CPU 85.5 -> 71.8 ms (-16 %) with S = 2; S = 4 is no faster (more pairs to distance-check). Epot identical to the 6 printed decimals over the 30 steps; polymer on 1 CPU thread over 200 steps: run averages agree to 2e-15 relative (exact up to summation order). Only active above `nb_cell_list_min_atoms` (800); below that the list holds every pair and is never rebuilt. Set 0 to get the old every-step rebuild (e.g. to reproduce earlier trajectories). (Before Sep 27, 2026 the verbosity-2 report added the one-time workspace upload, ~500 ms here, to every call - that was never a per-step cost.) |
| `-gfnff.gpu_block_size N` | 0 (auto) | Kernel launch block size. |
| `-gfnff.gpu_cn_pair_regen*`, `gpu_cn_pair_cutoff_factor` | see `curcuma -help gfnff` | CN-derivative pair list refresh during MD/opt. |

### GPU MD performance switches (CUDA only, Oct 2026)

Measurements and derivations: [GFNFF-PERFORMANCE-ROADMAP.md](GFNFF-PERFORMANCE-ROADMAP.md) (polymer_2x, 7320 atoms, RTX 5080, 16 threads). On ROCm they are ignored unless the row names ROCm. Machine-tested only; ROCm and the multi-GPU split not checked.

| Option | Default | Effect |
|---|---|---|
| `-gfnff.gpu_cn_pair_skin_bohr S` | 1.0 | Verlet skin of the CN pair list: rebuilt on the device only when an atom moved further than the skin criterion, not every MD step; negative = old behaviour. Applies to CUDA and, since Oct 13, 2026, ROCm (CN pair kernels 12.1 -> 0.95 ms per step at 7320 atoms; 4 rebuilds in 51 steps). CUDA roadmap: 94 rebuilds in 1000 steps instead of 397 in 400, `ff_total` median 24.0 -> 19.3 ms. |
| `-gfnff.gpu_eeq_overlap true\|false` | true | The device projected-PCG solver waits only for its inputs and runs on a high-priority stream parallel to dc6dcn and the charge-independent kernels. |
| `-gfnff.gpu_dispersion_overlap true\|false` | true | The D4 dispersion gradient starts right after dc6dcn, parallel to the EEQ solve (changes only the summation order of atomic adds). Roadmap, both overlaps together: wall clock 34.1 -> 30.9 ms/step (n=3); without the solver priority no gain. |
| `-gfnff.gpu_staged_list_upload true\|false` | true | Rebuilt D4/repulsion pair lists and (since Oct 14, 2026) the H-bond list go through a pinned double buffer filled by `-threads` host threads, copied asynchronously; CUDA polymer_2x, 9.4 M D4 pairs, 16 threads: 118 -> 48 ms per rebuild. Roadmap, together with the parallel host rebuilds: 33.7 -> 21.7 ms/step (1000 vs 250 fs, n=3). ROCm (Oct 14, 2026): repulsion lists and their gather CSR only, one pinned block kept across rebuilds, upload 25.4 -> 4.6 ms at 1.47 M pairs ([GPU_STEP_PAUSES_2026-10-14.md](GPU_STEP_PAUSES_2026-10-14.md)). |

Environment variables:

| Variable | Effect |
|---|---|
| `CURCUMA_GFNFF_PROFILE=1` | Setup phase timings summed over both q-loop topology passes (the verbosity-2 report shows only the last pass) and the GPU upload phases. |
| `CURCUMA_GFNFF_DISP_SPLIT=0` | Default on. With the Coulomb split active (several devices, N >= split threshold, no PBC) the D4 dispersion pairs (`k_dc6dcn_per_pair` + `k_dispersion`, ~16 ms per step on one A4500 at polymer_2x) are divided into contiguous ranges over the devices; 0 keeps them all on the first device. Results differ from one device only by FP64 summation order (SP gradient 1 vs 4 A4500: max 1e-15 Eh/Angstrom). Turn off for A/B timing. |
| `CURCUMA_GFNFF_STEPPHASE=1` | stderr, per GFN-FF GPU calculator call: `STEPPHASE` (ms for cn, dlogdcn, hb_cn, hbxb with `[F]` = refresh fired, launch, eeq, coulomb, total); `GAP <section> <ms>` for host sections no phase timer covers (> 0.5 ms: repulsion/dispersion list rebuilds, CN list); `HBXBSUB`/`HBXBDET`/`DHB` for the parts of a fired HB/XB refresh. |

Threading note: GFN-FF topology and parameter generation use OpenMP with the GFN-FF thread budget (`-threads N`); inside batch workers they use the batch's per-worker share. Before Sep 2026 they silently ran on one thread.

Reference numbers, polymer_2x GFN-FF single point (7320 atoms, 1502 fragments), one A4500: 6.8 s (5.6 s with `gpu_disp_pairs_on_device`), 4.8 GB host. CPU with 36 threads: ~6 s.

## 4. CPU knobs of the native GFN1/GFN2 SCF

These apply with and without a GPU (the eigensolve and the reduction run on the host
whenever the device-resident path is not active). Defaults are the measured optimum on the
36-core development box; the optimum is machine-dependent, which is what section 5 is for.

| Option | Default | Effect / when to change |
|---|---|---|
| `-threads N` | 1 | Intra-molecule threads, including the eigensolve. The dense divide-and-conquer eigensolve is memory-bandwidth-bound, so more threads than memory channels can lose: polymer (nao 3222) on 36 cores measured 8 threads 27.1 s, **16 threads 23.8 s**, 24 threads 26.5 s, 36 threads 28.6 s. |
| `-eigensolver_max_threads N` | 0 (follow `-threads`) | Caps the eigensolve alone. Use it when the integrals and the gradient want all cores but the eigensolve peaks earlier. Replaces `CURCUMA_EIG_MAX_THREADS`, which still works and still wins. |
| `-scf_reduce auto\|sygst\|trsm` | auto | Reduction of F to standard form with the cached Cholesky factor of S. `sygst` does half the flops but threads poorly, `trsm` does twice the flops as two BLAS3 solves and keeps scaling; `auto` takes `trsm` from `-scf_reduce_threads` up. Measured n=3222 on 36 cores (OpenMP OpenBLAS), dsygst vs 2x dtrsm: 620/900 ms at 1 thread, 195/192 at 8, 190/**123** at 16, 308/**156** at 36. The two agree to 8e-15 elementwise. |
| `-scf_reduce_threads N` | 8 | Where `auto` switches. Measure with `sygst` vs `trsm` at your thread count before changing it. |
| `-scf_mixed_precision true\|false` | true on CPU and on consumer GPUs, **false on full-rate-FP64 GPUs** (compute capability 6.0/7.0/8.0/9.x/10.x, e.g. H100/H200/Blackwell) | FP32 eigensolves far from convergence, FP64 near it. On an RTX A4500 turning it off costs a factor 1.9 (polymer, 9.2 -> 17.7 s); on an H200 it *is* the slower path (5.29 s vs 3.36 s per iteration) and can converge to a false fixed point. An explicit flag always wins over the device default. |
| `-scf_fp32_threshold X` | 1e-3 (CPU), 1e-5 (GPU) | Where the FP32 phase ends. Smaller = longer in FP32 = faster and less safe. Must stay above `-scf_threshold`. |
| `-scf_fp32_stall_patience N` | 3 | FP32 iterations without real progress (max\|dq\| not improving by 30 %) before the rest of the SCF runs in FP64. FP32 eigenvectors carry ~1e-7 noise, which puts a floor under dq on a large system; without the guard the SCF hovers there. 0 disables. |
| `-scf_fp32_false_fixpoint_factor X` | 10 | If an FP64 iteration reports a residual more than X times what FP32 last claimed, FP32 converged to a false fixed point and the rest runs in FP64. Measured on an H200: FP32 claimed 7.1e-6 at an energy 1.1 kcal/mol off while the truth was 9.4e-3. 0 disables. |
| `-eigensolver mkl\|native\|purify\|lobpcg` | mkl | `mkl` (LAPACK dsyevd) is the fastest on CPU; the others are the GPU-portable / research paths. |

What is deliberately NOT a flag, and why:

| Variable | Why it stays an environment variable |
|---|---|
| `CURCUMA_EIG_PROFILE`, `CURCUMA_GPU_PROFILE`, `CURCUMA_XTB_REDUCE_PROBE`, `CURCUMA_GPU_EIG_VERIFY_ALWAYS`, `CURCUMA_GPU_EIG_CORRUPT` | Diagnostics. They print or inject, they do not tune. |
| `CURCUMA_EIG_TRED2=eigen\|blocked\|scalar` | Only reachable with `-eigensolver native`, itself an opt-in research path. All three are bit-identical; see [SQM_EIGENSOLVE_GPU.md](SQM_EIGENSOLVE_GPU.md). |
| `CURCUMA_GPU_EIG_MG_REUSE=1` | Restores a known-wrong behaviour (cusolverMg returns bad eigenvectors on a reused descriptor) for testing. Making it a documented flag would advertise a setting nobody should use. |
| `CURCUMA_GFNFF_GPU_RESIDENT_HBQ=1` | A deliberate deviation from the reference (live instead of frozen H-bond charges), i.e. a physics experiment, not a performance knob. |

## 5. Measuring it on your own machine

`scripts/tuning_sweep.py` scans these knobs on the hardware it runs on and prints a command
line. Every value it tries is a CLI flag, so nothing has to be rebuilt:

```bash
python scripts/tuning_sweep.py test_cases/molecules/larger/polymer.xyz --method gfn2
python scripts/tuning_sweep.py big.xyz --method gfn2 --gpu cuda --repeats 3 --json sweep.json
python scripts/tuning_sweep.py big.xyz --gpu cuda --knobs gpu_eigensolver_devices,threads
python scripts/tuning_sweep.py big.xyz --gpu cuda --dry-run     # what would be measured
```

`--dry-run` lists the settings and the number of runs without starting any (15 knobs / 92
runs for gfn2 on a 4-GPU 36-core box), and after the baseline the script prints a rough
total time, so an unattended cluster run can be sized before it is submitted.

It runs a baseline, then each knob on its own, then the winners together (a combination can
be slower than its parts, so that is measured rather than assumed), and it **checks every
run's energy against the baseline** - a setting that moves the energy by more than
`--energy-tol` kcal/mol is printed as SUSPECT and never recommended. `--json` keeps the full
table including a hardware note, so a result from the cluster can be read here.

### A full sweep, measured

polymer_2x (7320 atoms, nao 15444), gfn2, 4x RTX A4500, 45 runs (`--repeats 1`), about 3 h.
Baseline 143.95 s, 12 SCF iterations, E = -11784.87804452 Eh. **All 43 runs that completed
returned that same energy to the last printed digit**, including the distributed eigensolve
and the distributed density.

| Knob | values (wall time) |
|---|---|
| `-threads` | 1 / 2 / 4 / 8 / 16 / 32 / 36: **144 s throughout** |
| `-eigensolver_max_threads` | 0 / 4 / 8 / 16: 144 s |
| `-scf_reduce` | auto 143, sygst 144, trsm 144 |
| `-scf_mixed_precision` | true 145, **false 420 (0.34x)** |
| `-scf_fp32_threshold` | 1e-3 **223**, 1e-4 **166**, 1e-5 **144** |
| `-scf_guess` | eeq 144, h0 157 (14 instead of 12 iterations) |
| `-gpu_eigensolver_devices` | all 144, **none 220 (1.53x for the split)**, 0,1 207 |
| `-gpu_eigensolver_backend` | auto 144, mp 144, **mg 2226 (0.06x), 37 iterations** |
| `-gpu_eigensolver_block` | 64 183, 128 145, 256 **140**, 512 148 |
| `-gpu_eigensolver_fp32` | true 145, false 176 |
| `-gpu_eigensolver_verify` | true 144.26, false 144.32 - **free** |
| `-gpu_density_devices` | all 144, **none 193 (1.34x for the split)** |
| `-gpu_multipole_otf` | auto 144, on 144, off: does not fit the card (see below) |
| `-gpu_sparse_integrals` | auto 144, on 144, off: no result within 40 min |

What this says beyond the individual numbers:

- **Both multi-GPU defaults earn their keep at this size** (1.53x eigensolve, 1.34x density)
  and only at this size - at 1410 atoms the same knobs measured flat, which is what the
  4000-basis-function gate is for.
- **`-threads` is irrelevant on the GPU path here.** The integrals are built on the device
  and the host multipole integrals are never built at all ("deferred, 17.2 GB"), so the whole
  setup phase measures 20.63 s at 1 thread and 20.53 s at 36, inside a run whose SCF alone is
  113 s on the device. This says nothing against threading on the CPU path, where the same method on
  polymer measured 44.2 -> 18.8 s.
- **Mixed precision is the one large lever on a consumer/workstation card** (2.9x), and
  `-scf_fp32_threshold` is a real dial rather than a safety switch: 1e-3 costs 55 % over the
  GPU default of 1e-5.
- **`mg` (cusolverMg) is a memory fallback, never a performance option**: 15x slower and 37
  instead of 12 iterations, because the per-solve verification keeps rejecting its
  eigenvectors. The energy is still exact, which is the verification doing its job.
- **Two `off` settings simply cannot run this system**, which is the empirical case for both
  `auto` defaults. `-gpu_multipole_otf off` is instructive about where the limit sits: the 7.2 GB
  of stored matrices still fit the 20 GB card, and the resident SCF loop's own buffers then do
  not, so the failure arrives one step later. It now says so ("the stored GFN2 multipole
  interaction matrices hold 7.2 GB of device memory - rebuild them per iteration instead")
  instead of the bare "GPU resident SCF step failed at iteration 0".

Two things to know before trusting a sweep:
- Run it the way the production job runs (same allocation, same `CUDA_VISIBLE_DEVICES`). The
  multi-GPU knobs in particular only mean something for what the process actually sees.
- The result is for that structure, method and mode. A knob that wins on a 7000-atom single
  point can lose on a 200-atom optimisation - the multi-GPU paths for instance do nothing
  below 4000 basis functions by design, and a sweep on a small molecule will say so.

## 6. Build

### Did this build get cuSOLVERMp, or only the fallback?

The distributed eigensolver is **optional**, and a build without it is easy to miss - that is
how the Sep 2026 H200 run ended up measuring the cusolverMg fallback (15x slower than a single
GPU on polymer_2x) while looking like "multi-GPU does not help". Three ways to check, in the
order they are cheap:

```bash
bash scripts/find_mgpu_libs.sh     # BEFORE configuring: what is on this machine, and the -D flags
cmake .. ... 2>&1 | grep "=== curcuma multi-GPU eigensolver"   # one grep-able summary line
ldd release/libcurcuma_cuda_mgpu.so | grep -Ei "cusolvermp|cublasmp|nccl"   # the built artifact
curcuma -methods            # "multi-GPU eigensolver: mp (cuSOLVERMp)" vs "mg (cusolverMg fallback)"
```

### Where the three libraries come from

They are NOT Python-only, and only one of them is in the CUDA toolkit:

| | in the CUDA toolkit? | other channels |
|---|---|---|
| cusolverMg (the fallback) | **yes** - a plain `module load cuda` has it | - |
| cuSOLVERMp | no | NVIDIA HPC SDK (`math_libs`), standalone tarball, conda `nvidia::libcusolvermp-dev`, PyPI `nvidia-cusolvermp-cu13` |
| cuBLASMp | no | HPC SDK, tarball, conda `nvidia::libcublasmp-dev`, PyPI `nvidia-cublasmp-cu13` |
| NCCL | no, but many clusters put it in the CUDA module or offer its own | HPC SDK (`comm_libs`), distribution package (this box: `/usr/lib/libnccl.so`), conda `nvidia::nccl`, PyPI `nvidia-nccl-cu13` |

On a cluster the least work is usually `module load nvhpc`, which carries all three.
`scripts/find_mgpu_libs.sh` searches `CUDA_HOME`/`NVHPC_ROOT`/conda/site-packages/`ldconfig` and
prints the matching `-DCUSOLVERMP_ROOT=... -DCUBLASMP_ROOT=... -DNCCL_ROOT=...`.

Two PyPI traps, checked against the index in Sep 2026: **`nvidia-cusolver` is the single-GPU
cuSOLVER** (its `-cu13` spelling is deprecated in favour of the unified name) and has nothing to
do with cuSOLVER**Mp**, whose package `nvidia-cusolvermp-cu13` is current and keeps the suffix;
and the suffix-less `nvidia-nccl` is a placeholder package that exists only to warn you.

The configure summary is one of
`cuSOLVERMp + NCCL`, `cusolverMg only - optional cuSOLVERMp missing`, or
`none - one molecule stays on one GPU (everything else works)`, and the "not found" branch names
each missing variable individually plus the pip/HPC-SDK command that supplies it.

**What is lost without them**: only the split of ONE molecule's eigensolve. Single-GPU runs,
batch distribution over several GPUs (`-gpu_devices`) and the distributed density
(`-gpu_density_devices`) are unaffected - and where only `mg` is available,
`-gpu_eigensolver_devices none` is faster than using it.

`-DCURCUMA_REQUIRE_MULTI_GPU_EIGENSOLVER=ON` turns the "not found" case into a configure error,
for machines where the distributed eigensolve is the point of the build. Default is OFF.

- Multi-GPU eigensolver (`CURCUMA_MULTI_GPU_EIGENSOLVER`, default ON): cusolverMg is found in the CUDA toolkit. cuSOLVERMp, cuBLASMp and NCCL are not part of the toolkit; pass their locations, e.g. from the pip wheels `nvidia-cusolvermp-cu13`, `nvidia-cublasmp-cu13`, `nvidia-nccl-cu13`:
  `cmake .. -DCUSOLVERMP_ROOT=<site-packages>/nvidia/cu13 -DCUBLASMP_ROOT=<site-packages>/nvidia/cublasmp/cu13 -DNCCL_ROOT=<nccl prefix>`
  (or an NVIDIA HPC SDK `math_libs` / `comm_libs` directory). The library directories are baked into the RPATH of `libcurcuma_cuda_mgpu.so`, so the run host must see the same paths (shared file system) or have them on `LD_LIBRARY_PATH`. `cmake` prints which backends were found.

- `CMAKE_CUDA_ARCHITECTURES` defaults to `75;80;86;89;90;100;120` with CUDA >= 13 (Turing to Blackwell, incl. H100/H200 sm_90). For fast local builds pass only your card, e.g. `-DCMAKE_CUDA_ARCHITECTURES=86`. A binary without the card's architecture falls back to JIT or fails the kernels.

### CUDA break-even for GFN-FF (measured Oct 17, 2026, AI-run)

RTX A4500 against 16 CPU threads (36 cores visible), `-md`, `-method gfnff`, 1 run per cell, ms per step from the
difference of two runs of different length (start-up cancels). `-gpu auto` and its ROCm atom threshold (`gfnff.gpu_auto_min_atoms`, 900) were
removed on Oct 26, 2026; the numbers below stay as the break-even data (CUDA break-even about 66 atoms, ROCm about 900).

| Atoms | CPU ms/step | GPU ms/step | GPU/CPU |
|---|---|---|---|
| 24 (caffeine) | 0.311 | 0.424 | 1.36 |
| 66 (trisaccharide) | 1.339 | 0.899 | 0.67 |
| 117 | 1.987 | 1.331 | 0.67 |
| 231 | 2.925 | 2.193 | 0.75 |
| 309 | 4.609 | 2.806 | 0.61 |
| 1410 (peo201-chain) | 32.494 | 3.728 | 0.11 |
| 2910 (peo201-water500) | 67.393 | 9.204 | 0.14 |

Per-step break-even lies between 24 and 66 atoms. Each GPU process pays a fixed start-up: caffeine `-sp` 0.78 to 0.95 s
with `-gpu cuda` against 0.02 to 0.05 s on the CPU (3 runs each). At 66 to 300 atoms the step gain is 0.4 to 1.8 ms, so the
GPU needs roughly 500 to 2000 MD steps to repay it; for single points and short optimisations of small molecules it
does not. Not measured: water boxes, charged or multi-fragment systems, other GPUs, CPU thread counts other than 16.

### GFN-FF split path: matrix build, fragment sums (Oct 17, 2026, AI-run, RTX A4500 x4, polymer_2x)

nsys of one MD step (13.3 ms of kernel timeline; wall 20 ms with the host stages) showed the critical path on device 0:
`k_eeq_build_cols` 2.1 ms, then ~18 PCG iterations of 0.29 ms of which the matvec is 0.1 ms, then Coulomb 2.2 ms.

| Change | Measured |
|---|---|
| `k_eeq_build_cols`: 2-D grid (no 64-bit `tid % N`, `tid / N`) | 2.10 -> 1.89 ms per call (400 calls, nsys) |
| same, far field with `rsqrt` instead of FP64 `sqrt` + division (1 ulp, as `k_coulomb_tiles`) | -> 1.33 ms per call |
| fragment sums for many small fragments: warp per fragment (`k_frag_sum_warp`), block per fragment otherwise | 21 -> 9.5 us per call (4095 calls) |
| MD step, `-gpu_split_devices all`, 3 runs of 600 steps (difference of 700 and 100 steps) | 20.33 21.30 20.61 ms before; 18.90 19.22 19.88 and 19.13 19.38 19.29 ms after (6 runs, median 19.25 vs 20.61 ms) |

Values: the fragment sums change their summation order and the far-field column elements differ by 1 ulp, so
trajectories are not bit-identical to the build before, but they stay run-to-run and GPU-count independent
(see GPU_REPRODUCIBILITY.md): 9 of 9 single-GPU MD final states equal, 4, 3 and 2 GPUs equal. Energy and gradient against
the CPU: polymer 4.6e-12 Eh / 4.4e-14 Eh/A. The remaining cost of the matrix build is FP64 throughput (1/64 rate on this GPU).

#### Tried and dropped (Oct 17, 2026, same setup, 3 runs of 600 steps per row, reverted)

| Idea | Result |
|---|---|
| Smaller share of columns/tiles/pairs for the own device (device 0 busy 10.7 of 13.3 ms, helpers 8 ms): weight 0.4 to 1.0 for Coulomb, matrix build and dispersion | Weight 0.4: device-0 Coulomb 2.28 -> 1.30 ms and build 1.40 -> 0.74 ms, but the helpers slow to 2.1-2.45 and 1.5 ms; step period 17.97 -> 18.30 ms (nsys), wall 19.46 (1.0), 19.34 (0.8), 19.53 (0.6), 19.41 (0.4) ms: no gain |
| PCG iteration: single-block dot (one launch instead of cublasDdot's two), projection sum + delta in one kernel | PCG span on device 0 4.59 -> 4.47 ms (nsys, 40 steps); wall 19.25 -> 19.10 and 19.44 ms: not measurable |
| Own block launched first in the split matvec; larger own share of the EEQ columns (weight 1.5 / 2 / 3) | wall 19.98 / 20.49 / 21.82 ms (weight 1.0: 19.79): worse; weights 0.7 and 0.5: 19.14 / 19.16 ms, same as before |

What the trace shows instead: one PCG iteration takes 0.26 ms, of which ~100 us is the own-device matvec kernel and ~100 us the gap
before it, which is the CPU launching the three helper blocks (about 30 us each: set device, wait event, peer copy,
kernel, peer copy, event). Removing that needs the helper chain as a CUDA graph per device or launches from a
helper thread; neither was tried. The PCG span (4.5 ms) is a quarter of the 18 ms step, the rest is the other kernels and the host stages
(integrator gap, HB/XB refresh), so a 0.1 ms change is below the 1 ms run-to-run spread of the step time.

### HB/XB refresh, where the time goes (Oct 17, 2026, AI-run, polymer_2x, `-threads 8`, 31 refreshes in 300 steps)

`CURCUMA_GFNFF_STEPPHASE=1`, median per refresh: detection 20.2 ms (`detectHB` 20.3, of which the pair loop 16.0 and the A-B pair list 4.0
at 8 threads; at 16 threads 9.7 and 2.6), `upload_hb_xb` 4.1 ms, `rebuildBondHBData` 1.95 ms, rest below 0.5 ms: 26.5 ms per refresh,
2.7 ms per step on average (`hb_update_force_every` 10), during which the GPU is idle. The detection reads the geometry of the current
step and Phase 1 (bonds with `nr_hb`, HB kernels) needs its result, so it cannot run ahead of or beside the GPU work of the same step
without changing which step uses the list; only the EEQ solve is independent of the lists, and moving it before Phase 1
means reordering the input events of the K3 overlap. Not done.

Done: `HBondSoA` upload through the staged pinned upload (`stagedHBondUpload`, same device bytes: the 300-step trajectory and final
state are identical to the build before) 4.1 -> 2.6 ms per refresh (HB and XB together), i.e. about 0.15 ms per step.
Remaining reducible parts: the pair loop of the detection (threads), `rebuildBondHBData` (serial host).

#### Helper chain of the split matvec as a CUDA graph (Oct 17, 2026, tried and dropped)

Each helper device's chain (copy of x, column dots, copy of y) was captured once per split setup and launched with
`cudaGraphLaunch`. `cudaMemcpyPeerAsync` is refused during capture ("operation not permitted when stream is capturing"); `cudaMemcpyAsync` with
`cudaMemcpyDefault` on the UVA addresses captures. Result (nsys, 100 steps, polymer_2x, 4x A4500): 5241 graph launches at 13.3 us average CPU time,
against ~14 us for the three plain calls (`cudaMemcpyPeerAsync` 5.4 us, kernel launch 3.3 us, second copy 5.4 us) they replace.
PCG span 4.80 (graphs off) vs 5.02 ms (on), step period 18.78 vs 19.21 ms, gap between two own-device column dots 161 vs 165 us; wall 19.03 / 19.31 ms
(off, 2 x 3 runs) against 19.15 / 19.29 ms (on). The earlier reading that the CPU launch of the three helpers causes the ~100 us gap before the
own kernel was therefore wrong: launching costs ~45 us in total. The 160 us between two matvecs is the ~20 small dependent kernels of the PCG
(dots, updates, two projections, the preconditioner, a host round trip for the convergence test), each 1-15 us long with ~1 us gaps.
What would shorten it is fewer dependent steps per iteration (one persistent kernel with grid synchronisation, or s-step/pipelined CG
with one reduction per iteration), not cheaper launches.

#### Pipelined CG for the split EEQ solve (Oct 17, 2026, tried and dropped)

Chronopoulos-Gear CG with the projected local-inverse preconditioner: `s = A p` by recurrence, one fused reduction
(`r.u`, `w.u`, `r.r`, one block of 1024 threads) per iteration instead of three dependent dots, the projection applied once to `w`; mixed-precision
rounds and restarts as before. It converged on every solve (polymer_2x, 60 MD steps).

| | classic loop | pipelined CG |
|---|---|---|
| matvecs per solve (nsys, 40 steps) | 17.4 | 19.4 |
| time between two own-device column dots | 157 us | 141 us |
| PCG span on device 0 | 4.63 ms | 4.93 ms |
| step period (nsys) / wall, 2 x 3 runs of 600 steps | 17.97 ms / 18.97 and 19.47 ms | 19.13 ms / 19.49 and 19.38 ms |

The iteration got 10 % shorter and needs 2 more matvecs (the recurrence for `s` lags one step), so the solve is not faster. A second attempt that lets
the helper kernel write its result straight into the own device's vector (no copy back) gave 19.06 and 19.19 ms: no change.

What the iteration timeline shows (polymer_2x, 4 GPUs): the critical path of a matvec is the last helper. The CPU issues the three helper chains
one after the other (events, peer copy, kernel, peer copy: 13-17 us each), the third helper's kernel starts ~50 us after the first, its result
arrives ~25 us after the kernel ends, and only then do the dots follow. The own-device kernel (97 us) and the rest of the iteration (~25 us of
small kernels, one host round trip) are in series with that. Remaining ideas, none tried: uneven column shares so the last-launched helper gets
fewer columns (also changes the matrix build per device), one host thread per helper device to issue the chains in parallel, an in-kernel
(persistent) PCG with grid synchronisation. Each is bounded by the PCG span of 4.5 of 18 ms.

### Sweeping the MD step (Oct 18, 2026)

`scripts/tuning_sweep.py ... --method gfnff --gpu cuda --md` scans the knobs of this document for the GFN-FF MD step time (ms per step from runs of
two lengths, warm-up run discarded because the first run writes the topology cache and understated the baseline: 14.0 against 19.2 ms/step). Run on the
A4500 box for validation (polymer_2x, 2 repeats, subset of knobs): `auto` 19.13, `none` 41.35, `all` 18.58, `0,1` 25.17, `0,1,2` 20.24 ms/step;
`gpu_eeq_mixed_rho` 1e-5 / 1e-4 / 1e-6: 19.05 / 19.13 / 19.20 ms/step (noise floor 0.6 ms); the energies of all settings agree with the baseline within the tolerance.
The H200 procedure is item 8 of docs/MULTI_GPU.md.

