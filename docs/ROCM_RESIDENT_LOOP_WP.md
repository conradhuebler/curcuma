# ROCm: device-resident GFN1/GFN2 SCF loop, work packages (from Oct 24, 2026)

Goal: close table A of [GPU_PARITY_CUDA_ROCM.md](GPU_PARITY_CUDA_ROCM.md) for the SCF loop, i.e. the ROCm backend sets
`supportsDevicePotential`, `supportsStoredMultipoleOnDevice`, `supportsResidentLoop` (+ `Gfn1`, open shell, preconditioner)
like the CUDA class `XtbGpuScfBackend` (`xtb_gpu_method.cpp`). Reference implementation: `cuda/xtb_gpu_context.cu`
(Stage 5/6: `beginPotential*`, `residentSolvePotential`, `buildDevicePotentialOnly*`, `convergedPotentials`,
`beginResidentLoop*`, `residentScfStep`, `residentLoopCharges`, `runBroydenUpdate`, `k_occupations*`, `sccEnergy`).
Host side (unchanged contract): `xtb_native.cpp` lines ~1090-1760 and 2160-2260.

Rules for every package:
- Port to `rocm/xtb_hip_context.{h,hip}` and `xtb_hip_method.cpp`; CUDA files stay untouched (they cannot be compiled here).
- Same numerics as CUDA and CPU (FP64 unless the CUDA code uses FP32 behind an existing switch); no new defaults that change
  results; each new path has a switch to turn it off and falls back to the existing ROCm host-loop path with a counted notice.
- Wave32 (gfx1200) is the only measured target; no new warp-size assumptions without a comment (Known Issue F-17).
- Acceptance per package (reference = CPU `-gpu none`, baseline ROCm = `/home/conrad/src/bench_wt/rocm/release/curcuma`):
  energy `|dE| <= 1e-9 Eh` against CPU (baseline ROCm vs CPU is 5.5e-12), gradient `|dG| <= 1e-5 Eh/A`, on caffeine, caffeine +1
  (doublet), 12 MOR41 structures (AcCl AcI Bz CO ED01 ED07 ED16a ED24 PR07 PR16 PR25 PR41), peo201-chain (1410 atoms),
  plus `ctest -R "gpu|rocm|sqm|xtb|grad"` (known failures: `xtb_cpscf`, `cli_sqm_12`).
- Timing per package: wall of peo201-chain `-sp` and of PR24 `-md` 100 fs, median of 3, before and after, same binary state.

| WP | Content | State |
|---|---|---|
| 1 | Device potential in the host loop: `supportsDevicePotential`, `beginPotential`/`beginPotentialStored`, `solvePotential`, `convergedPotentials` (Coulomb gemv, third order, stored multipole potential, D4 dE/dq, expanded v_ao on the device) | done Oct 24 (AI, machine-tested; `convergedPotentials` moved to WP2, see notes) |
| 2 | Resident GFN2 loop, closed shell: occupations, charges from populations, atomic moments, SCC energy, device Broyden, D4 reference weights, `beginResidentLoop`/`residentScfStep`/`residentLoopCharges` | done Oct 24 (AI, machine-tested; gradient potentials stay on the host, see notes) |
| 3 | Resident GFN1 loop (`beginResidentLoopGfn1`, atom-resolved third order) | done Oct 24 (AI, machine-tested; see notes) |
| 4 | Open shell on the resident loop (`k_occupations_uhf`, `setResidentOpenShell`) | done Oct 24 (AI, machine-tested; see notes) |
| 5 | Resident Broyden preconditioner hook (`supportsResidentPreconditioner`) | done Oct 24 (AI, machine-tested; GFN1 and GFN2, see notes) |
| 6 | Device solvation (`supportsDeviceSolvation`) | done Oct 24 (AI, machine-tested; GFN1 and GFN2, see notes) |
| 7 | Pseudo-diagonalisation on the device (`k_pseudo_*`; the eigensolve is 905 of 1200 ms per iteration on peo201-chain) | done Oct 24 (AI, machine-tested; ROCm default off, see notes) |
| 8a | Post-SCF on the device: GFN2 D4 (weights, 2-body, ATM, CN fold), GFN1 D3(BJ), D4 covalent CN, Hamiltonian CN download (`d4PostScfDevice`, `d3DispersionDevice`, `d4CovalentCnDevice`, `downloadCn`) | done Oct 24 (AI, machine-tested; see notes) |
| 8b | Gradient tail on the device, `convergedPotentials`/`prepareConvergedPotentials`/`gradientResident`, skipped host build of the multipole interaction matrices (`defer_amat`), repulsion energy on the device | done Oct 25 (AI, machine-tested; see notes) |
| 9 | GFN-FF: (A) device rebuild of the non-bonded repulsion list, (B) tabulated dc6dcn, (C) device Schur step of the dense EEQ with several fragments | done Oct 25 (AI, machine-tested; B default off by measurement; see notes) |
| 10 | GFN-FF reproducibility audit on ROCm and fixed-point accumulation (`fix_add`, `k_fix_*` of CUDA); same for the GFN1/GFN2 gradient kernels | done Oct 25 (AI, machine-tested; default on for both, see notes) |
| 11 | GFN2 multipole interaction on the fly (`supportsOnTheFlyMultipole`, `beginPotentialOnTheFly`, CUDA `k_multipole_potential_otf` / `k_energy_multipole_otf`) for systems above the 4 GB stored-matrix limit | done Oct 25 (AI, machine-tested; auto above the limit only, see notes) |

Notes per package are appended below by whoever works on it (what was tested, what not, what differs from CUDA).

### WP1 (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switch `-xtb.rocm_device_potential`, default true (PARAM in `native_xtb_method.h`, fallback default in
  `xtb_hip_method.cpp`). `HipScfBackend` (`xtb_hip_method.cpp`) overrides `supportsDevicePotential`,
  `supportsStoredMultipoleOnDevice`, `beginPotential`, `beginPotentialStored`, `solvePotential`, `lastError`,
  `devicePotentialFallback`; the context (`rocm/xtb_hip_context.hip`) has `potBegin`, `potBeginStored`,
  `resSolvePotential` and the kernels `k_mp_build_amat15`, `k_vsh_third`, `k_vat_add_d4`, `k_expand_vao`; the multipole
  potential reuses `k_mp_potential_stored` and the `dMpAmat` buffer of `-xtb.rocm_device_multipole` (15-block layout,
  so the device build writes the 6 distinct dd blocks instead of CUDA's 9). gamma q_sh is `rocblas_dgemv` on the resident
  `dgamma`. q_at is summed on the host from q_sh in the order of `unpackSCC`. D4 reference weights stay host-built.
- Not ported: `convergedPotentials`/`prepareConvergedPotentials`/`gradientResident`. The host calls them only on the
  resident loop (`xtb_native.cpp`, `use_resident_loop`); in the host loop the post-SCF potential is rebuilt on the host
  and the existing `gradient()` path runs unchanged. They belong to WP2.
- Host seam change (outside the ROCm files): `GpuScfBackend::devicePotentialFallback` (`xtb_native.h`, default false, so
  CUDA keeps the hard error) and in `xtb_native.cpp` a failed `solvePotential` continues the SCF with the host potential
  when the backend allows it. A failed begin already led to the host potential build. Either case prints one warning per
  method object (`reportSolverFallbacks`), not counted in `gpu_fallback` (no CPU fallback). Fault injection
  `CURCUMA_ROCM_FAIL_POT=begin|solve` (`solve` from the third call): caffeine energy equal to the switch off to 1e-14 Eh.
- Host work left per iteration: D4 weights (`buildRefWFlat`), occupations, Broyden, SCC energies incl. the host
  multipole energy, and the host interaction matrices are still built (`defer_amat` needs `supportsResidentLoop`).
- Acceptance (gfn2 `-sp -threads 10`, n = 1 each, CPU = `-gpu none`): caffeine, caffeine +1 doublet, 12 MOR41: dE on vs CPU
  max 5.5e-12 Eh (ED01; off vs CPU 5.5e-12 too), dG on vs CPU max 2.7e-6 Eh/A (ED01), on vs off max 5.7e-14 Eh,
  1.5e-14 Eh/A. peo201-chain: dE 8.2e-12 Eh (same off), dG 3.6e-8 Eh/A, on vs off 0 Eh, 1.3e-14 Eh/A.
  `ctest -R "gpu|rocm|sqm|xtb|grad"`: 206 of 208 passed, failures `xtb_cpscf`, `cli_sqm_12_gfn2_scf_precond` (known).
- Timing (median of 3, alternating, same binary, `-threads 10`): peo201-chain `-sp` 14.33 -> 14.18 s (potential build
  10.9 ms -> device per iteration, eigensolve 905 ms per iteration dominates); PR24 `-md` 100 fs 54.76 -> 54.60 s.
- `-opt` caffeine (33 steps): with mixed precision on, on and off differ from step 1 in the 6th digit of the gradient
  norm; with `-scf_mixed_precision false` both agree in every printed digit. Read: the FP32 SCF steps amplify the
  last-bit differences of the device sums. Not measured further.
- Not run: GFN2 with solvent (the host path stays, `supportsDeviceSolvation` false), d-shell and metal systems beyond
  MOR41, systems above 4 GB of interaction matrices (about 5770 atoms, then the host path), multi-GPU, peak VRAM, ROCm
  cards other than gfx1200.

### WP2 (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switch `-xtb.rocm_resident_loop`, default true (PARAM in `native_xtb_method.h`, fallback default in `xtb_hip_method.cpp`);
  needs `-xtb.rocm_device_potential`. `HipScfBackend` overrides `supportsResidentLoop`, `beginDispersionWeights`,
  `beginResidentLoop`, `residentScfStep`, `residentLoopCharges`, `residentLoopFallback`; the context has `resBeginWeights`,
  `resBeginLoop`, `resStep`, `resBroyden`, `resCharges`, kernels `k_res_*`, and `potStepCore` (the device half of WP1's
  `resSolvePotential`, shared). Open shell (`setResidentOpenShell` stays false), GFN1, solvent, `-scf_precond` other
  than none, `CURCUMA_SCF_PRECOND_FILE` and `-keep_diis true` keep the host-driven loop.
- Deliberate differences from the CUDA loop, so that the ROCm loop follows the host-driven loop: the SCC vector is packed
  in the host `packSCC` order; the D4 reference weights are built from the mixed input q_at (step 0: the guess q_at, as
  the host) where CUDA uses the previous output q_at; charges, moments, band and energy sums are fixed-order gathers in
  the host summation order (no atomics); the weight kernel is the current host `buildAtomRefW` incl. its underflow rule
  (CUDA `k_d4_build_refw` sets gwk = 0 there). No pseudo-diagonalisation (WP7): full eigensolve every step.
- Not ported: `convergedPotentials`/`prepareConvergedPotentials`/`gradientResident` (base false), so after the loop the
  host downloads the charges, rebuilds the potentials and runs the existing `calculateGradientGpu`; P and C are
  downloaded as before (`m_wfn_on_device` stays false because ROCm downloads the multipole integrals).
- Host seam change (outside the ROCm files): `GpuScfBackend::residentLoopFallback` (`xtb_native.h`, default false, CUDA
  keeps the hard error); in `xtb_native.cpp` a failed `residentScfStep` restarts that SCF from the unchanged initial
  guess on the host-driven loop. A failed begin already led to the host-driven loop. Either case turns the loop off for
  the rest of the run with one warning per method object, not counted in `gpu_fallback`. Fault injection
  `CURCUMA_ROCM_FAIL_RES=begin|step` (`step` from the third call): caffeine energy and gradient equal to the switch off
  (dE 0, dG 2.4e-16 Eh/A), 13 iterations both.
- Acceptance (gfn2 `-sp -threads 10 -dump_gradient`, n = 1 each, CPU = `-gpu none`): caffeine and 12 MOR41: dE on vs CPU
  max 5.5e-12 Eh (ED01; off 5.5e-12), dG on vs CPU max 2.7e-6 Eh/A (ED01), on vs off max 2.1e-13 Eh, 4.5e-14 Eh/A,
  SCF iterations on = off for all 13 (ED01 18 against CPU 19, as off). Caffeine +1 doublet took the host-driven loop:
  on = off exactly, dE vs CPU 2.1e-12 Eh. peo201-chain: dE on vs CPU 2.3e-12 Eh (off 8.2e-12), on vs off 1.0e-11 Eh,
  dG 3.6e-8 Eh/A, 12 iterations each. `ctest -R "gpu|rocm|sqm|xtb|grad"` with default true: 206 of 208 passed,
  failures `xtb_cpscf`, `cli_sqm_12_gfn2_scf_precond` (known).
- `-opt` caffeine 30 steps: with `-scf_mixed_precision false` on and off agree in every printed digit for all 30 steps;
  with mixed precision on they differ from step 1 in the 6th digit of the gradient norm (same pattern as WP1).
  `-md` PR24 100 fs `-seed 42`: energies, temperatures and the trajectory file agree in every printed digit.
- Timing (median of 3, alternating, same binary, `-threads 10`): peo201-chain `-sp` 14.18 -> 14.06 s; PR24 `-md` 100 fs
  54.50 -> 54.01 s. The eigensolve still dominates (WP7).
- Not run: d-shell and metal systems beyond MOR41, Tele = 0 (integer filling branch), Broyden history lengths other than the default,
  systems above 4 GB of interaction matrices, multi-GPU, peak VRAM, ROCm cards other than gfx1200.

### WP3 (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switch: the WP2 `-xtb.rocm_resident_loop` (same exclusions: `-scf_precond` other than none, `CURCUMA_SCF_PRECOND_FILE`,
  `-keep_diis true`) plus the generic `-scf_gpu_gfn1_resident` / `CURCUMA_GFN1_RESIDENT=0` that the host already reads;
  GFN1 does not need `-xtb.rocm_device_potential`. `HipScfBackend` overrides `supportsResidentLoopGfn1` and
  `beginResidentLoopGfn1`; the context has `resBeginLoopGfn1`, `potStepCoreGfn1` (gamma q_sh by `rocblas_dgemv`, v_at by
  `k_vsh_third` launched over the atoms) and `fockAndSolve` (Fock + eigensolve, now shared with WP1's `potStepCore`);
  `resStep`, `resBroyden`, `resCharges` and the `k_res_*` kernels serve both loops (`k_res_unpack_atoms` got a moments
  flag, the loop buffers moved to `resAllocLoop`). Fault injection `CURCUMA_ROCM_FAIL_RES=begin|step` covers GFN1.
- Follows the host-driven GFN1 loop: SCC vector = q_sh; the third order of step 0 uses the guess q_at (as the host,
  which reads `m_wfn.q_at` before the first unpack), later steps the q_at summed from the mixed q_sh; output q_at from
  the shell populations; E3 = sum q_A^3 Gamma_A / 3 in atom order. CUDA uses the q_at summed from q_sh in step 0 too.
- No host change. After the loop the host rebuilds the potentials and runs the existing gradient (`gfn1_resident`
  branch); halogen bond and D3(BJ) stay host post-SCF (WP8). Open shell: the begin succeeds, `setResidentOpenShell`
  returns false, the host-driven loop runs (caffeine +1: on = off exactly).
- Acceptance (gfn1 `-sp -threads 10 -dump_gradient`, n = 1 each, CPU = `-gpu none`): caffeine, fi-pyridine (halogen
  bond, I), 12 MOR41: dE on vs CPU max 3.4e-12 Eh (PR07; off 3.7e-12), dG on vs CPU max 2.6e-6 Eh/A (AcCl), on vs off
  max 2.7e-13 Eh, 1.4e-14 Eh/A; SCF iterations on = off for all (AcCl 9 against CPU 10, as off). peo201-chain: dE on vs
  CPU 1.1e-10 Eh (off 1.8e-12), on vs off 1.1e-10 Eh, also with `-scf_mixed_precision false` (dE 1.1e-10, dG on vs off
  2.1e-14 Eh/A, 10 iterations each): the charges agree, the energy difference is in the summation of the device energy
  sums (band Tr(P H0) and q gamma q as fixed-order device sums against the host Eigen sums); which of the two was not
  isolated. GFN2 after the change: caffeine, AcCl, Bz, PR07 on and off equal to before (dE 0, dG <= 1.9e-14 Eh/A).
  Fault injection begin and step: caffeine energy equal to off (dE 0), one warning each.
  `ctest -R "gpu|rocm|sqm|xtb|grad"`: 206 of 208 passed, failures `xtb_cpscf`, `cli_sqm_12_gfn2_scf_precond` (known).
- `-opt` caffeine 30 steps: with `-scf_mixed_precision false` on and off agree in every printed digit; with mixed
  precision they differ from step 2 in the 6th digit of the energy change (same pattern as WP1/WP2). `-md` PR24 100 fs
  `-seed 42`: with `-scf_mixed_precision false` the tables and the trajectory file are identical; with mixed precision
  the temperatures differ from the first print (10 steps) in the 7th significant digit.
- Timing (median of 3, alternating, same binary, `-threads 10`): peo201-chain `-sp` 18.64 -> 18.51 s; PR24 (120 atoms)
  `-md` 100 fs 37.16 -> 37.66 s, i.e. slower: 35.0 -> 37.4 ms per SCF iteration at 120 atoms (cause not measured; the
  GFN1 host loop has no D4 weight work to save, so the per-step device overhead shows).
- Not run: GFN1 above nao ~5270, where the host keeps P and C on the device (`m_wfn_on_device`) after the GFN1 loop;
  GFN1 with solvent (host-driven loop by the host condition), Tele = 0, d-shell metals beyond MOR41, multi-GPU, other cards.

### WP4 (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switch `-xtb.rocm_resident_open_shell`, default true (PARAM in `native_xtb_method.h`, fallback default in
  `xtb_hip_method.cpp`); acts only with `-xtb.rocm_resident_loop`. `HipScfBackend::setResidentOpenShell` (the host calls it
  after a successful begin for n_alpha != n_beta, `xtb_native.cpp` unchanged); the context has `setResidentOpenShell` /
  `resSetOpenShell` and the kernel `k_res_occupations_uhf` with `res_fermi_channel` (port of `XTB::fermiFillingChannel`
  and the CUDA `occupation_channel`); `resStep` switches the occupation kernel; a new begin resets the open-shell state.
  Density, charges, moments, energies and Broyden are the closed-shell code (columns up to the last occupied column
  over both channels). After the loop the host rebuilds `focc` and the entropy term from the downloaded eps
  (`occupationsFromEps`), and the device gradient gets `eps_eff = f eps / 2`, so W = sum f eps C C^T on ROCm as on the CPU.
- Difference from CUDA: the Newton number and derivative sums run in ascending orbital order on one thread (the host
  order) instead of a block tree.
- Acceptance (`-sp -threads 10 -dump_gradient`, n = 1 each, CPU = `-gpu none`). Open shell: caffeine +1, GMTKN55
  RSE43/P45, MB16-43/07, HEAVYSB11/pbme3 (Pb), HEAVYSB11/teme (Te), RC21/8e (+1), G21EA/EA_11 (-1), SIE4x4/nh32+_1.0 (+1)
  (doublets), MB16-43/O2, W4-11/ch2-trip, G21EA/EA_14 (-1), BH76/RKT19, YBDE18/ch2 (triplets), W4-11/n, W4-11/p,
  G21EA/EA_c- (-1), G21IP/o+ (+1) (quartets), MOR41 ED01 as triplet (Cr), PR07 +1 (W), PR24 +1 (120 atoms); GFN1 and GFN2,
  plus caffeine +1, ED01 triplet, MB16-43/07, W4-11/ch2-trip at `-xtb.electronic_temperature` 0, 300, 1000, 5000. The
  resident loop ran for every open-shell case (verbosity 2 line). GFN2 (35 runs): dE on vs CPU max 1.6e-11 Eh, dG max
  2.8e-6 Eh/A, on vs off max 2.0e-13 Eh, 4.2e-13 Eh/A. GFN1 (36 runs): dE on vs CPU max 8.7e-11 Eh, dG 8.9e-6 Eh/A (ED01
  triplet at T = 0, 97 iterations; off vs CPU the same), on vs off max 2.0e-13 Eh, 2.0e-11 Eh/A. SCF iterations on = off
  in every run; ROCm (on and off alike) differs from the CPU by up to 2 iterations for the ED01 triplet. ED01 triplet with
  GFN2 at T = 0 does not converge in 151 iterations on CPU, off and on alike.
- Closed shell with smearing (caffeine, Bz, PR07, ED01 at T = 0, 300, 1000, 5000; T = 5000 lowers E by up to 4e-2 Eh, so
  the fractional occupations matter): dE on vs CPU max 1.8e-11 Eh, dG 3.1e-6 Eh/A, on vs off 2.4e-13 Eh, iterations
  on = off. WP2/WP3 closed shell after the change: caffeine, AcCl, Bz, PR07, GFN1 and GFN2, energy equal to before
  (dE 0), same iteration counts.
- Switch off: caffeine +1 takes the host-driven loop, energy equal. Fault injection `CURCUMA_ROCM_FAIL_RES=begin|step` on
  caffeine +1 (GFN1, GFN2): energy equal to off (dE 0), one notice at the end (plus the verbosity-2 restart line for step).
- `-opt` caffeine +1, 30 steps, and `-md` 50 fs `-seed 42`, GFN1 and GFN2: with `-scf_mixed_precision false` on and off
  agree in every printed digit and in the coordinates of the final structure and trajectory; with mixed precision the opt
  differs in the last printed digit of the energy change from step 2 (GFN1) / step 7 (GFN2), the GFN1 MD temperatures in
  the 6th significant digit, the GFN2 MD not in any printed digit (same pattern as WP1-WP3).
- Timing (`-sp`, median of 3, alternating, same binary, `-threads 10`): PR24 +1 GFN2 1.14 -> 1.14 s (28 iterations),
  GFN1 0.96 -> 1.02 s (19 iterations; the per-step GFN1 overhead seen in WP3), peo201-chain +1 GFN2 43.07 -> 42.81 s
  (30 iterations, energy equal to 1e-8 Eh as printed).
- `ctest -R "gpu|rocm|sqm|xtb|grad"`: 206 of 208 passed, failures `xtb_cpscf`, `cli_sqm_12_gfn2_scf_precond` (known).
- Not run: transition-metal open shells beyond the artificial ED01 triplet and PR07 cation, fractional n_alpha / n_beta,
  open shell with solvent (host-driven loop), multi-GPU, other cards.

### WP5 (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switch `-xtb.rocm_resident_precond`, default true (PARAM in `native_xtb_method.h`, fallback default in `xtb_hip_method.cpp`);
  acts with `-scf_precond frag` on the resident loop; false restores the old rule (`-scf_precond` other than none keeps the
  host-driven loop, where frag is ignored with a warning). `CURCUMA_SCF_PRECOND_FILE` and `-keep_diis true` still keep the
  host-driven loop (CUDA also restarts the device history at every begin and leaves the file matrix to the CPU loop).
- `HipScfBackend` overrides `supportsResidentPreconditioner`, `supportsResidentPreconditionerGfn1` (new base virtual, default
  false), `residentPreconditionerTag`, `residentSetPreconditioner`; it maps the CUDA atom-major indices that `xtb_native.cpp`
  produces back to the host `packSCC` order of the ROCm loop (GFN1: identity). The context has `broydenSetPreconditioner`
  (`pcSet`), `pcApply`, `pcDrop`, `broydenPreconditionerTag/Active`, `broydenGuardFiredAt`, kernels `k_res_pc_diag`,
  `k_res_pc_blocks` (one thread per row, row sum in column order, no warp-size assumption), `k_res_broyden_u_pc`,
  `k_res_broyden_u_unprecond`; `resBroyden` follows CUDA `runBroydenUpdate` (guard on `rocblas_dnrm2` of F, x_best/F_best
  device copies). Without M the WP2 kernel sequence is unchanged.
- Host changes outside the ROCm files: `GpuScfBackend::supportsResidentPreconditionerGfn1` (`xtb_native.h`), its use in
  `XTB::selectScfPreconditioner` (`xtb_precond.cpp`, GFN1 device loop only where the backend says so), the text of the
  "ignored" reason and a comment in `xtb_native.cpp`. CUDA keeps the base false, so its GFN1 behaviour is unchanged.
- Acceptance (`-sp -threads 10 -dump_gradient`, n = 1 each, CPU = `-gpu none` with the same switches): urea3-water8, water8,
  unrelaxed urea-water 197, peo8 segmented (`-scf_precond_segment_atoms 3 -scf_precond_max_fragment_atoms 20`), caffeine,
  AcCl, Bz, PR07, GFN1 and GFN2: dE ROCm frag vs CPU frag max 1.1e-11 Eh (GFN2) / 6.6e-12 Eh (GFN1), dG max 3.5e-6 / 5.2e-6
  Eh/A; the device M was uploaded and the resident loop ran in every run (verbosity-2 lines). SCF iterations ROCm = CPU in
  15 of 16 pairs (19 of 20 with the charged and open-shell runs below); PR07 GFN2 8 against 10 (CPU at `-threads 1` 9, ROCm with `-scf_mixed_precision false` 8). PR07 is sensitive
  before any mixing: the energy of iteration 0 differs CPU vs ROCm by 5.5e-6 Eh and CPU 1 vs 10 threads by 3.3e-6 Eh without
  the preconditioner; not isolated further. Iterations none -> frag (CPU and ROCm alike): GFN2 urea3-water8 12/11 -> 6, water8
  8 -> 4, peo8 11 -> 6, caffeine 13 -> 4, AcCl 14 -> 5, Bz 8 -> 4, PR07 22 -> 10 (CPU) / 8 (ROCm), ureawater197 unrelaxed
  24 -> 33 (guard fired at 4 on both); GFN1 10 -> 7, 5 -> 3, 9 -> 7, 11 -> 4, 10/9 -> 4, 5 -> 3, 15 -> 5, 13 -> 8.
- peo201-chain (1410 atoms, segmented as above, M of 202 blocks, 20.4 MB GFN2 / 0.76 MB GFN1): dE ROCm vs CPU 1.1e-11 Eh
  (GFN2) / 1.0e-10 Eh (GFN1), dG 7.3e-8 / 4.6e-7 Eh/A, iterations 12 -> 8 (GFN2) and 10 -> 9 (GFN1) on both.
- Guard (unrelaxed urea-water 197, GFN2): `g2:4:best`, `g2:4:keep`, `g2:4:keep:best`, `g2:4`, `g1:3:best`, `off`: ROCm and CPU
  fire at the same iteration (4, 4, 4, 4, 2, none) with the same counts (33, 27, 27, 33, 23, 31) and dE <= 5e-13 Eh. GFN1 on
  this structure: no guard hit, 8 iterations each.
- Charged (hydronium-water2 +1) and open shell (caffeine +1 doublet): the host skips the preconditioner with its warning, the
  resident loop runs with the scalar seed, iterations as without the switch, GFN1 and GFN2.
- Without `-scf_precond`: caffeine, AcCl, Bz, PR07 and the other systems above, GFN1 and GFN2, plus peo201-chain: energy and
  iterations equal to the binary before WP5 (dE 0); gradients differ by <= 7.5e-16 Eh/A, the run-to-run scatter of one binary
  (2.1e-16 Eh/A on caffeine, atomic adds in the gradient kernels).
- Fault injection with frag (urea3-water8, unrelaxed urea-water 197, GFN1 and GFN2): `CURCUMA_ROCM_FAIL_RES=step` restarts
  the SCF on the host-driven loop with the scalar seed (frag is not applied there), energy and iterations equal to
  `-xtb.rocm_resident_loop false`, one notice; `begin` takes the host-driven loop with one notice and the "ignored" warning.
- `-md` urea3-water8 5 fs `-seed 42`, GFN2: SCF iterations per step 6 5 5 4 5 5 on CPU and ROCm, M reused 5 times on the
  device (tag path), MD energies equal in every printed digit.
- Timing (`-sp`, median of 3, alternating, same binary, `-threads 10`, segmented as above): peo201-chain GFN2 none 14.50 s
  -> frag 11.21 s (SCF 10.9 -> 7.0 s; M built in 0.6 s on the host), GFN1 19.16 -> 15.26 s (SCF 15.1 -> 10.9 s).
- `ctest -R "gpu|rocm|sqm|xtb|grad|scf_precond"`: 206 of 208 passed, failures `xtb_cpscf`, `cli_sqm_12_gfn2_scf_precond` (known).
- Not run: blocks above 6144 entries (the setup then fails and the scalar seed is used, as on CUDA), Hessian on ROCm with frag,
  MD longer than 5 fs, `-scf_precond_jacobian response` on ROCm, `-scf_precond_rebuild_rms` (re-upload on rebuild), solvent,
  multi-GPU, other cards; the frag path after a step fallback is not restored (scalar seed for the rest of the run).

### WP6 (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switch `-xtb.rocm_device_solvation`, default true (PARAM in `native_xtb_method.h`, fallback default in `xtb_hip_method.cpp`).
  `HipScfBackend` overrides `supportsDeviceSolvation` (needs the WP1 device potential), `beginSolvation`, and the new base
  virtuals `supportsDeviceSolvationGfn1` / `beginSolvationGfn1` (needs the WP3 GFN1 loop). The context has `solvBegin` (upload of
  the Born matrix B, limit 2 GB, `CURCUMA_ROCM_SOLV_MAX_MB` lowers it) and `solvApply` (v_at += B q_at by `rocblas_dgemv`,
  plus `rocblas_daxpy` of a constant shift), called in `potStepCore` after the D4 term (q_at = the multipole input q_at) and in
  `potStepCoreGfn1` after the third order (q_at = the third-order input q_at). A new potential begin or GFN1 loop begin turns
  the field off until the next upload.
- GFN1 uses q_solute = q_at + cm5 (`ALPBSolvation::soluteCharges`, private): the host passes B (`bornMatrix()`) and the shift
  B cm5 (`addPotential` at q = 0); the device sum B q + B cm5 differs from the host B (q + cm5) only by rounding. CUDA has no
  GFN1 path (base virtuals false), so CUDA GFN1 + solvent keeps its host-driven loop.
- Solvation energy and gradient stay on the host, as on CUDA and the CPU: the SCC energy of the loop does not contain the
  solvation energy on any path, `energySolvation()` runs post-SCF on the downloaded q_at, the post-SCF host potential rebuild
  includes `addSolvationPotential` (for the Pulay/W terms), and `ALPBSolvation::addGradient` adds dE_solv/dR.
- Host changes outside the ROCm files: the two base virtuals (`xtb_native.h`); in `xtb_native.cpp` the GFN1 resident-loop condition
  accepts a solvent when `supportsDeviceSolvationGfn1()` and B are available, and calls `beginSolvationGfn1` after the begin
  (failure keeps the host-driven loop). The GFN2 host path (`solv_blocks_device`, `beginSolvation`) is unchanged.
- Fallback: a failed upload takes the host-driven loop with the host solvation; a failed step returns that SCF to the host
  potential (WP1 fallback, or WP2/WP3 restart and then the WP1 fallback, the solvation failure is sticky per geometry); later
  solvated calculations stay on the host solvation; gas-phase calculations keep the device paths. One notice per method
  object, not counted in `gpu_fallback`. Fault injection `CURCUMA_ROCM_FAIL_SOLV=begin|step` (`step` from the third call) on
  caffeine and caffeine +1, GFN1 and GFN2: energy equal to the switch off (dE 0), dG <= 1.8e-16 Eh/A, equal iterations, one notice.
- Acceptance (`-sp -threads 10 -dump_gradient`, n = 1 each, CPU = `-gpu none`): caffeine, caffeine +1 doublet, acetate -1, AcCl,
  Bz, CO, ED01, PR07, ED24 x {water ALPB, water GBSA, methanol ALPB, toluene GBSA, hexane ALPB} x GFN1/GFN2 (90 triples):
  dE on vs CPU max 1.2e-11 Eh (GFN2, ED01) / 2.7e-11 Eh (GFN1, PR07; off vs CPU the same), dG max 4.6e-6 / 4.8e-6 Eh/A (ED01),
  on vs off max 2.3e-13 Eh, 3.9e-14 Eh/A; SCF iterations on = off in all 90, |on - CPU| <= 1; the device path (verbosity-2
  line) and the resident loop ran in every on run. peo201-chain water ALPB: dE on vs CPU 3.1e-11 Eh (GFN2) / 7.3e-11 Eh (GFN1),
  on vs off 2.6e-11 / 6.7e-11 Eh (the device energy sums of WP2/WP3), dG 7.1e-8 / 4.9e-8 Eh/A, 15 / 9 iterations on all
  paths; with a 10 MB limit (B is 16 MB) the upload fails, one notice, energy equal to off.
- Gas phase against the binary before WP6 (caffeine, AcCl, Bz, PR07, GFN1 and GFN2): energies and iterations equal, gradients
  within 4.7e-16 Eh/A (the run-to-run scatter of the atomic adds, see WP5).
- `-opt` caffeine water 20 steps (lbfgs) and `-md` 50 fs `-seed 42`, GFN1 and GFN2, `-scf_mixed_precision false`: on and off agree
  in every printed digit except the wall-time columns; final structure and trajectory identical up to the sign of zero.
- Timing (`-sp` water, median of 3, alternating, same binary, `-threads 10`): peo201-chain GFN2 17.39 -> 16.96 s, GFN1
  16.55 -> 16.42 s (the eigensolve still dominates, WP7).
- `ctest -R "gpu|rocm|sqm|xtb|grad|solv"` with default true: 238 of 240 passed, failures `xtb_cpscf`, `cli_sqm_12_gfn2_scf_precond`
  (known). No ctest covers ROCm with a solvent.
- Not run: mixed precision on with a solvent in opt/MD, solvents beyond the five listed (dmso, acetone, chloroform not run on ROCm),
  `-xtb.solvent_epsilon` without a named solvent, CPCM (not native), solvated `-scf_precond frag`, metal open shells in solvent,
  systems above the 2 GB limit, multi-GPU, other cards.

### WP7 (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switches: the generic `-scf_pseudo_diag`, `-scf_pseudo_diag_max_dq` (0.3), `-scf_pseudo_diag_max_steps` (12),
  `-scf_pseudo_diag_fp64` (false), read in `xtb_hip_method.cpp` and passed to `XtbHipContext::setPseudoDiagonalisation`.
  **ROCm default: off** (fallback `kRocmPseudoDiagDefault = false`; the PARAM default true is the CUDA one), see the
  acceptance below. Verbosity 2 names the path and prints the counts of the last resident SCF (`steps taken`, `rejected
  by the gate`). Trace `CURCUMA_PSEUDO_TRACE=1`, fault injection `CURCUMA_ROCM_FAIL_PSEUDO=begin|step`.
- Context (`rocm/xtb_hip_context.hip`): kernels `k_pseudo_coldot`, `k_pseudo_gap`, `k_pseudo_theta`, `k_pseudo_gram_diag`
  (block trees, no warp-size assumption), `pseudo_blas` overloads, `Impl::pseudoRotateOrthonormalise`, `pseudoDiagStep`
  (FP32), `pseudoDiagStep64`, `eigensolveStep` (gate, then full solve and storage of the eigenvectors; called by
  `fockAndSolve`), `psResetLoop` at both loop begins, `ps_last_dq` set at the end of `resStep`, the gradient builds W =
  C_o (2 G_oo) C_o^T from a rebuilt F after a pseudo step (CUDA `w_from_goo`). Same gates as CUDA: inside a resident step,
  closed shell, previous max|dq| < max_dq, at most max_steps in a row, gap > max(30 kT, 1e-3 Eh), all |theta| < 0.1.
- Differences from CUDA: eigenvectors are stored only by full solves inside the resident loop; the FP64 step keeps F
  (Tb = dS, M = dCw), so a failure after the rotation falls back to the full solve instead of failing the step; any failure
  or failed allocation turns the path off for the rest of the run with one notice. L is the probed dL of the device
  integral build (else a private factor), independent of `-xtb.rocm_cached_cholesky`. Memory: about 2 n^2 floats more
  (previous C and FP32 L; 83 MB at nao 3222), FP64 steps one n^2 double more.
- Acceptance (`-sp -threads 10 -dump_gradient`, n = 1 each; on = `-scf_pseudo_diag true`, off = `false`, CPU `-gpu none`):
  caffeine, caffeine +1, 13 MOR41 (12 + PR24), urea3-water8, the 16 open-shell GMTKN55 structures of WP4, GFN1 and GFN2
  (64 triples): dE on vs off max 1.9e-11 Eh, on vs CPU 1.9e-11 Eh, dG on vs off max 2.9e-6 Eh/A; iterations on vs off
  differ by 1 in 3 cases (ED01, urea3-water8 GFN2: on = CPU; PR16 GFN1 12 vs 11). Open shell: no pseudo step, on = off
  exactly (dE 0, dG <= 4e-16 Eh/A). Pseudo steps on closed shells 0 to 9 per SCF (PR24, PR25, ED07 GFN2: 0).
- peo201-chain: GFN2 dE on vs off 2.3e-12 Eh, dG 1.6e-7 Eh/A, 12 = 12 iterations, 5 of 12 steps pseudo; GFN1 2.7e-12 Eh,
  1.4e-7 Eh/A, 10 = 10, 3 of 10. Only FP32 steps qualify, and the FP32 phase ends at `scf_fp32_threshold`, so the
  pseudo steps are a minority of the iterations on this system.
- Water ALPB (caffeine, caffeine +1, PR07, peo201-chain, GFN1/GFN2): dE on vs off max 3.5e-11 Eh, dG 2.5e-6 Eh/A,
  iterations equal. `-scf_precond frag` (segmented as WP5; urea3-water8, caffeine, peo201-chain): dE max 4.7e-11 Eh,
  dG 5.7e-6 Eh/A, but iterations +2 for caffeine (4 -> 6, GFN1 and GFN2) and +1 for urea3-water8 GFN2 and peo201-chain
  GFN1. Caffeine GFN2 trace: the pseudo step in iteration 2 gives max|dq| 3.6e-4 (full solve 4.1e-4), the next step
  1.7e-4 against 3.4e-6: the one-rotation map evaluation is inexact at O(theta^2) and the fast preconditioned Broyden
  takes that error into its secant history. Not measured on CUDA.
- `-scf_pseudo_diag_fp64 true -scf_mixed_precision false` (caffeine, AcCl, PR07, ED24, urea3-water8, peo201-chain GFN2;
  the first five GFN1): dE on vs off max 7.6e-12 Eh, dG 1.9e-6 Eh/A (W from G_oo), peo201-chain 10 of 12 steps pseudo,
  iterations equal except AcCl GFN1 (10 vs 9). Not timed.
- `-opt` caffeine 30 steps and `-md` PR24 100 fs `-seed 42`, GFN1/GFN2: with `-scf_mixed_precision false` on and off
  agree in every printed digit (no pseudo step runs there); with mixed precision the opt differs from step 1 in the 6th
  digit of the energy change and gradient norm, the MD temperature in the 6th (GFN1) / 7th (GFN2) significant digit, the
  pattern of the FP32 path seen in WP1-WP4.
- Fault injection (caffeine, GFN1/GFN2): `begin` gives energy equal to off (dE 0, same iterations), `step` (from the
  third attempt) 2 pseudo steps then full solves, dE <= 3e-12 Eh against off; one notice each, no SCF restart.
- Timing (median of 3, alternating, same binary, `-threads 10`): peo201-chain `-sp` GFN2 14.11 -> 12.35 s (SCF 11.0 ->
  9.2 s), GFN1 18.60 -> 16.90 s; PR24 `-md` 100 fs GFN2 54.10 -> 54.08 s (no pseudo step on PR24).
- `ctest -R "gpu|rocm|sqm|xtb|grad"` with the ROCm default off: 206 of 208 passed, failures `xtb_cpscf`,
  `cli_sqm_12_gfn2_scf_precond` (known); no ctest turns the ROCm pseudo step on.
- Not run: d-shell metals beyond MOR41, Tele = 0 with the pseudo step, systems above peo201-chain (VRAM budget beyond
  nao 3222 not measured, allocation failure only by injection), FP64 pseudo timing, opt/MD with the FP64 pseudo step,
  multi-GPU, other cards.


### WP8a (Oct 24, 2026; AI-generated, machine-tested on gfx1200 only)

- Switches (PARAM in `native_xtb_method.h`, fallback defaults in `xtb_hip_method.cpp`), default true since the acceptance
  below: `-xtb.rocm_d4_postscf_device` (GFN2 D4), `-xtb.rocm_d3_device` (GFN1 D3(BJ)), `-xtb.rocm_cn_device` (Hamiltonian CN
  download and D4 covalent CN). `HipScfBackend` overrides `d4PostScfDevice`, `d3DispersionDevice`, `d4CovalentCnDevice`,
  `downloadCn` with the CUDA contracts; the host call sites (`calcDispersionEnergy`, `primeD4CnFromDevice`, the integral
  block of `xtb_native.cpp`) are unchanged. Context (`rocm/xtb_hip_context.hip`): `d4PostScf`, `d3Dispersion`,
  `d4CovalentCn`, `dlCn`, kernels `k_d4_cn`, `k_d4_cn_fold`, `k_d4_atm_c6mat`, `k_d4_nb_count`, `k_d4_nb_fill`,
  `k_d4_atm_nl`, `k_d4_atm_nl_b`, `k_d4_add2`, `k_d3_cn`, `k_d3_weights`, `k_d3_pairs`, `k_d3_cn_fold`; `k_res_d4_refw`
  (WP2) got an optional dW/dCN output (host order `expd/norm - gwk (dnorm/norm)`, non-finite -> 0).
- The post-SCF D4 reuses the `beginDispersion` data (element C6 block, sqrt r4r2, geometry), the 2-body kernel `k_d4_grad`
  and the WP2 reference tables of the same geometry (`rw_geom_ok`: set by `beginDispersionWeights`, cleared by
  `beginDispersion`). After a host-driven SCF (no tables for that geometry) the call declines and the host code runs, with a
  verbosity-2 line and no notice. A device failure turns that path off for the rest of the run with one warning per
  method object (not counted in `gpu_fallback`); fault injection `CURCUMA_ROCM_FAIL_POSTSCF=d4|d3|cn` (every call of that
  path). Verbosity 2 names the path of every call (`ROCm: ... on the device|host`), verbosity 3 the ATM kernel time.
- Differences from CUDA: weights follow the host underflow rule (gwk on the reference with the largest CN; CUDA
  `k_d4_build_refw` sets 0); ATM uses `pow` (CUDA: `q^5 cbrt(q)` for alp = 16); ATM runs one block of 64 threads per atom
  with a fixed shared-memory tree (`k_d4_atm_nl_b`, CUDA: warp per atom or thread per atom); `CURCUMA_ROCM_ATM_THREAD=1`
  selects the thread-per-atom kernel with the summation order of the host loop. The neighbour-list offsets are scanned on the
  host (nat ints down and up). The ATM C6 matrices (2 nat^2 doubles) stay allocated between calls.
- Acceptance (`-sp -threads 10 -dump_gradient`, n = 1 each; off = all three false, CPU = `-gpu none`): caffeine, caffeine +1,
  fi-pyridine, 13 MOR41 (12 + PR24), the 16 open-shell GMTKN55 structures of WP4, GFN1 and GFN2 (64 structure/method
  pairs), configurations off, each switch alone, all on: on vs off max 8.5e-14 Eh (GFN2, PR25; GFN1 0), dG on vs off
  1.2e-14 Eh/A; D4 alone and D3 alone give dE 0 (below the resolution of the total energy, about 1e-14 Eh) and dG <= 6.2e-16
  Eh/A, so the dispersion energy and gradient agree with the host code to that level; the residue of "all on" comes from
  the CN switch (device CN in the multipole damping radii). On vs CPU max 8.5e-12 Eh, 2.7e-6 Eh/A (as off). SCF iterations
  equal in all ROCm configurations; the device path ran in every on run (verbosity-2 lines). peo201-chain: GFN2 on vs off
  0 Eh, 1.3e-14 Eh/A, on vs CPU 2.3e-12 Eh; GFN1 on vs off 4.6e-13 Eh, 4.4e-16 Eh/A, on vs CPU 1.1e-10 Eh (off 1.1e-10).
  Water ALPB (caffeine, caffeine +1, PR07, fi-pyridine, peo201-chain, GFN1/GFN2): on vs off max 1.4e-12 Eh (peo201-chain
  GFN2), dG 1.6e-14 Eh/A, on vs CPU <= 7.2e-11 Eh. Off against the binary before WP8a (caffeine, caffeine +1, AcCl, Bz, PR07,
  GFN1/GFN2): dE 0, dG <= 4.1e-16 Eh/A (the WP2 weight kernel is unchanged for the loop).
- Fault injection d4, d3, cn (caffeine, GFN1/GFN2): one warning, energy equal to off (dE 0).
- `-opt` caffeine 30 steps and `-md` PR24 100 fs `-seed 42`, GFN1/GFN2, `-scf_mixed_precision false`: on and off agree in
  every printed digit of the optimisation steps and the MD table (except the wall-time column), the trajectories are
  identical, final structures differ only in the sign of zero.
- Timing (median of 3, alternating, same binary, `-threads 10`, all on vs all off): peo201-chain `-sp` GFN2 14.33 -> 13.41 s
  (post-SCF energy block 1.98 -> 1.11 s; D4 incl. ATM and CN fold 1.24 s -> 0.37 s, of which the ATM kernel 0.36 s against
  1.22 s for the host-fed `k_d4_atm`; the thread-per-atom list kernel took 1.34 s), GFN1 19.33 -> 18.64 s (post-SCF energy
  block 2.15 -> 1.39 s); PR24 `-md` 100 fs GFN2 54.15 -> 51.29 s, GFN1 37.62 -> 37.67 s (no change within the scatter).
- `ctest -R "gpu|rocm|sqm|xtb|grad|d4|d3|disp"` with the defaults true: 207 of 209 passed, failures `xtb_cpscf`,
  `cli_sqm_12_gfn2_scf_precond` (known).
- Not run: CUDA comparison (no CUDA toolchain here), d-shell metals beyond MOR41 and PR07, ATM above peo201-chain (the C6
  matrices are 2 nat^2 doubles: 1.6 GB at 10000 atoms, allocation failure then falls back), the host-driven SCF with the
  switches on beyond caffeine (`-xtb.rocm_resident_loop false`: D4 declines with the verbosity-2 line, CN on the device,
  energy equal to all off, dG 2.1e-16 Eh/A), multi-GPU, other cards. The CPU reference and the GFN2 D4 `m_d4_variational_qresp`
  path were not varied.

### WP8b (Oct 25, 2026; AI-generated, machine-tested on gfx1200 only)

- Switches (PARAM in `native_xtb_method.h`, fallback defaults `kRocm*Default` in `xtb_hip_method.cpp`), default true since the
  acceptance below: `-xtb.rocm_grad_resident` (A: `prepareConvergedPotentials`/`gradientResident`),
  `-xtb.rocm_grad_tail_device` (B: `gradientTailDevice`, `repulsionEnergyDevice`), `-xtb.rocm_defer_host_integrals` (C: lazy
  host copies). `HipScfBackend` overrides the four CUDA virtuals with the CUDA contracts plus `downloadOverlapInto`/`downloadH0Into`
  and the new base virtual `deferHostCopies` (default false). Context: `convPotentials` (output charges of the loop -> D4 weights
  by `k_res_d4_refw` -> `potBuildCore` / `potBuildGfn1`, the potential halves split off `potStepCore` / `potStepCoreGfn1` with an
  unchanged kernel sequence), `gradient(..., resident)`, `gradTail`, `repEnergy`, kernels `k_rep_energy`, `d_mp_grad_pair`,
  `k_grad_tail_mp`, `k_grad_tail_dcn`, `k_grad_tail_cn` (ports of CUDA, one thread per atom).
- Host changes outside the ROCm files (`xtb_native.h/.cpp`, CUDA path unchanged because `deferHostCopies` is false there): the
  download of the GFN2 multipole integrals is skipped where the existing deferral (`m_mp_ints_deferred`) applies, and the post-SCF
  host potential rebuild is skipped on the resident loop while S/H0 are still deferred (it only fed the debug `m_F`), not with
  `d4_charge_source cpscf`. `defer_amat` was already active on ROCm since WP1/WP2; the post-SCF rebuild and the gradient's host
  rebuild built the 18 nat^2 host matrices anyway, so C and A are what keep them off the host.
- Differences from the host: repulsion summed per atom over j < i, then in atom order; the gradient tail kernels are gathers that
  evaluate each pair from both sides (CUDA the same). `convPotentials` declines with an active device solvent (the host gate
  `!m_solvation` never asks then); a failure inside `gradientResident` takes the host gradient of that geometry (CPU loops).
- Lazy host consumers checked, all on against all off (dE <= 3.6e-14 Eh, dG <= 1.0e-14 Eh/A, iterations equal): host-driven loop after `CURCUMA_ROCM_FAIL_RES=step|begin`
  (closed and open shell), `CURCUMA_ROCM_FAIL_POT=begin|solve`, `CURCUMA_ROCM_FAIL_SOLV=step`, `-scf_precond frag` (urea3-water8),
  water ALPB, `-xtb.rocm_resident_loop false`, `d4_charge_source eeq|cpscf`, `-hessian` AcCl (15 frequencies equal to off, 0.1 cm-1
  from CPU, the printed digit). Fault injection `CURCUMA_ROCM_FAIL_GRAD=resident|tail|defer` (caffeine, GFN1/GFN2): one warning,
  dE 0 and dG <= 3.6e-16 Eh/A against off (defer seen only with A off, where the gradient's host rebuild needs S/H0).
- Acceptance (`-sp -threads 10 -dump_gradient`, n = 1 each; configurations off, A, B, C, all on; CPU `-gpu none`): caffeine,
  caffeine +1, fi-pyridine, 13 MOR41, 16 open-shell GMTKN55, GFN1 and GFN2 (64 pairs): on vs off max 1.4e-14 Eh (B, ED07), dG
  8.2e-16 Eh/A; on vs CPU 8.5e-12 Eh, 2.7e-6 Eh/A (as off); SCF iterations equal in all ROCm configurations; path lines in every on
  run. peo201-chain: GFN2 on vs off 1.4e-12 Eh (from B, the repulsion summation; A and C give 0), dG 1.1e-15 Eh/A, on vs CPU 9.1e-13
  Eh; GFN1 0 Eh, 5.6e-16 Eh/A, on vs CPU 1.1e-10 Eh (off the same). peo201-water500 (2910 atoms): GFN2 3.6e-12 Eh, GFN1 9.1e-13 Eh,
  dG <= 1.2e-15 Eh/A (no CPU reference run).
- `-opt` caffeine 30 steps and `-md` PR24 100 fs `-seed 42`, `-scf_mixed_precision false`, GFN1/GFN2: all printed digits of the
  optimisation and MD tables equal except the time column; trajectories identical; final structures differ in the sign of zero.
- Peak host RSS (`-sp` with gradient, off -> C -> all on): peo201-chain GFN2 2172 -> 1376 -> 914 MB, GFN1 1602 -> 1428 -> 1101 MB;
  peo201-water500 GFN2 6700 -> 2854 -> 681 MB, GFN1 3786 -> 2561 -> 827 MB (above nao ~5270 the converged P/C stay on the device
  once the multipole integrals are deferred, as on CUDA). The caffeine baseline of the process is about 450 MB.
- Timing (median of 3, alternating, all on vs all off, `-threads 10`): peo201-chain `-sp` GFN2 13.61 -> 12.42 s, GFN1 18.67 -> 17.30 s;
  with gradient GFN2 14.04 -> 12.86 s, GFN1 19.25 -> 17.87 s; PR24 `-md` 100 fs GFN2 51.27 -> 49.65 s. peo201-chain GFN2 phases: setup
  1117 -> 793 ms (no integral download), post-SCF energy 1111 -> 422 ms (no host rebuild, Fock, S/H0 host build).
- Not run: CUDA comparison, the host gradient fallback after a real device-gradient failure (only by the injection in
  `convPotentials`), solvated runs with A (host gate keeps the host rebuild), Vulkan, multi-GPU, other cards.

### WP9 (Oct 25, 2026; AI-generated, machine-tested on gfx1200 only)

GFN-FF on ROCm, rules as above but reference = CPU `-method gfnff -gpu none` and ROCm with the three switches off.
- A, device repulsion list rebuild: trait `has_device_repulsion_build = true` (`gfnff_hip_method.h`), so the wrapper honours
  `-gfnff.gpu_rep_pairs_on_device` on ROCm. `FFWorkspaceHip::generateNonbondedRepulsionOnGPU` (`rocm/gfnff_rocm.hip`) with
  `k_rep_rows_count`, `k_rep_rows_build`, `k_rep_pairs_static` and `k_disp_csr_lower` for the lower CSR half: pairs in
  ascending (i, j), gather CSR built on the device (CUDA: arrival order, no CSR). Fallback: host build from the same
  tables and the host coordinate mirror (`repulsionHostBuildFromTables`), counted, host build for the rest of the run;
  fault injection `CURCUMA_ROCM_FAIL_REPDEV=1`; test dump `CURCUMA_ROCM_REP_DUMP=<file>` (device and host lists, sorted).
  Only MD/opt can see it: the setup list is the host build, `-sp` never rebuilds.
  MD 200 steps `-dt 1 -seed 42` with dumps, on vs off: peo201-chain 15 rebuilds, urea400-water1000 27 rebuilds, every
  rebuilt list identical to the host list of the same step (pair set and alpha/repab/r_cut bitwise); energies per step
  on vs off max 2.3e-13 Eh (peo201-chain) and 6.0e-12 Eh (urea400-water1000, from step 102; the gather sums each row in
  ascending partner order, the host list in cell-list order).
- B, tabulated dc6dcn: switch `-gfnff.rocm_dc6dcn_tables` (method options `RocmGfnffOptions`, read by
  `GFNFFHipComputationalMethod` and published to the workspace constructor through `rocmGfnffOptionsSlot()`), the CUDA
  environment name `CURCUMA_DC6DCN_TABLES=0|1` overrides it. Kernels `k_d4_partner_tables`, `k_dc6dcn_per_pair_tab`
  (port, compact types `ctype_elem`/`atom_ctype`); with tables the fused gather is off. Fault injection
  `CURCUMA_ROCM_FAIL_DC6TAB=1` (table allocation fails, counted, uv path kept). peo201-chain `-sp` gradient: dE tab vs off
  6e-14 Eh, dG 1.0e-16 Eh/A. MD 200 steps, median ms/step of 3 runs (fused uv = current default / uv per pair
  (`CURCUMA_ROCM_DISP_FUSED=0`) / tab): peo201-chain 6.11 / 6.14 / 6.18, urea400-water1000 43.0 / 47.0 / 47.7, PR24
  0.999 / 0.921 / 0.912 (bimodal 0.91 or 1.0 in every configuration). Choice: the fused uv gather stays the default
  (`rocm_dc6dcn_tables` false); the tables are slower than the uv per-pair kernel at 6200 atoms.
- C, device Schur step of the dense EEQ: switch `-gfnff.rocm_device_schur` (method options as B), default true since the
  MD checks below. `eeqDenseSchurHip` (`rocm/eeq_solver_hip.hiph`) behind `solveWithDeviceRHSAndGPUSchur` (one fragment)
  and `...SchurGeneral` (several): the unchanged dense core (`eeqBuildFactorSolveDevice`, split off `eeqBuildFactorSolve`),
  then `k_eeq_schur_sums1_hip` / `k_eeq_schur_sums_hip` (one thread per sum, atoms ascending, the host order),
  nfrag + nfrag^2 numbers to the host, the wrapper's Gaussian elimination (`eeqSchurLambdaHost`), `k_eeq_schur_apply_hip`
  (multiply and subtract rounded separately, as the host build without FMA). The small-system kernel (N <= 128) keeps
  priority. Failure (or `CURCUMA_ROCM_FAIL_SCHUR=1`): counted once, host Schur step for the rest of the run.
  The trait `has_device_schur` stays false: the wrapper calls both entry points whatever the trait says; the trait
  only adds the CUDA-only lu/ldlt/pcg/batched branches (no ROCm functions, `EEQSolver::gpuSupports` rejects those names
  for rocm) and would stop the wrapper reading `eeq_rocm_cpu_fragment_threshold`. The wrapper's verbosity-1 path text
  still says "WP2 dense Cholesky + CPU Schur" (template string); the ROCm solver prints its own verbosity-2 line.
  `-sp` never reaches the device EEQ (the first call reuses the CPU charges), so C was checked in MD.
- Found and fixed while testing C (affects the dense path with the switch off as well): the dense matrix build ran on the
  solver's own blocking stream without an order against the workspace streams that upload the coordinates. MD of the
  water dimer with `CURCUMA_ROCM_EEQ_SMALL=0`: 1 of 5 runs of the binary before WP9 differed by 1e-4 Eh from the first
  step on; after the fix (`ev_order` in `eeqBuildFactorSolveDevice`: event on the single-stream main stream or on the
  null stream, then `hipStreamWaitEvent`) 10 of 10 runs agree for off and for on.
- C acceptance (MD 30 steps `-dt 1 -seed 42`, on vs off, energies per step): 15 cluster structures (charged ones with
  `-gfnff.frag_charge_model reference`), urea-water-cluster197 with `-gfnff.eeq_rocm_cpu_fragment_threshold 0` (39
  fragments), peo201-chain with `-eeq_solver.eeq_ppcg_min_nfrag 0` (dense, 1410 atoms): with the small-system kernel off
  max |dE| 1.7e-13 Eh (peo201-chain; clusters <= 4.3e-14 Eh), device Schur in every step; with the small kernel on the
  device Schur runs for N > 128 or nfrag > 8, max 1.4e-13 Eh. Fault injection: one counted notice, energies as off.
- C timing (MD 200 steps, median ms/step of 3): urea-water-cluster197 CPU fragment route (default, threshold 16) 0.98,
  dense device route off 2.99, on 2.90; macrocycle-bgal (231 atoms, 2 fragments) 2.94 / 2.96; peo201-chain exact
  19.5 / 19.3. The threshold 16 stays: at 197 atoms and 39 fragments the CPU route is 3x faster than the dense device route.
- Acceptance of all three (`-sp -gradient -dump_gradient -threads 10`, n = 1 each; CPU `-gpu none`; base = ROCm binary
  before WP9; off = the three switches false; on = all three true; def = defaults A on, B off, C on): all 95 MOR41 and
  every 8th GMTKN55 `struc.xyz` by sorted path (308). MOR41: dE vs CPU max 1.3e-10 Eh (PR21) for base, off, on and def
  alike, 0 of 95 differ at 8 decimals, dG vs CPU max 2.1e-9 Eh/A; on vs off and def vs off max 1.8e-15 Eh, 4.2e-17 Eh/A;
  off vs base 1.8e-15 Eh. GMTKN55 sample: dE vs CPU max 9.9e-14 Eh (all four), 0 differ at 8 decimals, dG 2.4e-13 Eh/A;
  on vs off 8.9e-16 Eh, 2.8e-17 Eh/A. 9 charged multi-fragment structures run on the CPU method (existing counted fallback,
  base the same). `-sp` reaches B only (A and C need geometry steps), hence the MD checks above.
- MD 200 steps `-dt 1 -seed 42`, energies per step against off (all three switches): peo201-chain def 2.8e-13, on 2.8e-13 Eh;
  urea400-water1000 def 3.6e-12, on 4.7e-12 Eh; PR24 2.5e-14 Eh. Rebuilds happened (diagnostics `gpu_rep_update`), 15 and 27.
  A fault injection in MD (peo201-chain, 60 steps): one notice, 4 host builds, energies as off within 1.7e-13 Eh.
- `ctest -R "gfnff|gpu|rocm|grad"` with the defaults: 115 of 115 passed (`xtb_cpscf` and `cli_sqm_12` are not in this set).
- Timing (median of 3, alternating; MD median ms/step from the diagnostics; off / def / on): peo201-chain 6.19 / 6.17 / 6.24,
  urea400-water1000 43.2 / 43.1 / 47.7 (on includes the tables), PR24 1.00 / 0.98 / 0.93 (bimodal 0.91 or 1.0 in every
  configuration, no rebuild below `nb_cell_list_min_atoms`); repulsion refresh per rebuild step: peo201-chain 3.66 -> 0.76 ms,
  urea400-water1000 21-31 -> 7.6 ms. `-sp` wall peo201-water500 (2910 atoms) 0.57 s in all three (no switch acts in `-sp`).
- No edit to the shared template (`gfnff_gpu_method_impl.h`) or the CUDA traits. ROCm class `GFNFFHipComputationalMethod`
  now has its own constructor and overrides `setMolecule`, `calculateEnergy`, `setParameters` (options scope).
- Not run: CUDA (no toolchain), multi-GPU split and its mixed-precision EEQ (one GPU here), fixed-point accumulation (WP10),
  implicit solvent on the device, PBC with the device repulsion rebuild, `-opt` runs, MD longer than 200 steps, systems
  above 6200 atoms (the row build is O(N^2) per rebuild), other cards.

### WP10 (Oct 25, 2026; AI-generated, machine-tested on gfx1200 only)

Run-to-run reproducibility of the ROCm GFN-FF step and of the ROCm GFN1/GFN2 gradient, port of the CUDA fixed-point
scheme ([GPU_REPRODUCIBILITY.md](GPU_REPRODUCIBILITY.md)). Reference = CPU `-gpu none`; before = ROCm binary before
WP10; off/on = the new binary with the switch false/true. Runs `-threads 10`, a fresh directory per run; sp = `-sp
-gradient -dump_gradient` (17 digits, energy 15 decimals), MD = `-dt 1 -seed 42`, fingerprint = md5 of the trajectory
plus the full-precision averages, velocities, geometry and gradient of `<basename>.final.json`.
- Audit, GFN-FF (`rocm/gfnff_rocm.hip`): order-dependent FP64 atomics in `add_grad` (every term kernel, the gathers and
  `k_slot_gather` add per atom into the shared gradient), `blockReduceAddEnergy` (one atomic per block per energy
  slot), dE/dCN (`k_dispersion*`, `k_bonds`), `zz_hb` (`k_bonds`), CN sums (`k_cn_compute_pairs`), HB CN
  (`k_hb_cn_per_atom`), the shared partner gradients of `k_coulomb_tiles_hip`, `k_coulomb_self` (no launch site).
  The gathers of table B row 8 sum each atom in a fixed order, but their per-atom result still went through an FP64
  atomic. `eeq_solver_hip.hiph` has no floating-point atomics (segmented fragment sums, integer LDS counters).
- Before, n = 10 each, all runs differ: sp distinct dumps 10/10 for caffeine, PR24, peo201-chain, urea400-water1000
  (max |dE| 0 / 0 / 5.7e-14 / 2.3e-13 Eh, max |dG| 3.5e-17 / 2.1e-17 / 7.6e-17 / 2.2e-16 Eh/A); MD PR24 100 fs 5
  distinct of 10, peo201-chain and urea400-water1000 200 fs 10 of 10 (spread of <Epot> 0 / 2.3e-13 / 1.0e-12 Eh).
- Change: switch `-gfnff.rocm_fixed_point` (PARAM in `gfnff.h`, `RocmGfnffOptions::fixed_point`, read by
  `GFNFFHipComputationalMethod`), default true, device constant `d_gfnff_fix` (`setGfnffFixModeHip`, uploaded when it
  changes). `fix_add`, `acc_add`, `blockReduceAddEnergyFix`, `k_fix_to_double`, `k_fix_copy` with the CUDA scales
  (forces, dE/dCN, `zz_hb` and the dispersion energy 2^-52, other energy slots 2^-48 per thread value, CN 2^-56).
  Conversions as CUDA: `zz_hb` before `k_hb_alpha_chainrule`, dE/dCN before `k_coulomb_postprocess` /
  `k_cn_chainrule`, gradient after the chain rule, the energy slots before the download, CN and HB CN before their
  first reader; the verbosity-3 snapshots through `k_fix_copy`. Two differences from CUDA: the chunk loop of
  `fix_add` for |v| >= 64 is one integer (same accumulator bits modulo 2^64, bounded cost), and a non-finite or
  out-of-range contribution sets `d_gfnff_fix_bad`, the step then returns NaN with one warning (CUDA wraps silently).
  `k_coulomb_tiles_hip` with TI == TJ (N >= 2000) walks the j tile in lock step (one barrier per k, every thread on
  its own slot) and adds the partner forces in k order without atomics; the smaller tiles use fixed-point shared
  atomics. With the switch off every kernel takes the old FP64 atomic path.
- After, switch on, n = 10 each: one distinct fingerprint in all seven cases (the four sp, the three MD). Switch off,
  same binary: 10 distinct of 10 (caffeine, peo201-chain, urea400-water1000 sp), spread as before. Three concurrent
  processes (`-threads 4`): peo201-chain MD 12 of 12 and urea400-water1000 MD 12 of 12 equal to the solo run (off: 11
  distinct of 12); the CPU thread count did not change the MD result (4 and 10 threads, same md5). Other paths, MD
  100 fs, n = 5 each, one fingerprint per case: peo201-water500 default, with `CURCUMA_ROCM_REP_GATHER=0
  CURCUMA_ROCM_HB_GATHER=0 CURCUMA_ROCM_DISP_FUSED=0`, with `-gfnff.rocm_dc6dcn_tables true -gfnff.gpu_coulomb_implicit
  false` (off: 5 of 5 distinct); urea-water-cluster197 (CPU fragment EEQ), water dimer (single stream, small EEQ
  kernel), macrocycle-bgal (these three with the binary before the lock-step change, which acts only at N >= 2000).
- Accuracy against the CPU, switch on, n = 1 per structure: caffeine dE 8.0e-15 Eh, dG 1.5e-13 Eh/A; 8 MOR41 (AcCl AcI
  Bz CO ED01 ED07 PR07 PR24) max dE 1.9e-13 Eh (PR07), max dG 8.0e-14; peo201-chain 1.3e-12 / 3.9e-14; urea400-water1000
  7.8e-12 / 5.3e-14 (off: 1.8e-11 / 5.1e-14 over 10 runs). MD averages on vs before: <Epot> 1.9e-13 (PR24), 1.6e-12
  (peo201-chain), 1.0e-11 Eh (urea400-water1000). The energy carries one rounding of up to 2^-49 Eh per thread value.
- Timing (median of 3 alternating, MD ms/step = median `step_total` of `-md_diagnostics_timing`, steps 10-200; before /
  off / on): urea400-water1000 43.35 / 43.38 / 45.60 with fixed-point shared atomics in the 128 x 128 Coulomb tile,
  43.21 (off) / 44.00 (on, +1.8 %) after the lock-step walk; peo201-chain 6.19 / 6.20 / 6.31 (+1.8 %); PR24 0.987 / 0.998
  / 0.874 (integer atomics are native on gfx1200, FP64 atomics a compare-and-swap loop); peo201-chain `-sp` wall 1.18 /
  1.19 / 1.19 s. Default on by the rule of the work package (slowdown on urea400-water1000 MD <= 5 %).
- GFN1/GFN2 (`rocm/xtb_hip_context.hip`): the SCF loop has no floating-point atomics (the resident loop gathers charges
  and moments, `k_res_pc_blocks` adds at most one non-zero value per element). FP64 atomics were left in `k_grad_repulsion`,
  `k_grad_cn_onsite`, `k_grad_h0_pulay`, `k_grad_coulomb` (gradient, dE/dCN) and `k_multipole_moments` (host-loop
  moments). Switch `-xtb.rocm_fixed_point` (PARAM in `native_xtb_method.h`, default true, `XtbHipContext::setFixedPoint`):
  `xfix_add` (2^-52) and `k_xfix_to_double` before the download; an out-of-range or non-finite contribution
  (`xfix_bad`) fails the device call, so the host path takes over. Before, n = 10: gfn2 caffeine and PR24, gfn1
  caffeine sp 10 distinct dumps each (energies equal at 15 decimals, dG spread <= 4.6e-16 Eh/A); gfn2 PR24 and gfn1
  caffeine MD 20 fs 10 distinct of 10. After: one distinct fingerprint in all five cases (n = 10). On vs CPU (gfn2,
  caffeine, the 8 MOR41, peo201-chain; gfn1 caffeine): max dE 5.5e-12 Eh, max dG 2.7e-6 Eh/A (ED01), identical in the
  shown digits for on and off (SCF convergence, not the accumulation); on vs before dG <= 2.0e-14 Eh/A. Timing (median of
  3, off / on): peo201-chain gfn2 `-sp -gradient` 13.82 / 13.81 s, PR24 gfn2 MD 100 fs 50.65 / 50.64 s.
- `ctest -j1 -R "gpu|rocm|sqm|xtb|grad|gfnff"`: 271 of 273 passed, failed `xtb_cpscf` and `cli_sqm_12` (known), 2
  disabled. The required set without `gfnff` (208 tests) gave 206 of 208 with the GFN-FF part alone.
- Still not order-independent: library reductions (rocBLAS `ddot`, `dsymv`/`ssymv`, `dgemm`; rocSOLVER `potrf`,
  `syevd`/`sygvd`) were not audited and the rocBLAS atomics mode is left at its default; repeated runs were identical
  in every case above. Integer atomics give list slots (`k_disp_pairs_build`, `k_generate_cn_pairs_write`, HB detection
  counters); every sum over these lists is now integer or ordered. `k_eeq_reduce_sums` and `k_eeq_reduce_fragment_sums`
  in `gfnff_rocm.hip` keep FP64 atomics but have no launch site.
- One unexplained event: in the first concurrent urea400-water1000 test (3 processes, binary with the CUDA chunk loop)
  one run hung the GPU and the reset killed the two others (3 of 9 runs); not reproduced in 33 concurrent runs after
  the bounded `fix_add` (21 before and 12 after the lock-step change), and 9 concurrent runs of the binary before WP10
  completed. Cause not determined.
- Not run: CUDA (no toolchain), PBC, `-opt`, MD longer than 200 steps, systems above 6200 atoms, GFN-FF with implicit
  solvent (CPU fallback anyway), the multi-GPU split, the range flags by fault injection, other cards.

### WP11 (Oct 25, 2026; AI-generated, machine-tested on gfx1200 only)

- Before: the WP1 device potential stores 15 nat^2 doubles of multipole interaction blocks (`dMpAmat`) with a limit of 4 GB
  (`HipScfBackend::kMpMaxBytes`, nat up to about 5770). Above it `beginPotentialStored` and the host upload `beginPotential` both
  failed, the host built the 18 nat^2 matrices and the SCF potential for the rest of the run (one warning), and without the
  device potential the resident loop and the WP8b gradient potentials were off. The host rule of `-gpu_multipole_otf auto`
  (`xtb_native.cpp`, 1 GB, nat > ~2635) asked for the on-the-fly path from 2635 atoms on, but ROCm did not offer one.
- Change: switch `-xtb.rocm_otf_multipole` (PARAM in `native_xtb_method.h`, fallback `kRocmOtfMultipoleDefault` in
  `xtb_hip_method.cpp`), default true. `HipScfBackend` overrides `supportsOnTheFlyMultipole` (needs the WP1 device potential),
  `beginPotentialOnTheFly` and the new base virtual `preferStoredMultipole(nat)` (`xtb_native.h`, default false), which keeps the
  stored path in the automatic choice while the 15 nat^2 doubles fit 4 GB. So with `auto` only systems above the ROCm limit change
  path; `-gpu_multipole_otf on` or `CURCUMA_GPU_MP_OTF=1` force the on-the-fly path for any size (the test route), `off` keeps the
  stored path. Host change in `xtb_native.cpp`: the automatic choice asks `preferStoredMultipole` (CUDA unchanged, base false).
- Context (`rocm/xtb_hip_context.hip`): `potBeginOtf` (shared begin with an `otf` flag: no `dMpAmat`, a resident stored block of
  the same context is released, geometry and damping radii in `dMpGeo`), kernels `k_mp_potential_otf`, `k_mp_energy_otf` and the
  element function `d_mp_elems15`, dispatch in `launchMpPotential` / `launchMpEnergy`, used by `potBuildCore` (host-driven device
  potential, resident step, WP8b `convPotentials`) and the energy of `resStep`. The gradient tail (`k_grad_tail_mp`, WP8b) already
  evaluated the pairs on the fly and needs no matrices. No atomics in the new kernels, so `-xtb.rocm_fixed_point` is not involved.
- Difference from CUDA: the ROCm kernels are the stored kernels with each read replaced by the element formula of
  `k_mp_build_amat15` (block per atom, same tree), so the result equals the stored path bit for bit in every run below; CUDA
  uses a thread or a warp per atom with its own summation order. CUDA switches at 1 GB (`auto`), ROCm at 4 GB.
- Acceptance (`-sp -method gfn2 -threads 10 -dump_gradient`, n = 1 each; cpu = `-gpu none`, base = binary before WP11, def = new
  default (stored), otf = `-gpu_multipole_otf on`): caffeine, caffeine +1 doublet, AcCl AcI Bz CO ED01 ED07 PR07 PR24, peo201-chain
  (11 structures) in three configurations: gas phase, `-scf_precond frag` (peo201-chain segmented as WP5), water ALPB (33 per
  variant). otf vs def and def vs base: dE 0, dG 0, equal iterations in all 33; the on-the-fly line and the resident loop in every
  otf run. Against the CPU (all three variants alike): gas dE max 5.5e-12 Eh, dG 2.7e-6 Eh/A; frag 1.2e-11 Eh, 5.9e-6 Eh/A;
  water 3.4e-11 Eh, 4.6e-6 Eh/A. Caffeine +1 with frag: the host skips the preconditioner (open shell), as before.
- `-opt` caffeine (`-scf_mixed_precision false`, 40 energy evaluations, each with the on-the-fly begin): final structure identical
  to def, optimisation table equal except the time column. `-md` PR24 100 fs `-seed 42`: trajectory fingerprint (WP10 md5) of otf
  equal to def and to base. Fault injection `CURCUMA_ROCM_FAIL_POT=begin` with otf on caffeine: one warning, energy equal to
  `-xtb.rocm_device_potential false` at 15 decimals.
- Above the limit: lattice of 2940 randomly oriented H2 (5880 atoms, 3.6 A spacing, nao 5880, generator in the session scratchpad,
  not in the repository), stored blocks would need 4.15 GB. auto (on the fly) vs `-xtb.rocm_otf_multipole false` (host potential,
  the old path) vs CPU: dE 1.2e-10 / 5.9e-11 Eh against the CPU, 6.0e-11 Eh between the two, dG 5.5e-6 Eh/A against the CPU for
  both (convergence), 7.3e-15 Eh/A between them; switch off equal to base bit for bit; 6 iterations each. Wall median of 3
  alternating 28.9 s (on the fly) against 40.1 s (host potential), CPU 111.7 s (n = 1); SCF 24.7 s on both GPU paths (eigensolve),
  setup 0.53 against 7.65 s, post-SCF 0.60 against 3.92 s; peak host RSS 817 MB against 11.2 GB, peak VRAM 7.76 against 7.65 GB.
- Device memory limit: peo201-chain (nao 3222) and the H2 lattice need about 27 nao^2 doubles of VRAM (7.76 GB at nao 5880), so on
  this 16 GB card a system above 5770 atoms fits only at about 1 AO per atom; a methane lattice of 6050 atoms (nao 9680) would
  need about 20 GB (estimate from that ratio, not run), and the bulk structures above the limit in `test_cases/structures/bulk`
  (urea400-water1000, nao about 14000; peo201x2-water1500, nao about 15400) do not fit either. The on-the-fly path removes the
  4.15 GB of stored blocks and the host matrices, not the nao^2 buffers.
- Timing below the limit (median of 3 alternating, `-threads 10`): peo201-chain `-sp` with gradient base 13.2 s, def 13.4 s, otf
  13.6 s (SCF 10.96 against 11.05 s for 12 iterations); PR24 `-md` 100 fs base 50.63 s, def 50.65 s. otf saves 0.24 GB of VRAM on
  peo201-chain (2.84 -> 2.60 GB).
- `ctest -j1 -R "gpu|rocm|sqm|xtb|grad|gfnff"`: 271 of 273 passed, failures `xtb_cpscf`, `cli_sqm_12_gfn2_scf_precond` (known),
  2 disabled.
- Not run: CUDA (no toolchain), open-shell and solvated systems above the limit, `-md`/`-opt` above the limit, systems above
  5880 atoms, multi-GPU, other cards. The OTF kernels were not timed separately (inside the SCF time above).
