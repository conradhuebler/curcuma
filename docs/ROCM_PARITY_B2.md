# ROCm parity B2: pressure and NPT for the periodic GFN-FF on ROCm (Nov 24, 2026)

Claude Generated. 🤖 AI-generated, ⚙️ machine-tested on one RX 9060 XT (gfx1200); human production testing pending.
Item B2 of [ROCM_PARITY_PLAN_2026-11.md](ROCM_PARITY_PLAN_2026-11.md). CUDA reference: A1-GPU (commit `277765fc`, device
virial) and HB-X stage 2 (commit `2b91d73e`, crossing kernels) of [WP-PERIODIC-NONBONDED.md](WP-PERIODIC-NONBONDED.md).
Base: `feature/rocm-parity` `d4dbdfc5` plus the coordinator patches `base_a2.patch` and `b1_only.patch` (A2, B1).
Nothing committed.

## What was ported

`src/core/energy_calculators/ff_methods/rocm/gfnff_rocm.hip`:

- Helpers `vir_outer`, `vir_pair`, `blockReduceAddVirial`, `vir_fix_round` (CUDA text; 64-bit shuffle masks, wave32
  lane indexing as `blockReduceAddEnergy`). The virial slots are fixed point (2^-48 Eh) in both accumulation modes
  (`-gfnff.rocm_fixed_point` true or false), as on CUDA. `VirArr<VIR, N>`: per-thread storage, empty for VIR false.
- Multi-atom terms: `k_angles_body`, `k_dihedrals_body`, `k_inversions_body`, `k_storsions_body`, `k_xbonds_body`,
  `k_hbonds_body` (`template <bool PBC, bool VIR>`) and `k_batm_body<VIR>` with the VIR blocks of the CUDA bodies,
  used only by the new kernels `k_*_vir`. HB: the virial rows are taken from the values that go to `add_grad` or to the
  per-triple slots (slot mode, `k_slot_gather`), so both write paths carry the virial.
- Deviation from CUDA (measured reason): CUDA routes `k_*<PBC>` through the shared body. On gfx1200 that changed the
  ISA of 5 of the 7 former kernels (`k_dihedrals`, `k_inversions`, `k_storsions`, `k_xbonds`, `k_hbonds`, both PBC
  instantiations; register and load scheduling, e.g. `k_dihedrals<false>` 5957 -> 5921 instructions), also with the
  empty `VirArr`. So `k_*<PBC>` keep the former source verbatim and the body is a copy with the VIR blocks; any later
  change to a local-term kernel has to be made in both places.
- Pair terms and chain rules as separate kernels after the step's gradient kernels: `k_bonds_virial`,
  `k_hb_alpha_virial`, `k_repulsion_virial` (bonded and non-bonded list), `k_cn_chainrule_virial` (CN range
  `m_cn_r2_max`), and `k_dispersion_virial_gather` instead of the CUDA `k_dispersion_virial`: the ROCm gradient-mode
  dispersion is the per-atom gather, whose C6 on the fused path comes from `d_disp_uv` (the per-pair `disp.C6` is not
  refreshed then), so the virial kernel walks the same CSR with the same C6 source (half of each pair per direction).
- D4 tail: the C8-only `k_d4_tail_atoms` pass into `GFNFF_VIR_TAIL_E8` (V-LRC); that slot follows the energy
  accumulation mode, the tensor slots are converted with `fixToDoubleHip` always.
- HB-X: `vir_cross_add`, `k_hb_crossing`, `k_xb_crossing` (CUDA text); the per-triple HB energies are recorded while
  the crossing term is on (`ensureHBTripleEnergyBuffer`, `hbTripleEnergyTarget(..., m_record_hb_energies || want_cross)`).
- Host: `d_virial` buffer (`GFNFF_VIR_BUFFER`, `cuda/gfnff_virial_slots.h`), zeroed per gradient step with the virial;
  launch macro `GFNFF_HIP_PBC_VIR_LAUNCH` (VIR kernel only with a cell and a requested virial); BATM takes the FP64
  kernel with the virial; no HIP graphs while the virial is on (phase 1 and phase 2, as CUDA);
  `FFWorkspaceHip::setWantVirial` / `downloadVirial` (declared in the shared header since A1-GPU, now defined).

Wrapper: `qm_methods/gfnff_hip_method.h` `has_periodic_virial = true`. That enables the shared paths of
`gfnff_gpu_method_impl.h` unchanged (`setWantVirial`/`setHBCrossing` per step, `assembleHybridVirial()`,
`followHostCellChanges()` after `GFNFF::updateCell`) and removes the B1 refusal (`setBackendVirialUnavailable`).

Tests: `test_cases/unit/test_gfnff_gpu_pbc_virial.cpp` builds a ROCm target `test_gfnff_rocm_pbc_virial`
(`GFNFF_PBC_VIRIAL_ROCM`, HIP method class, `-gpu rocm` in (e), 77 without a HIP device), ctest `gfnff_rocm_pbc_virial`
(labels `gfnff_gpu;gpu;pbc;simplemd;rocm`, lock `gpu`). `test_gfnff_gpu_pbc_ewald.cpp`: `kHybridVirial` on ROCm =
`GFNFFRocmBackend::has_periodic_virial` (model virial in (a), virial FD in (b), virial in (c) now run on ROCm).

## Measurements

One RX 9060 XT (gfx1200), ROCm 7.2.4, build `-DUSE_ROCM=ON -DUSE_MKL=ON -DUSE_PCH=ON -DUSE_AVX2=ON`, every GPU run under
`flock /tmp/claude-1000/rocm-gpu.lock`. Base = this worktree after the two patches, built before the B2 edits (binary,
plugin and `gfnff_rocm.hip.o` saved).

| gate | result |
|---|---|
| `gfnff_rocm_pbc_virial` (a): GPU vs CPU per term, 5 structures (20 H2O, 60 H2O, H3O+(H2O)2 + 18 H2O Q = +1, caffeine, acetic dimer + F-I...pyridine; L 44 Bohr) | 100 checks pass; worst term component 5.1e-13 Eh, pressure virial 5.5e-13 Eh (limit 1e-7); energy <= 5.6e-13 Eh, gradient <= 1.6e-13 Eh/A; non-zero and checked: bond, angle, dihedral, inversion, D4 pair, repulsion, HB, XB, BATM, Ewald pair and background, D4 tail, V-LRC, HB crossing, CN chain rule |
| (b) FD under cell scaling at fixed fractional coordinates (isotropic + 3 axial, eps 1e-5/2e-5, Richardson), 20 H2O and the composite | 8/8 pass, worst 4.5e-9 Eh (limit 1e-7) |
| (c) second GPU instance, 60 H2O | energy, gradient and every tensor bit for bit; energy and gradient with the virial == without, bitwise |
| (d) `GFNFF::updateCell` (scale 1.001) on the GPU instance vs fresh instance, 20 H2O | energy, gradient, virial difference 0 |
| (e) SimpleMD NPT, Berendsen every 2 steps, 20 steps, CPU vs `-gpu rocm` | both finish, 10 couplings each, cell edge difference 0, pressure 647.0869 bar both (difference 1.6e-10 bar), no GPU fallback |
| `gfnff_rocm_pbc_ewald` with the virial checks on | 172 checks pass; model virial device vs CPU / host-Ewald <= 5.0e-12 Eh (20 checks, limit 1e-7), virial FD <= 5.3e-9 Eh (16 checks, limit 1e-7), 970 H2O two steps bit for bit incl. the virial |
| CLI NPT, 1728 H2O (5184 atoms, L 37.25 A), 100 steps dt 0.5 fs, CSVR 300 K, Berendsen tau_p 100 fs every 2 steps, 4 threads: CPU, ROCm twice | status columns 1-14 of `in.md.csv` (22 rows) identical in all three, trajectory md5 identical in all three, final state md5 identical between the two ROCm runs; final edge 36.79215733249 A, CPU vs ROCm 8.4e-12 A; `-gpu_strict true`, no fallback |
| no cell, 6 structures (water, caffeine, acetic acid dimer, trisaccharide, macrocycle-bgal, peo201-chain; `-sp -dump_gradient -gpu rocm -gpu_strict true -threads 1`), base vs B2 | 6/6 md5-identical; second B2 run 6/6 identical |
| ISA (`a2/isa.py`, gfx1200 code object of `gfnff_rocm.hip.o`), base vs B2 | all 118 base kernels instruction-identical (both PBC instantiations of the local terms included); 14 new kernels (7 `*_vir`, 5 pair/chain-rule virial kernels, 2 crossing) |
| ctest `gfnff_rocm_pbc_(refuse|terms|hybrid|virial|ewald)`, `cli_gfnff_gpu_0*_rocm` | 8/8 pass |
| `scripts/gpu_parity_check.py` | `has_periodic_virial` equal (true/true), 0 unexplained differences; info: `k_dispersion_virial` CUDA only, `k_dispersion_virial_gather` ROCm only |

Timing (n = 1, machine shared with other agents): the 100-step NPT run of 1728 H2O took 18.0 / 17.5 s wall on ROCm and
39.6 s on the CPU (4 threads), start-up included. Not a performance measurement.

## Not tested

- CUDA build of the changed shared files (no CUDA toolkit here); the CUDA sources and the CUDA trait were not touched,
  `test_gfnff_gpu_pbc_virial.cpp` changed only by the ROCm defines and the `-gpu` name variable.
- XB crossing term: zero in every structure run (no XB triple inside the shell); only HB crossing was compared non-zero.
  sTors under a cell (no structure has sTors), case-4 H-bonds, `-gfnff.rocm_fixed_point false` (FP64 accumulation; the
  virial slots stay fixed point there, but an out-of-range contribution is then not reported: `d_gfnff_fix_bad` is read
  only in fixed-point mode), mixed precision, `CURCUMA_ROCM_DISP_FUSED=0` (per-pair C6 path of the dispersion virial),
  `CURCUMA_ROCM_HB_GATHER=0`, multi-stream mode with the virial in the unit test boxes (all below the 500-atom single-stream
  threshold; only the 5184-atom CLI run is above it, and the stream mode of that run was not checked in the log), C-rescale barostat, long NPT runs, wave64 hardware, gfx1150, Windows.
- `ctest -L structinv`, full `ctest -j1`.

## Not implemented relative to CUDA

- L1 fast list rebuilds (`list_fast_rebuild`, `k_hb_detect_pbc`, `k_hb_hnb`, row-block builds, `CURCUMA_GFNFF_LISTHASH`):
  not ported. ROCm keeps its row builds for D4/repulsion and the HB detection that queues every hydrogen with a cell
  (A2): correct, slower list rebuilds after a barostat step.
- Multi-GPU with a cell: off on both backends.

## Proposed central lines (coordinator)

- TODO.md: replace "ROCm periodic virial and NPT (B2)" by "ROCm L1 fast list rebuilds with a cell (k_hb_detect_pbc,
  list_fast_rebuild): docs/ROCM_PARITY_B2.md 'Not implemented'".
- AIChangelog.md: `2026-11-24 ROCm B2: device virial of the periodic GFN-FF hybrid on ROCm (*_vir kernels, pair/chain-rule
  virial kernels, HB/XB crossing), -compute_pressure and the barostat with -gpu rocm; has_periodic_virial true (docs/ROCM_PARITY_B2.md)`.
- GPU_PERIODIC_GAPS.md: virial/NPT on ROCm done Nov 24, 2026 (gfx1200); L1 open on ROCm.
- GPU_PARITY_CUDA_ROCM.md / gpu_parity_known.json: `has_periodic_virial` equal; `k_dispersion_virial` (CUDA) vs
  `k_dispersion_virial_gather` (ROCm): same term, ROCm follows its gather gradient kernel.
- ff_methods/CLAUDE.md "Device virial of the CUDA hybrid": add "ROCm since Nov 24, 2026 (B2): `k_*<PBC>` keep the former
  source, `*_body` is a copy with the VIR blocks (a shared body changed the gfx1200 ISA); change both".
- docs/NPT.md: the barostat also runs with `-gpu rocm -gfnff.gpu_periodic true`.
- The refusal texts in `gfnff_gpu_method_impl.h` ("use the CPU method or CUDA") are now unreachable on both backends.

## Plugin ABI

No change needed: `FFWorkspaceGPU`/`FFWorkspaceHip` declarations are unchanged (the virial members and methods were
declared since A1-GPU), the new buffer lives in the plugin-private `FFWorkspaceHipImpl`, no `ComputationalMethod`
member or entry point changed; the trait is compiled into the plugin.
