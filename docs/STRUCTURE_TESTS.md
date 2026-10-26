# Generated structure property tests

AI-generated (Oct 2026), machine-tested only. Property tests that need no stored reference value, generated from the
structure library (`test_cases/structures/manifest.json`) and a policy (`test_cases/struct_test_policy.json`).

## What runs

`scripts/gen_struct_tests.py` selects (structure, method) pairs at configure time and writes one CTest entry per pair
(`structinv_<id>_<method>`, label `structinv`; with CUDA also `structinv_gpu_*`, with ROCm `structinv_rocm_*`, labels `gpu`, `RESOURCE_LOCK gpu`).
`scripts/struct_invariants.py` runs the checks through the `curcuma` command line (`-sp`, `-dump_gradient`, Eh/Angstrom):

| test | passes when |
|---|---|
| translate / rotate | energy unchanged, gradient unchanged / rotated with the molecule |
| permute | energy unchanged by a fixed shuffle of the atom order |
| sumforce / torque | sum of the forces and of r x g is zero |
| fd | analytic gradient equals the central finite difference of curcuma's own energy (Richardson, h = 0.005 A, first / middle / last atom, tight `-scf_threshold 1e-10` for gfn1/gfn2) |
| repeat / threads | two identical runs agree; `-threads 1` against `-threads 4` |
| compressed / compressed_gpu | geometry scaled to 0.30 about its centroid (every D4 reference far below the CN): finite-difference gradient (h = 0.001 A) and force sum on the CPU; CPU against the GPU gradient. Written for the D4 weight underflow, which gave wrong gradients on the host and on CUDA |
| gpu | CPU against `-gpu cuda` or `-gpu rocm` (energy and gradient, `-gpu_strict true`; skipped with exit 77 when the backend is not available); the ROCm variant was generated but has not run on ROCm yet |

Current selection, family `invariants` (prefix `structinv`): 12 structures with at most 12 atoms (10 assumed neutral closed shell, `dinitrogen`,
`water-dimer.roo2.8`) x gfn1, gfn2, gfnff = 36 CPU tests, 36 CUDA tests. Family `compressed` (prefix `structcomp`): 14 structures with at most 16 atoms x gfn1, gfn2 = 28 CPU tests, 28 CUDA tests; with CUDA 128 tests in total, 128 of 128 pass in 73 s wall (`ctest -L structinv -j8`, one RTX A4500).
`python3 scripts/gen_struct_tests.py --plan` lists the selection and the reason for every skipped structure.

## Tolerances

Set from the worst case of the pilot run times a margin of at least 20 (`measured` block of the policy), not from a wish.
Energy 1e-8 Eh, gradient 1e-8 Eh/A, sum of forces 1e-10, torque 1e-7, fd 1e-6 (gfnff 1e-5), repeat 1e-12, threads 1e-10,
CPU against GPU 1e-8 Eh and 1e-5 Eh/A; compressed: fd 1e-5, force sum 1e-8, CPU against GPU gradient 1e-6 (worst measured 3.4e-7 and 1.4e-9). A limit is loosened only through `known_deviations` (structure + method + test + reason, printed with
the result) or the method-wide `permute_energy` of gfnff (1e-4 Eh, Known Issue #34: HCN 6.8e-5 Eh measured).

## Findings of the first run

- N2: the plain central difference left 1.2e-4 to 1.8e-4 Eh/A for all three methods (truncation, anharmonic bond); Richardson extrapolation brought it to 4e-8.
- gfnff permutation: hydrogen-cyanide differs by 6.8e-5 Eh under renumbering (Known Issue #34, now covered by a bounded limit).
- gfnff `ammonia.ideal-c3v`: analytic and finite-difference gradient differ by 6e-5 Eh/A independent of the step. The inversion energy has a kink at exactly this
  geometry (Known Issue #42, a known deviation with a 2e-4 limit).

**Negative control** (Oct 26, 2026): with the CUDA plugin from before the D4 weight fix, `compressed_gpu` fails on `acetic-acid-dimer` gfn2 with max|dg| 0.435 Eh/A (limit 1e-6); with the current plugin it passes at 7e-11.

## Charge and spin

The manifest has charge/spin for 23 of 80 structures. The generator uses a structure without them only if the policy lists it under
`assumed_neutral_closed_shell` and the electron count is even. This list is an AI proposal that the operator accepted "for now" on Oct 26, 2026; charge/spin are still
not written into the manifest. Other structures are skipped, the plan says why.

## Not covered

Only structures up to 12 atoms, CPU single process; no size-consistency test of clusters, no open-shell or charged structure (none has a confirmed
charge/spin), no periodic systems, no uff/qmdff, no ROCm/Vulkan, no energy against an external program (no xtb/tblite/orca here). Tolerances rest on one
run per pair (n = 1), not on a distribution.

## Adding a structure or method

Add the structure to the library with charge and spin (`structlib.py add`), or confirm it in the policy list, and re-run `cmake .` in the build directory.
A new method goes into `families.invariants.methods` of the policy.
