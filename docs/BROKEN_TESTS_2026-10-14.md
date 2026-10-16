# Failing ctest entries after merging origin/master (Oct 14, 2026)

State: branch `feature/rocm-launch-latency` at `1dc6d366` (origin/master `e0ab0c23` merged in), `release/` build with
`USE_BLAS=ON`, `USE_GFNFF_FAST_EXP=ON`, `USE_ROCM=ON`. `ctest -j1`: 11 of 353 fail. A separate build of
`origin/master` (no ROCm) fails the same CLI tests, so none of them comes from the GPU work on this branch chain.
Each statement carries its sample size; "not checked" means no run exists.

| Test | Cause (status) | Evidence |
|---|---|---|
| `cli_confscan_01` to `07`, `cli_gfnff_06`, `cli_gfnff_07` (9 tests) | FIXED Oct 16, 2026. `CurcumaLogger::set_verbosity` changed the global level only; `get_verbosity()` reads the thread override of `CurcumaMethod` first, so every `set_verbosity` call inside a method scope on that thread had no effect (confscan success line, EEQ solver-path lines in GPU MD). `set_verbosity` now also updates the thread's override when one is active (`src/core/curcuma_logger.cpp`). | After the change: the 9 tests pass (n=1 run each); `ctest -R '^cli_'` 69 of 70 pass, 3 disabled, the remaining failure is `cli_sqm_12` (row below). Not checked: other scopes that relied on the old no-effect behaviour (no test shows one). |
| `cli_analysis_02_scattering_options` | Level mismatch: "Using logarithmic q-spacing" is `CurcumaLogger::info_fmt` (verbosity >= 2, `src/capabilities/analysis.cpp:1617`), the test calls without `-verbosity`. Either the call becomes `result_fmt` or the test passes `-verbosity 2`. | Line count at `-verbosity` 1/2/3: 0/1/1 (n=1 each). |
| `cli_simplemd_14_gfnff_react_h_recombination` | Open. The MD runs to the end but prints no `REACT` line at verbosity 1 to 3. Muting through the thread override was tried at the print site and changed nothing. Not checked: whether `topology_mode react` reaches the GFN-FF object and whether `updateReactiveTopologyIfNeeded()` (`gfnff_method.cpp:2596`) is called and scans. | 3 runs (verbosity 1, 2, 3) and 1 run with the print site patched, all 0 `REACT` lines. |
| `cli_sqm_12_gfn2_scf_precond` | Native GFN2 SCF for `test_cases/cli/sqm/12_gfn2_scf_precond/nacl.xyz` (4 isolated ions) does not converge (last max\|dq\| 0.92), with and without `-scf_precond frag`. Converges in a build with `USE_BLAS=OFF`, `USE_GFNFF_FAST_EXP=OFF`; fails with both ON, also on plain `origin/master`. Probably BLAS (fast exp only touches GFN-FF), not separated. Not a thread effect (`OPENBLAS_NUM_THREADS=1` unchanged). | 1 build per setting, 3 variants each of the same command. |
| `xtb_cpscf` | Two tolerance checks: (C) H2 response 8.2e-5 against 1e-6, H2O isolation 1.70e-4 against 1.5e-4. Fails before this merge too (also on a HEAD build without the GPU changes, reported by a subagent, not re-run by the author). Cause not examined. | 1 run. |

Scope of the logger finding: about 600 calls of `CurcumaLogger::set_verbosity` / `get_verbosity` remain outside the logger
(for example `gfnff_method.cpp`, `forcefield*.cpp`, the capabilities, the xTB methods). Each is ineffective whenever a
`CurcumaMethod` on the same thread has set the override. Which of them hide output a test needs is not checked beyond the
tests above.

Known and unrelated: `test_orca_interface` (needs ORCA; recorded in `qm_methods/CLAUDE.md`) failed before the merge and
does not appear in the list above after it (not checked why).
