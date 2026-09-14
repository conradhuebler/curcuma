# rev_over_preset: finishing the working-tree change (Sep 12, 2026)

Status: **done, uncommitted**. Only `gfnff.h` + `gfnff_method.cpp` still differ from HEAD;
`ff_workspace.{h,cpp}` and `ff_workspace_gfnff.cpp` are back at HEAD content (they carried only
`CURCUMA_REVTRACE` prints - all removed, zero matches left in `src/`).

## What the change is
- New PARAM `rev_over_preset` (gfnff.h:344), **default `stage1a`**: the fit is opt-in (operator
  decision 2026-09-12, until stage 3 re-measures the barriers).
- `stage1a` sets no table - cleared maps let `revOverP()` fall back to `rev_over_p` 0.3 and
  `revValence()` to the nominal valences, `rev_over_shift` stays 0.5. Equals HEAD behaviour.
- `fit2026-09-12` = p_over H 0.0 / C 0.209059 / N 0.339041 / O 0.994136, shift 0.869859,
  valence N 2.536341 / O 2.530618. Mirrors `params/rev_over_fit_2026-09-12.json`.
- Documented limitation: the preset's shift is applied only while `rev_over_shift` still reads the
  registry default 0.5, so an explicit `0.5` cannot be told from "untouched". Measured in (c).
- Unknown preset -> warn + stage1a, with `m_parameters["rev_over_preset"]` normalised so the log
  line and `-export_run` report what ran.
- `InitialiseMolecule()` react log also prints shift / k / p(H,C,N,O) / valence(N,O) / preset.
## Acceptance (all in `release/`, ASCII output only)
a) `make -j4` exit status **0**. Warnings of the 3 recompiled TUs, HEAD vs working tree, same
   flags (`-Wno-inline`): **67 vs 67, normalized sets byte-identical** -> no new warnings.
b) `ctest -R gfnff`: **63 passed / 1 failed of 64**; the failure is pre-existing (last section).

c) Preset equivalence. Single points with `-gfnff.rev_enabled true`; energies from the 12-digit
   `# energy` line of `-dump_gradient`; fresh dir + `-gfnff.cache_topology false` per run.
   Geometries: linear H3 at 0.85 A, CH4+H at 1.25 A (over-coordinated, so E_over bites; a relaxed
   CH3OH frame separates the presets by only 6.1e-07 Eh).

   | run | H3 [Eh] | CH5 [Eh] |
   |---|---|---|
   | default (no flag) | -0.091839394148 | -0.630376167896 |
   | `-gfnff.rev_over_preset stage1a` | -0.091839394148 | -0.630376167896 |
   | `-gfnff.param_file rev_over_stage1a.json` | -0.091839394148 | -0.630376167896 |
   | `-gfnff.rev_over_preset fit2026-09-12` | -0.094531584764 | -0.679200442685 |
   | `-gfnff.param_file rev_over_fit_2026-09-12.json` | -0.094531584764 | -0.679200442685 |
   | fit + `-gfnff.rev_over_shift 0.7` | -0.094531584764 | -0.670756567434 |
   | fit + `-gfnff.rev_over_shift 0.5` | -0.094531584764 | -0.679200442685 |
   | unknown preset `nonsense` | -0.091839394148 | -0.630376167896 |
   | rev off (`-method gfnff`) | -0.094540060389 | -0.680919186371 |

   preset == param_file in both cases to **< 1e-12 Eh** (all 12 digits). stage1a vs fit differ by
   **2.69e-03 Eh (H3)** / **4.88e-02 Eh (CH5)**. shift 0.7 moves CH5 by **8.44e-03 Eh** (H3 is
   insensitive, the fit has p(H)=0); shift 0.5 == fit default = the documented limitation. Log
   line verified: default `preset stage1a` shift 0.500 p 0.300 x4; with the flag `fit2026-09-12`
   shift 0.870 p 0.000 0.209 0.339 0.994; with the override shift 0.700.

d) Compared against a binary built from **HEAD content** (files swapped in, built, swapped back;
   md5 verified). Identical to all 12 digits, rev off **and** rev default: CH3OH -0.755670961012 /
   -0.755670272794, H3 -0.094540060389 / -0.091839394148, CH5 -0.680919186371 / -0.630376167896 -
   the new `stage1a` default reproduces HEAD exactly.

## The one ctest failure is not from this change
`cli_simplemd_18_gfnff_rev_nve_vs_gfnff` reports `ratio(dt=0.25) 101.3333 outside [0.5,1.5]`. The
**HEAD binary gives the identical 101.3333**, and `-gfnff.rev_over_preset fit2026-09-12` gives
102.6667 - the preset is irrelevant. The test divides two endpoint Etot differences, and its own
header calibration (1.005 / 1.012) is reproduced exactly from the **second-to-last** table row
(gfnff 2.060e-04 / 8.300e-05, rev 2.070e-04 / 8.400e-05) - that is the stale final-row Etot which
commit 4ef8d30c fixed. The test is calibrated pre-4ef8d30c; the fix belongs in the test (fit a
slope instead of using endpoints, Known Issue #32), not in `src/`.
