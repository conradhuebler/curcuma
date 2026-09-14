# Batch topology reuse ignored the geometry - FIXED (2026-09-13, AI/machine-tested)

Worktree `curcuma-topo` (branch `fix/topo-cache-reuse`), binary md5 **61fe5fe8**, baseline
(pre-change) binary built from the same tree via `git stash`, md5 **fa88d019**. All batch runs
`-no_bmt -threads 1 -gfnff.cache_topology false` in a fresh directory per point.

## A. Cause

`-batch_reuse_topology true` reuses ONE `GFNFF` object; frames 1..N only get `updateGeometry()`.
The existing invalidation (`needsFullTopologyUpdate()` + `GeometryChangeDetector`, >0.5 Bohr)
refreshed only **`TopologyInfo`** (`getCachedTopology()`). The force-field interaction lists -
bonds, angles, torsions, the bonded/non-bonded repulsion partition, the EEQ fragments - are built
once in `initializeForceField()` -> `generateGFNFFParameterSet()` and were **never rebuilt**, so
every frame was evaluated on frame 0's bond graph (at hcn 1.4 r_eq reuse keeps the C-H bond term
-0.58983 Eh with no RepulsionNonbonded, where a fresh perception has -0.48380 Eh plus +0.02742).

Fix (`ff_methods/gfnff.h` / `gfnff_method.cpp` / `main.cpp`): `perceiveGeometricBonds()` (the
geometric perception factored out of `getCachedBondList()`) is compared each `Calculation()`
against `m_ff_bond_graph`, recorded in `generateGFNFFParameterSet()`; on mismatch
`rebuildForceFieldForCurrentGeometry()` re-perceives, regenerates the lists and re-partitions the
workspace (the init sequence), and drops a captured static CN/charge state (gfnff-fast).
Per-frame, not per-0.5-Bohr: H2 gains and loses a bond between two ticks of that trigger, which is
what made the H2 falsifier survive the threshold-only version. New PARAM
**`-gfnff.reuse_topology_check`** (Bool, default false = old behaviour), enabled automatically by
`-batch_reuse_topology true`; `-gfnff.reuse_topology_check false` is the "trust frame 0" opt-out
(verified to reproduce the pre-change number exactly, -0.574817462182 Eh). `topology_mode=react`
and `constant` are excluded: react's bond set is owned by the hysteresis scan and is
path-dependent by design.

## B. Falsifier - frame 0 = r_eq vs frame 0 = the 1.4 r_eq frame, E read at the 1.4 r_eq frame

| bond | before, gfnff+reuse | before, gfnff-fast | after, both |
|---|---:|---:|---:|
| hcn_HC-H | 1.343e-01 Eh (84.31 kcal) | 1.745e-01 (109.48) | **0.000e+00** |
| h2_H-H | 1.541e-01 Eh (96.69 kcal) | 1.861e-01 (116.79) | **0.000e+00** |
| ch4_C-H | 1.054e-01 Eh (66.13 kcal) | 1.297e-01 (81.39) | **0.000e+00** |
| h2o_O-H | 2.718e-03 Eh (1.71 kcal) | 7.103e-02 (44.57) | **0.000e+00** |
| react (all four) | -15.69 ... +1.01 kcal, unchanged | - | - (excluded) |

Residual exactly 0, i.e. bit-identical between the two seeding orders. h2/ch4/h2o now equal the
fresh single point exactly; hcn reuse stays 8.3e-4 Eh (0.52 kcal/mol) high because the reuse
protocol generates its parameters at the frame where the topology last changed (1.30 r_eq) while
a fresh point generates at 1.40 (Coulomb +8.96e-4, dispersion -6.5e-5; bond/angle/torsion/
repulsion exact) - a protocol difference, no longer an order dependence.

## C. Regressions checked

- Homogeneous batch (caffeine MD trajectory 12 frames x20 = 240, max per-atom displacement 2.54 A,
  0 topology rebuilds): max |dE| = **0.000e+00 Eh over 240 frames**, 0 frames differing. Wall
  time min of 3: base 0.0809 s, fixed 0.0802 s (no measurable change). Water (3 atoms, 240
  frames) costs +0.03 ms/frame for the extra perception; its trajectory had exploded (O-H
  19-34 A in the last 4 of 12 frames), and there the fix changes those 80 frames by <=2.0e-8 Eh -
  the intended correction (no bond perceived at 19 A).
- gfnff single point + gradient, exact: caffeine E -4.672737068614 (dE 0.000e+00), benzene
  -2.362725526194 (dE 0.000e+00), max |dG| component 0.000e+00 both.
- Build: exit 0; warning sets of `gfnff.h`, `gfnff_method.cpp`, `main.cpp` identical to the
  baseline full build (45 / 12 / 157 messages, no additions).
- `ctest -R gfnff`: 60/65 pass, 5 fail (98, 99, 106, 108, 109), 90.2 s. `ctest -R cli_simplemd_`:
  17/22 pass, same 5 fail, 85.9 s. All 5 reproduce identically on the baseline binary
  (98/99 "Total energy drift 0.446525 Eh"; 106 "T_max 7149.634619"; 108 rev slopes
  4.2266e-03/6.1795e-03; 109 "mean Epot -0.3995 kcal/mol") -> pre-existing. Scoping the check to
  the reuse path restored the plain-gfnff NVE control in test 108 to its baseline slope
  (-9.4780e-06 / 2.7042e-05; the unscoped version had made it 9.8e-4, i.e. it altered MD).
  No `release_tblite/dumps` tree in this worktree (configure said so); none of these tests need it.

## D. Part B - the react bond drop (measured, NOT fixed; 7 X-H bonds, `revgfnff` react, ascending grid)

| bond | r_eq [A] | r_drop [A] | r/r_eq | deciding switch | dE_jump at the drop | with rev transitions off |
|---|---:|---:|---:|---|---:|---:|
| hcn_HC-H | 1.0691 | 1.9245 | 1.80 | rev transition coord | < 5e-7 Eh | 2.75 r_eq, +0.009012 Eh (+5.65 kcal) |
| h2_H-H | 0.7415 | 1.1123 | 1.50 | rev transition coord | < 5e-7 Eh | 2.50 r_eq, +0.011059 Eh (+6.94) |
| hf_H-F | 0.9233 | 1.6620 | 1.80 | rev transition coord | < 5e-7 Eh | 3.00 r_eq, +0.011201 Eh (+7.03) |
| hcl_H-Cl | 1.2788 | 2.3018 | 1.80 | rev transition coord | < 5e-7 Eh | 2.75 r_eq, +0.002920 Eh (+1.83) |
| ch4_C-H | 1.0914 | 1.9645 | 1.80 | rev transition coord | -4.61e-4 Eh (-0.29 kcal) | 2.75 r_eq, -0.002959 Eh (-1.86) |
| h2o_O-H | 0.9618 | 1.7312 | 1.80 | rev transition coord | < 5e-7 Eh | 2.75 r_eq, -0.066128 Eh (-41.5) |
| nh3_N-H | 1.0143 | 1.8257 | 1.80 | rev transition coord | < 5e-7 Eh | 2.75 r_eq, -0.045312 Eh (-28.4) |

Attribution (bisection of the frame count until `REACT bond broken` appears, 5 %-of-r_eq grid, so
the crossing carries ~0.05 r_eq of grid error): the drop follows `rev_bo3_center` exactly
(1.6 -> 2.0 moves every drop 1.80 -> 2.25 r_eq) and moves in when `rev_tr_prebreak` is raised
(0.5 -> 0.9 gives 1.80 -> 1.50), while `rev_bo_break` (0.02 -> 0.9) and
`react_bond_break_factor` (2.6 -> 1.6) change **nothing**: the stage-1b transition coordinate
removes the bond before either the bond-order weight switch or the react hysteresis can. With the
rev machinery off (`-gfnff.rev_enabled false`, break factor 1.6) the drop moves out to
2.50-3.00 r_eq and becomes a real step (last column).

Reading for the orchestrator: the 1.8 r_eq topology drop is **energy-neutral** (dE_jump < 5e-7 Eh
in 6/7 cases), because the transition blend has already carried the bonded terms to ~0 there. The
well is not truncated by the removal; it is truncated by the bond-weight decay before it - at the
drop the model sits only +9.8...+27.9 kcal/mol above its own minimum where the reference is
+30.3...+94.5 (h2o +18.7 vs +89.9 at 1.73 A). So this is the stage-3a(iii) join radius, which the
design review already expects to move out (~2.6x), not extra hysteresis: raising the break factor
or the weight threshold moves the drop *outward* but leaves the well shallow.