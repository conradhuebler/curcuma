# STALE_CN_STATUS - GFN-FF on a reused calculator read a stale CN (plain GFN-FF, `-method gfnff`)

Sep 24, 2026. Opus agent. No `git commit`. Plain-GFN-FF bug, first seen in package 23
(`P2P3_STATUS.md` 8.1). Scratch: `<scratchpad>/stale/` (section 9). AI-generated and
machine-tested only; human production testing pending.

## Recommendation

**Two stale-CN defects, not one. Both fixes are correctness fixes with no measured regression
in single-point numbers. I recommend making both unconditional, with no flag.**

- **Fix A (the reported bug): IN THE WORKING TREE.** The Coulomb self-energy term chi(CN) read a
  CN that only gradient calls refreshed. Cl2- falsifier: 1.77e-2 -> **2.41e-5 Eh/A**, exactly
  as P2P3_STATUS.md predicted. Effect on per-structure single points: GMTKN55 2462 / MOR41 95 /
  S30L-CI 90 all bit-identical at 8 decimals (max 1.5e-10 Eh at full precision). MD is
  unchanged, because every MD step is a gradient call. What does change: every energy-only call
  on a reused calculator. User-visible cases: the native `-opt.optimizer lbfgs` never converged
  on caffeine (5000 iterations) and now converges in 44 steps; the energy-based Hessian
  (`-hess 2`) of H2O was 81 cm-1 off and is now 0.3 cm-1 off.
- **Fix B (found while verifying A): READY AS A PATCH, NOT APPLIED to the shared tree.** The
  2.41e-5 left after A is not the true residual either. On fresh calculators the FD matches the
  analytic gradient to 2e-12. The rest comes from the D4 pair C6: it was baked once at topology
  build and never refreshed, in both energy AND gradient calls, while the analytic gradient
  already contains dC6/dCN. With A+B, every FD residual in `gfnff_sqe` falls to 1e-10..1e-11
  (Cl2- 1.8e-11; the documented CH4+H react "known residual" 2.36e-4 -> 5.3e-11). The energy
  that `-opt` reports now equals a fresh single point at the optimised geometry. Before, it was
  off by up to 0.16 kcal/mol (UPU23/2h), and `cli_curcumaopt_07`'s 17 optimised energies had a
  2.8e-5 Eh spread; with B they fall within 4e-6 Eh of each other. Patch:
  `test_cases/revgfnff/_log/STALE_CN_fixB.patch` (`git apply --check` clean). I left it out of
  the shared tree because a concurrent agent is building from it (section 8), and because it
  changes MD/opt numbers, which is the operator's decision.
- **The costs of B, stated plainly.** On the 1410-atom polymer, single thread, it adds 5-11 ms
  to an energy-only call of ~25 ms and stays within noise on a gradient call of ~65-80 ms. The
  machine was loaded (load average ~22), so these timings are rough. B makes the already-failing
  `cli_curcumaopt_07` fail 2 more frames, by 1.2e-5 and 1.9e-5 Eh against a 1e-5 tolerance,
  because its golden values encode the old scatter. On UPU23/2h, `-opt` lands in a different
  local minimum, 0.29 kcal/mol higher (n=1, flexible RNA backbone). I found nothing where the
  stale behaviour compensated for a real error. The one partial cancellation I saw was between
  the stale Coulomb term and a frozen-topology difference (WATER27/H2O4 perturbed, section 6);
  that is not a physics compensation.
- **Not fixed, and not stale CN:** (i) kept-vs-fresh differences of up to ~30 kcal/mol on 16/151
  perturbed structures remain, in the Angle/Inversion/Bond/Repulsion/HB terms. This is the
  designed frozen topology (hybridisation, angle/inversion setup, HB and fragment perception),
  which the reference freezes in MD as well. (ii) Caffeine NVE has a dt-independent drift of
  ~+12 uEh/ps that A and B leave unchanged and that is not caused by COM/rotation removal. Its
  cause is not identified. (iii) The GPU path (`gpu_only`) has its own C6/self-energy handling;
  not checked.

## 0. Baseline and method

- Handed-over binary `build_rev/curcuma` md5 `d9908823` = `curcuma_base`.
- A concurrent agent (soft mu q0 rule, `MU_CUSP_STATUS.md`) edited `gfnff.h`, `eeq_solver.*`
  and `gfnff_method.cpp` and rebuilt `build_rev` during this work (section 8). All A/B
  comparisons therefore come from ONE private source snapshot built three ways in
  `<scratchpad>/stale/bld` (`variant.py` toggles only my call sites):
  `curcuma_sbase` 8857fe63 (no fix), `curcuma_sA` 3f702f64 (A), `curcuma_sAB3` 82213ebf (A+B).
  `curcuma_sbase` is bit-identical to `d9908823` on all 2462 GMTKN55 + 185 MOR41/S30L-CI
  energies (full precision), so the snapshot's foreign changes do not touch plain GFN-FF.
- Full `ctest` on the unmodified `build_rev` (before anyone rebuilt it): 13 failures, the
  documented set (confscan_dtemplate, test_orca_interface, xtb_cpscf,
  cli_curcumaopt_07_opt_multixyz, cli_confscan_01..07, cli_simplemd_18/20).

## 1. Root cause A

`FFWorkspace::postProcess()` evaluates `-sum_i q_i (chi_base_i + cnf_i sqrt(CN_i))` from
`FFWorkspace::m_cn`. Exactly one function, `setCNDerivatives()`, wrote `m_cn`, and
`GFNFF::prepareCNAndEEQ()` called it only in its gradient branch. The energy-only branch
recomputes `m_last_cn` and solves the EEQ with it, so the charges are right, but it never
handed the CN over. An energy-only call therefore read the CN of the last gradient geometry. With
no earlier gradient call it read `chi_static`, built from the topology-build CN. Per-structure
single points never showed the bug: their first call reads `chi_static` at their own geometry.
`NumGradFixedCharges()` had the same gap. Static-CN mode (WP-S1) is unaffected: it freezes
`m_last_cn` itself (`reuse_cn`) and hands the frozen value over either way; `gfnff-fast` MD is
bit-identical base vs A vs A+B.

**Fix A** (in the tree): new `FFWorkspace::setCN()` (`ff_workspace.h`), called next to
`setD3CN(m_last_cn)` on every CPU call in `prepareCNAndEEQ()` and in the
`NumGradFixedCharges()` energy lambda (`gfnff_method.cpp`, ~l.1497 and ~l.12382).

## 2. Root cause B (patch only)

`GFNFFDispersion::C6` in the workspace pair list is set once in `GenerateDispersionPairsNative()`
(topology build) and never updated. The Fortran reference evaluates C6(CN) every call
(`gfnff_engrad.F90:323`, `d3_gradient(..., cn, dcn, ...)`), and curcuma's own gradient applies
dC6/dCN via `m_dc6dcn_ptr`. On a reused calculator (MD, opt, FD, batch reuse) the energy was
therefore a frozen-C6 energy. The force was not the derivative of that energy: the dr-part used
the frozen C6, the CN-part used the current dC6/dCN.

**Fix B** (`STALE_CN_fixB.patch`, 4 files, additions only): `D4ParameterGenerator::
refreshC6WeightsForCN()` recomputes the Gaussian CN weights and the C6 half-contraction at the
current CN, with no P1a threshold. The dc6dcn derivative and its P1a cache are untouched.
`FFWorkspace::forEachD4PairList()` visits the installed list and every stored rev-gfnff corner's
list. In `prepareCNAndEEQ()` every pair's C6 is refreshed, threaded over the pool above 4096
pairs; the refresh is skipped in static-CN mode and on the GPU-only path.

## 3. The falsifier (Cl2-, q=-1, 2.73 A), max |g - g_FD| in Eh/A

`fdcheck.py`: FD energies on one reused calculator, energy-only (reuse-E) or with a gradient call
at every point (reuse-G), or on fresh calculators (fresh-E).

| binary | reuse-E | reuse-G | fresh-E |
|---|---:|---:|---:|
| base | 1.766e-02 (h = 1e-3/1e-4/1e-5, identical) | 2.413e-05 | 2.0e-12 |
| A | **2.413e-05** | 2.413e-05 | 2.0e-12 |
| A+B | **2.6e-12** | 2.6e-12 | 2.0e-12 |

Kept-vs-fresh per term (topology built at 2.73 A), kcal/mol, base: Coulomb +1.1e-3 / -0.11 /
-0.33 / -0.85 at 2.7301 / 2.72 / 2.70 / 2.65 A, plus Dispersion 1.5e-6 .. 1.5e-3. A removes
the Coulomb part, B the dispersion part: A+B agrees with fresh to 1.4e-13 kcal/mol.

`test_gfnff_sqe` FD lines (plain-gfnff reference residual in brackets are the same numbers):

| case | base | A | A+B |
|---|---:|---:|---:|
| 2a Cl2- 2.73 A | 1.77e-2 | 2.41e-5 | 1.8e-11 |
| 2b HCOO-...HF | 9.34e-4 | 9.42e-6 | 2.5e-10 |
| 2c CH4+H react in flight ("known residual") | 2.36e-4 | 4.64e-5 | 5.3e-11 |
| 4d F2- 1.60 A (FD via gradient calls) | 6.40e-5 | 6.40e-5 | 1.9e-10 |

So the comments in `test_gfnff_sqe.cpp` (l.~162: "inherently hard ... fragment/EEQ corner";
l.~215-235: "KNOWN RESIDUAL ... not FD truncation") describe this bug, not a method limitation.
I did not edit them: the file is in a concurrently edited area. The orchestrator should update
them and can tighten `tol_corner` 3e-4 -> ~1e-6 once B is in.

## 4. Per-structure single points (fresh calculator each, full precision via `-batch`)

| set | n | A vs base | A+B vs base |
|---|---:|---|---|
| GMTKN55 | 2462 | 2 differ, max 3.2e-14 Eh | 2 differ, max 3.2e-14 Eh |
| MOR41 + S30L-CI | 185 | 3 differ, max 1.5e-10 Eh (MOR41/PR21) | same |

At the 8-decimal CLI format of all earlier campaigns: 0/2462 and 0/185 changed. MAD/max vs xtb,
per-subset reaction MADs and WTMAD-2 are therefore unchanged (no structure moves by 1e-8 Eh);
`gmtkn55_reactions.py` was not re-run for that reason.

## 5. Reused-calculator paths

**Kept-vs-fresh sweep** (`kept/kept.py`, deterministic): 154 structures (3 per GMTKN55 subset +
every 6th MOR41), 11 perturbed frames each (sigma 0.03/0.06 A), topology kept from frame 0,
max |E_kept - E_fresh| per structure, kcal/mol:

| binary | mode | median | p90 | >0.01 | >0.1 | >1 |
|---|---|---:|---:|---:|---:|---:|
| base | energy-only | 0.40 | 2.8 | 128 | 98 | 47 |
| base | gradient | 0.032 | 0.99 | 99 | 46 | 15 |
| A | either (identical) | 0.031 | 1.0 | 98 | 45 | 15 |
| A+B | either | **0.00056** | 1.1 | 66 | 39 | 16 |

Remainder after A+B: term diffs of the worst cases are Angle/Inversion (MOR41 ED32/PR34,
FH51, BHPERI) or Bond/Repulsion/HBond/Coulomb (WCPT18/ts1h2o, WATER27 clusters). This is the
frozen topology by design, not CN. 3-5 structures per binary are missing because of a
pre-existing crash (section 7).

**MD** (`md/runmd.sh`: acetic-acid-dimer NVE, AHB21/21 NVT q=-1, react 4 H at 6000 K = the
`cli_simplemd_14` scenario, gfnff-fast): base vs A identical except the wall-time column. A+B
changes trajectories (forces now use the current C6); gfnff-fast is bit-identical in all three
variants.

NVE, 10 ps, dt 0.25 fs, 4 perturbed replicas each (`nve/rep.sh`), Etot fit, uEh:

| system | variant | slope (uEh/ps) | residual sd | max per-print step |
|---|---|---|---:|---:|
| caffeine | base | 13.1 / 14.5 / 11.6 / 10.1 | 8.4 / 9.1 / 7.6 / 8.1 | 41 |
| caffeine | A+B | 13.3 / 14.6 / 13.5 / 11.7 | 4.3 / 6.2 / 4.3 / 4.1 | 16 |
| acetic-acid dimer | base | -9 / 0 / -19 / -19 | 26 / 24 / 33 / 42 | 107-129 |
| acetic-acid dimer | A+B | -14 / -15 / -41 / +68 | 21 / 24 / 34 / 223 | 79-776 |

B halves caffeine's Etot fluctuation. The ~12 uEh/ps drift is the same at dt 0.5 and 0.25, and
the same without COM/rotation removal: an unrelated, unidentified inconsistency. The
acetic-acid rare mEh jump events occur in both variants (base: 1 of 5 runs, A+B: 1 of 5),
consistent with chaos, n=5 per arm. `cli_simplemd_18` plain-gfnff slope at dt 0.125:
1.44e-6 (base) -> 5.65e-7 (A+B), n=1.

**Optimisers** (caffeine, AHB21/21): auto/lbfgspp/diis/ancopt give the same md5 base vs A. Native
`lbfgs`: base stalls at -4.677370 for 5000 iterations with no output, A converges in 44
steps, A+B in 42. Reported minimum energy vs a fresh single point at that minimum
(`optset/run.sh`, 11 structures, kcal/mol):

| structure | base | A+B |
|---|---:|---:|
| UPU23/2h | -0.159 | +0.0004 |
| IDISP/F22l | -0.024 | +0.0001 |
| HAL59/FI_pyr | -0.016 | -0.0001 |
| Amino20x4/ARG_xak | +0.0066 | -0.0002 |
| caffeine | +0.0044 | 0.0000 |
| MOR41/PR22 | -0.059 | -0.013 |
| MOR41/ED17 | -1.28 | -1.21 (frozen-topology angle term, not CN) |
| other 4 | <= 0.003 | <= 0.002 |

**Energy-based Hessian** (`-hess 2`) vs gradient-based (`-hess 1`), same binary:

| system | base | A | A+B |
|---|---:|---:|---:|
| H2O | 81.0 / 52.0 cm-1 | 0.30 / 0.23 | 0.20 / 0.17 |
| caffeine (>200 cm-1) | 26.8 / 4.2 | 21.3 / 2.2 | 21.3 / 2.2 |

The entries are max / mean deviation. Caffeine's remaining 21 cm-1 is not stale CN.

Paths not measured individually but fixed by A: `GFNFF::NumGrad()`/`NumGradFixedCharges()`,
`EnergyCalculator` numerical gradient, `-batch_reuse_topology true` (`revgfnff_classa.py`'s kept
protocol), C++ FD tests using `CalculateEnergy(false)`. Not affected: Casino (re-initialises
every step), ConfScan (a fresh calculator per conformer).

## 6. Where a number gets "worse"

- `cli_curcumaopt_07` (already failing on 2 frames, golden drift): with B, 4 frames fail. The
  new ones are frames 11 and 14 at 1.2e-5 and 1.9e-5 Eh against the 1e-5 tolerance. The 17
  optimised energies: base -9.295164 .. -9.295192 (spread 2.8e-5); A+B -9.295179 .. -9.295183
  (spread 4e-6). The golden file encodes the frozen-C6 scatter and should be regenerated if B
  goes in.
- UPU23/2h `-opt` reaches a minimum 0.29 kcal/mol higher (fresh energies of the two minima).
  This is a different local minimum on a floppy backbone, n=1, not a systematic effect.
- WATER27/H2O4 perturbed frame 7: kept-vs-fresh goes -16.5 -> -28.0 kcal/mol with A, because
  the stale Coulomb term (+11.5) partly cancelled a frozen-topology difference. A makes
  energy-only equal to gradient calls (-27.9 before too); nothing physical was being
  compensated.
- `cli_simplemd_14` (react, H recombination, 6000 K): the scenario gives 1 formed / 1 broken
  (base) vs 1 formed / 0 broken (A+B), which is chaos. The test's criteria pass in all variants.

## 7. ctest, full suite (snapshot builds, `CURCUMA=<bld>/curcuma ctest -j16`, 305 run)

**Identical 21 failures for base, A and A+B**: the 13 documented ones, plus `gfnff_sqe`. Its only
failing line is B2/3d, which the concurrent mu-q0 edits in the snapshot cause (it fails in the
snapshot base too), not this work. The last 7 are `parameter_io_tests` and `cli_errors_01..06`,
which run `<source>/release/curcuma`, a path that does not exist in a copied tree (same artefact
as in FRAG_CHARGE_STATUS.md 12). Inside the failing tests, only `cli_curcumaopt_07` changes
(section 6). `cli_simplemd_18/20` fail for the same reasons in all three.

New permanent test **`cli_gfnff_06_stale_cn_energy_only`**
(`test_cases/cli/gfnff/06_stale_cn_energy_only/`, registered in `test_cases/cli/CMakeLists.txt`,
label `gfnff;stale_cn`). It checks FD on a reused calculator, energy-only == gradient-call
energy, and reused Coulomb == fresh. It **fails 4/4 checks on the pre-fix binary** (1.77e-2,
1.36e-3, 1.36e-3, 5.0e-4) and passes 4/4 with A and with A+B. It was not run through ctest
here; `cmake .` in the build dir is needed to register it.

Side finding (pre-existing, all three binaries): a multi-frame `-batch true` run WITHOUT topology
reuse can SIGSEGV in `GFNFF::classifyBondType` (`calculateTopologyInfoOnce`). Seen on 12
perturbed frames of MOR41/PR22, MB16-43/32 and AL2X6/al2me4. Every frame alone, and every
0+k pair, runs fine, so the crash depends on history: state carried between successive fresh
calculators in one process. Not investigated.

## 8. Process incidents

- **Concurrent agent on the same files.** gfnff.h/eeq_solver.*/gfnff_method.cpp changed under me
  from 21:35, and `build_rev/curcuma` was rebuilt at 21:40 (md5 df729c0c). Fix A was already in
  the tree then, so that agent's binaries since 21:40 contain A (energy-only/FD numbers on reused
  calculators differ from d9908823). B was in the tree only from ~21:44 to ~21:48, then
  removed. **But `build_rev/curcuma` was rebuilt at 21:49:21 and its md5 `287e2f11` equals my
  snapshot A+B (first-version B) build exactly: the current `build_rev` binary CONTAINS fix B,
  while the source tree no longer does.** Any numbers the mu agent took with it include B
  (kept-calculator dispersion and MD/opt change slightly, and so does its FD). Its next `make`
  drops B again, so results would silently shift between two of its builds. Both A-only and
  A+B were tested here and both are sound, but its before/after pairs must come from the same
  build.
  My first `ctest` after A ran against that foreign rebuild; I discarded it and replaced it with
  the snapshot runs.
- **/tmp (tmpfs) ran full** (~22:20, 94G/94G). My 11 GB snapshot build was part of it, together
  with `frag/` 17G, `mu/` 12G and `x2scope/` 5G. I deleted my build dir (11 GB freed). Other
  agents' runs around 22:20 may have failed with "No space left on device", so their outputs
  from that window should be checked. My own runs affected: one caffeine Hessian, rerun. Free
  space afterwards: 23 GB. My snapshot `src/` (2.7 GB) is still there, for reproduction.

## 9. Reproduction (`<scratchpad>/stale/`)

`fdcheck.py BIN XYZ Q h` (falsifier), `keptterms.py BIN MULTI.xyz Q` (per-term kept vs fresh),
`g55full.py`/`setfull.py` (full-precision GMTKN55 / MOR41+S30L-CI), `diff.py`, `kept/kept.py`,
`md/runmd.sh`, `nve/run.sh`+`rep.sh`+`ana.py`, `opt/reopt.sh`, `optset/run.sh`, `hess/fr.py`,
`variant.py SRC base|A|AB`, `src/` (snapshot), `fixB_block.txt`. Binaries: `curcuma_base`
(d9908823), `curcuma_A` (ad2988d3, shared-tree build, contaminated by foreign edits, used only in
sections 3-5 where it matches `curcuma_sA`), `curcuma_sbase`, `curcuma_sA`, `curcuma_sAB3`.

## 10. For the orchestrator

- Apply B: `git apply test_cases/revgfnff/_log/STALE_CN_fixB.patch`, once the concurrent agent
  is done, if the operator agrees. Then regenerate the `cli_curcumaopt_07` golden values and
  tighten the `test_gfnff_sqe` 2c tolerance.
- Docs not written here: a Known-Issues entry in CLAUDE.md / `docs/GFNFF_STATUS.md` (the Cl2-
  "free-ion gradient residual" and the CH4+H "known residual" were this bug), and an
  AIChangelog line.
- The `cn_cutoff_bohr`/topology CN vs per-step CN agree to ~1e-12 on all 2647 structures, which
  is why A and B are invisible in single points.

## 11. Fix B applied and verified (Sep 25, 2026)

Applied `git apply test_cases/revgfnff/_log/STALE_CN_fixB.patch` cleanly (4 files:
`d4param_generator.cpp/h`, `ff_workspace.h`, `gfnff_method.cpp`). No other agent was touching
`build_rev` this time (checked `ps`/`git status` first); did a full clean rebuild anyway
(`rm -rf build_rev`, fresh `cmake`+`make -j24`, 0 errors) rather than trust the handed-over
binary. **`build_rev/curcuma` md5 `999b8f902e43564ab0c4c81a23464923`.**

**`cli_curcumaopt_07_opt_multixyz` golden values regenerated.** Ran `-opt helicen.xyz -method
gfnff -threads 1 -no_bmt` with the Fix-B binary, then a FRESH `-sp` single point at each of the
17 optimised geometries (not the `-opt`-reported number, per the task's instruction). Confirmed
myself: the 17 fresh energies span only **9.6e-7 Eh** (tighter even than the ~4e-6 estimate
above) and all 17 land in the SAME minimum. This is qualitatively different from the OLD golden
file, where frames 02 and 12 were -9.126721 / -9.141227 Eh - **105.7 / 96.6 kcal/mol higher** -
i.e. the pre-Fix-B binary's native `-opt` misconverged on 2 of the 17 starting frames (same
failure class as caffeine's native-lbfgs stall in section 5, just landing somewhere instead of
stalling). Verified `-threads 4` gives the same 17 minima as `-threads 1`. New 6-decimal golden
file written; re-ran ctest with `CURCUMA=<build_rev>/curcuma` (needed - without it, ctest
silently falls back to `release/curcuma`, per the project's own documented gotcha) after
`cmake .` to refresh the build-tree copy: **20/20 PASS** (was 18/20 before regeneration, with
frames 11 and 14 the ones over tolerance).

**UPU23/2h**: not applicable to this golden-value set. `helicen.xyz`/`cli_curcumaopt_07` is a
17-conformer helicene molecule, unrelated to the GMTKN55 `UPU23` structure. UPU23/2h only
appears in section 5's separate `optset/run.sh` 11-structure sample and is not encoded as a
golden value anywhere in the test suite; there is nothing to change for it here. Stating this
explicitly rather than silently skipping the instruction.

**Regression check, isolating Fix B specifically.** Rather than diff against the shared
`_run/energies.json` caches (last populated Sep 9, before most of this branch's other WIP
work - which would conflate Fix B with everything else since), built a second binary from the
SAME tree with only Fix B reverted (`git apply -R`, incremental rebuild, copy binary aside,
`git apply` again, rebuild again - restored binary md5 verified identical to the first clean
build). This holds every other in-flight change on `reactff2-llm` fixed and isolates Fix B's
own effect. Ran `gmtkn55_compare.py --method gfnff`, `mor41_validation.py --method gfnff`,
`s30lci_gfnff_compare.py` against BOTH binaries (patched copies pointing at each binary + a
private scratch RUNDIR; `cur|*|gfnff` keys dropped before each run so "cur" is always freshly
computed by that binary; cached `xtb` entries reused unchanged). Only GFN-FF is touched by Fix
B (`forEachD4PairList`/`refreshC6WeightsForCN` are called only from `GFNFF::prepareCNAndEEQ()`;
gfn1/gfn2 use the native xTB SCF's own D4 path and never reach this code), so gfn1/gfn2 were not
re-checked.

| set | n | no-Fix-B vs xtb | with-Fix-B vs xtb | cur energies: no-B vs with-B |
|---|---:|---|---|---|
| GMTKN55 (gfnff) | 2462 (2460 scored) | MD=-0.474 MAD=0.860 max=131.608 RMSD=8.216 kcal/mol | identical to 3 decimals | 0/2462 differ at the 8-decimal (1e-8 Eh) CLI print resolution |
| MOR41 (gfnff, reaction-level) | 41 reactions / 95 structures | MD=-5.035 MAD=15.990 max=139.970 RMSD=35.733 kcal/mol | identical to 3 decimals | 0/95 differ at 1e-8 Eh |
| S30L-CI (gfnff) | 30 complexes / 90 fragments | MAD=0.387 max=9.934 RMSD=1.818 kcal/mol | identical; `results.csv` byte-for-byte identical | n/a (byte-identical CSV) |

Fix B changes nothing measurable in any of the three sets, at the resolution the CLI exposes
(~1e-8 Eh per structure; 3 printed decimals on the aggregate kcal/mol stats) - an independent
confirmation of section 4's full-precision finding, though not a reproduction of the exact
1.5e-10 Eh digit (these CLI-driven scripts only print 8 decimals, so 0/2557 differing at all is
compatible with, not a re-measurement of, that sub-1e-8 number).

**Aside, not caused by Fix B, out of scope here**: the absolute GMTKN55/MOR41-vs-xtb MAD above
(0.86 / 16.0 kcal/mol) is markedly higher than the ~0.26 kcal/mol documented in CLAUDE.md Known
Issue #25 for a clean GFN-FF port. Since the no-Fix-B binary (which already contains everything
else currently in this WIP tree) shows the identical elevated numbers, this predates Fix B
entirely and is presumably this branch's other in-flight rev-gfnff work (e.g. the mg3/well-form
default change in the recent commit log) trading vanilla-xtb port fidelity for physical accuracy
elsewhere. Not investigated further; it does not change the Fix-B verdict (identical before and
after).

**Full `ctest`** (`CURCUMA=<build_rev>/curcuma ctest -j8`): **294/306 passed, 12 failing** (was
13). The 12: `confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`, `cli_confscan_01..07` (7
tests), `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `cli_simplemd_20_gfnff_rev_h_budget` - exactly
the pre-existing documented set MINUS `cli_curcumaopt_07_opt_multixyz`, which now passes (20/20
sub-checks) after the golden-value regeneration above. `gfnff_sqe` and the new
`cli_gfnff_06_stale_cn_energy_only` both pass. No new failures introduced anywhere in the suite.

**Fix B is fully in**: patch applied to the working tree (not reverted), `build_rev/curcuma` md5
`999b8f902e43564ab0c4c81a23464923` (rebuilt clean, and separately confirmed to reproduce
byte-for-byte after a revert+reapply+rebuild round trip). Nothing left unresolved from this
task; the two section-10 suggestions this task did not cover (tightening `test_gfnff_sqe`'s 2c
tolerance; the Known-Issues/AIChangelog doc entries) remain open for the orchestrator.
