# ConfSearch — Roadmap, Open TODOs & Known Issues

Status of the ConfSearch efficiency/robustness work (branch `confsearch`, Jun 2026).
All features below are 🤖 AI-generated, ⚙️ machine-tested only — **human production testing pending**.

## Done (committed)

- **Phase A** (`bd10da3`): temperature-triggered RATTLE (`rattle_threshold_temp`/`rattle_hot_mode`),
  opt-in topology abort (`topo_check`), opt-in Epot abort (`epot_abort`/`epot_abort_window`),
  seed energy window vs. global minimum + funnel (`seed_energy_window`/`seed_window_schedule`),
  optimised minima fed back into the bias pool (`opt_feedback_bias`, persistent/prune-exempt).
- **Phase B** (`0067b26`): symmetry-permutation-aware RMSD-MTD bias — ConfScan's reorder rules are
  cached and applied as a **smooth sum-over-images** (one Gaussian per symmetry image, NOT a hard
  min → no discontinuous force). Weighted-Kabsch foundation (`RMSDDriver` per-atom weights).
- **Phase C "couple"** (`6d6fccd`): `bias_calibration=couple` sets the MTD hill width from the
  dedup RMSD scale (`alpha = ln2/(bias_couple_factor*rmsd)^2`).
- **Phase C "cluster"/"weighted"** (pending commit): learn the basin radius and per-atom RMSF from
  opt+filter clustering; guarded (≥2 distinct minima, d_inter clamp) so they degrade gracefully.
- **RATTLE mode-2 fix**: the activation gate was `if (rattle == 1)`, so `-rattle 2` (X-H only) ran
  plain Verlet with 0 constraints despite being documented and handled by `InitConstrainedBonds`.
  Now `rattle == 1 || rattle == 2`. Cleaned up the constraint printout (summary + element-labelled
  per-bond list).

## Open TODOs / Roadmap

### 1. Verbosity ownership rework (BIG, cross-cutting)
**Problem:** the global `CurcumaLogger` verbosity is mutable shared state. Every `CurcumaMethod`
sub-object (`ConfScan`, `SimpleMD`, the Opt drivers) sets it to its own level in its constructor
and never restores it; `OptimizerDriver::Optimize` and the energy-method setup clamp it to 0. Net
effect: after the first sub-call the parent's level is gone, so ConfSearch's per-cycle logs were
invisible at default verbosity, and `SimpleMD::InitConstrainedBonds` had to fall back to `std::cout`
(its `CurcumaLogger` calls were silently dropped even at `-verbosity 3`).
**Patched locally** (not architecturally): RAII restore in `OptimizerDriver::Optimize`; ConfSearch
re-asserts `m_verbosity` after each sub-phase. **Proper fix:** make verbosity scoped — e.g. a RAII
guard in the `CurcumaMethod` base (save in ctor, restore in dtor) so a child restores the parent's
level on destruction, and stop the energy-setup path from leaking 0. Once done: move the RATTLE
report (and the ConfSearch re-asserts) back to `CurcumaLogger`, gate per-bond lists at verbosity ≥ 3,
and drop the local patches. **Anchor:** `FIXME` in `SimpleMD::InitConstrainedBonds` (simplemd.cpp).

### 2. FIXME: `CitationRegistry::cite` thread race (crash at threads>1)
`GFNFFComputationalMethod::calculateEnergy` calls `CitationRegistry::cite` on every energy
evaluation; under ConfSearch `threads>1` several MD workers call it concurrently → data race in the
global registry → SIGSEGV (in `cite`) / SIGABRT (malloc), triggered during the NaN/instability
cleanup. Reproduced: triose `threads=4 startT 600`; `threads=1/2` ok. **Workaround:** run ConfSearch
MD with `threads=1` and use the 16 cores via parallel processes. **Fix options:** (a) make
`CitationRegistry::cite` thread-safe (`std::call_once`/mutex); (b) cite once at method init, not in
the hot `calculateEnergy` path.

### 3. Phase C validation (experimental)
- `cluster`/`weighted` are bounded/guarded but **unvalidated**: the basin-radius and RMSF estimates
  rest on a heuristic assignment (RMSD + energy tol) that lumps everything into one cluster when only
  one minimum is found (then calibration is skipped). Needs runs with several distinct minima.
- The dedup (ConfScan) still uses the **unweighted** RMSD; for full consistency the RMSF weights
  should also feed ConfScan's metric. Not yet plumbed.
- `bias_couple_factor`/`bias_energy_tol` defaults are guesses — tune against real benchmarks.

### 4. MD stability with wide hills — PARTIALLY ADDRESSED (Jun 2026)
Wide `couple`/`cluster` hills + a dense shared pool sum many simultaneous Gaussians → large
cumulative bias force → possible NaN blow-up at high T (seen at `threads=4`, masked by issue #2).
**Root cause for the cross-run heating** (`-startT 500`, "MD temperature climbs with every MD
run"): the shared pool persists across all runs/cycles, the hill height `W_i = k·counter_i` grows
on every visit (never capped/reset — `pruneByCounter(1)` removes nothing), and `wtmtd` damps only
the COLVAR energy, not the force. So each run inherits a taller, non-conservative bias that
out-paces the thermostat → `<T>` climbs run-by-run until a NaN blow-up (baseline triose: 552→596
→654→explode).

**Implemented (opt-in, default off → behaviour-neutral):**
- `rmsd_mtd_max_height` (Int, 0=unbounded): cap the counter in the force, `W_i = k·min(counter_i,
  cap)`. Strongest bound (triose cap=3 holds `<T>`~500–525 K).
- `rmsd_mtd_freeze_inherited` (Bool): freeze the heights of structures present at a run's start
  (snapshot in `prepareRun()`); only this run's own deposits grow (and inherited ones are not
  bumped). Bounds the cross-run escalation while keeping geometry sharing.
- `temp_abort` + `temp_abort_factor`(1.5) + `temp_abort_delta`(300 K): abort an MD run when the
  running-mean `<T>` runs away from the target (mirrors `epot_abort`; either threshold trips,
  `<=0` disables one). Catches the runaway before the NaN.

**Also fixed here (pre-existing latent bug):** `topo_check`/`epot_abort`/`temp_abort`/the new MTD
bounds are SimpleMD-owned PARAMs, so a flat CLI flag (`-topo_check`, `-temp_abort`, …) is
auto-routed to `controller["simplemd"]`, which `ConfSearch` never consumed (it only forwarded its
own `controller["confsearch"]` members, all left at `false`) → the flags were silently ignored.
`ConfSearch::start()` now reads each gate from `controller["simplemd"]` (fallback = ConfSearch
default). The three gate aborts also print via `fmt::print` (visible despite issue #3) instead of
the swallowed `warn_fmt`.

**Still open:** the *intra-run* blow-up with very wide `couple`/`cluster` hills at high T
(coupling `k`↔`alpha`, total-force capping, or a true well-tempered force damping) is not yet
addressed — the caps/abort above bound the cross-run heating and provide a safety net, but a
single run with an over-wide hill can still spike. 🤖 AI-generated, ⚙️ machine-tested only.

## Verification anchors
- Bit-identity of the Phase-B MTD refactor (perms off): byte-identical `.mtd.xyz`/`.bias.xyz` vs the
  pre-B8 build on a fixed-seed ConfSearch (`-seed 42 -threads 1`).
- Regression: `ctest -R "cli_simplemd_|cli_confscan_|cli_curcumaopt_|cli_rmsd_"`.
