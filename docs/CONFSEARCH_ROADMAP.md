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

### 1. Verbosity ownership rework — LARGELY RESOLVED (Jun 2026)
**Problem:** the global `CurcumaLogger` verbosity is mutable shared state. Every `CurcumaMethod`
sub-object set it to its own level in its ctor and never restored it; the optimizer + energy-method
setup clamp it to 0. After the first sub-call the parent's level was gone, so ConfSearch's per-cycle
logs were invisible and `SimpleMD::InitConstrainedBonds` had to use `std::cout`.

**What was done:**
- **Base RAII** in `CurcumaMethod` (ctor captures the parent's global level, dtor restores it). This
  scopes the level for all standalone commands and same-thread nesting (e.g. ConfSearch's main-thread
  `ConfScan` filter restores the ConfSearch level on its own).
- **Thread-pool boundary:** the global verbosity is a non-thread-safe shared static; the MD and
  batch-opt phases run/destroy their sub-objects on `CxxThreadPool` worker threads, where the base
  RAII restore happens on the *worker* thread and the static is left clamped (verified: 1 worker →
  level intact, 2 workers → left at 0). So the pool-owning helpers `ConfSearch::PerformMolecularDynamics`
  and `PerformOptimisation` re-assert `m_verbosity` on the main thread before returning. The 7
  per-phase `set_verbosity` re-asserts in `ConfSearch::start()` were removed (base RAII + these two
  helper-boundary restores cover them).
- **Energy-setup boundary:** `EnergyCalculator` construction + `setMolecule` (GFN-FF param-gen) still
  leaves the global clamped to 0, so `SimpleMD::Initialise` re-asserts `m_verbosity` right after the
  energy-method setup; the `InitConstrainedBonds` RATTLE report is now on `CurcumaLogger`
  (summary ≥ 1, per-bond/angle list ≥ 3, and correctly silent at 0 — the old `std::cout` leaked even
  in silent mode).

**Verified:** ConfSearch per-cycle logs survive both temperature cycles and multi-worker MD without
the re-asserts; the RATTLE report shows at `-v 3` (summary at `-v 1`), hidden at `-v 0`.

**Remaining caveat (NOT fixed, out of scope):** the verbosity is still a single shared static, so a
*truly parallel* ConfSearch (`activeThreadCount > 1`) can transiently race the level across workers
(cosmetic — workers may briefly log at the wrong level mid-run; the main-thread state is restored at
the pool boundary). A fully clean fix would make the verbosity `thread_local`, but that breaks the
main→worker level propagation, so it is deferred. The two helper-boundary + one energy-setup
re-asserts are the residual "local patches" — kept deliberately, with comments, because the shared
static cannot be cleanly scoped across threads.

### 2. ~~`CitationRegistry::cite` thread race (crash at threads>1)~~ — RESOLVED (Jun 2026)
Original report: `GFNFFComputationalMethod::calculateEnergy` cites on every energy eval; under
ConfSearch `threads>1` several MD workers call it concurrently → SIGSEGV/SIGABRT during the
NaN/instability cleanup (repro: triose `threads=4 startT 600`).
**Resolution:** (a) `CitationRegistry::cite` was already fully `std::mutex`-guarded on every path
(since `9721f55`) and `Citations::database()` is a thread-safe function-local `static`, so the
registry race itself was gone; (b) a `thread_local` fast path now skips the global lock after a
thread's first sighting of a key, removing the per-step lock traffic on the hot path; (c) the
concurrent crash dumps got per-instance filenames (`<basename>.unstable.json`/`.final.json`/
`_step_*.json`) so workers no longer clobber a shared `unstable_curcuma.json`. The residual
`threads=4 startT 600` blow-up was the **bias-heating runaway** (issue #4 above / TODO #4), now
bounded. **Verified:** triose `threads=4 startT 600` runs to completion (exit 0) even with 4
concurrent instabilities (4 distinct `*.unstable.json`, no crash); with `-temp_abort`/
`-rmsd_mtd_max_height` it is fully clean and prints one citation summary. ConfSearch MD may now use
`threads>1` with gfnff when the bias is bounded. (`confsearch.cpp:111` still collapses the inner MD
to 1 thread when ConfSearch parallelises *cycles* — that is to avoid nested thread pools, unrelated.)
The shared base `curcuma_restart.json` (`CurcumaMethod::TriggerWriteRestart`) is still a single name
across instances, but it is a benign last-writer-wins dump (no crash; ConfSearch does not use MD
restart). 🤖 AI-generated, ⚙️ machine-tested only.

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

**Implemented controls:**
- `rmsd_mtd_freeze_inherited` (Bool) — **ON by default for ConfSearch**: freeze the heights of
  structures present at a run's start (snapshot in `prepareRun()`); only this run's own deposits grow
  (and inherited ones are not bumped). Bounds the cross-run escalation while keeping geometry sharing.
- `temp_abort` + `temp_abort_factor`(1.5) + `temp_abort_delta`(300 K) — **ON by default for
  ConfSearch**: abort an MD run when the running-mean `<T>` runs away from the target (mirrors
  `epot_abort`; either threshold trips, `<=0` disables one). Safety net that catches a runaway run.
- `rmsd_mtd_max_height` (Int, 0=unbounded) — opt-in: cap the counter in the force,
  `W_i = k·min(counter_i, cap)`. Tightest temperature (triose cap=3 holds `<T>`~500–525 K) but
  caps exploration, so it stays off by default.

**Default choice (measured on triose `startT 500→400, repeat 2`):** A=current/none → max⟨T⟩ 2e21 K,
**4 blow-ups**, 20 conformers; B=`freeze`+`temp_abort` → max⟨T⟩ 984 K, **0 blow-ups**, **42
conformers**; C=`cap20`+`temp_abort` → 616 K, 0 blow-ups, 25 conformers; D=`temp_abort` only →
4e10 K, **2 blow-ups**, 5 conformers. `temp_abort` alone is insufficient (sudden single-step NaNs
still blow up + constant aborting starves exploration); **freeze + temp_abort wins** (no blow-ups,
best yield), so both are ON by default. The exact `-startT 500` default command now runs with
0 blow-ups (was 4+ NaN explosions to ~1e63 K).

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
