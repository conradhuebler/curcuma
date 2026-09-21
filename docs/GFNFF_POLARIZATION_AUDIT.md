# GFN-FF Polarization / Geometry-Response Audit (June 2026)

**Question:** Where does curcuma's native GFN-FF silently stop being a *polarizable*
force field — i.e., where are EEQ charges (or other geometry-dependent quantities:
CN, dCN/dx, distances, D4 C6, pi-bond-orders, topology) computed once at the start
and then frozen/reused, instead of being recomputed so they respond to the current
geometry?

**Bottom line:** The normal default path (CPU and GPU) is correctly polarizing —
charges, CN, dCN, distances and D4 are re-solved/recomputed at every geometry. The
caches in the hot path are true *warm-starts* that converge to the current-geometry
solution. There is exactly **one default-active silent freeze**:
`m_phase2_historically_implausible`. Everything else that freezes is either opt-in +
documented (`static_charges`/`static_cn`/`static_all`/`skip_phase2`) or
correct-by-design (HB/XB charges held at Phase-1, matching the Fortran reference).

Verdict scheme: **POLARIZING** = re-solved every geometry · **WARM-START** = re-solved
every geometry, cache only speeds it up (converges to the same current solution) ·
**FROZEN** = computed once / reused stale across geometry changes.

## Verdict table — every cache/reuse site

| # | Site (file:line) | Quantity | Default | Verdict |
|---|---|---|---|---|
| 1 | `gfnff_method.cpp:1735` computeSharedDistances | distance matrix | always | POLARIZING |
| 2 | `gfnff_method.cpp:1151-1207` prepareCNAndEEQ | CN (`m_last_cn`) | on | POLARIZING |
| 3 | `gfnff_method.cpp:1300-1325` CN derivatives | dCN/dx | on | POLARIZING |
| 4 | `gfnff_method.cpp:1344/1410` calculateFinalCharges | Phase-2 charges | on | POLARIZING |
| 5 | `eeq_solver.cpp:3559-3577` RHS build | EEQ RHS (chi+cnf·sqrt(CN)) | always | POLARIZING |
| 6 | `eeq_solver.cpp:2089-2103` Cholesky cache hit | A-factor (re-solves current RHS) | `eeq_refactor_eps_bohr=0.05` | WARM-START |
| 7 | `eeq_solver.cpp:3420-3539` matrix cache | Coulomb off-diagonal block | `eeq_matrix_rebuild_eps_bohr=0.0` (off) | WARM-START |
| 8 | `eeq_solver.cpp:1657/1862` PCG warm-start | PCG initial guess | on (PCG) | WARM-START |
| 9 | `d4param_generator.cpp` updateCNValuesForGradient | D4 C6 + dc6dcn | `d4_cn_cache_threshold=0.01` | WARM-START |
| 10 | `gfnff_gpu_method.cpp:645-710` GPU EEQ | GPU factor (re-solves current RHS) | `eeq_rmsd_threshold=0.0` | WARM-START |
| 11 | `gfnff_method.cpp:8781-8838` topo.json cache | Phase-1 topology data only | on if `.topo.json` | POLARIZING (P1 frozen by design) |
| 12 | `gfnff_method.cpp:879-932` Tier-1 topology cache | bonds/angles/torsions, hybridization, **pi_bond_orders**, rings | on (>0.5 Bohr rebuild) | WARM-START (see risk #3) |
| 13 | `gfnff_method.cpp:1033-1048` Tier-2 dynamic | CN in topo struct | on | POLARIZING |
| 14 | `forcefieldthread.cpp:2371-2382` Coulomb q read | per-pair q_i,q_j | on | POLARIZING |
| 15 | `forcefieldthread.cpp:2525+` HB/XB q read | H/X-bond charges (Phase-1) | on | FROZEN BY DESIGN (correct) |
| 16 | `eeq_solver.cpp:993/3103` historically_implausible | all Phase-2 charges | trigger-based | **FROZEN (genuine trap)** |
| 17 | `gfnff_method.cpp:569/1138` static_charges/cn/all | charges and/or CN+dCN+D4 | false (opt-in) | FROZEN (opt-in) |
| 18 | `gfnff_method.cpp:1226` m_skip_eeq_recalc | Phase-2 charges | false (no runtime setter) | FROZEN (dead/diagnostic) |
| 19 | `eeq_solver.cpp:817/984` skip_phase2 PARAM | Phase-2 charges | false (opt-in) | FROZEN (opt-in) |

## FROZEN (non-polarizing) cases — ranked by risk to a normal MD/opt run

### 1. `m_phase2_historically_implausible` — HIGHEST RISK, the genuine silent trap
- **Set at** `eeq_solver.cpp:3692, 3725, 3796, 3811`; **consumed at** `eeq_solver.cpp:993-998`
  and `3103-3108`.
- **Trigger (any one):** (a) all solvers fail (empty return) with no valid cached charges;
  (b) Phase-2 returns NaN/Inf; (c) `max|q2| > max(10, 5·max|q1|)`; (d) Coulomb-energy
  ratio `|E2-E1|/(|E1|+1) > 5`.
- **Behavior:** sets the flag + `m_phase2_implausible_natoms=natoms`. From then on **every**
  call with that atom count returns Phase-1 topology charges and skips the solve. The flag
  is **never reset** in-process (no reset code exists); the only escape is the atom count
  changing.
- **Scope:** per-`EEQSolver` instance (per process), not persisted to disk.
- **Consequence:** a single bad geometry mid-trajectory (transient close contact / near-singular
  EEQ matrix during an aggressive opt step or hot MD step) **permanently** demotes the run to
  fixed Phase-1 charges. The FF stops polarizing for the rest of the run; Coulomb forces become
  geometry-frozen-charge forces. The only signal is a one-time `warn()` (verbosity >=1) plus an
  info line (>=2). Activates on its own (no opt-in), silent after the first warning, permanent.
  Thresholds are loose so it rarely trips on well-behaved systems — but when it does it is sticky
  and invisible.
- **Recommended fix:** make the demotion *transient* — fall back to Phase-1 for the offending
  step only and re-attempt Phase-2 on the next geometry (the bad geometry usually passes), or
  reset the flag when a later solve succeeds. Connects directly to the EEQ solver-robustness
  work (a more robust solver produces fewer implausible Phase-2 results in the first place).

### 2. `static_charges` / `static_all` (#17) and `skip_phase2` (#19) — opt-in foot-gun
- User-requested (`-static_charges`, `-static_all`, `-skip_phase2`). `static_all` also freezes
  CN/dCN/D4. Documented "equilibrium dynamics only / invalid for charge-transfer or ionic
  dynamics." `static_charges` emits a startup `warn()`. Risk only if misapplied outside the
  valid regime. Not a silent default trap.

### 3. `pi_bond_orders` in the Tier-1 topology cache (#12) — low but non-zero
- `pi_bond_orders` (`gfnff_method.cpp:9045`) is geometry-dependent (Hückel solver on the current
  geometry) but lives in Tier-1, so it is refreshed **only when a full topology rebuild fires
  (>0.5 Bohr max displacement)**. Between rebuilds the conjugation pattern feeding the N-angle
  `f2` factor and torsion `sumppi` is held constant. For normal MD/opt (no atom >0.5 Bohr from
  the rebuild reference) the error is small/second-order, but it is the one geometry-dependent
  quantity not refreshed every step (unlike CN). Worth flagging as the "parametrize-then-static"
  pattern even though the impact is minor.

### Not a bug: HB/XB charges at Phase-1 (#15)
Both the Fortran reference and the curcuma CPU path evaluate the H/X-bond term with charges
frozen at topology build (Phase-1), `forcefieldthread.cpp:2525+`. The GPU matches this. The
opt-in `CURCUMA_GFNFF_GPU_RESIDENT_HBQ=1` path that refreshes HB charges to live values actually
*deviates* from the reference (`ff_workspace_gpu.cu:2165-2172`) and is OFF by default.

## Geometry-dependent quantities — recompute status

- **CN / dCN/dx** (#2/#3): recomputed every `Calculation()` (the Feb-2026 "CN derivative computed
  once at init" bug is fixed). GPU `computeCN` runs every step unless `m_frozen_cn` (only set by
  `static_cn`). POLARIZING.
- **Distance matrix** (#1): recomputed every call. POLARIZING.
- **D4 C6 / dc6dcn** (#9): recomputed from current CN, 0.01-CN-delta warm-start skip. WARM-START.
- **pi_bond_orders** (#12): geometry-dependent but Tier-1 (>0.5 Bohr rebuild) — risk #3 above.
- bonds/angles/torsions equilibrium params and hybridization: legitimately fixed-at-init topology.

## Fixes applied (June 2026)

1. **`m_phase2_historically_implausible` made transient** (`eeq_solver.cpp`): removed the two
   consume-site "skip forever" returns; the flag is reset on the next plausible solve. Phase-2
   is re-attempted every step; a bad step still falls back to Phase-1 for that step only.
2. **Reaction-field cache-staleness bug fixed** (exposed by #1): the SchurCholesky factor cache
   (`m_chol_cache`, `eeq_refactor_eps_bohr`) keyed on geometry+CN but the matrix is
   `A_nn + B(reaction field)`; under ALPB solvation B changes while geometry/CN do not → stale
   cache → ~3 mEh-wrong solvated charges. The freeze had masked it. Fix: `!m_reaction_field`
   added to `cache_size_ok` and `can_persist` so solvation bypasses the cache. Verified:
   `gfnff_solv_*` 28/28 pass, gas-phase bit-identical, `lu`/`pcg`/`ldlt` unaffected.

## Fix applied (July 2026) — same cache, the general case

3. **Phase-1/Phase-2 cache-key collision.** Fix #2 above patched one *symptom* of a deeper
   defect: `solveWithSchurCholesky` has no explicit notion of which phase called it. It infers
   "I am Phase 2" from `m_pending_geometry`/`m_pending_cn` being non-empty (`eeq_solver.cpp`
   `cache_size_ok` ~:2113-2119, `can_persist` ~:2166-2169). Those buffers are written **only**
   by `calculateFinalCharges` (Phase 2) and were **never cleared**, so the documented contract
   ("Phase 1 arrives with empty pending buffers") held only for the *first* solve on a given
   solver instance. From the second topology build onward Phase 1 looked like Phase 2.

   The damaging variant is the one that follows an `invalidateCholeskyCache()` — exactly what
   `GFNFF::getCachedTopology()` does before every full topology rebuild: the factor is dropped
   but the *key* survives, so Phase 1 passes `can_persist`, stores its **topological-distance**
   factor, and Phase 2 of the same call then hits it (`max_dr == 0`) and solves the geometric
   system with it.

   Fix: clear both buffers at the top of `calculateTopologyCharges()` **and** inside
   `invalidateCholeskyCache()`, so the factor and its key are always dropped together. This
   subsumes the `!m_reaction_field` special case in #2 (which is retained).

   Measured on `test_gfnff_topology_reentrancy` (new regression guard): before the fix, a second
   topology build across a cache invalidation shifted H2O Phase-1 `qa` by **4.5e-2 e**, Phase-2
   charges by 4.5e-2 e and `dgam` by 6.7e-3 (CH4: 6.6e-3 / 7.6e-3 / 1.8e-3); after, all
   differences are exactly 0.

   **Scope, honestly stated:** no production-path impact could be demonstrated. Single points,
   `-opt` and a 300-step caffeine `-md` (154 full topology rebuilds) are **bit-identical**
   before vs after, because on a rebuild the inline topology cache restores the Phase-1 charges
   and skips the Phase-1 EEQ solve entirely, so the corrupted path is not reached. The bug is
   real and reproducible at the API level; the fix removes a latent trap and restores the stated
   invariant, but it is not known to have changed any published number.

---
*Audit by background investigation agent, June 2026. Source-of-truth = the verdict table above;
re-verify file:line before acting (line numbers drift).*
