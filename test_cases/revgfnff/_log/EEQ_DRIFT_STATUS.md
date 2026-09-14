# curcuma EEQ charge drift along the class-A scans (EEQ half of the Hirshfeld falsifier)

AI-generated (measurement, 2026-09-13), machine-evaluated. No `src/` changes, no build.
Binary: frozen copy of `build_rev/curcuma`, **md5 071361107efba0bc346f40dc1a276725**.
Scripts + raw data: `scratchpad/eeqdrift/{drive_nc.py,fast.py,analyze.py,final.py}` with
`charges.json`, `fast.json`, `final.json`.

## Method

- **Protocol** = `CLASSA_FROZENCN.md` mode 1 (react / kept-r_eq topology): frames =
  `ref/_geom/<mol>/ref.xyz` + the 20 class-A grid points **ascending r** (`r/r_eq` 0.75..3.5,
  index 5 == r_eq); `-method revgfnff -gfnff.topology_mode react -batch true
  -batch_reuse_topology true -gfnff.cache_topology false -verbosity 4 -threads 1`.
- **Charges**: `gfnff_diag_charges.json` (`topo_qa` = Phase 1, `energy_charges` = Phase 2) read
  from a **prefix run** `[ref]+points[:k+1]` (fresh temp dir, cache off) = charges of point k.
  Every prefix run's batch energies equal the full run's for frames 0..k: **all identical**.
  Cache-free reruns of everything: bit-identical to the cached runs.
- **Reference**: `ref/H/<sys>_rks/energies.json` `hirshfeld`, same point index (rks).
- **Protocol proof**: my ΔCoulomb (r_eq -> 3.5 r_eq, react) = C=C -35.5, O-O -35.7, N-O -26.1,
  C=N -22.3, N-N -20.1, N=N -19.3, O-Cl -15.6, C-O -11.3 kcal/mol = the
  `FABLE_ROADMAP_REVIEW.md` Sec 2.1 list (-35/-34/-26/-22/-20/-19/-16/-11) to <= 1 kcal/mol.
- Charges are **step functions** (jumps 0.05-0.28 e in one grid step) at react topology-corner
  events (h2o2 log: O1-H4 formed r_scan 2.972 A, blend revert, then symmetric) - **not** at the
  scanned bond's break.

## Table 1 - dq (r_eq -> last bonded grid point), rks Hirshfeld vs curcuma Phase-1 EEQ.

| bond | r_eq A | far A | r/r_eq | break r A (r/r_eq) | dqEEQ_i | dqREF_i | dqEEQ_j | dqREF_j | max abs diff | 0.1 e | abs ratio |
|---|---:|---:|---:|---|---:|---:|---:|---:|---:|---|---:|
| C=C  | 1.3269 | 4.644 | 3.50 | 5.015 (3.78) | -0.284 | -0.016 | -0.284 | -0.016 | **0.267** | FAIL | 17.7x |
| N-O  | 1.4523 | 4.357 | 3.00 | 4.391 (3.02) | -0.341 | +0.048 | +0.009 | -0.046 | **0.388** | FAIL | 7.2x |
| C=N  | 1.2672 | 4.435 | 3.50 | 4.789 (3.78) | -0.316 | -0.054 | -0.003 | +0.078 | **0.262** | FAIL | 5.8x |
| N-N  | 1.4923 | 4.477 | 3.00 | 4.512 (3.02) | -0.188 | +0.011 | -0.188 | +0.011 | **0.199** | FAIL | 16.8x |
| N=N  | 1.2382 | 4.334 | 3.50 | 4.680 (3.78) | -0.103 | +0.015 | -0.103 | +0.015 | **0.118** | FAIL | 6.7x |
| O-O  | 1.4694 | 4.041 | 2.75 | 4.165 (2.83) | -0.108 | +0.011 | -0.108 | +0.011 | **0.120** | FAIL | 9.6x |
| O-Cl | 1.7236 | 5.171 | 3.00 | 5.212 (3.02) | -0.074 | +0.041 | -0.000 | -0.033 | **0.114** | FAIL | 1.8x |
| C-O  | 1.4303 | 4.291 | 3.00 | 4.325 (3.02) | -0.131 | +0.022 | +0.047 | -0.068 | **0.152** | FAIL | 6.1x |
| O-H scale | 0.9618 | 2.885 | 3.00 | 3.272 (3.40) | +0.185 | +0.024 | -0.285 | +0.003 | **0.288** | FAIL | 7.7x |
| C-H scale | 1.0914 | 3.274 | 3.00 | 3.712 (3.40) | -0.009 | +0.056 | -0.018 | -0.084 | 0.066 | OK | 0.2x |

- **Falsifier fails for all 8 drifting bonds** (0.114-0.388 e >> 0.1 e); only C-H passes. EEQ
  moves 1.8-17.7x the physical charge; O-O / N-O / N-N / C=C even move the **wrong way**.
- Bond presence: break at r/rcov ~ 1.66 for every bond. O-O first (2.83 r_eq); C=C / C=N / N=N
  reach the grid end (3.5 r_eq); the rest break at 3.02-3.40 r_eq (far point 3.00). None unreadable.
- Phase 2 differs materially on heteronuclear bonds (N-O j +0.009 vs -0.044, C=N j -0.003 vs
  +0.062, C-O j +0.047 vs +0.007, H2O O-H i +0.185 vs +0.141); symmetric bonds agree <= 0.05 e.
  Phase-1 numbers above are primary (as instructed).

## Table 2 - is the drift the charge movement? (kcal/mol; `gfnff-fast` = CN+charges frozen at r_eq)

| bond | dCoul obs | frozen charges | charge part | charge/obs | pair est q_iq_j d(1/r) | abs factor |
|---|---:|---:|---:|---:|---:|---:|
| C=C  | -35.5 | -0.2 | -35.2 | 99% | +5.2 | 6.7 |
| N-O  | -26.2 | +0.2 | -26.3 | 101% | -4.2 | 6.3 |
| C=N  | -22.3 | +5.1 | -27.4 | 123% | +12.3 | 2.2 |
| N-N  | -20.0 | -1.6 | -18.5 | 92% | -6.1 | 3.0 |
| N=N  | -19.3 | -0.5 | -18.9 | 98% | -5.6 | 3.4 |
| O-O  | -35.9 | -0.3 | -35.7 | 99% | -11.6 | 3.1 |
| O-Cl | -15.6 | +3.1 | -18.6 | 120% | +4.3 | 4.3 |
| C-O  | -11.2 | +5.2 | -16.4 | 146% | +9.2 | 1.8 |
| O-H  | +7.0 | +31.7 | -24.8 | net 0.3x | +66.3 | 0.4 |
| C-H  | -0.2 | +0.3 | -0.4 | both ~0 | +0.6 | 0.7 |

- **Yes: 92-146% of the Coulomb drift is the charge movement**; the geometry-only (frozen-charge)
  part is <= 5.6 kcal/mol everywhere. The bonded pair's own term is only 1.8-6.7x the total and
  changes sign, so the drift is the whole charge redistribution incl. the EEQ per-atom self-energy.
- Stage-2 decision input: the drift is the charge model's response, not bond physics (reference
  moves <= 0.078 e where EEQ moves 0.10-0.34 e). Fitting bond D to it would fit a charge-model error.

## Class-S fix - `scripts/revgfnff_fit.py` (not `revgfnff_ref.py`); 2 one-line edits, verified

- `_topology_index`: `cls == "C"` -> `cls in ("C", "S")`. All four `ref/S/*` now select **idx 19 =
  largest scan value d=6.00 A** (separated geometry); before: smallest (d=2.30, contact geometry).
- `_order_points`: `reverse=(cls == "C")` -> `("C", "S")`. Required because the default mode is
  `react`, where the sort (not `_topology_index`) decides frame 0; without it frame 0 stayed
  d=2.30 even though `_topology_index` named the d=6.00 frame. Now d=6.00 first in react and static.
- Regression check over 230 systems of classes A/C/E/H/S, old vs new: exactly 8 differ, all class S
  (4 systems x 2 modes), no other class changes.
- The trailing null `monomer_A/B` entries still in `ref/S/*/energies.json` are dropped by the
  loader (no gradient) so they never reach `_topology_index`; if they did, `all(v is not None)`
  would fall back to index 0 (contact geometry). Not fixed here - remove them or guard the caller.

## Open / not measured

- Only the 8 drifting bonds + C-H / O-H; the other class-A bonds are not in this table.
- dq is read at the last bonded grid point, never extrapolated; grid spacing (0.2-0.5 r_eq) is
  coarser than the charge jumps, so a jump position is resolved to one grid step.
- The step-like behaviour is observed (react log events), not root-caused to a specific corner
  rule; per-corner EEQ (stage 1b) is the likely owner.
