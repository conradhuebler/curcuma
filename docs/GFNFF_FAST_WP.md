# Fast GFN-FF — Work Package (Draft, June 2026)

**Goal:** a deliberately approximated, faster GFN-FF variant for large systems / long MD,
trading the per-step polarization response (and some long-range accuracy) for speed —
exposed as a named, opt-in method, NOT a silent demotion.

This is the *controlled* version of the one accidental freeze trap documented in
[GFNFF_POLARIZATION_AUDIT.md](GFNFF_POLARIZATION_AUDIT.md) (`m_phase2_historically_implausible`):
instead of a run silently losing polarization after a bad step, the user opts into a
clearly-labelled fast mode whose approximations are known and warned about.

## Why
The default GFN-FF re-solves Phase-2 EEQ charges and recomputes CN/dCN/D4 every geometry
step. On large systems the O(N^2) D4 pair build dominates (mixture2: 6200 atoms ~90 s,
~13.7M D4 pairs). For pre-equilibrated large systems where the electrostatic response to
small geometry changes is second-order, freezing the charge/CN parametrization + cutting
the O(N^2) terms is a large speedup at bounded accuracy cost.

## Building blocks (already exist — see polarization audit)
- `static_charges` (`gfnff_method.cpp:569`) — freeze Phase-2 EEQ charges after the first
  solve; skip the per-step EEQ solve.
- `static_cn` — freeze CN + dCN/dx + D4 (C6 / dc6dcn).
- `static_all` — both of the above.

`gfnff-fast` is essentially `static_all` + distance cutoffs + relaxed caches, packaged as
one method so it round-trips through method resolution and BMT naming.

## Proposed interface
- Method string **`gfnff-fast`** (mirrors `gfnff-d3`, `uff-d3`), resolved in MethodFactory
  to native GFN-FF with the fast preset. (Alternative: a `-gfnff.fast true` flag; a named
  method is preferred for round-trip/BMT.)

## Bundle (preset)
1. **Charges**: parametrize EEQ once at the initial geometry, reuse for all steps
   (`static_charges`). The dominant correctness sacrifice (non-polarizing).
2. **CN / D4**: freeze CN/dCN and the D4 C6/dc6dcn (`static_cn`) — skips the per-step
   13.7M-pair D4 build, the measured perf killer.
3. **Cutoffs** (ranking pending the perf investigation, see
   [GFNFF_PERFORMANCE_LEVERS.md](GFNFF_PERFORMANCE_LEVERS.md)): D4 dispersion + Coulomb pair
   cutoffs with neighbor lists, turning the O(N^2) walls into O(N).
4. **Relaxed caches**: larger `eeq_refactor_eps_bohr` / topology-rebuild thresholds.

## Accuracy / validity caveats (MUST document + warn at runtime)
- **NON-POLARIZING**: invalid for charge-transfer, ionic dynamics, large conformational
  change, or reactions. Intended for equilibrium dynamics / relaxation of pre-equilibrated
  large systems only.
- Distance cutoffs violate Hellmann-Feynman vs the full energy and degrade MD energy
  conservation (already documented for `eeq_distance_cutoff`) — quantify the drift before
  defaulting any cutoff on.
- Emit a startup `warn()` (as `static_charges` already does).

## Validation plan
- Energy + gradient of `gfnff-fast` vs full `gfnff` on `test_cases/eeq_mixture_fractions/`
  (mixed/water/urea, 50–3007 atoms): tabulate speedup AND accuracy loss per approximation.
- NVE MD energy-conservation: full vs fast on a fraction (drift/ns).

## Depends on / open
- `docs/GFNFF_PERFORMANCE_LEVERS.md` (perf investigation) — defines which cutoffs to bundle
  and their measured impact.
- Operator decision: method name `gfnff-fast` vs flag; which approximations are default-on.

## Status — IMPLEMENTED (June 2026, AI/machine-tested)
`-method gfnff-fast` is live: `method_factory.cpp` routes it to native GFN-FF with the
`static_charges=true` + `static_cn=true` preset and emits the non-polarizing warning. SP is
bit-identical to `gfnff` (mixed_0051, mixture2 — the capture happens after the single solve);
the speedup is in MD/opt where the per-step EEQ + CN/D4 are frozen after step 1.

Measured (mixed_0502, 40-step NVE, 4 threads): `gfnff` 5.9 s vs `gfnff-fast` 5.0 s (~15%).
The win grows on large many-fragment systems where the EEQ Phase-2 solve is a larger share.
GPU paths (`-gpu cuda|rocm`) honor the preset too (routed through with the static flags).

**Note:** the HB-list cost (Lever 1) is now reduced for ALL gfnff runs (exact cell-list
generation, ~1.5x on mixture2 SP) — see `docs/GFNFF_PERFORMANCE_LEVERS.md` — so `gfnff-fast`
adds the charge/CN-freezing win on top. The accuracy-changing distance cutoffs (D4/Coulomb)
remain an OPEN optional extra for `gfnff-fast` (not yet bundled; would need a smooth taper to
avoid the ~0.1 Eh energy jump documented in the perf levers doc).

Human production testing (long-MD stability, accuracy-vs-full on real systems) pending.
