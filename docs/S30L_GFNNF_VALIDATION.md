# S30L: Native GFN-FF vs xtb 6.6.1 GFN-FF — Validation + Fixes

Status: AI-generated, machine-tested (Jul 2026). NOT human production-tested.
The AI does not assign ✅ TESTED / ✅ APPROVED.

## Outcome (after F1/F2/CLI fixes; F3 diagnosed-deferred)

| group | n | MAD | max | RMSD | verdict |
|------|--:|----:|----:|-----:|---------|
| clean neutral (1,2,3-10,13,14,17-22) | 18 | 0.30 | 1.62 | 0.54 | reproduces xtb |
| F3 bond-term offset (11,12,15,16) | 4 | 9.7 | 12.6 | 10.2 | diagnosed, deferred |
| charged, fixed (24,25,26,27,28,29,30) | 7 | 1.5 | 3.0 | 1.6 | reproduces xtb |
| charged residual (23) | 1 | — | 27.2 | — | complex-term residual |
| **all 30** | 30 | **2.70** | 27.2 | 6.3 | was MAD 433 pre-fix |

- 22/30 systems within ~1.5 kcal/mol of xtb.
- Pre-fix (post-merge) MAD was 433 kcal/mol (charged systems off by 100s-1000s);
  the F2 + negative-charge CLI fix took charged MAD 1619 → 1.5 kcal/mol (7/8).
- xtb vs reference_s30l: MAD 14.8 (context only — GFN-FF is not meant to match the
  high-level S30L reference).

## Fixes applied

### F2 — charged host-guest complexes (FIXED)
Root cause: `GFNFF::calculateTopologyInfo` (gfnff_method.cpp) set `qfrag.assign(nfrag,0)`
and only assigned `m_charge` to fragment 0 when `nfrag==1`. For nfrag==2 (host-guest
complex) `qfrag` stayed `[0,0]`, so the per-fragment EEQ constraint forced total charge 0
regardless of `m_charge` — energy error sign/magnitude ∝ Q. Fix mirrors xtb 6.6.1
(`external/xtb/src/gfnff/gfnff_ini.f90:474-502`): at the Phase-1 EEQ solve site, for
nfrag==2 & m_charge!=0 try BOTH `[m_charge,0]` and `[0,m_charge]`, keep the assignment
with lower EEQ electrostatic energy (`calculateEEQEnergy`); nfrag>2 → charge on fragment 0.
Cache guard changed to check the SUM of qfrag (not qfrag[0], which is 0 for the
charge-on-fragment-1 winner). No change to Coulomb/D4/repulsion (they inherit charges).

### CLI negative-charge bug (FIXED, F2-enabling)
`-charge -1` was parsed as `-charge` (boolean true) + `-1` (separate flag), so
`controller["charge"].get<int>()` on boolean true returned 1 — every negatively-charged
species got +1. Root cause: the CLI value-consume condition treated any '-'-prefixed next
arg as a flag. Fix: `isFlag()` lambda (src/main.cpp) returns false for tokens that are
negative numbers ('-' followed by digit/'.'), so `-charge -1` flows into `stod("-1")`.
Verified: `-charge -1` -> "Charge -1". General fix — helps any numeric flag taking
negative values (spin, thresholds).

### F1 — NaN in angle-bending gradient (FIXED, defense-in-depth)
Root cause: `UFF::AngleBending` (forcefieldfunctions.h) had no guard on `1/sin(theta)`,
`1/|rij|`, `1/|rkj|`, and `costheta` was not clamped before `acos`. A near-collinear angle
underflowed `sin(theta)` to 0 -> `1/0=Inf` -> `0*Inf=NaN` in `derivate`, propagating to all
per-term gradients; `checkForNaN` then returned 0.0 for the energy (masking a finite,
correct energy). Fix: clamp `costheta` to [-1,1] before `acos`; epsilon-guard `sintheta`
(zero derivate when < 1e-14, mirroring `calculateDihedralAngle`); min-distance guard on
`rij.norm()`/`rkj.norm()`. Secondary (energycalculator.cpp): energy-NaN returns 0.0
(unchanged), but gradient-NaN-only now returns the finite energy + error flag instead of
substituting 0.0 (no longer lies to the optimizer). Note: post-merge system 2 no longer
triggers the NaN (origin's GFN-FF changes resolved that specific geometry); the guard is
defense-in-depth for any future near-collinear angle (UFF/QMDFF/GFN-FF).

### F3 — host bond-term offset (DIAGNOSED, DEFERRED)
Systems 11, 12, 15, 16: host energy ~0.05-0.63 Eh too negative vs xtb, mostly cancels in
ΔE, residual 4-13 kcal/mol. Term-by-term comparison for hosts 11/A and 15/A vs xtb:
- **Dispersion, Repulsion, BATM: EXACT match** (to ~1e-10 Eh) -> the bond graph
  (perception) and non-bonded exclusions are CORRECT, identical to xtb.
- **Bond stretching term: the discrepancy** (11/A: cur -19.3445 vs xtb -18.6174,
  -0.727 Eh; 15/A: cur -8.3741 vs xtb -8.2322, -0.142 Eh — curcuma too negative).
- Angle/Torsion/Coulomb: small secondary diffs.

So F3 is NOT a bond-perception bug (porting xtb's `gfnffrab`/`getnb` would be wrong and
risky — perception already matches). It is a **bond-stretching term parameter / bond-type
assignment** difference for macrocyclic hosts, most likely in the Hückel π-bond-order ->
bond-type -> (k, r0, alpha) assignment. Deferred: a targeted bond-term alignment is needed
and is risky for the 18 clean neutrals where the bond term already matches. Not fixed.

### F2 cache bug (FIXED, Jul 2026) — charged nfrag==2 complexes on cached re-runs
Found via term-by-term decomposition: charged host-guest complexes (S30L 25/26/27/28/30,
charge on fragment 0) gave the **correct** energy on a fresh run but a **wrong** Coulomb
(by ~Q²γ, i.e. hundreds-thousands of kcal) on any **cached re-run** (topology cache hit).
Fresh `25/AB` (q=+4): Coulomb +5.526 Eh (matches xtb). Cached re-run: Coulomb −0.068 Eh.

Root cause: the inline topology-cache load in `GFNFF::calculateTopologyInfo`
(gfnff_method.cpp, ~line 8947) restored `topology_charges`, `dxi`, `dgam`, `alpeeq` but
**not `qfrag`**. On a cache hit the F2 both-assignment trial is skipped (it sits inside
`if (!topology_from_cache)`), so `qfrag` stayed at the `[0,…,0]` default from line ~8829,
Phase-2 EEQ then enforced total charge 0, and the Coulomb term was computed with
neutralised charges. A fresh run sets `qfrag` via the F2 trial, so fresh ≠ cached.

Secondary: the cache **write** guard (~line 2951) validated `qfrag[0] ≈ m_charge`, which
rejected valid placements with the charge on fragment 1 (`qfrag=[0,m_charge]`, e.g. S30L
23/AB charge on the guest) — those were never cached, forcing a fresh recompute every call
and making the cached/fresh paths behave inconsistently.

Fix: (1) restore `topo_info.qfrag` and `eeq_topology_input.qfrag` from the cache on a
cache hit; (2) validate the write guard by the **sum** of `qfrag` (mirrors the load guard
at ~line 2698). Verified: fresh vs cached full-30 run are now byte-identical
(max |ΔΔE| = 0 kcal); 25/AB cached re-run gives Coulomb +5.5257 Eh (was −0.068). No
neutral-system change (qfrag=[0,0] from cache == default). The harness `results.md`
numbers are unchanged (they were from fresh runs), but any production re-use of a cached
charged complex (opt/MD) is now correct.

### System 23 (residual, diagnosed)
23/AB (q=+1): A and B match xtb **exactly** (to ~1e-7 Eh). The complex is −0.0433 Eh
(−27 kcal/mol) too negative. Term-by-term (correct xtb charges, see decomposition tool):
**Coulomb −13.8 kcal** (cur +1.1023 vs xtb +1.1242 — curcuma's Phase-2 EEQ charges
slightly too weak), **Bond −8.7 kcal** (F3-style, see below), Repulsion −2.3, BATM −0.2.
F2 placement is correct (charge on fragment 1 = guest, lower EEQ energy: E_a=+0.0041
vs E_b=−0.0632). xtb itself is 30 kcal/mol off the reference here, so 23 is a
disordered/hard system. The residual is a Phase-2 EEQ charge difference for this
charged complex, not a bookkeeping bug.

### Term-by-term decomposition (Jul 2026)
Tool: `scripts/s30l_term_decomp.py` (runs curcuma −verbosity 2 + xtb --gfnff --sp in a
clean workdir with `.CHRG`, prints per-term deltas). Note: curcuma's `Torsion + Inversion
+ sTors` == xtb's `Torsion` (xtb has no separate inversion line) — the script folds them.

Outlier ΔE drivers (cur − xtb, kcal/mol), after the F2-cache-fix:
- **23 (−27.2)**: AB Coulomb −13.8 + Bond −8.7 + Rep −2.3 (A,B exact). Charged complex.
- **11/12 (−11.8/−12.6)**: same host. Bond −456 kcal host-intrinsic (F3, Hückel
  π-bond-order → bond-type → (k,r0,α) too stiff) + Torsion +47; mostly cancels in ΔE,
  residual from imperfect cancellation. Guest 12/B matches when Torsion+Inversion folded.
- **15/16 (−9.6/−4.9)**: same host. Bond −89 + Coulomb −14 host-intrinsic (F3);
  mostly cancels.
- **7/8 (+1.3/+1.6, borderline ok)**: torsion deficit −9.7/−6.9 kcal **intrinsic to
  specific monomers** (7/A, 8/B; cur ~60-80 % of xtb), cancels in ΔE. Not a ΔE driver.
- **25/26/27/28 (±1.4)**: large Coulomb offset in AB (±800 to ±3500 kcal) that **fully
  cancels** with the charged monomer A; net ΔE small. (Pre-fix these were wrong on cached
  re-runs — now fixed.)
- **30 (−3.0)**: small charged residual.

**Common cause for the 4 biggest neutral outliers (11/12/15/16) is the BOND stretching
term** (F3): the Hückel π-bond-order → bond-type → (k, r0, α) assignment makes curcuma's
bonds ~0.05-0.73 Eh too stiff for the two macrocyclic hosts. Dispersion, Repulsion,
Coulomb (neutral), BATM match xtb to ~1e-10 Eh → bond perception/graph is correct; it is
the bond PARAMETERS. Deferred: a targeted bond-term alignment is needed and is risky for
the 18 clean neutrals where the bond term already matches. Not fixed.

## Setup / reproduce

- 30 S30L host-guest complexes, A=host, B=guest, AB=complex, Turbomole coord (Bohr).
  Reference association energies (kcal/mol) in `reference_s30l`.
- Both engines run on identical xyz (coord->Å). Charge from `.CHRG` -> curcuma `-charge`,
  xtb `.CHRG` file. ΔE = E(AB)-E(A)-E(B) [kcal/mol].
- Harness: `scripts/s30l_gfnff_compare.py` (subset: `python ... 1 5 25`).
- Term decomposition: `scripts/s30l_term_decomp.py` (e.g. `python ... 25 AB 23 AB`).
- Results: `test_cases/s30l_test_set/_run/{results.csv,results.md}` + per-component logs.

```bash
cd release && make -j4
python scripts/s30l_gfnff_compare.py          # all 30
```

**Cache note**: curcuma writes `_run/{i}/{A,B,AB}.topo.json` topology caches. After
changing GFN-FF topology logic, delete them before re-running
(`find test_cases/s30l_test_set/_run -name '*.topo.json' -delete`); a stale cache can
mask a fix (the cache guard catches total-charge mismatches but not placement).

## What was tested / not tested

- Tested: 30 S30L triples (90 SP x 2 engines), gas phase, default GFN-FF, D4, no solvent.
- Not tested: solvent, optimisation, MD, gradient correctness vs xtb (energies only),
  systems outside S30L, F3 bond-term root cause (Hückel π-orders), system 23 complex term.
- Pre-existing (from the origin/master merge, NOT these fixes): `cli_curcumaopt_07_opt_multixyz`
  fails (gfnff multi-frame energy tolerance ~0.002 Eh > 1e-5; stale golden values from
  origin's GFN-FF changes). Confirmed failing without these fixes too.

Human production testing pending until the operator removes this note.