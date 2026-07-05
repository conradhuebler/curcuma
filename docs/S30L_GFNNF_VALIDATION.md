# S30L: Native GFN-FF vs xtb 6.6.1 GFN-FF — Validation + Fixes

Status: AI-generated, machine-tested (Jul 2026). NOT human production-tested.
The AI does not assign ✅ TESTED / ✅ APPROVED.

## Outcome (after F1/F2/CLI + F3 bond/Hückel + metal_type fixes)

| group | n | MAD | max | RMSD | verdict |
|------|--:|----:|----:|-----:|---------|
| clean neutral (1,2,3-10,13,14,17-22) | 18 | 0.30 | 1.62 | 0.54 | reproduces xtb |
| F3 fixed (11,12,15,16) | 4 | 0.43 | 0.94 | 0.50 | reproduces xtb |
| charged, fixed (24,25,26,27,28,29,30) | 7 | 1.5 | 3.0 | 1.6 | reproduces xtb |
| charged residual (23) | 1 | — | 27.2 | — | complex-term residual |
| **all 30** | 30 | **1.44** | 27.2 | 5.0 | was MAD 433 pre-fix, 2.70 before F3 |

- 28/30 systems within ~1.5 kcal/mol of xtb (only 23 and the borderline 7/8 outside).
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

### F3 — host bond-term + Hückel + Coulomb offsets (FIXED, Jul 2026)
Systems 11, 12, 15, 16: host energy too negative vs xtb, residual 4-13 kcal/mol.
Term-by-term (`scripts/s30l_term_decomp.py`, `scripts/s30l_bond_compare.py`) showed:
- **Dispersion, Repulsion, BATM: EXACT match** (to ~1e-10 Eh) -> the bond graph
  (perception) and non-bonded exclusions are CORRECT, identical to xtb.
- **Bond stretching term: the discrepancy** — driven by the Hückel π-bond-order (pibo)
  feeding fpi/pi_shift/bstrength. Three root causes, all fixed:

**F3a — S hybridization** (`gfnff_method.cpp` `determineHybridization`): S(Z=16) with CN=2
took the geometry-based branch (C-S-C bent -> sp2), but xtb `gfnff_ini2.f90:285-298` sets
group-16 (O,S) CN=2 -> sp3 unconditionally. sp2 S gave bsmat[2][2]=1.24 (too stiff) and
nel=1 in the Hückel (sp2) instead of nel=2 (sp3) -> over-delocalised pibo on C-S -> bond
energy too negative by ~0.7 Eh on the S30L 11/12 hosts (cucurbituril-type). Fixed 11/12
(-11.8/-12.6 -> -0.33/-0.31 kcal/mol).

**F3b — halogen hybridization + π-system membership** (`gfnff_method.cpp`,
`huckel_solver.cpp`): (1) Br(35)/I(53) with CN=1 fell through to hyb=1 (sp) because the
"halogen" range was z 9-17; xtb sets group-7 halogens CN=1 -> hyb=0. sp atoms are
π-candidates, so I wrongly entered the π-system. (2) `hoffdiag` default for Z>=18 was 1.0
(the C reference) instead of xtb's 0 (only B,C,N,O,F,S,Cl are parameterised) -> heavy
halogens wrongly coupled in the Hückel. (3) π-candidate detection included ANY sp/sp2 atom
and missed F/Cl (hyb=0, hoffdiag>0). Rewrote: π-candidate iff element has hoffdiag>0
(B,C,N,O,F,S,Cl); Case-2 π-bond only to a heteroatom (not sp3 C); halogens CN=1->hyb=0,
CN=2->hyb=1. Verified by `scripts/s30l_hyb_audit.py`: 0 π-membership mismatches vs xtb
across all 30 (was 60 F + 12 I + 190 C mismatches).

**F3c — `metal_type` array misalignment** (`gfnff_par.h`): the array had only 83 entries
— the K-Kr row missed Kr(36)=0 and the Rb-Xe row missed I(53)=0 and Xe(54)=0 (plus
Ag(47) was 1 not 2). The 3-element shortfall shifted every Z>=36 index by 1-3, so I(Z=53)
read index 52 (=2, transition metal) instead of 0. `calculateDgam` then used the TM
ff=-0.9 for I instead of the halogen ff=-0.07 -> dgam(I)=-0.176 (should be -0.0137) ->
gameeq(I)=-0.136 (negative hardness) -> iodine over-polarised by ~0.19 e -> S30L 15/16
hosts (C24F15I3) Coulomb -14 kcal/mol. Array rewritten to match xtb `metal(86)`
(`gfnff_param.f90:318-325`) exactly. Fixed 15/16 (-9.6/-4.9 -> +0.13/+0.94 kcal/mol);
15/A Coulomb now -0.19647438 Eh = xtb -0.19647438 to 1e-9.

No regression on the 18 clean neutrals (all unchanged or improved). MAD 2.70 -> 1.44.

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

**The 4 neutral outliers (11/12/15/16) are FIXED** (F3a/b/c above): bond term now matches
xtb (pibo to ~1e-4, bond E to ~1e-6 Eh), Coulomb matches for the I-host. The remaining
outliers are charged (23, 25-28, 30) and the borderline 7/8 (torsion, cancels in ΔE).

## Setup / reproduce

- 30 S30L host-guest complexes, A=host, B=guest, AB=complex, Turbomole coord (Bohr).
  Reference association energies (kcal/mol) in `reference_s30l`.
- Both engines run on identical xyz (coord->Å). Charge from `.CHRG` -> curcuma `-charge`,
  xtb `.CHRG` file. ΔE = E(AB)-E(A)-E(B) [kcal/mol].
- Harness: `scripts/s30l_gfnff_compare.py` (subset: `python ... 1 5 25`).
- Term decomposition: `scripts/s30l_term_decomp.py` (e.g. `python ... 25 AB 23 AB`).
- Per-bond params: `scripts/s30l_bond_compare.py` (curcuma vs xtb pibo/fc/r0/alpha).
- Hybridization audit: `scripts/s30l_hyb_audit.py` (per-atom hyb + pi membership vs xtb).
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
  systems outside S30L, F3 bond-term root cause (Hückel π-orders) on non-S30L molecules,
  system 23 complex term.
- Pre-existing (from the origin/master merge, NOT these fixes): `cli_curcumaopt_07_opt_multixyz`
  fails (gfnff multi-frame energy tolerance ~0.002 Eh > 1e-5; stale golden values from
  origin's GFN-FF changes). Confirmed failing without these fixes too.

Human production testing pending until the operator removes this note.