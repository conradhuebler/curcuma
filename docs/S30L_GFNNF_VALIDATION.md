# S30L: Native GFN-FF vs xtb 6.6.1 GFN-FF — Validation + Fixes

Status: AI-generated, machine-tested (Jul 2026). NOT human production-tested.
The AI does not assign ✅ TESTED / ✅ APPROVED.

## Outcome (after F1/F2/CLI + F3 bond/Hückel + metal_type + ipis fixes)

| group | n | MAD | max | RMSD | verdict |
|------|--:|----:|----:|-----:|---------|
| clean neutral (1,2,3-10,13,14,17-22) | 18 | 0.30 | 1.62 | 0.54 | reproduces xtb |
| F3 fixed (11,12,15,16) | 4 | 0.43 | 0.94 | 0.50 | reproduces xtb |
| charged host, fixed (24,25,26,29) | 4 | 0.18 | 0.58 | 0.25 | reproduces xtb (25/26 = 0.0) |
| charged guest residual (23,27,28,30) | 4 | 8.5 | 27.2 | 11.0 | .CHRG quirk + torsion |
| **all 30** | 30 | **1.38** | 27.2 | 5.0 | was MAD 433 pre-fix, 1.44 before ipis |

- 28/30 systems within ~1.5 kcal/mol of xtb (23 and 30 outside; 7/8 borderline torsion).
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

**F3d — Hückel electron count missing the π-system charge (ipis)** (`gfnff_method.cpp`,
`huckel_solver.cpp`): xtb subtracts the π-system charge from the Hückel electron count
(`gfnff_ini.f90:874`: `nelpi = nelpi - ipis(pis)`), where `ipis(pis)` is the charge
localised on π-system `pis` (computed at `gfnff_ini.f90:517-562` via qheavy +
neutralize-fragment + re-EEQ + dqa·1.1). curcuma counted only the per-atom π-electrons
(26 for S30L 25/A π-system 1) and never subtracted ipis, so for charged hosts the
Hückel had 2 electrons too many → wrong occupations → pibo too high (0.785 vs xtb 0.703)
→ bond term +157 kcal/mol on 25/A. Fixed by computing ipis in `calculateTopologyInfo`
(qheavy + neutralize each π-system's EEQ fragment + re-run topology EEQ + dqa·1.1,
rounded) and subtracting it from nelpi in `HuckelSolver::calculatePiBondOrders`.
Verified via a python Hückel re-solve: nel=24 gives pibo(18,20)=0.70305 = xtb exactly
(was 0.7852 with nel=26). Neutral molecules have ipis=0 (unchanged) — which is why the
18 clean neutrals matched all along. Fixed 25/26 (−1.4 → 0.0004/−0.0000 kcal/mol) and
30/B (now exact). MAD 1.44 -> 1.38. No neutral regression (ipis=0 for neutral).

### Remaining residuals (charged guest + torsion)
- **23 (−27)**: xtb `.CHRG` quirk — curcuma correctly places +1 on the guest, xtb puts
  it on the host (one-line `.CHRG` fallback); with the charge on the guest curcuma == xtb
  exactly (verified). Not a curcuma bug.
- **27/28 (1.45/1.39)** and **7/8 (1.27/1.62)**: **torsion term** — curcuma computes
  ~80% of xtb's torsion energy (e.g. 27/A cur 0.076 vs xtb 0.096 Eh, −12 kcal on the
  host, partly cancels in ΔE). Separate torsion-generation/parameter issue (open).
- **30/AB (+3.98)**: a bond pibo difference with ipis=0 (the −1 does not localise on
  the host π-systems); previously masked by 30/B's error (now fixed). Open.

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

### System 23 (RESOLVED as harness/xtb `.CHRG` quirk — curcuma is correct)
23/AB (q=+1): A and B match xtb **exactly**. The complex ΔE was −27 kcal/mol off the
harness xtb. Root cause is NOT curcuma:

- S30L 23 has A(host)=q0, B(guest)=q+1, so the +1 belongs on the **guest**.
- curcuma's F2 both-assignment trial places it on the guest (the lower-EEQ-energy
  fragment), which is physically correct.
- xtb's `.CHRG` reader (`gfnff_setup.f90:208-222`) reads line 1 = total charge, line 2 =
  per-fragment charges. The S30L `.CHRG` has only ONE line, so line 2 falls back to line 1
  → `qfrag=[charge, 0]` → the WHOLE charge on **fragment 1 (the host)**. The F2 trial in
  `gfnff_ini.f90:478` is dead code (it requires `qfrag(2:nfrag) > 999`, but `qfrag` is
  always pre-initialised to 0 by `gfnff_setup.f90:122/198/224`), so xtb NEVER auto-detects
  the placement — it always charges fragment 1.
- **Verified**: xtb run with a two-line `.CHRG` (`+1` / `0 1` → charge on guest) gives
  23/AB total −18.93093998 Eh, Coulomb +1.10227010 Eh — **identical to curcuma**
  (−18.93093729, +1.10227011) to ~1e-5 Eh. So with the charge on the correct fragment,
  curcuma == xtb; the −27 kcal "gap" is entirely xtb mis-placing the charge on the host.

This is a harness/reference artefact, not a curcuma bug, so it is left as-is (curcuma is
the more correct of the two). A harness that passed per-fragment `.CHRG` to xtb would
close it, but the S30L fragment order in the AB xyz is not always [A, B], so a blanket
`qA qB` line 2 breaks 24/27/28/29/30 (verified: MAD 1.44→12.3). xtb itself is also
~30 kcal/mol off the S30L reference for 23. The 25/26/27/28/30 residuals (±1.4–3 kcal)
are the same class of small charged-fragment effect.

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