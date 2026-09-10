# rev-gfnff — where curcuma could be better than the reference

⚠️ AI-generated, machine-tested only — not human production tested.

A running list of places where curcuma's GFN-FF **deliberately does, or could do, something
different from the reference** because the reference is demonstrably wrong on the physics.
Everything here is separate from port fidelity: the default is to reproduce
`pprcht/gfnff` (and, where they agree, `xtb`), because that is what the GFN-FF parameters
were fitted against. This file collects the cases where that costs accuracy, so a future
revised force field can pick them up deliberately rather than rediscover them.

**The distinction that decides each case**: were the GFN-FF parameters fitted *against* the
quirk? If yes, "fixing" it invalidates the parametrisation and must stay opt-in. If the
quantity is a hard-coded ab-initio constant or a plain coding slip, curcuma may lead.

**Scope, going forward**: this file is also where the parametric adaptation towards
**react-gfnff** ([GFNFF_REACT_TOPOLOGY.md](GFNFF_REACT_TOPOLOGY.md)) is collected. Port
fidelity gets curcuma to reproduce GFN-FF; the entries here are what makes it describe the
chemistry better than GFN-FF does — including refitting parameters where the functional
form itself is the limit. Bond-breaking and transition states (entries 4, 6, 9, 10) are the
obvious targets, and the `reactff2` branch feeds in through entries 11-15, which say what a
refit has to supply before the reactive mode can drop its empirical filters.

## Already deviating from the reference by default

### 1. `nh_linear_fix` — GEODEP angle rule creates artefact minima at N-H centres
Started on branch `confsearch` (`51830efa`), ported and then refined here (see
[GFNFF_STATUS.md](GFNFF_STATUS.md)). The reference's "input angle > 160 deg -> sp" fallback
makes a thermally stretched =N-H its own equilibrium, so an optimisation walks into the
artefact. **External verdict** on the same two geometries (formamidine, only the imine
C-N-H angle varied), ΔE(179° − 119°):

| method | kcal/mol |
|---|---:|
| r²SCAN-3c | +24.1 |
| GFN2 | +24.0 |
| GFN-FF, guard on (curcuma default) | +14.0 |
| GFN-FF, guard off (= xtb/pprcht) | **−66.4** |

The reference inverts the sign; curcuma keeps it. `-gfnff.nh_linear_fix false` restores
bit-faithfulness. Status: **done, default on.**

### 2. `storsion_reference_loop_bug` — sTors evaluates only its last entry
`gfnff_engrad.F90:494-503` loops `do i=1,m` but calls `sTors_eg(m, ...)` with the array
SIZE, so pprcht *and* xtb evaluate one entry m times and drop the rest. Its `erefhalf` is a
DLPNO-CCSD(T) diphenylacetylene value, **not a fitted parameter**, so summing correctly
cannot invalidate anything. curcuma sums; `-gfnff.storsion_reference_loop_bug true`
reproduces the reference. Status: **done, default correct.**

## Open — physics known to be wrong in BOTH implementations

### 3. EEQ over-rewards charge delocalisation (charged radicals, 2c-3e bonds)
Measured on GMTKN55 `G21EA/EA_25`, the dichlorine radical anion Cl₂⁻ at 2.73 Å.
E(Cl₂⁻) − E(Cl⁻) − E(Cl):

| method | kcal/mol |
|---|---:|
| r²SCAN-3c (ORCA 6, UHF) | **−41.5** |
| GFN2 (xtb, UHF) | −33.9 |
| GFN-FF, nfrag = 2 (reference behaviour, charge pinned on one Cl) | −6.6 |
| GFN-FF, nfrag = 1 (charge delocalised over both) | −106.3 |

Experiment puts the dissociation near −30 kcal/mol. **Neither GFN-FF setting is right**, and
the fragment count is the only lever — a discrete switch that brackets the truth without
hitting it. Term decomposition (curcuma):

| | bond | Coulomb | total |
|---|---:|---:|---:|
| nfrag = 2 | −11.5 | +5.3 | −6.6 |
| nfrag = 1 | −11.2 | **−94.8** | −106.3 |

The bond term is essentially unchanged; the entire difference is the EEQ. Going from
q = (−1, 0) to (−0.5, −0.5) is worth **−94.8 kcal/mol** to GFN-FF's EEQ self-energy
(which scales as q², so splitting one unit charge over two centres halves the penalty),
where the true stabilisation is ~−35. A revised model needs a delocalisation term that is
not just the EEQ self-energy — a **method change**, not a parameter or perception fix.

Curcuma currently reproduces the reference here (Known Issue #17: the second q-loop pass
keeps pass 1's fragmentation, `gfnff_ini.f90:467`), because that is both faithful and the
smaller error. Affects the whole charged-complex family: AHB21, BH76 SN2 transition states,
G21EA, SIE4x4, CHB6.

### 3b. Carbene angle: GFN-FF bends CH2- 30 deg wrong either way
Bending scan of the methylene anion CH2- (C-H fixed at 1.129 A, H-C-H varied), relative
energies in kcal/mol:

| H-C-H | GFN-FF carbene (θ0=145) | GFN-FF reference (θ0=120) | GFN2 | r²SCAN-3c |
|---:|---:|---:|---:|---:|
| 90 | 13.77 | 10.42 | 2.56 | 1.58 |
| **100** | 8.87 | 5.72 | **0.00** | **0.00** |
| 110 | 5.27 | 2.62 | 1.66 | 1.14 |
| 120 | 2.77 | 0.81 | 6.65 | 4.57 |
| 130 | 1.17 | **0.00** | 14.27 | 9.72 |
| 145 | **0.00** | 0.02 | 28.81 | 19.02 |

GFN2 and r²SCAN-3c agree on a minimum at **100°** (experiment ~102°). The reference-faithful
treatment minimises at 130°, the carbene one at 145° — so restoring port fidelity (Known
Issue #18) also moved GFN-FF *towards* the truth, unlike case 3. But 130° vs 100° is still a
30° error: GFN-FF has no equilibrium angle that fits a bent carbanion. Unlike case 3 there is
no conflict here, only a ceiling.

### 4. GFN-FF cannot do MOR41 reaction thermochemistry at all
Against DLPNO-CCSD(T) (Table S1, 41 reactions): pprcht MAD 62.6, curcuma 63.1, xtb 71.7
kcal/mol, where GFN2 reaches ~12. A method limitation, not a port issue — recorded so nobody
re-opens it as a bug. See [MOR41_VALIDATION.md](MOR41_VALIDATION.md).

### 5. Bond perception in charged species is chemically right but energetically worse
Same root as (3), stated separately because it may deserve its own treatment: curcuma's
charge-shrunk radii correctly perceive Cl₂⁻ as **bonded** (it is — a real radical anion at
~2.7 Å), and the reference's own second pass agrees (`#bonds : 1`). But GFN-FF's bond
parameters are fitted for ordinary 2-electron bonds and overbind a 2-centre-3-electron bond.
Getting the perception right therefore does not help until the bond term can express a
half bond.

### 6. The angle term is ~3x too soft at an SN2 transition state
Umbrella scan of the F...CH3...F- transition state (C-F and C-H fixed, the three F-C-H
angles varied together), energy relative to the D3h point:

| F-C-H | GFN-FF | GFN2 | r²SCAN-3c |
|---:|---:|---:|---:|
| 80° | **2.13** | 6.03 | **6.43** |
| 85° | 0.48 | 1.49 | 1.61 |
| 90° | 0.00 | 0.00 | 0.00 |

GFN2 reproduces r²SCAN-3c almost exactly; GFN-FF is a factor of three too soft. Known Issue
#19 restored port fidelity by removing a non-reference `|qa| < 1` guard around the angle
`fqq`, which **softens** the F-C-H constants further (0.2631 -> 0.224) — i.e. away from the
physics, like case 3. It was still the right call: the accidental stiffening only ever
applied to angles at a full-unit-charge fragment, so it was not a principled correction.
A revised model would need a stiffer bend at a hypercoordinate centre generally.

*Caveat on the numbers*: the scan geometry is symmetric, so its perception differs slightly
from the GMTKN55 structure (there one fluoride is its own fragment). The GFN-FF column is
therefore the reference-faithful `fqq` regime throughout — the comparison measures the term,
not the guard.

### 7. MB16-43 is out of reach for GFN-FF, so the pprcht-vs-xtb split there is moot
`MB16-43` is the last large GMTKN55 residual (curcuma-vs-xtb MAD 10.4, max 105) and the one
cluster where **pprcht and xtb disagree with each other** by 30-170 kcal/mol per structure.
Arbitrating that split needs a yardstick, so the set's own published decomposition reactions
(`MB16-43/.res`, e.g. `-2·E(35) -20·E(H2) +2·E(LiH) + ... = 228.1748 kcal/mol`) were computed
with all three engines:

| engine | MAD vs published reference | max |
|---|---:|---:|
| pprcht | **383.4** kcal/mol | 1135.6 |
| curcuma | 384.0 | 1152.7 |
| xtb | 398.3 | 1135.6 |

(42 of 43 reactions; one engine failed on reaction 03.)

**All three are wrong by two orders of magnitude more than they differ from each other.** The
split is irrelevant for accuracy here; only port fidelity is at stake, and the reference is
itself ambiguous, so there is nothing to chase.

Pipeline validated against r²SCAN-3c on two reactions, including the worst case:

| reaction | published | r²SCAN-3c | GFN-FF |
|---|---:|---:|---:|
| 35 | 228.2 | **228.1** | 300.1 (ppr) / 314.7 (cur) / 181.8 (xtb) |
| 22 | 706.2 | **704.2** | −150.0 (all three) |

Two contributing reasons, both structural: MB16-43 molecules are "mindless" random main-group
species whose bonding GFN-FF was never fitted for, and several are **open-shell radicals**
(structure 22 has 75 electrons, a doublet) which GFN-FF has no concept of at all. Same class
as case 4 (MOR41 reaction thermochemistry), one order of magnitude worse.

### 8. AL2X6 dimerisation comes out with the wrong SIGN in every GFN-FF
`AL2X6` is the second cluster where pprcht and xtb disagree (`al2me6`, 26.3 kcal/mol in the
repulsion term; curcuma matches pprcht there exactly). The set's own published reference is
the dimerisation energy `2 AlX3 -> Al2X6`, so the split can be arbitrated directly:

| reaction | published | r²SCAN-3c | pprcht | xtb | curcuma |
|---|---:|---:|---:|---:|---:|
| al2h6 <- 2 alh3 | 38.5 | **40.6** | −32.7 | −32.7 | −32.7 |
| al2f6 <- 2 alf3 | 51.6 | **52.0** | −2.4 | −2.4 | −5.2 |
| al2cl6 <- 2 alcl3 | 32.5 | **32.4** | −23.3 | −23.3 | −23.3 |
| al2me6 <- 2 alme3 | 23.1 | **24.7** | −8.9 | **−35.2** | −8.9 |

(positive = dimer bound; r²SCAN-3c reproduces the published values to 1-2 kcal/mol, which
validates the reaction pipeline.)

**Every GFN-FF variant gets the sign wrong**: it predicts the dimers unbound by 3-35 kcal/mol
where they are bound by 23-52, an error of 40-70 kcal/mol. The 3-centre-2-electron Al-X-Al
bridge is simply not in the model. The pprcht-vs-xtb split is therefore moot here too, though
pprcht (and with it curcuma) happens to be the closer of the two on `al2me6`: 32 kcal off
against xtb's 58.

Same conclusion as case 7 (MB16-43) and case 4 (MOR41): where the two references disagree,
they are both far enough from reality that the disagreement carries no information.

### 9. The pprcht-vs-xtb split sits almost entirely outside what GFN-FF can do

The whole GMTKN55 sweep was re-evaluated as REACTIONS rather than single points: each
subset's `.res` file is a `tmer2++` script that carries both the stoichiometry and the
published high-level reference value, so every reaction can be scored from the cached
single points. Pipeline checks: it reproduces the AL2X6 dimerisation numbers of entry 8
exactly, and r2SCAN-3c reproduces the published value on the three reactions run below
(0.3 / 3.2 / 16 kcal/mol).

**Where GFN-FF is usable, the two references are indistinguishable:**

| reaction class | n | MAD curcuma | MAD xtb | max split | split > 1 kcal |
|---|---:|---:|---:|---:|---:|
| pure conformers | 285 | 1.49 | 1.47 | 5.13 | 1 |
| non-covalent complexes | 239 | 9.41 | 9.41 | 0.01 | 0 |
| isomerisations | 102 | 28.10 | 28.12 | 8.38 | 2 |
| charged non-covalent | 63 | 33.81 | 33.78 | 1.03 | 4 |
| bond breaking / atomisation / open shell | 397 | 118.96 | 120.58 | ~170 | 39 |

GFN-FF is genuinely accurate only for conformers (1.5 kcal/mol), and there the two
implementations differ by more than 1 kcal/mol on exactly ONE of 285 reactions. For
non-covalent complexes they agree to 0.01 kcal/mol. The split lives in the last row.

**Where the split IS large, pprcht is usually the better one** — over all 541 evaluable
reactions, restricted to those where the two implementations disagree:

| split threshold | n | curcuma closer | xtb closer | mean split | mean own error |
|---|---:|---:|---:|---:|---:|
| > 1 kcal/mol | 42 | 21 | 21 | 24.0 | 228.5 |
| > 5 kcal/mol | 22 | 15 | 7 | 43.7 | 290.4 |
| > 20 kcal/mol | 11 | 10 | 1 | 75.4 | 404.2 |

**The decisive single case** is the largest split in the set, MB16-43 reaction 41
(`-2*[41] -16 H2 +4 BH3 +6 CH4 +F2 +2 NaH +2 AlH3 +2 S2`; structure 41 is a doublet and
S2 a triplet, so GFN-FF is doubly out of its depth):

| method | kcal/mol |
|---|---:|
| published (GMTKN55) | +160.29 |
| r2SCAN-3c | **+176.33** |
| curcuma = pprcht | +146.04 |
| xtb | **-23.04** |

pprcht is in the right range with the right sign; xtb has the wrong sign and is ~190
kcal/mol out. Two smaller checks, for contrast:

| reaction | published | r2SCAN-3c | curcuma | xtb |
|---|---:|---:|---:|---:|
| ISOL24 i11 (isomerisation) | +36.90 | +33.71 | +89.27 | +97.65 |
| ICONF N3P3H12 (conformer) | +12.16 | +11.88 | -8.00 | -2.87 |

The ICONF case is the one conformer where the split matters (5.1 kcal/mol) — and there
BOTH implementations get the sign wrong, xtb merely less badly. So the pattern is not
uniform; it is a tendency that only becomes clear at the extreme end.

**Verdict**: tracking pprcht rather than xtb is supported by the external references, most
clearly at the biggest divergences. But this is not where GFN-FF's error lives. Chasing
the remaining split further would buy nothing for any application the method is suited to,
and the numbers above are the reason MB16-43 and AL2X6 are treated as closed
(entries 7 and 8).

### 10. The bond charge factor `fqq` overflows to NaN for strongly ionic bonds

**curcuma deviates by default.** The reference writes the logistic naively
(`gfnff_ini.f90`, non-metal bond branch):

    qafac = qa(i)*qa(j)*70
    fqq   = 1 + qfacbm0 * exp(-15*qafac) / (1 + exp(-15*qafac))

For an ionic bond that overflows. `|qa_i * qa_j| > 0.675` already puts the exponent past
709, `exp` returns +Inf, and `Inf/(1+Inf)` is NaN — so the force constant, the bond energy
and the entire single point become NaN. GMTKN55 `PX13/hf_4_ts` and `hf_6_ts` (cyclic (HF)n
proton-transfer transition states with qa = +-1.007 on every atom) return **NaN from
pprcht AND from xtb**; the value the GMTKN55 harness records from xtb for them is the
angle term alone, i.e. junk.

curcuma clamps the logistic to its analytic limits outside `|t| <= 300`. The branch is
**bit-identical** wherever the reference produces a finite number — verified over all 2460
GMTKN55 structures that pprcht can evaluate: MAD 0.01230, max 26.9978, and the same
per-threshold counts before and after the change, with only the two NaN structures moving
from "curcuma has no energy" to "pprcht has no energy". No parametrisation was ever done
against a NaN, so this is a numerical accident and not a fitted quirk, which is why it is
the default rather than opt-in.

**It buys correctness, not accuracy.** The two structures now yield finite energies, and
the PX13 barriers they complete show what those energies are worth:

| (HF)n proton transfer | published | curcuma | xtb |
|---|---:|---:|---:|
| hf_2 | 42.3 | 66.8 | 66.8 |
| hf_3 | 20.7 | 137.6 | 137.6 |
| hf_4 | 14.7 | **223.1** | 404.3 (junk) |
| hf_5 | 14.6 | 220.0 | 220.0 |
| hf_6 | 16.6 | **293.2** | 609.7 (junk) |

**curcuma covers a case both reference implementations drop.** That is the point of the
guard: a strongly ionic bond no longer takes the whole single point down, so an optimiser,
an MD or a conformer search walks through such a geometry instead of aborting. The value it
returns is not accurate — the PX13 barriers are 3-18x too high, because GFN-FF cannot
describe a transition state where the bonding topology changes. Making those numbers right
is a parametrisation question and belongs to the react-gfnff work, not here.

**When does it actually fire?** The trigger is a BOND whose two topology charges multiply
to `qa_i * qa_j < -0.675`. `CURCUMA_BONDDUMP=1` now prints `qaprod=` per bond, so the
headroom of any system is measurable. Scanned over the reference sets:

| set | n | largest negative qa product | over 0.675 |
|---|---:|---:|---:|
| S30L-CI (host-guest, multi-fragment, explicit counterions) | 86 | 0.140 | 0 |
| ALKBDE10 (bonded ion pairs: LiF, NaCl, KF, CaO, ...) | 10 | 0.157 | 0 |
| IL16 / AHB21 / CHB6 / DIPCS10 / SIE4x4 / G21EA / PX13 | 176 | 1.014 | 2 |

So it is NOT "ionic system" or "multi-fragment" that does it — GFN-FF's topology charges
stay modest there (LiF reaches only qa = +-0.32). It takes a strongly polar bond in a
SYMMETRIC BRIDGING arrangement: the (HF)n proton-transfer rings put qa = +-1.007 on every
atom because each hydrogen sits midway between two fluorines. `hf_5_ts` (0.578) and
`AHB21/3` (0.544) are the closest non-overflowing cases, at about 80 % of the threshold.

That geometry — a hydrogen or an ion halfway between two acceptors — is exactly what a
**reactive** run walks through, which is why this guard matters more for react-gfnff than
for any static benchmark.

**A harness bug found on the way.** `xtb` prints `TOTAL ENERGY NaN Eh` for these two
structures but still prints finite values for the individual terms above it, and
`scripts/gmtkn55_compare.py`'s last-resort parser scanned every line containing "energy"
and took the first number followed by "Eh" — the ANGLE energy. So the comparison recorded
`PX13/hf_4_ts` from xtb as -0.000479853103 Eh, a fabricated number, and the two "xtb"
columns of the table above (404.3 and 609.7) came from it. The parser now detects NaN and
reports a failure instead. Anything that reads an external code's output needs to fail
loudly when that code fails; a fallback that keeps looking until it finds *a* number will
eventually find the wrong one.

**What the whole port campaign bought, honestly.** The reaction-level numbers of entry 9
are the answer: over 285 conformer reactions curcuma sits at MAD 1.49 kcal/mol and xtb at
1.47; over 239 non-covalent complexes both are at 9.41. curcuma, pprcht and xtb are
indistinguishable in accuracy. The port fixes moved curcuma from *a wrong implementation*
of GFN-FF (up to 500 kcal/mol on a single structure) to *a correct* one, which is worth
doing because nothing downstream — MD, conformer search, optimisation — can be reasoned
about otherwise. It did not, and could not, make GFN-FF more accurate. Any real accuracy
gain has to come from the deliberate deviations collected in this file.

## What react-gfnff needs from a refit

The reactive topology mode ([GFNFF_REACT_TOPOLOGY.md](GFNFF_REACT_TOPOLOGY.md)) works, and
`reactff2` is where it is developed. It reaches a running bond-forming, bond-breaking MD
through **empirical filters standing in for missing physics**, each of them switchable and
each of them documented there as such. Those filters are the shopping list for a refit: an
entry below is not a bug to fix in the port, it is a place where the functional form or its
parameters have to change before the stopgap can go.

These entries are **not** externally arbitrated the way entries 1-10 are. They are known
gaps rather than measured deviations, and they are marked as such so nobody quotes a number
from here that has not been measured.

### 11. A topology change is a step in the potential energy

Rebuilding the bond list regenerates every bonded term, the repulsion partition and the EEQ
constraints at once, so the energy jumps: measured H + H → H₂, about −450 kJ/mol arriving in
a single step; a break at the default factor returns +21 to +34 kJ/mol. `dE_jump` records it
per event. The consequence is that the mode is **NVT-only** — a thermostat has to absorb the
jumps, and NVE drifts at every event.

Everything else on this list follows from this one. A refit target that removes it: bonded
terms whose contribution goes smoothly to zero over the formation/break window, so the two
topologies agree in energy where they are exchanged, and the 3-/4-body damping that already
switches semi-smoothly is joined by the repulsion re-partition. That is a change of
functional form, and the parameters were fitted for a fixed topology, so it cannot be done
by re-tuning alone.

### 12. The formation and break radii are geometric factors, not energetics

Formation is optimistic (1.6 × the covalent sum) because the non-bonded repulsion wall
limits the capture radius; retention is conservative (2.6) because the Gaussian well of the
bond term decays slowly, so that breaking removes only ~20-35 kJ/mol of residual well
instead of ~480. Both numbers are chosen for the behaviour they produce, not derived.

The refit target is the per-pair, energy-based criterion already named as an open refinement
there: remove a bond where **its own** well has decayed past a threshold. That needs the
bond term's depth and width to be trustworthy far from equilibrium, which is exactly where
GFN-FF was never fitted — see entry 6 for the same problem at a transition state.

### 13. The valence cap is a rule where an energy belongs

A new bond forms only while both partners' used valence (Σ bond orders, σ = 1 plus the
Hückel π order) stays within the element valence plus one exchange slot. Without it a hot,
confined system over-bonds into a cluster whose energy eventually turns NaN. The caps
(H/F/halogens 1+1, O 2+1, N 3+1, B 3+1, C 4+1, hypervalence-capable and metals 6) are
chemical common sense, not a fit.

A refit replaces the rule with an **over-coordination energy** rising beyond the element
valence (ReaxFF-style). Its shape and offset are what the QM reference has to supply:
deliberately hyper-coordinated species, computed at a level that resolves the penalty. That
is a genuine parametrisation task and the clearest candidate for "QM-learned" in this file.

### 14. The refractory period hides an energy pump

A pair whose bond broke may not re-form for 10 scans (~25 fs). Without it the form/break
cycle pumps the recombination energy through the thermostat repeatedly: measured for
N₂ + 3 H₂ at 3500 K in a 3.5 Å wall over 20 ps, 301 events without cap and refractory
against 35 with, and a NaN in the first case.

This is a timer compensating for entry 11. With a continuous energy at the exchange point
there is nothing to pump and the refractory period has no work left to do. It is listed
separately because it is the cheapest test of whether a refit succeeded: if the refractory
period can be set to zero without the event count exploding, the discontinuity is gone.

### 15. There is no reference data for the reactions this is aimed at

GFN-FF is not parametrised for transition states, and entry 6 measures the angle term at ~⅓
of the correct stiffness at an SN2 saddle. A reactive force field is judged on barriers and
paths, so a refit needs reference paths, not just reference minima: NEB or comparable at
GFN2-xTB or higher for the target reactions. The current target is ammonia synthesis
(N₂ + 3 H₂); Miller-Urey-type chemistry is the candidate after it.

Until such a set exists, the react mode demonstrates machinery rather than energetics, and
the documentation says so in those words. Producing that set is the prerequisite for every
entry above, because none of them can be shown to have improved anything without it.

## How to add to this list

An entry belongs here when an **external** reference (r²SCAN-3c, GFN2, DLPNO, experiment)
has been consulted and shows the reference implementation to be wrong — not merely different
from curcuma. Record the numbers, not the impression.
