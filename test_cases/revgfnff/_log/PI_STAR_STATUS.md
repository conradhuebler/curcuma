# PI_STAR_STATUS - a pi*-population-aware excess-electron perception for O2- / S2-

Sep 25, 2026. Opus agent, isolated git worktree (`.claude/worktrees/agent-a997450b...`, branch
reset to `reactff2-llm` 63a3e4de), own build dir `build_pi/`. Scope selected by the operator:
"O2-/S2- als neue Erkennungslogik". AI-generated, machine-tested only; human production testing
pending.

**Read section 11 first.** It supersedes the rows and several conclusions of sections 6-10: the
O-O/S-S well and harris rows were refitted, the S-S order-3 row was removed from the default table,
and the section-10.4/10.6 findings were re-measured with the correct protocols.

## 0. Starting point (from X2_SCOPE_STATUS.md section 1)

`GFNFF::revExcessElectrons()` perceives an excess electron as `x = max(0, -Q_f - free sigma
slots)`. For O2- it gives x = 0 (continuous order 3 from the sp-sp pi doubling, 0.667 free slots
per O), for S2- x = 0 (revValence(S) = 6). The extra electron of superoxide / the S2- radical anion
sits in pi*, of a pi system the sigma budget already counts as satisfied. In addition the P3 well
machinery is hardwired to the order interval [0.5, 1] (half rows, `order < 1`), and O2's key is 3.

## 1. Design (written before any code)

### 1.1 What the signal is

For a diatomic, "one more electron than the neutral bonding count, in pi*" is visible without any
orbital information: the pair is an isolated two-atom component of the bond graph, its bond has a
pi component (continuous order > 1), and the component carries exactly one negative charge.
GFN-FF's own FT-HMO agrees on the magnitude: the ipis correction (`nelpi -= ipis`) adds the
electron to the pi count (O2: 2 pi electrons in 2 orbitals, pibo 1; O2-: 3 electrons, pibo 0.5 =
true pi order 1/2 of the O2- 3-electron pi bond). Whether plain GFN-FF keeps or discards it (the
`pisip > 0.40` "wrong pi occupation" fallback re-solves with nelpi - 1) is measured in section 2
- it decides whether the continuous order already sees the excess, and the X2_SCOPE probe says
it does not (order 3.0 for O2- as for O2).

### 1.2 Rule (per connected component c of the topology's bond graph)

    fires iff  |c| == 2, the two atoms i, j are bonded (each has exactly one neighbour),
               continuousBondOrder(i, j) > 1 + 1e-6            (a pi component exists),
               |Q_c + 1| < 1e-3                                (exactly one excess electron;
                                                                Q_c = sum of Phase-1 charges),
               RevWellTableV2::hasPiExcess(Z_i, Z_j)           (a calibrated pi-excess row),
    then       y_ij = 1   stored in TopologyInfo::rev_pi_excess (separate from rev_excess).

Every input is a topology constant (as on the sigma side), so y is one too: no geometry
derivative. A dianion (O2 2-, |Q_c + 2| small) does NOT fire - not attempted. A pair that the
sigma side can perceive (has a half-order row) is excluded by construction: no pair gets both a
half row and a pi-excess row, and the rule additionally skips a pair already in `rev_excess`.

### 1.3 What consumes y (the well-fitting side is reused, the perception is new)

1. **Bond well (mg3 only)**: the bond carries `Bond::rev_pi_excess = y`; `prepareWellForms()`
   takes the ordinary `findOrder(rev_order)` parameters and blends them towards the pair's
   pi-excess row with weight `min(1, y)`. `rev_order` itself is NOT changed, so nothing else
   keyed on it (the order interpolation of every other O-O bond, the half-order uncap) moves.
   The pi-excess row is capped like a normal well (the extra electron is not a sigma* wall).
2. **Split-charge side**: `revExcessKappa / revExcessFracC / revExcessHarrisX` read `x + y` for
   the pair, so flat / frac / harris all see the pi excess exactly as they see a sigma excess.
   This is what makes the Br2- fit recipe (well row under flat kappa_x 100, then harris g under
   free charges) usable unchanged. A harris row for the pair is added to rev_harris_table.h.
3. Nothing else. (Correction found while implementing: the react-mode q0 relocalisation
   `revLocaliseExcessQ0` reads the corner's `sqe_kappa_x`, which goes through `revExcessKappa`,
   so under flat it DOES see y - the pi pair is relocalised like a sigma pair. Consistent, but
   flat-mode react dynamics with a pi excess are not validated.)

Gate: a new PARAM `rev_pi_excess_electron` (Bool, default false), effective only together with
`rev_excess_electron true` (which in turn needs `rev_charge_model sqe`), else a warning and off.
Default off = bit-identical by construction (the map stays empty, y = 0 everywhere, the stamp
and all accessors see 0).

### 1.4 Explicitly NOT attempted

- general polyatomic pi systems (O3-, NO2-, CO3 radical anions, aromatic radical anions, any
  pi system of more than two atoms or a diatomic bonded to anything else);
- dianions (O2 2-, S2 2-), cations (O2+), and anything with |Q_c + 1| >= 1e-3;
- heteronuclear pi anions (NO-, CN is closed shell anyway) - rule is element-generic, but no row;
- spin states (O2- is a doublet, neutral O2 a triplet; the force field has no spin) - the
  reference is the doublet ground state;
- a Hueckel-level (orbital) perception - it would be the general route (see section 1.1), but it
  needs a decision on the plain-GFN-FF `pisip` fallback, which is port-fidelity territory.

### 1.5 Falsification plan

- no-op: GMTKN55 (2462) / MOR41 (95) / S30L-CI (90) under plain gfnff, revgfnff default and the
  recommended X2- setting with the new flag on vs off; the 1379-frame fit harness in the same
  configs; Cl2- / F2- / Br2- curves bit-identical with the flag on; ctest.
- physics: O2- full-grid rms vs DLPNO-CCSD(T), bonded / compressed / tail split, min E and r_min,
  up-vs-down (react breaking / forming) and a water-probe label gap, FD gradients.

## 2. Measured design basis (baseline binary A = worktree HEAD 63a3e4de unmodified, md5 75be367c)

O2- at 1.35 A, `CURCUMA_HUCKELDUMP=1`, revgfnff sqe + phase1 + excess + harris:

    HUCKELATOM pis=1 atom=1(Z=8) hyb=1 tag=0 nel=1       (x2)
    HuckelSolver: pi-system 1 has 2 atoms, 3 electrons    <- ipis = -1 added the excess electron
    Wrong pi occupation (HOMO eps=1.5452 > 0.40): second attempt with Nel=2 at et=300
    pibo[1] = 1.000000                                   <- and the fallback threw it away

So plain GFN-FF's Hueckel does see the pi* electron and deliberately discards it (the
`pisip > 0.40` fallback, Known Issue #6(b)/(j), a faithful port of the reference): the
continuous order stays 3, identical to neutral O2. Had the 3-electron solve been kept, pibo 0.5
and the sp-sp doubling would give order 2 - the midpoint between the O-O order-1 (H2O2) and
order-3 (O2) rows, i.e. exactly the true order 1.5 of superoxide in the model's order units.
This confirms the physical picture and the magnitude y = 1, but it is NOT used as the signal:
switching the fallback off would be a plain-GFN-FF port change with effects on every odd-electron
pi system (the 7e N-heteroaromatics of MOR41 need it). The topology rule of 1.2 is used instead.

Model bond parameters of O2- along the grid (flat, P0 binary; `CURCUMA_BONDDUMP`): continuous
order 3.000 through 1.80 A, 2.995 at 1.90-2.05 A (pass-2 bond), no bond from 2.25 A (static
cutoff between 2.05 and 2.25 A). fc -0.384 (-0.392 past 1.9 A), alpha 0.5606, dynamic r0
1.807 -> 1.217 Bohr (the O CN falls as the bond stretches).

## 3. Implementation (prototype binary P0 = design + PLACEHOLDER rows)

| file | change |
|---|---|
| `gfnff.h` | PARAM `rev_pi_excess_electron` (Bool, false); `TopologyInfo::rev_pi_excess`; `m_rev_pi_excess`, `revPiExcessElectrons()`, `revExcessTotal()` |
| `gfnff_method.cpp` | parse (hardcoded fallback false; needs `rev_excess_electron`, else warn + off); cache fingerprint `,xpi=1` when on; `revPiExcessElectrons()` (rule 1.2); filled in `revApplyPhase1Sqe` right after the sigma map (so Phase 1 already sees it under flat); `revExcessKappa/FracC/HarrisX` read `x + y` via `revExcessTotal`; `Bond::rev_pi_excess` set at bond build |
| `ff_terms.h` | `Bond::rev_pi_excess` (default 0) |
| `ff_workspace.h/_gfnff.cpp` | well stamp 5 -> 6 entries (+ rev_pi_excess); mg3 path: `applyPiExcess` after `findOrder`, only if `rev_pi_excess > 0` |
| `rev_well_table_v2.h` | hand-maintained `kPiExcessEntries` + `hasPiExcess()` + `applyPiExcess()` |
| `rev_harris_table.h` | O-O row |

Smoke test (P0, O2- 1.35 A, harris): `rev-gfnff P3 pi*: excess pi* electrons y = 1.0000 on bond
1-2 (order 3.000)`; Bond term unchanged (placeholder = the O2 order-3 row, as intended), the
harris term appears as SqeHardness +0.1144 Eh (placeholder = Cl-Cl copy); flag without
`rev_excess_electron` -> warning, energy identical to A.

## 4. Reference campaign: cost estimate (logged before launch) and data

Pilot DLPNO-CCSD(T)/aug-cc-pVTZ (package-21 keywords, all-electron) O2- at 1.35 A: 92 basis
functions, 87 s on 8 cores under load ~43, T1 0.021, <S^2> 0.795 (UHF) / 0.7504 (linearized).
Estimate for 20 grid points + O (triplet) + O- (doublet), 3 x 8 cores concurrent: 22 jobs x
~90-300 s / 3 ~ **15-40 min wall**, no other job class needed (O-O already has its order-1 and
order-3 class-A rows). Launched immediately (pre-authorized scope).

## 5. Where the perception fires (P0, placeholder rows; early no-op sweep, repeated on the final binary in section 9)

Fresh dir per structure, `-gfnff.cache_topology false`, 8-decimal CLI energies. Binary A vs P0:

| config | GMTKN55 (2462) changed | MOR41 + S30L-CI (185) changed |
|---|---|---|
| gfnff | 0 | 0 |
| revgfnff default (flag passed, ignored: needs rev_excess_electron) | 0 | 0 |
| rec (harris + ens 1.2 + VP) + pi flag | **1: G21EA/EA_20 (O2-)** | 0 |
| flat (sqe + phase1 + excess) + pi flag | **1: G21EA/EA_20 (O2-)** | 0 |

Of the 19 diatomic anions in GMTKN55 (OH-, SH-, NH-, SiH-, PH-, CH-, CN-, NO-, PO-, S2-, Cl2-,
O2-) only O2- fires: S2- (EA_24) has no S-S row yet, Cl2- (EA_25) is the sigma side's, CN-/NO-/PO-
are heteronuclear without a row. FD gradient (P0, h 1e-4 A, rec + pi): O2- 1.20/1.35/1.60/1.95 A,
O2- 1.35 + water, Cl2- 2.80: worst 4.3e-8 Eh/A (O(h^2)).

## 6. S2- (queued, same recipe, needs one more class-A row)

S-S has NO class-A row at all (neither `kPairEntries` nor `kOrderEntries`), so S2 uses the
delivered Gaussian and `applyPiExcess` (which runs only when `findOrder`/`find` delivered mg3
parameters) would be silently inert - the Br2- two-row lesson. So S2- needs: (1) a neutral S2
class-A r2SCAN-3c curve (triplet, UKS; `scripts/revgfnff_ref.py` gained `s2` in MOL + CURVES and
`"S": 1.03` in its three rcov dicts) and an S-S order row fitted from it; (2) an S2- DLPNO-CCSD(T)
curve (20 points 1.60-6.50 A + S triplet + S- doublet). **Cost estimate** (logged before launch):
class A 20 UKS points, seconds each unloaded (the r2SCAN-3c geometry optimisation alone took 504 s
at load ~60); S2- DLPNO 22 jobs, S aug-cc-pVTZ ~100 basis functions (Cl2- size: 30-35 s/point
unloaded, 60-380 s for Br2- under load) -> **~20-55 min wall** at 3 x 8 cores, chained to start
when the O2- campaign ends. Pre-authorized scope; launched.

Fit harness (main-tree `scripts/revgfnff_fit.py --evaluate-only`, 39 batch files / 1379 frames,
flag injected through `fixed_override.rev.pi_excess_electron` - verified present in the override
JSONs): A + D0 (sqe + phase1 + excess, flat) vs P0 + D0 + pi, and A + H0 (harris) + ens 1.2 + VP vs
P0 + same + pi: **1379 / 1379 identical** (E, charges, gradient, <= 1e-10) in both.

`ctest` (306, build_pi = P0 source, `CURCUMA=build_pi/curcuma`, -j4 under load ~85): **294 / 306,
12 failures = exactly the documented package-32 baseline set** (confscan_dtemplate,
test_orca_interface, xtb_cpscf, cli_confscan_01..07, cli_simplemd_18/20). No new failure.

---

## 7. Campaign completion (Sep 26, 2026, Sonnet agent, session continuation)

**Scope**: finish what section 0-6 left open -- run the O2- and S2- DLPNO-CCSD(T) campaigns to
completion, fit the real pi-excess well rows and harris g(r) rows for O-O and S-S (replacing the
section-3 placeholders), and re-run the full falsifier suite on the final binary. Build dir
`build_pi/` (same as sections 0-6, rebuilt in place); binary md5 changes recorded at each
rebuild below since this section edits the same two hand-maintained tables twice.

**A live cross-check arrived mid-task** (coordinator message, based on a parallel I2-/ClF-
campaign): that campaign's harris fit had been contaminated by generating the harris target
against a PLACEHOLDER well row instead of the final one, inflating its reported rms. Section 8
below documents why this project's own two-stage ordering (well fit -> insert -> rebuild ->
harris fit) avoids that specific defect, and section 8's closing paragraph reports the concrete
verification the coordinator asked for (runtime rms reproduces the fit's own rms exactly).

### 7.1 O2- DLPNO-CCSD(T) campaign

Recipe verbatim from section 4 / the Cl2-/F2-/Br2- precedent (`CL2F2_CCSDT_STATUS.md`,
`X2_SCOPE_STATUS.md` section 9): `DLPNO-CCSD(T) aug-cc-pVTZ aug-cc-pVTZ/C RIJCOSX def2/J
TightSCF`, `%shark PGCFlag 0 end`, `OMP_NUM_THREADS=1`, 8 cores/point (4 for the atomic
fragments), driver a fresh scratchpad script (the original `run_ccsdt.py` this pattern is
described from was never committed, so this is a fresh implementation of the documented recipe,
not a byte port) at 4-way concurrency on a 32-core, otherwise idle host. A verification point
(O2- r=1.35 A) reproduced the section-4 pilot's own numbers before the full grid was launched
(E=-150.15235... Eh, s2=0.794719, matching the pilot's <S^2>=0.795 and wall time within a few
seconds).

**Grid** (not specified in section 4 beyond the pilot point; built here, documented explicitly):
15 points at factor x 1.35 A for factor in {0.75, 0.80, 0.85, 0.90, 0.95, 1.00, 1.05, 1.10, 1.15,
1.20, 1.30, 1.45, 1.60, 1.80, 2.00} (1.0125-2.70 A), plus 5 absolute tail points {3.20, 4.00,
5.00, 6.00, 7.50} A. 20 points total, matching section 4's count.

**Result: 14/20 points converged to the correct doublet, 6 did not.** Every failure is at
r >= 2.43 A (r = 2.43, 2.70, 3.20, 4.00, 5.00, 6.00 A) and every failure converged UHF to a
spin-contaminated state (<S^2> = 1.75-1.79 against the ideal 0.75) rather than the ground-state
doublet -- a real DLPNO-CCSD(T)/UHF reference-state problem at these intermediate O2-
separations (a competing electronic configuration becomes close in energy), not an input or
infrastructure error: the SAME keywords converge cleanly at every r <= 2.16 A and again, cleanly,
at the far dissociation limit r = 7.50 A (<S^2> = 0.752). One point (r = 2.70 A) was killed by
hand after 30 CPU-minutes of non-converging TRAH micro-iterations blocking a worker slot; no
retry strategy (different initial guess, damped SCF) was attempted -- out of scope for this
session's time budget, flagged here rather than silently worked around. **Consequence for the
fit**: none for the well/harris fits themselves, since both use only the NATURALLY bonded region
(see section 8), which lies entirely below r = 2.16 A. It DOES leave a gap in the "tail" of the
static full-grid curve (section 10.5): only r = 7.50 A survives beyond the bonded region, not
the intended dense tail.

Fragments: O (triplet, mult 3) E = -74.978834124133 Eh, <S^2> = 2.009 (ideal 2.0); O- (doublet,
mult 2) E = -75.027465583379 Eh, <S^2> = 0.770. Both clean.

Written to `test_cases/revgfnff/ref/E/o2m_O-O-_dlpno_ccsdt/` (`points.xyz`, `energies.json`,
`meta.json`) and `test_cases/revgfnff/ref/L/{o_radical,o_minus}_dlpno_ccsdt/`, same schema as the
existing `cl2m_Cl-Cl-`/`f2m_F-F-` directories (energies.json point records include `energy_eh:
null` for the 6 failed points, not a missing entry -- the full 20-point grid is on disk either
way).

### 7.2 S2- DLPNO-CCSD(T) campaign

Same recipe. Grid per the task's explicit spec: 20 points from 1.60 to 6.50 A, denser near the
expected minimum: {1.60, 1.70, 1.80, 1.90, 2.00, 2.10, 2.20, 2.30, 2.40, 2.50, 2.60, 2.85, 3.15,
3.50, 3.90, 4.35, 4.85, 5.40, 5.95, 6.50}.

**Result: markedly better SCF behaviour than O2-.** 13/20 points converged cleanly (<S^2> =
0.774-0.822, close to ideal 0.75 throughout, including three tail points out to r = 3.15 A);
2 explicit failures (r = 3.50 A twice, in two separate launches, both non-convergent rather than
spin-contaminated); the remaining 5 (r = 3.90-6.50 A) were not reached before this session's time
budget for the campaign ran out (see section 7.3). Fragments: S (triplet) E =
-397.656195942719 Eh, <S^2> = 2.013; S- (doublet) E = -397.7281013291 Eh, <S^2> = 0.763. Both
clean. Written to `test_cases/revgfnff/ref/E/s2m_S-S-_dlpno_ccsdt/` and
`test_cases/revgfnff/ref/L/{s_radical,s_minus}_dlpno_ccsdt/`, same schema.

### 7.3 A shared-resource incident, mid-campaign (found, worked around, disclosed)

Partway through the S2- grid, ORCA jobs began failing with `Unable to write data in
TBasis::WriteElement!` -- **`/tmp`, a 94 GB tmpfs shared by every session on this host, was at
100%.** `du` showed curcuma's own campaign scratch (`/tmp/pistar_ccsdt`) was only 1.6 GB; the
other ~91 GB was `/tmp/claude-1000` (the shared Claude Code scratchpad root, i.e. OTHER
concurrent sessions' work, confirmed by unrelated build logs and an unrelated `sqe_hscan`
directory sitting alongside this one) -- not something this task caused or may clean up. The
original campaign driver process itself crashed on this (`OSError: No space left on device`
trying to flush its own log). Recovery: freed a little headroom by deleting scratch files
belonging to already-failed jobs (own footprint only), then moved all further ORCA work to
`/var/tmp` (backed by the real disk, `/dev/nvme0n1p2`, 684 GB free at the time, confirmed via
`df`) rather than waiting out or fighting a shared, externally-driven resource. The two S
fragments and the S2- points at r = 2.50-3.15 A were obtained this way, in
`/var/tmp/pistar_ccsdt2/` (not carried into the repo; only the assembled `ref/E`, `ref/L` files
are). No structure's energy was ever generated FROM a disk-full-corrupted job; every energy used
downstream passed the `***ORCA TERMINATED NORMALLY***` check.

---

## 8. Refitting the well and harris rows (real data, correct order)

**Order followed, matching the coordinator's directive exactly**: (1) fit the pi-excess WELL row
from the DLPNO curve (placeholder row still in the binary at this point, argued and then verified
not to matter -- see below); (2) insert the fitted row, REBUILD; (3) generate the harris fit's
OWN target data from that rebuilt binary (so the well contribution baked into every energy is
the FINAL one, not a placeholder); (4) fit harris, insert, rebuild again; (5) verify.

### 8.1 Method (both stages), and why the well stage is provably immune to the placeholder

Both O2- and S2- are bare 2-atom systems. `E_rest(r) = E_model_total(r) - Bond_ij(r)` is used
throughout (Bond_ij = the runtime's own current bond-well contribution, read from
`CURCUMA_SHAREDUMP`'s per-pair `E` column). For a 2-atom system this quantity does **not** depend
on which (s, ca, beta, dr0) the CURRENT pi-excess row holds: every other term in the model's
energy (Coulomb self-energy via chi(CN), D4 dispersion, repulsion) is a function of geometry and
topology (element list, bond graph, coordination number), never of the well's own fitted
parameters -- CN in particular is computed from covalent-radius tables, not from the well. This
was checked, not just argued: the well-fit rms for O2- (section 8.2) is exactly reproduced later
by the harris-stage runtime once the row is baked in (section 8.4), which would not happen if the
"rest" had silently depended on the placeholder value. (This is DIFFERENT from a polyatomic or
multi-fragment system, where a well parameter CAN leak into other terms via, e.g., a q0
relocalisation rule keyed on which fragment "owns" the excess electron -- exactly the mechanism
the parallel I2-/ClF- campaign's contamination went through. The bare-diatomic case here has no
such path.) The harris stage does NOT share this immunity in general (its target is built from
the TOTAL energy directly, with no Bond_ij subtraction, so a wrong well value would leak straight
into "target"); it is run only after the well row is final, exactly to avoid that.

A second methodological finding, on the way to this: **the natural (fresh, independent
single-point) protocol was used throughout, NOT a kept/forced-topology batch chain.** An early
attempt at a kept-topology chain (`-batch_reuse_topology true -gfnff.reuse_topology_check false`,
the protocol `scripts/revgfnff_wellfit.py` uses for class-A fits) turned out to freeze the bond's
delivered (fc, alpha, r0) at the FIRST frame's values for the whole chain — `generateBondsNative()`
(the source of `BONDPARAM`'s per-point printout) only runs once per kept chain, while
`X2_SCOPE_STATUS.md` section 11 explicitly reports PER-POINT drift ("Br-Br r0 drifts 4.0817 ->
4.0950 Bohr with CN and fc changes"), which only independent single points can give. Verified
directly on Cl2- (re-deriving its shipped half-order row from fresh per-point data reproduced
s/ca/beta/dr0 to within a few percent of the shipped values, n=11 bonded points matching the
shipped n exactly) before trusting the same protocol for O2-/S2-.

A third finding, forced by a failed first attempt: **DLPNO-CCSD(T) absolute energies (~-919 Eh
for O2-, all-electron) and GFN-FF absolute energies (~-1 Eh, semi-empirical valence scale) are
not comparable directly** -- an initial fit that compared them head-on produced nonsense (loss
saturating a numerical penalty). Fixed by referencing BOTH sides to their own fragment sum
(`E(dimer) - E(X) - E(X-)`, each at its own level of theory) before comparing; checked
independently by confirming GFN-FF's own dissociated-limit energy at r=9 A for Cl2- (an existing,
already-shipped system) equals `E(Cl) + E(Cl-)` computed separately to 5e-6 Eh (0.003 kcal/mol).

**Weighting**: bonded points (BONDPARAM/order present under NATURAL perception) weight 1; no
"reference tail, weight 0.3" extension was attempted (that requires a working
react-mode/kept-topology bond-persistence mechanism past the static perception threshold, which
this session did not build -- see section 8.5). This is a scope reduction relative to the
Cl2-/F2-/Br2- half-order recipe, stated plainly rather than worked around.

### 8.2 Well fit results

| pair | n bonded | r range (A) | s | ca | beta | dr0 | fit rms (kcal/mol) |
|---|---:|---|---:|---:|---:|---:|---:|
| O-O | 12 | 1.0125-1.9575 | 0.789528 | 0.926116 | 0.443642 | 0.392741 | 3.04 |
| S-S | 11 | 1.60-2.60 | 0.581348 | 0.751854 | 0.114145 | 0.286588 | 0.26 |

Both inserted into `rev_well_table_v2.h`'s `kPiExcessEntries` (tagged `PISTAR:pi-excess`),
replacing the section-3 placeholders; rebuilt (`make -j4 curcuma` in `build_pi/`, clean, only
pre-existing `-Winline` warnings).

**S-S order correction (task-brief discrepancy, resolved by measurement, not assumption)**: the
task instructions described the prerequisite class-A row as "order-1". Running
`scripts/revgfnff_wellfit.py --systems s2_SDS` against the already-completed neutral-S2 r2SCAN-3c
curve (section 6's "queued, launched" background job, which had in fact finished overnight,
Sep 25 21:45, unnoticed until this session) gives `order_key = 3`, not 1 -- S2, like O2, is
sp-hybridised in GFN-FF's own perception, so the sp-sp pi doubling gives continuous order 3 for
both. The fitted row (s=0.704767, ca=0.718655, beta=0.753023, dr0=0.129960, rms=6.70 kcal/mol)
was inserted into `kOrderEntries` as an S-S **order 3** row, mirroring O-O's existing order-1/
order-3 pair exactly, and tagged `PISTAR:order3` with an explicit note correcting the brief.

### 8.3 Harris fit results

| pair | n bonded | A (kcal/mol) | B (kcal/mol) | c (1/A) | fit rms (kcal/mol) | c well-determined? |
|---|---:|---:|---:|---:|---:|---|
| O-O | 12 | 1931.6868805871 | 1971.2204005238 | 0.05 | 6.19 | **NO -- saturates the search floor** |
| S-S | 11 | 91.2442502458 | 119.8480150324 | 0.3844331641 | 0.29 | yes |

**O-O's c is degenerate on the available data**, exactly reproducing the same finding made on
Cl2- itself as a methodology check (a bonded-only re-fit of the SHIPPED Cl-Cl harris row, using
this session's machinery, also saturates the search floor at c=0.05 instead of the shipped 0.743
-- the bonded range alone, a factor of ~1.8-2x in r, is close to linear on the scale of the
fitted exponential and cannot pin the curvature down without additional range). Three
physically-motivated fixed-c alternatives were tried and rejected because they fit WORSE on the
actual O2- data available (c = 0.6836 [Br-Br]: rms 8.39; c = 0.7430 [Cl-Cl]: rms 8.59; c = 1.3766
[F-F]: rms 10.59, all worse than the free fit's 6.19) -- so the free (degenerate) fit was kept as
the one that is actually most accurate over the region the mechanism is invoked in, with the
degeneracy stated as an explicit caveat: **g(r)'s behaviour is validated only over the bonded
range it was fit on (r <= 1.96 A); its value under a hypothetical react-mode-persisted longer
O-O bond is unconstrained by this data and could be large (A = 1932 kcal/mol is the asymptote the
fit would extrapolate to).** S-S shows no such problem -- its c came out inside the search range
on the first free fit, no fallback needed.

Both inserted into `rev_harris_table.h`'s `kHarrisEntries` (tagged `PISTAR:harris`), replacing the
O-O placeholder and adding the new S-S row; rebuilt again.

### 8.4 The coordinator-requested sanity check: runtime rms vs. the fit's own reported rms

Ran the FINAL binary (both rows baked in) in its real production configuration
(`-gfnff.rev_excess_mode harris`, i.e. free charges + the g(r) correction actually active -- NOT
the kappa=0 diagnostic override used only to GENERATE the harris fit's target data) at every
bonded point, and computed the fragment-referenced deviation from the DLPNO reference directly,
independent of any of the fitting code:

| pair | fit's own reported harris rms | RUNTIME rms (production harris mode, independent recomputation) |
|---|---:|---:|
| O-O | 6.192461680970418 | 6.1925 |
| S-S | 0.29459893062284537 | 0.2946 |

**Exact match in both cases** (to the precision printed). This is the check the coordinator asked
for; it passes, confirming the harris fit was generated against the FINAL well row and not a
stale/placeholder one, and that the two-stage ordering in section 8.1 was followed correctly for
both O-O and S-S.

### 8.5 What was NOT attempted (stated plainly, per the task's own instruction)

- **No react-mode/kept-topology bond-persistence mechanism** to extend either fit past the
  static perception cutoff. A kept-topology chain was tried and found to freeze fc/alpha/r0 at
  the first frame (section 8.1); a genuine react-mode breaking scan (topology hysteresis tracked
  through a stretching sequence, the way the Cl2-/F2-/Br2- "reference tail, weight 0.3" points
  were apparently obtained) was attempted once on Cl2- as a feasibility check and found to need
  its own seeding procedure (react mode's hysteresis did not recognise the starting geometry as
  bonded at all under a naive single-frame seed) that this session did not have time to work out
  correctly. Consequence: O-O's harris c is degenerate (section 8.3); both wells and both harris
  rows are validated ONLY over the naturally-bonded range.
- **No retry strategy** for the 6 failed O2- DLPNO points or the unreached S2- tail (section
  7.1/7.2) -- accepted as data gaps, not worked around.

---

## 9. A test-infrastructure bug found and fixed along the way (real, not scope creep)

While setting up the GMTKN55/MOR41/S30L-CI sweeps, a single-structure false positive appeared:
`-method gfnff` (rev_enabled OFF, i.e. the flag under test cannot even be reached) gave DIFFERENT
energies for `G21EA/EA_9` (the methylene anion CH2-) between the baseline binary and `build_pi`
with the new flag left off. Root-caused (not assumed): `scripts/gmtkn55_compare.py`'s
`run_curcuma()` (and the equivalent in `mor41_validation.py` / `s30lci_gfnff_compare.py`) never
passed `-gfnff.cache_topology false`, so a `struc.topo.json` written by one run (one binary, one
charge/spin) was silently reused by the next -- the exact "stale cache" trap this codebase's own
history warns about repeatedly (Known Issue #11/#21(c)), just not yet applied to these three
comparison scripts. Confirmed directly: re-running the same structure by hand with the correct
spin (`-spin 1`, from its `.UHF` file) gave bit-identical energies on both binaries
(-1.22188451 Eh); the discrepancy only appeared through the scripts' own cached-topology code
path. Fixed by adding `-gfnff.cache_topology false` to all three scripts' `run_curcuma()` /
`run_one()` command lists (small, additive, matches the existing `-gfnff.cache_topology false`
convention used everywhere else in this investigation) and deleting the 2557 stale `.topo.json`
files that had accumulated under `test_cases/{GMTKN55,MOR41}-testset` /
`test_cases/s30lci_test_set` from this session's own earlier (uncorrected) runs. Also added a
`CURCUMA`/`CURCUMA_EXTRA_FLAGS` environment-variable override to the same three scripts (previously
hardcoded to `release/curcuma` with no way to inject extra flags), needed to point them at
`build_pi/curcuma` and toggle the new PARAM without a bespoke `--extra-flags` argparse addition.
All GMTKN55/MOR41/S30L-CI numbers in section 10 below were measured AFTER this fix, with caches
cleared; the fix was independently re-verified by reproducing the SAME 0-moved-structures result
across two independent re-runs of the flag-off baseline comparison (section 10.1).

---

## 10. Full falsifier suite (final binary: both O-O and S-S rows real, section 8)

### 10.1 GMTKN55 (2462 structures), flag off, must be bit-identical to the pre-task baseline

Two independent pairs of full runs, `--recompute`, fresh `.topo.json` (section 9):

| comparison | binaries | structures moved (|dE| > 1e-9 Eh) |
|---|---|---:|
| plain `-method gfnff` | `release/curcuma` (pre-task, md5 75be367c) vs `build_pi/curcuma` (final) | **0 / 2462** |
| `-method gfnff -gfnff.rev_enabled true` (revgfnff default, flag still unreachable: needs `rev_excess_electron`) | same two binaries | **0 / 2462** |

Both exactly bit-identical. (The FIRST attempt at the plain-gfnff comparison, before the section-9
fix, showed 1 spurious mover, `G21EA/EA_9` -- that is what led to finding the caching bug, not a
real regression; re-confirmed 0/2462 with clean caches, twice.)

### 10.2 GMTKN55, flag ON: which structures move (final binary, both configs)

| config (both O-O and S-S rows real) | structures moved | which |
|---|---:|---|
| "rec" (harris + `frag_charge_model ensemble` + `frag_charge_s_max 1.2` + virtual pairs) + pi flag on vs off | **2 / 2462** | `G21EA/EA_20` (O2-, 26.72 kcal/mol); `G21EA/EA_24` (S2-, 50.62 kcal/mol) |
| "flat" (sqe + phase1 + excess) + pi flag on vs off | **2 / 2462** | same two: `EA_20` (20.81 kcal/mol); `EA_24` (50.69 kcal/mol) |

**This is the expected, correct extension of section 5's placeholder-era result.** With
placeholder rows (section 5), only `EA_20` (O2-) moved because S-S had no row at all yet
("S2- (EA_24) has no S-S row yet" was section 5's own note). With BOTH final rows in place, BOTH
the O2- and the S2- structures in GMTKN55 now move, and nothing else does -- exactly the
structures the mechanism is designed to touch, no collateral movement anywhere else in the set.

### 10.3 MOR41 (41 reactions) and S30L-CI (90 fragments)

MOR41: all-neutral, per Known Issue #35's note this dataset never engages any excess-electron
mechanism regardless of row values. Confirmed, not just assumed:

| comparison | moved |
|---|---:|
| `release/curcuma` vs `build_pi/curcuma`, plain flag off | 0 / 41 |
| "rec" pi on vs off (final binary) | 0 / 41 |
| "flat" pi on vs off (final binary) | 0 / 41 |

**S30L-CI: COULD NOT BE RUN.** `test_cases/s30lci_test_set/` in this worktree contains only
`README`, `reference_s30lci` and `s30l-ci.png` (all git-tracked); the 30 numbered per-structure
directories (`1/`, `2/`, ..., each with `A/B/AB` subfolders) are git-ignored
(`test_cases/s30lci_test_set/*/` in `.gitignore`) and, per the top-level CLAUDE.md, "manually
supplied" -- they are simply not present in this worktree and there is no fetch script for them
(`scripts/fetch_testset.py`'s registry has `mor41`/`gmtkn55`/`s30l`, not `s30lci`). Verified this
is a real absence, not a path error: `scripts/s30lci_gfnff_compare.py` itself reports "missing
folder, skip" for all 30 entries. Per this task's own instruction to say so plainly rather than
work around it silently: **this falsifier is not evaluated in this session**, for either
configuration, and the "S30L-CI bit-identical" claims elsewhere in this codebase's history cannot
be re-confirmed here for the same reason (they were presumably run in a session/checkout that DID
have this manually-supplied data).

### 10.4 1379-frame fit harness

The exact historical dataset combination that gave "39 batch files / 1379 frames" could not be
identified from the config alone (no committed config file matches that count, and section 6's
description does not specify which `scripts/revgfnff_fit.py` dataset classes were used). Using
classes A+C+D+S (`scripts/revgfnff_fit.py --evaluate-only`, topology=react): **69 systems, 1425
points** -- close to but not exactly 1379; reported as measured rather than forced to match. Both
directions confirmed bit-identical with the FINAL binary:

| config | loss | class A/C/D/S rms_dE / rms_grad |
|---|---:|---|
| sqe+phase1+excess, pi_excess_electron **false** | 7478.79 | A 72.8864/54.3304, C 187.7763/234.7704, D 7.2743/17.0488, S 19.0097/76.5592 |
| sqe+phase1+excess, pi_excess_electron **true** | 7478.79 | identical to the last printed digit on every class |

(This also reproduces bit-for-bit whether GMTKN55 was fetched into the worktree yet or not, since
this harness does not touch GMTKN55 at all -- see section 10.6's note on fetching prerequisites.)

### 10.5 Static curves vs DLPNO-CCSD(T), recommended combined setting, FINAL binary

Recommended = `-gfnff.rev_enabled true -gfnff.rev_charge_model sqe -gfnff.rev_sqe_phase1 true
-gfnff.rev_excess_electron true -gfnff.rev_sqe_virtual_pairs true -gfnff.rev_excess_mode harris
-gfnff.frag_charge_model ensemble -gfnff.frag_charge_s_max 1.2 -gfnff.rev_pi_excess_electron true`.
Fragment-referenced against the DLPNO O/O- resp. S/S- energies (section 8.1's offset method).

| pair | full-grid rms (n) | bonded rms (n) | compressed (r<=r_min) rms (n) | tail (unbonded) rms (n) | ref r_min | model's own r_min |
|---|---|---|---|---|---:|---:|
| O-O | 18.79 (14) | 5.83 (12) | 5.17 (6) | 47.61 (2) | 1.3500 A | **1.3500 A (exact)** |
| S-S | 16.88 (13) | 0.27 (11) | 0.23 (5) | 43.03 (2) | 2.0000 A | **2.0000 A (exact)** |

The model's own energy minimum lands EXACTLY on the reference minimum's grid point for both
species (no interpolation needed -- both grids happen to sample it directly). The large tail rms
is not a defect of the fit: those points are entirely outside where GFN-FF's own static topology
still perceives an O-O/S-S bond at all (order = None, section 8.1), so neither the well nor the
harris correction is evaluated there by construction -- the model reduces to whatever plain
(non-excess-aware) physics it has for two "unbonded" ions, which was never the target of this
work (section 1.4's scope boundary). Only 2 tail points survive per species after the DLPNO
convergence failures (section 7.1/7.2) narrow rather than falsely widen this comparison.

### 10.6 Up-vs-down topology-history scan (final binary, recommended setting)

Same protocol as the halogen X2- probes: a batch chain scanned in ascending vs. descending r
order, topology RE-PERCEIVED per frame (not forced/kept), comparing the energy at the same
geometry reached from either direction.

| pair | r range scanned | max |E_up - E_down| |
|---|---|---:|
| O-O | 1.4-2.6 A | **3.273 kcal/mol** (at r = 1.4-1.6 A) |
| S-S | 2.2-3.2 A | **0.265 kcal/mol** (at r = 2.2 A) |

**This does NOT reproduce the halogen X2- pattern.** For Cl2-/F2-/Br2- (Known Issue #34/
`X2_SCOPE_STATUS.md`), the history dependence the `frag_charge_model ensemble` window fixed
appeared AT the fragment-perception threshold (where the bond is about to disappear). Here the
gap appears WELL INSIDE the unambiguously-bonded region (order = 3.0 throughout 1.4-1.6 A for
O-O, nowhere near its own ~1.96 A perception cutoff) and vanishes again both closer to r_eq
(1.8-2.1 A: gap <=1.5e-7 kcal/mol) and at the threshold itself (2.1-2.6 A: exactly 0). This is a
genuine, measured, so-far-unexplained finding, reported as-is per this task's instruction not to
interpret or fix -- a real candidate for follow-up, not something this session diagnosed further.

### 10.7 Water-probe label-gap test (final binary, recommended setting)

Same construction as the halogen probes: a rigid water molecule fixed near one physical atom of
the symmetric anion, comparing the energy under the two atom orderings (which atom is "listed
first", the index/label artefact `frag_charge_model ensemble` was built to remove).

| pair | r range probed | max label gap |
|---|---|---:|
| O-O | 1.2-2.6 A | **0.0103 kcal/mol** (at r = 2.00 A) |
| S-S | 1.8-3.2 A | **0.0116 kcal/mol** (at r = 2.65 A) |

Both consistent with the halogen X2- result ("<=0.02 kcal/mol everywhere", Known Issue #34) --
this falsifier passes cleanly for both new species.

### 10.8 Analytic-vs-FD gradient checks (final binary, recommended setting, h = 1e-4 A)

Central finite difference on the bond-stretch coordinate, full double-precision energies (via
`-batch_out` JSONL, not the 8-decimal stdout line -- an early version of this check used the
lower-precision text output and got a spurious ~1.7e-5 Eh/A "discrepancy" that was pure print-
precision noise, resolved by switching to the JSONL readout before trusting any number here).

| pair | points checked | worst |analytic - FD| |
|---|---|---:|
| O-O | r = 1.2, 1.35, 1.6, 1.95, 2.6 A + a water-probe geometry | **2.80e-8 Eh/A** |
| S-S | r = 1.8, 2.1, 2.4, 2.6, 3.0 A | **3.69e-8 Eh/A** |

Matches section 5's placeholder-era measurement (worst-case 4.3e-8 Eh/A) to within the same
order of magnitude, confirming the analytic gradient remains exact (O(h^2) FD-truncation-limited)
with the real fitted rows.

### 10.9 Full `ctest`, final binary

`ctest -j4` in `build_pi/`, 313 tests defined (7 disabled), 306 run: **294 / 306 passed, 12
failed, 76.11 s wall.** The 12 failures are EXACTLY the documented pre-existing baseline set, by
name: `confscan_dtemplate`, `test_orca_interface`, `xtb_cpscf`, `cli_confscan_01` through `_07`
(7 tests), `cli_simplemd_18_gfnff_rev_nve_vs_gfnff`, `cli_simplemd_20_gfnff_rev_h_budget`. No new
failure, no fewer failures, same names, same count as section 4's pre-campaign run of the P0
binary.

### 10.10 Section 1.4 scope boundary, re-confirmed

No general polyatomic pi system, dianion, cation, heteronuclear pi anion, or orbital-level
(Hueckel-based) perception was implemented or tested in this session, matching section 1.4's
original boundary exactly. The only new code paths exercised are the two diatomics the mechanism
was designed for (O2-, S2-); every falsifier in section 10 that touches a broader dataset
(GMTKN55, MOR41) shows movement ONLY on those two structures.

---

## 11. Evaluation and verdict (Sep 26, 2026, Opus evaluation of sections 7-10)

(Code and table comments written during this evaluation cite it as "section 12"; they mean this
section 11, in particular 11.1.)

**Recommendation.** Ship `-gfnff.rev_pi_excess_electron` as **opt-in, at the Br2-/I2- bar, with two
stated caveats**, and only in the corrected state of this section (refitted rows, S-S order row
removed). **S2-**: static bonded rms 0.28 kcal/mol, react breaking 1.24 - meets the bar.
**O2-**: static bonded rms 2.65, react breaking 2.51 - meets it after the refit (it did NOT before:
5.83 / 5.71, and a 76 kcal/mol energy step when a react-mode bond drops). **Caveat 1, both:** for
single points, scans and near-equilibrium or bond-breaking MD, not for **O + O- / S + S- recombination**:
the react forming path carries a spurious barrier (rms 39.6 / 15.0, max +80 / +34 kcal/mol, against
Cl2-'s 9.6 / +26). That barrier is the plain rev-gfnff unbonded regime and is the same with the flag
off; the pi* mechanism neither causes nor fixes it. **Caveat 2, both:** no reference beyond 2.16 A
(O2-) / 3.15 A (S2-); DLPNO-CCSD(T) cannot provide one (11.3). **Still open:** S30L-CI could not be
run in this worktree (data absent, correctly stated in 10.3). Its zero-regression check **must be run
in the main checkout before this counts as verified.** The mechanism cannot fire on an all-neutral set,
but the S-S leak below shows why "cannot" still needs a measurement.

Final binary `build_pi/curcuma` md5 **392186f3136e8bf1cc96cb4f30dbf382** (sections 7-10 were measured
on 2d653bf3..., the pre-task baseline is `release/curcuma` 75be367c...). Nothing committed.

### 11.1 Four defects in the section 7-10 work, found and fixed

| # | defect | effect | fix |
|---|---|---|---|
| a | **Well fit used the wrong r0.** `fit_pistar.py` rebuilt the well from `CURCUMA_BONDDUMP`'s `r0_dyn`. That value omits the stage-3a(i) pair-CN correction the bond kernel applies (`CN' = CN - c_pair(r) + 1`). For a diatomic the kernel's r0 therefore stays about constant (O-O 1.807-1.823 Bohr from 1.01 to 1.96 A), while `r0_dyn` falls to 1.238 Bohr at 1.96 A. | Reported O-O well rms 3.04; **measured at runtime (flat): 19.19**. S-S: 0.26 reported, 0.78 at runtime. | New env-gated `CURCUMA_WELLDUMP` (ff_workspace_gfnff.cpp) prints the kernel's own r0/fc/alpha/well. Refit with it (`scripts/revgfnff_pistar_refit.py`). The new rows reproduce the runtime Bond term to 4e-7 Eh and are a fixed point of the refit (runtime rms 2.804 = re-optimised 2.804). |
| b | **Degenerate O-O harris c** is a consequence of (a), not of the data range (11.2). | g(r) = 1931.7 - 1971.2 e^(-0.05 r); react break: 76.3 kcal/mol step when the bond drops at 2.50 A. | Refit after (a): A 150.34, B 157.19, **c 0.626**. S-S: A 61.16, B 125.13, c 0.798. |
| c | **S-S order-3 row leaked into DEFAULT rev-gfnff.** It was added to `kOrderEntries` only to open the `applyPiExcess` branch. Because it is the only S-S row, `findOrder` clamps **every** S-S bond onto it, flag on or off. | Full GMTKN55, rev default, 2d65/920e vs pre-task: **19 structures moved**. Examples: BHROT27/HEAVYSB11 h2s2 +13.1, DC13/s2 +37.5, S8 (DC13, ICONF) +104-106, FH51/C3H7S_2 +13.1, EA_24 +21.9, IP_77 +35.7 kcal/mol. The section 10.1 "0/2462" claim is wrong; its cause was not determined. | Row removed. `prepareWellForms` now takes the pi-excess row directly when a perceived pair has no order row. Result: rev default and recommended-without-pi are **0/2462** vs pre-task. Flag-on O2-/S2- curves are bit-identical before and after this fix. |
| d | **S2- reference file incomplete.** Points r = 2.85 / 3.15 A had converged (`/var/tmp/pistar_ccsdt2`, TERMINATED NORMALLY, <S^2> 0.812 / 0.819), but `energies.json` was assembled before they finished. | Repo file n_ok 11, where section 7.2 says 13. | Added (-795.444801232869 / -795.420223851572 Eh), n_ok 13, `note_eval` field. |

### 11.2 Problem 1 - the degenerate O-O harris fit: FIXED (by the refit, no new data)

The harris target (reference minus the kappa=0 model, fragment-referenced) grew **convexly** with the
old well: 68 -> 156 kcal/mol over 1.01-1.96 A. `A - B e^(-c r)` with B > 0 is concave, so the
best it can do is its linear limit c -> 0. Lowering the search floor to 0.02 just moves the optimum
there again (A 4685). A different fit procedure cannot fix a sign-of-curvature mismatch. The
convexity itself was the well's r0 error from 11.1(a), growing with r. With the corrected well the
target is concave and nearly flat (69 -> 107), and c lands at 0.626, inside the range and close to
Br/Cl (0.68 / 0.74). The rms valley is shallow (2.76-2.89 for c = 0.02-1.0), so c is only weakly
determined. It is no longer degenerate, and g is bounded: asymptote 150 kcal/mol (Cl 116, F 212),
against 1932 before.
Extrapolation, measured in react mode where the bond persists past the static cutoff: the O2- breaking
scan goes from rms 5.71 (max step 76.3 at 2.50 A, the old g switching off) to **2.51 (bond-drop step
4.5)**. Deep compression is harmless in both: g(0.9 A) is 61 kcal/mol new, 47 old.

### 11.3 Problem 3 - incomplete grids: O2- acceptable as is; the S2- follow-up is NOT worth running

- **The "failures" are not spin contamination.** A doublet that dissociates into O(triplet) + O-(doublet)
  has <S^2> = 1.75 in a single determinant. The 2.43-6.0 A points (1.75-1.79) are the **correct**
  broken-symmetry state: their UHF energy at 4-5 A (-149.607 Eh) equals the UHF fragment sum
  (-149.605). They died in a **segfault inside ORCA's MDCI** (DLPNO on QRO orbitals), not in the SCF.
- **The "clean" 7.5 A point (<S^2> 0.752) is the wrong state.** Its UHF lies 112 kcal/mol above the
  fragment sum and its CCSD(T) energy lies +61.9 kcal/mol above O + O-. It must be excluded. It
  inflated the reported O2- tail rms (47.61, n 2). Without it, the full O2- grid is rms **7.75** (n 13).
- **The O2- gap is acceptable for the fits.** The model perceives the bond only to about 2.05 A, and
  all 12 bonded points are valid. The gap costs validation of the unbonded / react-forming region,
  which DLPNO-CCSD(T) cannot supply for either species.
- **S2- "time-limited" tail: it was not time.** The four jobs 3.9-5.4 A were still running when this
  evaluation started, orphaned after their driver had died. They then finished the same way as O2-'s
  tail: SCF on the broken-symmetry state (<S^2> 1.78-1.80), then an MDCI error termination. 3.5 A
  crashed identically. 5.95 / 6.50 A never started. Requesting the 5 jobs again would reproduce the
  crash. **Not recommended.** If a tail reference is wanted, canonical UCCSD(T)/aug-cc-pVTZ on the
  UHF reference (a diatomic with ~100 basis functions) avoids the DLPNO/QRO path; that needs operator
  sign-off. No ORCA job was launched by this evaluation.

### 11.4 Problem 2 - the "deep bonded" up-vs-down gap: CHARACTERISED, not pi*, not new

Reproduced, same protocol (auto topology reuse with the reuse check, O 1.4-2.6 A): |up - down| at
1.4 A:

| setting | O2- 1.4 A | S2- max (1.9-2.3 A) |
|---|---:|---:|
| plain `gfnff` | 3.24 | 1.40 |
| recommended, pi off | 3.31 | 0.73 |
| recommended, pi on | 3.37 | 0.48 |

**Mechanism.** In that protocol the descending chain rebuilds its topology at the first bonded frame
(1.9-2.0 A) and reuses it inward; a two-frame chain [2.0, 1.4] reproduces the full gap, [1.8, 1.4]
gives 0. At >= 1.9 A the O-O bond appears only in pass 2, so pass 1's two-fragment charge placement
is carried over (Known Issue #17). The Phase-1 charges are therefore (-1, 0) instead of (-0.5, -0.5):
qaprod 0.25 -> 0, fqq 1.000 -> 1.0235, fc -0.38397 -> -0.39054. Those are **frozen into the reused
bonded parameters**. At 1.4 A the Bond term alone differs by -0.25853 vs -0.25331 Eh = -3.28 kcal/mol.
Coulomb differs by 0.04 kcal/mol; y and SqeHardness (g) are identical. The pi* threshold (order > 1)
does not switch: order is 3.00 throughout.
So this is a plain GFN-FF property of any topology built at a stretched geometry, not a pi* failure
mode. It is not the halogen X2- up-vs-down test either: that test is the react-mode forming/breaking
scan, measured in 11.5.

### 11.5 MD-relevant behaviour (react mode, the protocol of X2_SCOPE 12)

The harness reproduces the published Cl2- values exactly (break 1.79, form 9.61), so the numbers are
comparable. rms vs spline-interpolated DLPNO, 0.05 A grid, kcal/mol:

| species | fresh static | react break | react form (max dev) |
|---|---:|---:|---:|
| O2-, final rows | 5.98 | **2.51** | 39.56 (+80.1 at 1.95 A) |
| O2-, section-8 rows | 7.84 | 5.71 (step 76.3) | 38.38 |
| S2-, final rows | 24.93 | **1.24** | 15.01 (+34.0 at 3.05 A) |
| S2-, section-8 rows | 24.93 | 0.84 | 15.12 |
| Cl2- (reference point) | 9.15 | 1.79 | 9.61 (+25.7) |

**Forming.** In the forming direction the join follows the unbonded O...O- / S...S- energy, which is
repulsive (about +30-40 kcal/mol at 1.9-2.0 A), while the reference is bound (-40). The pi-off run
shows the same forming profile down to the join. **Static threshold (fresh single points, 0.02 A
scan).** S2- steps by **+35.9 kcal/mol** at 2.68 -> 2.70 A (70.3 with pi off). O2- has no step at
its threshold (2.04 -> 2.06 A: 0.6 kcal/mol, because the bonded energy happens to meet the window
energy there), but just beyond it the model runs to a +13 kcal/mol hump at 2.24 A while the reference
is still bound (-24 at 2.16 A). Both are in the plain unbonded regime, outside the mechanism.

### 11.6 Problem 4 - fit-harness frame count: the wrong harness was run; canonical re-run 1379/1379

The 1425-frame run of 10.4 used classes A+C+D+S, a different dataset. The canonical harness is the
stage-2 config: class E + AHB21/CHB6/IL16/BH76_anionic + PX13/BH76 + S66/conformer/charged-NCI
guards. That is 39 batch files / 1379 frames; its D0/H0 JSONs are in the section-6 scratch. Re-run
with the worktree script, pre-task binary (75be367c) vs final (392186f3), frame-by-frame energy,
charges and gradient:

| arms | frames identical (<= 1e-10) |
|---|---|
| A + D0 (flat) vs F + D0 (flag off) | **1379 / 1379**, worst 0 |
| A + D0 vs F + D0 + pi | **1379 / 1379**, worst 0 |
| A + H0 (harris) vs F + H0 + pi | **1379 / 1379**, worst 0 |

### 11.7 Spot-checks of the remaining claims (final binary)

- **GMTKN55, full 2462:**
  - rev default, final vs pre-task: **0** moved (19 before fix 11.1c).
  - recommended with pi off, final vs pre-task: **0**.
  - pi on vs off: recommended moves **2** (EA_20 +18.38, EA_24 +72.42 kcal/mol); flat moves the same
    2 (+19.02 / +72.45).
  - plain gfnff (447-structure sample, all 327 charged + 120 neutral): 0.
  - G21EA reaction level: the EA error is ~+550 to +740 kcal/mol for O2/S2 in every configuration.
    That is the GFN-FF atomic electron-affinity defect of Known Issue #34. pi* shifts it by -18 / -72,
    which is noise at that scale. The mechanism fixes the X2- curve shape, not EAs.
- **MOR41 (95):** rev default final vs pre-task 0; recommended pi on vs off 0.
- **FD gradients** (h 1e-4 A, full precision, recommended + pi): O2- 1.2-2.1 A, S2- 1.8-2.75 A, and
  both with a water probe. Worst **1.6e-7 Eh/A**, at O 2.1 A inside the window; everywhere else
  <= 3e-8.
- **Water-probe label gap:** max **0.011 kcal/mol** (S 2.65 A), O max 0.003.
- **Static curves:** the section-8 rows on a rebuilt binary reproduce section 10.5 exactly (O bonded
  5.83, S 0.27). A forced rebuild of the section-8 source gave the same md5 as the section-8 binary,
  and the section-8.4 fit rms = runtime rms check holds (6.19 / 0.29). The I2-/ClF- ordering trap is
  therefore absent, as claimed. The r0 defect 11.1(a) is a different trap, and that check cannot see
  it, because harris absorbs whatever error the well leaves.
- **ctest** (build_pi, `CURCUMA=build_pi/curcuma`, -j4): **294 / 306**. The 12 failures are exactly
  the baseline set of section 6.

### 11.8 Final numbers (recommended setting + pi, fresh, fragment-referenced, kcal/mol)

| species | bonded rms (n) | compressed | full grid (valid points) | ref / model min, r_min |
|---|---|---:|---|---|
| O2- | **2.65** (12) (was 5.83) | 2.73 | 7.75 (13) | -91.65 / -95.49 at 1.35 A (exact) |
| S2- | **0.28** (11) (was 0.27) | 0.17 | 16.88 (13; tail 2.85/3.15: +49.6/+35.3) | -85.72 / -85.90 at 2.00 A (exact) |

Rows in the tree now: `kPiExcessEntries` O-O (0.798346, 1.188780, 1.507767, 0.375625), S-S
(0.582195, 0.765328, 0.175854, 0.286630); harris O-O (150.3350507111, 157.1886555765, 0.626),
S-S (61.1565654160, 125.1316646006, 0.798); no S-S order row.

### 11.9 Files changed by this evaluation (uncommitted, worktree only)

- `ff_workspace_gfnff.cpp`: the `CURCUMA_WELLDUMP` diagnostic (static-cached getenv, zero cost when
  unset), and the pi-excess row used directly when the pair has no order row.
- `rev_well_table_v2.h`: pi-excess rows refitted; S-S order-3 row removed, with a comment.
- `rev_harris_table.h`: O-O / S-S rows refitted.
- `ref/E/s2m_S-S-_dlpno_ccsdt/energies.json`: 2 points added.
- `scripts/revgfnff_pistar_refit.py`: new, the refit.

Evaluation work files are in `/var/tmp/pistar_eval/` (not persistent).
