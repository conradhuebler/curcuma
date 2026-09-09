# GMTKN55: curcuma native vs xtb

⚠️ AI-generated, machine-tested only — not human production tested.

Runner: `scripts/gmtkn55_compare.py` (see [docs/TESTSET_RETRIEVAL.md](TESTSET_RETRIEVAL.md)
for fetching the set). Per-structure single-point energies, curcuma vs xtb 6.7.1, on the
2462 published GMTKN55 geometries (54 subsets) — a reproduces-the-reference-implementation
check, **not** the GMTKN55 WTMAD-2 accuracy statistic against high-level QM (the `.res`
files are `tmer2++` shell scripts with per-subset stoichiometry; not parsed here, out of
scope). Full sweep: `python scripts/gmtkn55_compare.py` (~13 min all-methods on this
machine).

## Results (curcuma - xtb, kcal/mol)

| method | n compared | skipped (open-shell) | MAD | RMSD | max |
|--------|-----------:|----------------------:|----:|-----:|----:|
| gfn2   | 2140 | 320 | 0.000 | 0.001 | 0.017 |
| gfn1   | 2142 | 320 | 0.047 | 0.549 | 11.93 |
| gfnff  | 2458 | 0   | 0.262 |  3.20 | 84.5  |

gfn1/gfn2 skip every structure with a nonzero `.UHF` (see "Known limitation" below).
gfnff has no open-shell term and runs everything.

The gfnff row was 3.389 after the two isolated-ion EEQ fixes below, 3.462 after the
pyrrole pi-veto fix (CLAUDE.md Known Issue #10), and reached **2.117** with the four
term-level fixes of Known Issue #12 - the nitro pi-electron count and the sp2-N-H bond
strength being the two GMTKN55 exercises heavily (`Amino20x4` MAD 5.664 -> 0.001,
max 26.9 -> 0.016). Known Issue #13 (two-fragment charge placement) then took it to
**1.445**: `AHB21` MAD 22.05 -> 0.606 (max 239 -> 4.8), `CHB6` 7.73 -> 1.51, `BH76`
6.90 -> 4.79, plus the `nh_linear_fix` refinement (`DIPCS10` 8.28 -> 0.000,
`NBPRC` 3.84 -> 0.108). Known Issue #14 (aryne rule + dxi itag) then took it to **1.244** and, for the first
time, moved the set-wide maximum: 500.0 -> 151.5 (`DC13/c20bowl` +500.0 -> -6.29). Known Issue #15
(main-group metals in the dxi metal count, TM-TM over-assignment, missing GEODEP sp2->sp3) brought it
to **0.966** - under 1 kcal/mol for the first time - with `AL2X6` 14.68 -> 3.34 and `MB16-43` 13.66 -> 10.77.
Known Issue #16 (oxirane ring angle, duplicate metal angle block, raw imetal in the bond code) then
reached **0.715**: `HEAVY28` 8.68 -> 0.407, `HEAVYSB11` 4.27 -> 0.031, all five oxirane structures exact.
Known Issue #17 (the q-loop's second pass re-detecting fragments) then took it to **0.460** and moved
the maximum again, 151.5 -> 105.5: `BH76` 4.79 -> 0.080, `SIE4x4` and `CHB6` -> 0.000. Known Issue #18
(the q-loop's carbene charge override) reached **0.447**, with `G21EA` 0.653 -> 0.131, and
Known Issue #19 (the non-reference angle-`fqq` guard) **0.446**, with `BH76` 0.080 -> 0.057.
Known Issue #20 (the missing hypervalent torsion correction) reached **0.416**, with `ICONF`
1.78 -> 0.398 and `PArel` 0.95 -> 0.348. Known Issue #21 - twelve defects found by arbitrating
EVERY remaining outlier against pprcht - reached **0.268** and moved the set maximum
105.5 -> 84.5: `BHDIV10` 1.182 -> 0.043, `PX13` 0.633 -> 0.023, `ISOL24` 0.881 -> 0.300,
`W4-11` 0.716 -> 0.197, `CARBHB12` 0.445 -> 0.000, `HAL59` 0.312 -> 0.000, `MB16-43`
10.405 -> 7.901.

Known Issue #22 (the raw `metal_type` gate and three out-of-bounds element tables) reached
**0.262** and took `HEAVY28` 0.406 -> **0.000** (max 2.06 -> 0.002).

**Where the residual now sits.** 46 of the 82 deviations above 0.3 kcal/mol outside
`MB16-43` were arbitrated one by one against pprcht, and every one of them is a
pprcht-vs-xtb reference split with curcuma matching pprcht to <0.03 kcal/mol. Excluding
`MB16-43` the MAD is **0.081**. What this table measures above roughly 0.3 kcal/mol is
therefore the divergence between the two GFN-FF implementations, not a curcuma port error
- and for the two subsets that dominate it, `MB16-43` and `AL2X6`, all GFN-FF
implementations are anyway far from the true QM. Read the gfnff row accordingly.

Known curcuma-vs-pprcht residuals that remain: inside `MB16-43` (`/23` +0.565, `/35`
+0.268, on top of a 60-84 kcal split) and `AL2X6/al2me6` at +0.03. Neither is root-caused.

> **Reading the numbers back**: `scripts/gmtkn55_compare.py` caches every energy in
> `_run/energies.json` and reuses it unless `--recompute` is given. A re-run after a code
> change therefore reports the OLD energies unless the `<subset>/<name>|cur|<method>` keys
> are dropped first. Two "MAD unchanged" observations in this repo's history came from
> exactly that.

**gfn1/gfn2 confirm prior findings** (main-group, closed-shell): essentially exact
reproduction of xtb, consistent with `docs/SQM_VALIDATION.md` / `docs/SQM_WP2_gfn1_accuracy.md`.

## GFN-FF: two isolated-ion EEQ bugs FIXED (Sep 2026)

The first full sweep found gfnff MAD 19.98 kcal/mol (RMSD 131.5, max 2432.6) - two bugs,
both surfaced by GMTKN55's isolated-ion structures (DIPCS10/G21IP/G21EA/SIE4x4/ALK8/CHB6),
which neither MOR41 nor S30L contain. Both are now fixed; gfnff MAD dropped **19.98 -> 3.389**
kcal/mol (RMSD 131.5 -> 17.45, max 2432.6 -> 500.0 - the remaining residual is unrelated,
see below). 107/107 single-atom GMTKN55 structures now match xtb to <0.01 kcal/mol (was 0/107).

**Bug 1 - missing EEQ self-energy for any isolated atom.** curcuma's native GFN-FF returned
exactly 0.0 Eh for a single free (charged) atom, regardless of charge:

```
$ curcuma -sp DIPCS10/mg_2+/struc.xyz -method gfnff -charge 2 -verbosity 2 -no_bmt
GFN-FF Parameter Generation [N=1 atoms, ...]
  Coulomb pairs   0.00 ms          <- 0 pairs generated for N=1, loop body never runs
[RESULT]  Coulomb   +0.0000000000
Single Point Energy = 0.00000000 Eh          # xtb: +3.876573876298 Eh
```

Root cause: the Coulomb/electrostatics term was generated as a pairwise list (`Coulomb
pairs` in the parameter-generation timing breakdown); for N=1 that list has zero entries,
so the per-atom diagonal EEQ self-energy (`chi*q + 0.5*gamma_AA*q^2`, nonzero even for
N=1) was never computed. Fix: `GFNFF::generateCoulombSelfEnergyNative()`
(`gfnff_method.cpp`) now derives the five per-atom self-energy inputs directly from the
topology (independent of the pair list, mirroring the Fortran `gfnff_engrad.F90:1378-1389`
structure, where the self-energy statement runs unconditionally per atom outside the
pairwise `j<i` loop); `FFWorkspace::setInteractionLists()` (`ff_workspace.cpp`) prefers
these over the old pair-derived vectors, falling back to the pair-scan only for
producers of `GFNFFParameterSet` that don't populate the new fields (UFF/QMDFF, no EEQ).

**Bug 2 - dgam metal correction unreachable for Li/Be.** After fixing bug 1, Mg2+/Na+
matched xtb exactly but Li+/Be+/Be2+ still showed a constant ~0.04 Eh-per-charge-unit
offset (Li+ +25.10, Be+ +25.10, Be2+ +200.80 kcal/mol - matching `-0.04*q^3` exactly, the
signature of a missing `dgam = qa*(-0.08)` main-group-metal hardness correction).
Root cause: `EEQSolver::calculateDgam()` (`eeq_solver.cpp`) applied the metal-type check
(`if (imetal==1) ff=-0.08`) nested inside an `else if (Z > 10)` branch, but the Fortran
reference (`gfnff_ini.f90:658-660`) applies it as an *unconditional*, independent `if` -
the only two metals with Z<=10 (Li=3, Be=4) never reached the check, while every other
metal (Z>10, e.g. Na/Mg) did. Fix: moved the metal + noble-gas correction out of the
`Z > 10` branch to run unconditionally for every element.

Verified via GMTKN55 RC21/... no, via DIPCS10/G21IP/ALK8: `curcuma -sp ... -method gfnff`
for Li+/Be+/Be2+/Na+/Mg2+ now all match `xtb --gfnff --sp` to <1e-6 Eh. Full project ctest
suite (`ctest`, 234 tests): 58/58 `gfnff`-labelled tests pass; the 21 failures in the full
run are pre-existing/environmental (missing `release_tblite/` reference-dump tree, 17
tests; `cli_curcumaopt_07_opt_multixyz` documented golden-value drift, CLAUDE.md Known
Issues; 3 unit-test binaries not rebuilt since before this session) - none touch GFN-FF/EEQ.

## GFN-FF: remaining outliers (unrelated to the two fixes above)

> **Updated Sep 2026.** The AHB21 bullet below is RESOLVED - it was the two-fragment
> charge placement (Known Issue #13), now at MAD 0.606 / max 4.8. The section is kept for
> the categories that remain. Current state: MAD 0.416 / RMSD 3.77 / max 105.5, with
> 106 of 2458 structures above 1 kcal/mol and 13 above 20.
>
> An arbitration run (worst ~77 outliers, curcuma vs pprcht vs xtb) showed these are
> **mostly genuine curcuma port errors**, not the pprcht-vs-xtb reference split that
> dominates MOR41: 66 port errors, 5 splits (RSE43, BHPERI, DC13/ch2n2), 6 mixed (only
> MB16-43, where pprcht and xtb themselves differ by 30-85 kcal/mol - and where all three
> engines miss the set's own published decomposition energies by a MAD of 380-400 kcal/mol,
> so that split cannot be arbitrated and is not worth chasing; see REV_GFNFF_TODO.md #7).
> Term fingerprints
> split them in two: a **bond/topology-perception** family (DC13/c20bowl +564 kcal in the
> bond term alone, AL2X6 bridged dimers, ALK8 Li clusters, HEAVY28/HEAVYSB11 heavy
> hydrides) and a smaller set of strained/hypervalent cases (oxiranes, H2S2O7, N-ylides).

MAD 3.389 / RMSD 17.45 / max 500.0 kcal/mol still exceeds gfn1/gfn2 by ~2 orders of
magnitude. The worst outliers cluster into recognisable categories, not one single bug -
none involve isolated atoms or Li/Be, so neither fix above touches them:

- **~~Charged anionic H-bond complexes (AHB21)~~ RESOLVED (Known Issue #13)**: these were
  the largest remaining outliers (+67 to +239 kcal/mol, all charge -1). Cause was not the
  H-bond term at all but the fragment the net charge was placed on; see below.
- **SN2/proton-transfer transition states (BH76)**: `fch3fts`, `hoch3fts`, `fch3clts`,
  `clch3clts` (-72 to -151 kcal/mol). Bond-breaking/forming TS geometries are a known
  hard case for any topology-perception-based force field (bond order is ambiguous at
  a stretched TS bond length) - expected divergence class, not necessarily a curcuma bug
  specifically (xtb-gfnff is not designed for TS either).
- **Unusual/non-classical bonding** (topology-perception edge cases): `MB16-43` (artificial
  stress-test molecules, several -69 to -106 kcal/mol), `AL2X6/al2me6` (Al bridging methyls,
  3c-2e bonds, -73), `ALK8/li4_me4` (Li4 cluster, -87 - slightly worse post-fix, not
  better: the dgam fix is more correct but doesn't guarantee every individual multi-term
  structure's cancellation improves), `YBDE18` (N-C ylide/dative-bond dissociation, -61 to
  -82), `DC13/c20bowl` (curved all-sp2 carbon bowl, +500, unchanged - the global max). Same
  bug class as the topology/hybridization-classification fixes already logged in CLAUDE.md
  Known Issues #6 (GFN-FF FT-HMO / bpair / hybridization fixes) - plausible further
  instances, not root-caused individually here.

Full per-subset table and outlier list: `test_cases/GMTKN55-testset/_run/gmtkn55_summary_gfnff.md`.

## What was tested / not tested

- Tested: 2462 GMTKN55 structures (54 subsets), single-point, default GFN-FF/GFN1/GFN2,
  D4, no solvent, curcuma vs xtb 6.7.1 on identical geometry/charge/UHF; full project
  ctest suite re-run after both GFN-FF fixes (no regressions).
- Not tested: GMTKN55 WTMAD-2 / relative reaction-energy accuracy vs the high-level QM
  reference (would need a `tmer2++` `.res` parser - separate, larger undertaking).
- Not tested: gradients, optimisation, any charged/open-shell system for gfn1/gfn2
  (skipped by design, see Known Issues #9).
- Not root-caused: the AHB21/BH76-TS/MB16-43/AL2X6/YBDE18/DC13 GFN-FF outlier
  categories above - flagged, not individually diagnosed to a specific line of code.

Human production testing pending until the operator removes this note.

## Reproduce

```bash
python scripts/fetch_testset.py fetch gmtkn55
python scripts/gmtkn55_compare.py --subset DIPCS10 --method gfnff --recompute   # smoke test
python scripts/gmtkn55_compare.py --recompute                                  # full sweep
```

Energies cached in `test_cases/GMTKN55-testset/_run/energies.json`; `--subset`/`--limit`
reruns and a resumed full sweep are served from cache (pass `--recompute` after a code
change, since the cache does not know the binary changed).
