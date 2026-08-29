# GFN-FF metal-reduced neighbour lists (Fortran four-list construction)

> 🤖 AI-generated, ⚙️ machine-tested. Human production testing pending.

Status: implemented July 2026. Ports the Fortran multi-list neighbour construction
(`external/gfnff/src/gfnff_ini2.f90:128-130, 197-202`) into native GFN-FF, replacing the
single full-connectivity adjacency list.

## Why

Native GFN-FF kept **one** bonded adjacency list including metal bonds. Fortran keeps
**four** and assigns hybridization from a metal-reduced, eta-aware mixture. Consequences of
the single-list approximation:

- Carbonyl O in Cr(CO)6 was assigned `hyb=2`; Fortran gives `hyb=1`. This shifted the C-O
  pibo and left PR01/PR02/PR05 at +0.5..0.6 kcal/mol.
- **Cyclopentadienyl and alkene carbons bonded to a metal were assigned `hyb=3` (sp3)**
  instead of sp2, because their metal bond inflated the coordination number. This was the
  larger error by far.
- `pi_fragments` was unusable in the angle `feta` correction, `btyp=6` (eta) was never
  assigned, and `itag` never reached the Hueckel solver.

## The four lists

All derive from the same bond list; they differ only in filtering.

| List | Field | Fortran | Filter |
|---|---|---|---|
| full | `nb_full` | `nbf`, getnb icase=1 | none |
| HC-reduced | `nb_hc` | `topo%nb` during the hyb loop, icase=2 | drops **all** bonds of atoms whose full CN exceeds `hc_crit` (4 if `periodic_group <= 2`, else 6) |
| metal-free | `nb_nometal` | `nbm`, icase=3 | drops metals (`mchar > 0.25` or `metal_type > 0`) and heavy atoms above `normcn` (`Z > 10`) |
| mixture | `adjacency_list` | `nbdum` -> final `topo%nb` | per atom: `nb_nometal` if eta-coordinated, else `nb_full` |

`periodic_group` is **negative** for the d-block, so every transition metal gets
`hc_crit = 4`. In Cr(CO)6 the Cr has CN=6 > 4, so all six Cr-C bonds vanish from `nb_hc`
**from both endpoints** — which drops each carbonyl C to `nb_hc == 1` and fires the CO->sp
rule on the oxygen. That is the entire mechanism behind the carbonyl fix.

### The ordering subtlety

Fortran spells three different lists `topo%nb` in this region. Inside the hybridization
loop `topo%nb` is still the **HC-filtered** list; it is only overwritten with the mixture
afterwards, at `gfnff_ini2.f90:335`. So in `determineHybridizationFortran`:

```
nb20i   = adjacency_list[i].size()                    // the mixture
nbdiff  = nb_full[i].size() - nb_hc[i].size()         // how much the HC filter removed
nbmdiff = nb_full[i].size() - nb_nometal[i].size()    // how much the metal filter removed
```

and the NO2 / B-N / N-SO2 rules, the CO->sp rule, and `nnNearestNoM` all index into
`nb_hc`. Feeding any of them the mixture silently changes the result.

## Consumer mapping

| List | Consumers |
|---|---|
| `nb_full` | `detectMolecularFragments` (Fortran `gfnff_ini.f90:468` fragments on `nbf`, so an eta ligand is never split off and `qfrag` stays correct); eta detection; `nbdiff`/`nbmdiff` |
| `nb_hc` | **only** inside `determineHybridizationFortran` |
| `nb_nometal` | eta detection; `findSmallestRings` (Fortran `gfnff_ini.f90:692` perceives rings on `nbm`) |
| `adjacency_list` (mixture) | everything else — pi systems, angles, torsions, inversions, dispersion ATM, BATM, EEQ, HB/XB; `neighbor_lists` is an alias |

## mchar

`metallic_character` ports `gfnff_ini.f90:249`:

```
mchar(i) = exp(-0.005 * en(Z_i)^8) * sum_j || d logCN_i / d R_j || / (cn(i) + 1)
```

Two things matter:

1. It uses Fortran `param%en` (added as `GFNFFParameters::gfnff_en`), **not** the `rab_en`
   table used by `computeRabEstimate`.
2. The CN cutoff is `cnthr = 100 - log10(accuracy)*50` = 100, and that is a **squared**
   threshold, so the radius is **10 Bohr** — not the 6 Bohr neighbour-list default. Using 6
   leaves mchar ~1e-4 low.

The N x N derivative tensor is not materialised: for `m != i` the norm reduces to
`|dlogdcn_i * derfCN(i,m)|`, and the self term collapses to one vector sum, so the whole
thing is one O(N*k) pass.

## Validation

**Against an instrumented Fortran build** (temporary per-atom print after
`gfnff_ini2.f90:130`, reached via `-method xtb-gfnff`), over all 95 MOR41 structures:

- `nb_full` counts: **94/95 exact**
- `nb_hc`, `nb_nometal` counts: **94/95 exact**
- `mchar`: **94/95 exact** to the 6-decimal print (the 95th differs by 1 ulp of the print format)

**Regression suites:**

- `ctest -R "^gfnff_val_"` — **18/18** (unchanged)
- `ctest -R gfnff` — **75/76**, the one failure being the pre-existing `gfnff_gpu_vs_cpu_h2`
  missing-fixture failure (`external/gfnff/test/h2.xyz` is absent from the repo)
- MOR41 gfnff vs xtb 6.6.1:

| metric | before | after |
|---|---|---|
| MAD | 57.317 | **7.295** |
| max | 457.32 | **37.80** |
| RMSD | 118.371 | **11.907** |
| MD (bias) | +46.870 | **-6.012** |
| within 1.0 kcal/mol | 33/95 | **39/95** |
| total \|error\| | 5445.2 | **693.1** |

Step by step: 57.32 (baseline) -> 16.68 (topology rewire) -> 9.15 (torsion enumeration)
-> 8.35 (getnb criterion) -> 7.30 (gated q-loop).

The seven previously-exact metal complexes (ED02, ED03, ED27, ED28, ED31, ED39, PR03) moved
by **0.000 kcal/mol**. The original carbonyl targets went to zero: PR01 +0.62 -> -0.00,
PR02 +0.51 -> -0.00, PR05 +0.51 -> -0.00.

Hybridization changes: 217 atoms across 48 structures. Dominant classes are C sp3->sp2
(118, Cp/alkene carbons — the eta fix) and O sp2->sp (73, carbonyls — the original target).

## Torsion enumeration

`generateTorsionsNative` and `generateSTorsionsNative` (`gfnff_torsions.cpp`) rebuilt their
own neighbour table from the raw `getCachedBondList()`, bypassing the topology adjacency
entirely. Fortran builds its torsion list from `topo%nb` (the nbdum mixture), where an
eta-coordinated carbon has its metal bond stripped - so enumerating on the full bond list
invents torsions **through the metal centre** that the reference never generates.

This was by far the largest remaining error after the topology rewire. Per-term decomposition
against `xtb --gfnff` over the 12 worst structures attributed **61% of all remaining error to
torsion** (e.g. ED29 Pt(C10H18): native torsion +0.2011 Eh vs xtb +0.0077 Eh, a 26x
over-count worth +121 kcal/mol). Both generators now consume `topo.adjacency_list`.

Effect: MOR41 MAD 16.679 -> 9.146, RMSD 32.7 -> 19.4, and the bias MD +10.85 -> -3.28.
Protected set unmoved (0.000), gfnff_val still 18/18.

Remaining error after this, by term (12 worst structures): bond ~27%, angle ~9%, everything
else <2% each. The bond share is dominated by ED07 and PR40, the two known deviations below;
the residual metal **bond** term (`btyp>=5`, `docs/GFNFF_METAL_BOND_ANALYSIS.md`) is still
unimplemented and is the main remaining lever.

## getnb criterion and the q-loop (Jul 2026)

Both originally deferred, both now done.

**getnb criterion.** The bond criterion is now the Fortran one
(`gfnff_ini2.f90:111-126` + `:361-419`):

```
rco  = gfnffrab(Z_i, Z_j, normcn_i, normcn_j)   // computeRabEstimate, R0 at "normal" CN
rco -= qa_i*f1 + qa_j*f2                        // charge shrink, fq=0.23 doubled for metals
rco *= fat(Z_i)*fat(Z_j)
bonded if r < fm * 1.25 * rco                   // fm: rthr2=1.00 for TM (no-op), 1.025 for main-group metals
```

replacing the older `1.3*(rcov_i+rcov_j)*fat_i*fat_j` heuristic. Native `nbf` (the getnb
`icase=1` list, `getCachedBondList()`) now reproduces Fortran on **95/95** MOR41 structures:
on ED07 it gives 68 bonds, exactly xtb's `#bonds: 68`, with the agostic C-H...W / C...W pairs
correctly excluded.

**Bond-term wiring fix (Jul 2026).** Getting `nbf` right was necessary but not sufficient: the
bond-term generators (`generateGFNFFBonds`, `generateBondsNative`) did **not** consume it. They
each re-derived the bond list from an independent `1.3*(rcov_i+rcov_j)*fat*fat` rule, so the
authoritative getnb list and the list actually feeding the bond energy could disagree. On ED07
they did (getnb 68, heuristic 70), reintroducing the two spurious metal contacts getnb had just
removed -- worth -0.0704 Eh. Both generators now enumerate `getCachedBondList()` directly.
Effect: ED07 -37.80 -> +6.34 kcal/mol; set MAD 7.295 -> 7.153, max 37.80 -> 35.08.

**q-loop.** `calculateTopologyInfo()` now runs `calculateTopologyInfoOnce()` up to twice,
mirroring `gfnff_ini.f90:258-263`. The charges only reach the model *through* the bond list,
so if the shrink changes no bond, pass 2 is a mathematical no-op - the second pass is
therefore gated on the bond list actually changing. Measured against Fortran: that happens on
exactly **1 of 95** MOR41 structures (PR40). The gate is also a safeguard against the
re-entrancy bug listed below.

Effect: ED07 +113.83 -> **-37.80**, PR40 +103.63 -> **-3.83**; MOR41 MAD 9.146 -> **7.295**,
max 113.83 -> **37.80**, RMSD 19.4 -> **11.9**. Protected set still 0.000, gfnff_val 18/18.

**ED07 was a neighbour-list problem after all** (corrected Jul 2026). The earlier reading --
"the -37.8 kcal is the bond term, i.e. the unimplemented `btyp>=5` branch" -- was wrong on both
counts. The `btyp>=5` metal branch is implemented (`53d6aeb`/`14fc648`), and per-bond comparison
(`scripts/s30l_bond_compare.py --xyz`) showed the residual was two **extra bonds** (W...C 2.92 A,
W...H 2.23 A) that xtb does not have, not a per-bond parameter error -- the bond term was reading
a 70-bond list while getnb said 68. Fixed by the bond-term wiring change above (ED07 -> +6.34).

The fix exposed a genuine metal bond-**parameter** error, previously masked: **PR07** (Kubas
ED07+H2) transiently regressed to -27.05 because its W-H bonds now used the correct list but got
`kbond -0.0585` where xtb has `-0.032` (1.83x). Root cause (fixed the same session): hydrogen
(group 1) was wrongly excluded from `mtyp=1`, mis-setting `fqq`/`fcn` on metal-H bonds, and the
`fsrb2` EN-scaling was keyed on `mtyp>0` rather than an actual metal bond. Fixing both -- against
the Fortran analyzer's per-bond fqq=1.0964 / fcn=0.25842 -- brought PR07 to +3.53 and a whole
class of hydride complexes (PR09, PR08, ED17, PR34/35/33/17, PR06, ED14, PR10) to within a few
kcal, with organic ligands byte-identical. Net MOR41 gfnff MAD 7.295 -> 5.288, within-1 39 ->
45/95, zero structures worse.

## What was NOT done / known deviations

- ~~The `getnb` distance criterion is not ported~~ - **DONE (Jul 2026)**, see "getnb criterion
  and the q-loop" above.
- ~~The two-pass q-loop is not ported~~ - **DONE (Jul 2026)**, gated on the bond list
  actually changing; see above.
- **`calculateTopologyInfoOnce()` is not re-entrant.** Calling it twice on unchanged input
  converges to a slightly different Coulomb energy (H2O: -0.0890528 -> -0.0881507 Eh) and
  then stays at that fixed point, so some one-shot state transitions on the first call. It is
  NOT accumulation (a 3rd and 4th call change nothing further) and it is NOT the EEQ Cholesky
  or A_nn cache (invalidating both changes the wrong value to a different wrong value); a
  persistent flag such as `m_phase2_historically_implausible` (`eeq_solver.h:958`) is the
  most likely candidate. This is a **pre-existing latent bug**, exposed rather than caused by
  the q-loop work. The q-loop sidesteps it by only re-running when the bond list changes, but
  it should be fixed on its own merits - any future code that rebuilds topology twice will
  hit it.
- **Not tested**: periodic systems, lanthanides/actinides, radicals, anything outside MOR41
  and the 18-molecule neutral main-group validation set. The validation set contains **no S
  and no P**, so the hypervalent rules (`nbdiff==0 && Z>10 -> hyb=5`) and the `normcn`
  hypervalent filter are exercised only incidentally via MOR41.
- **Unrelated open gap**: the GFN-FF transition-metal **bond** term (unimplemented Fortran
  `btyp>=5` branch, native metal bonds ~1.95x too strong) is untouched and remains the
  dominant contributor to the residual MAD — see `docs/GFNFF_METAL_BOND_ANALYSIS.md`.

## Debugging

Both are env-gated and free when unset:

```bash
CURCUMA_NBDIAG=1  curcuma -sp mol.xyz -method gfnff   # per-atom nbf/nbhc/nbm/nbdum + mchar
CURCUMA_HYBDIFF=1 curcuma -sp mol.xyz -method gfnff   # atoms where legacy and Fortran hyb disagree
```

Set `m_use_fortran_hyb = false` (`gfnff.h`) to fall back to the pre-July-2026 heuristic for
bisection.

**Cache note**: the topology cache format was bumped to **version 2**. A v1 `*.topo.json`
holds a hybridization from the old heuristic and an adjacency that is not the mixture, so it
is now rejected and regenerated.
