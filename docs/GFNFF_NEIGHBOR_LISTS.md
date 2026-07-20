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
| MAD | 57.317 | **9.146** |
| max | 457.32 | **113.83** |
| RMSD | 118.371 | **19.391** |
| MD (bias) | +46.870 | **-3.284** |
| within 1.0 kcal/mol | 33/95 | **39/95** |

(Intermediate: MAD 16.679 after the topology rewire alone; 9.146 after also pointing the
torsion generator at the topology adjacency - see "Torsion enumeration" below.)

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

## What was NOT done / known deviations

- **The `getnb` distance criterion is NOT ported.** Native keeps its own geometric bond
  criterion as `nbf`. This is deliberate: `gen%rthr2 = 1.00` (`gfnff_param.f90:720`) makes
  Fortran's metal radius enlargement an exact no-op for transition metals, and native `nbf`
  already reproduces Fortran on 94/95 MOR41 structures.
  - **ED07** is the exception: native bonds an agostic C-H...W contact that Fortran does not
    (`nbf(W)` 7 vs 5). `nb_hc` and `nb_nometal` absorb this completely, but `nbdum` and
    `nbdiff` do not. ED07 regressed -21.11 -> +113.83 kcal/mol. Porting the
    `gfnffrab`/`normcn`/charge-shrink criterion is the fix, tracked as future work.
- **The two-pass q-loop is NOT ported.** `gfnff_ini.f90:258` sets `qa = 0` before the loop,
  so pass 1 runs unshrunk and native corresponds to pass 1. Pass 1 differs from pass 2 on
  exactly **1 of 95** MOR41 structures: **PR40**, which regressed +61.68 -> +103.63 kcal/mol.
  Implementing the q-loop requires the `getnb` criterion first (a charge-shrink correction to
  a radius formula that never had one would be worse than nothing), so it is blocked on the
  item above.
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
