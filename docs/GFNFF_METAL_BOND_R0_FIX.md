# GFN-FF Metal Bond Equilibrium Distance (r0) Fix (Jul 2026)

Status: 🤖 AI-generated, ⚙️ machine-tested. Human production testing pending.

Follow-up to [GFNFF_METAL_BOND_ANALYSIS.md](GFNFF_METAL_BOND_ANALYSIS.md) (force-constant
branch, committed 53d6aeb) and [GFNFF_ETA_ANGLE_FIX.md](GFNFF_ETA_ANGLE_FIX.md) (angle,
committed 31d3a58). After those, the metal-carbonyl M-C **force constant** was correct but the
M-C **equilibrium distance r0** was wrong, over-binding Cr(CO)₆/Mn(CO)₅H by ~13/8 kcal.

## Two root causes in the bond r0

The bond r0 is `r0 = (ra+rb+shift)*ff` in curcuma; Fortran (`gfnff_ini.f90:1289`) is
`r0 = rtmp + vbond(1)` with `rtmp=(ra+rb)*ff` and the shift `vbond(1)` added **outside** ff.

### 1. Non-metal hyb/bbtyp shifts wrongly applied to metal bonds
The `[GFNFF_ETA_ANGLE_FIX]` metal-hyb rule made octahedral TM centers `hyb=0`. That exposed a
latent bug: curcuma's base-shift block (`getGFNFFBondParameters`, the bbtyp/F-F/sp3/sp
corrections) ran for **all** bonds, so the `(hyb1,hyb2)==(0,1)` "X-sp correction" `+0.14`
fired on the M-CO bond (metal hyb=0, carbonyl C hyb=1) — inflating r0 and over-binding. In
Fortran these shifts live in the **non-metal branch** (`gfnff_ini.f90:1129-1206`), which metal
bonds never enter (they take the mutually-exclusive metal branch). **Fix:** gate the
bbtyp/F-F/sp3/sp base shifts on `!bond_has_metal`. The gen_rabshift and the heavy-heavy shift
(both applied after the if/else in Fortran, line 1268/1276) stay for all bonds.

### 2. Metal shift mis-scaled by ff
For a metal bond curcuma's `rtmp=(ra+rb)*ff` is **bit-equal** to Fortran's `rtmp` (a metal bond
has `pibo=-99`, i.e. no Goedecker pi-shortening and no SK correction). Fortran then adds
`vbond(1)=rabshift+metal_shift` **outside** ff, but curcuma multiplied the whole shift by ff,
mis-scaling the large ~-0.41 metal shift by `(1-ff)~1.7%` → ~0.6 kcal on Cr(CO)₆. **Fix:** for
metal bonds compute `r0 = (ra+rb)*ff + (rabshift+metal_shift)`. Gated to metal bonds: the
**organic** path must keep the all-inside-ff form, because there curcuma's `rtmp != Fortran's
rtmp` (the Goedecker `rabd` folds pi/SK corrections that curcuma reconstructs via `pi_shift`
inside ff); moving the shift outside ff regresses organics (CO2 by 0.2 kcal — verified).

## Result (native vs xtb 6.6.1 `--gfnff`)

The M-C bond parameters now match the Fortran per-bond reference exactly. Cr(CO)₆ Cr-C:
r0 = 3.0585 Bohr, k_b = -0.06837, alpha = 0.4472 — all bit-equal to the analyzer.

| structure | before (HEAD 31d3a58) | after | note |
|-----------|----------------------:|------:|------|
| PR01 Cr(CO)₆ | -13.005 kcal | **+0.616** | residual now the ligand C-O (coordination) |
| PR05 Mn(CO)₅H | -8.5 kcal | **+0.090** | |
| PR02 / PR03 | (large) | -0.666 / -0.524 | |
| ED02 Fe(CO)₄ | — | **-0.000** | exact |
| ED03 Ni(CO)₃ | — | 0.202 | |
| ED39 RhCl(CO)₂ | — | -0.000 | exact |

Neutral regression: `ctest -R gfnff_val` 18/18; C2H4/Bz bit-identical; the fixes are metal-bond
gated. MOR41 gfnff reaction MAD vs xtb 74.2 → 73.4 (slight; the set is dominated by **other**
unfixed native-gfnff-on-TM failures — many TM structures are still 100s of kcal off, a separate
and larger deficiency than the bond term).

## Remaining (separate, smaller)
- **Ligand C-O bond under coordination** (~0.1 kcal/CO, the PR01 0.616 / ED03 0.202 residual):
  free CO is exact, but a metal-coordinated CO's internal C-O r0 differs by ~0.009 Å — curcuma's
  `pi_shift`-inside-ff reconstruction of Fortran's Goedecker `rabd` diverges when the ligand's
  ff shifts under coordination. Same family as the **pre-existing** CO2 0.2 kcal (organic-bond
  model), independent of this fix. Reworking the organic bond-distance model is out of scope.
- Native gfnff remains **not TM-ready** overall (many MOR41 structures 100s of kcal off from
  other unimplemented/incorrect TM terms). This fix closes only the metal **bond-stretch** term.
