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

### 2. All shifts mis-scaled by ff (generalized — fixes ligand C-O under coordination too)
Fortran's `r0 = rtmp + vbond(1)` adds the **entire** shift `vbond(1) = rabshift + shift`
**outside** the EN factor ff, where `shift` includes the pi-bond-order shift
(`gfnff_ini.f90:1174` `shift = hueckelp*(bzref-pibo)`) — Fortran does **not** fold pi into
`rtmp`. curcuma's `rtmp = (ra+rb)*ff` is bit-equal to Fortran's `rtmp` for **every** bond
checked (metal-C, and the ligand C-O: 2.13519 vs 2.13516 Bohr), and curcuma's `pi_shift` uses
the identical formula/constants (`hueckelp=0.340`, `bzref=0.370`). curcuma previously multiplied
**all** shifts by ff (`(ra+rb+shift)*ff`), which only matches when `ff~1` (free light-atom
bonds); it mis-scaled by `(1-ff)` whenever ff departs from 1 — the ~-0.41 metal-C shift
(`ff~0.983`, ~0.6 kcal) **and** the coordinated ligand C-O (`ff~0.949`, ~0.1 kcal/CO; free CO is
exact because its ff~1). **Fix (all bonds):** `r0 = (ra+rb)*ff + total_rabshift`
(`total_rabshift = rabshift + pi_shift + metal_shift`, everything outside ff). Neutral bonds
stay exact (their ff~1, so inside/outside ff agree to <1e-5 kcal); the coordinated carbonyls
now match — Fe(CO)₄, Ni(CO)₃, RhCl(CO)₂ exact. (An earlier metal-only version of this fix left
the organic path all-inside-ff; generalizing it fixed the ligand C-O under coordination.)

## Result (native vs xtb 6.6.1 `--gfnff`, all with fresh topology caches)

Per-bond parameters now match the Fortran per-bond reference: Cr(CO)₆ Cr-C r0=3.0585 Bohr,
k_b=-0.06837, alpha=0.4472; the ligand C-O r0 0.9679→**0.9592 Å** (Fortran 0.9590).

| structure | before (HEAD 31d3a58) | after | note |
|-----------|----------------------:|------:|------|
| ED02 Fe(CO)₄ | — | **0.00** | exact |
| ED03 Ni(CO)₃ | — | **0.00** | exact (C-O fix) |
| ED39 RhCl(CO)₂ | — | **0.00** | exact |
| PR03 | — | **0.00** | exact (C-O fix) |
| PR01 Cr(CO)₆ | -13.0 kcal | **+0.62** | residual = coordinated-Hückel pibo |
| PR02 | — | +0.51 | |
| PR05 Mn(CO)₅H | -8.5 kcal | **+0.51** | |

Neutral regression: `ctest -R gfnff_val` 18/18; CO/C2H4/C2H6/C3H8/MeOH/Bz/COD/PhOH/MeCN all
exact (<1e-4 kcal — their bond ff~1, so inside/outside-ff agree). MOR41 gfnff reaction MAD vs
xtb ~74 (unchanged; the set is dominated by **other** unfixed native-gfnff-on-TM failures — many
TM structures are still 100s of kcal off, e.g. ED04 438 kcal — a separate and far larger
deficiency than the bond term).

> **Measurement note:** curcuma writes a per-structure `mol.topo.json` topology cache (bond
> parameters incl. r0). It is keyed on geometry, **not** the code version, so a cache written by
> an older binary is served verbatim and masks parameter changes. Delete `*.topo.json` before
> comparing builds.

## Remaining (separate, smaller)
- **Coordinated-Hückel pibo** (PR01 0.62 / PR05 0.51 kcal): the r0 **formula** is now exact, but
  the C-O `pibo` fed into it differs (native 0.9952 vs Fortran ~0.9961) because the Hückel
  pi-bond-order solve in a metal complex diverges slightly (free CO is exact). Propagates through
  the identical `hueckelp*(bzref-pibo)` shift to ~0.0003 Bohr on C-O r0. A separate Hückel issue
  (same family as the S30L torsion-pibo residuals), not the bond r0 model.
- **CO2** (0.20 kcal, pre-existing, no metal): independent organic-model outlier, unchanged by
  this fix.
- Native gfnff remains **not TM-ready** overall (many MOR41 structures 100s of kcal off from
  other unimplemented/incorrect TM terms). This fix closes only the bond-stretch r0.
