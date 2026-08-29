# PR40 GFN2-D4 Dispersion Residual — Root Cause & Fix (2026-07)

## Summary

The last non-exact MOR41 structure for native GFN2 vs the tblite reference, **PR40**
(Ti/Al/Cl/C/H), over-bound the D4 dispersion by **+1.98e-3 Eh**. Root cause found and
fixed; PR40 now matches tblite to **~5e-9 Eh** with **zero regression** on any control.

## Root cause (one line)

`src/core/energy_calculators/dispersion/d4_corrections_data.cpp` — the `d4_sscale_data`
map was **missing the entry for reference system Na (`refsys=11`)** (and Ne, `refsys=10`),
so the D4 reference-polarizability inner-shell subtraction silently vanished for every
reference whose *reference system* is Na, making those reference C6 values ~2.4x too large.

## The mechanism (verified end-to-end, not inferred)

Native GFN2 builds the D4 reference C6 from Casimir-Polder integration of the reference
polarizabilities, applying dftd4's `set_refalpha_gfn2` correction
(`reference.f90:343-377`):

```
aiw          = sscale(refsys)·secaiw(:,refsys)·zeta(ga, hardness(refsys)·gc, ...)
alpha(:,ir)  = max( ascale·(alphaiw(:,ir) − hcount·aiw), 0 )
```

Native's `computeC6Reference` (`d4param_generator.cpp:1020`) looks up `sscale` by the
reference-system atomic number:
```cpp
double sscale_i = (d4_sscale_data.find(refsys_i) != end()) ? d4_sscale_data.at(refsys_i) : 0.0;
```
`d4_sscale_data` had keys `{1,2,6,7,8,9,17}` but **not 10 or 11**. dftd4's table has
`sscale(10)=sscale(11)=0.5` (`reference.inc`). For any reference with `refsys=11` the
lookup fell through to `0.0`, so `hcount·sscale·secaiw = 0` → **no inner-shell
subtraction** → the reference polarizability stayed at its raw (too-large) value → C6
too large.

`refsys=11` (Na) is used by exactly **one reference (the 4th, highest-CN) of Sc, Ti, V,
Zr, Nb, Hf, Ta** (`reference.inc`: `refsys(4, {21,22,23,40,41,72,73}) = 11`). No other
element is affected. That is why the bug is geometry/coordination-specific: it only bites
when one of those metals has a CN high enough to weight its 4th reference.

### Why PR40 specifically

PR40's Ti has D4 CN = **3.9113** (bit-identical to tblite — CN was never the issue). Its
reference weights `[0, 0, 0.586, 0.449]` (also bit-identical to tblite's `gwvec`) split
across Ti references 2 (cn 3.43) and **3 (cn 4.44, the refsys=Na one)**. Native's Ti-ref3
C6 was ~2.4x too large, e.g. c6(Ti_ref3, C_ref0) native **155.79** vs dftd4 **64.46**
(refs 0/1/2 matched to 6 digits). ED40a/ED40b also contain Ti but at lower CN that never
weights ref3 — which is why they were already exact and stay exact.

### Evidence chain (all built & run, no guessing)

1. Native D4 CN == tblite D4 CN to ~1e-6 (`dump_dftd4_atomic_c6`), so CN was fine.
2. Feeding native-EEQ vs tblite-Mulliken charges through tblite's own C6 machinery moved
   the 2-body energy by only 6e-5 Eh — charge source was not the cause.
3. 2-body/ATM split: native 2body = −0.0171409, ATM = +0.0001495; tblite 2body ≈ −0.01497
   (ATM ≈ −4.6e-5). The residual is entirely the 2-body C6.
4. Per-pair weighted C6 (native `weightedC6Gfn2`/`contractC6Gfn2` vs tblite
   `get_atomic_c6`): every Ti pair ~1.6x too large; all non-Ti pairs exact.
5. `buildAtomRefW` weights for Ti = `[0,0,0.58594,0.44910]` == tblite `gwvec` → weights
   correct.
6. Per-reference c6(Ti_a, C_b): refs 0/1/2 exact, **ref 3 ~2.4x high** → isolated to the
   reference polarizability of Ti ref3.
7. Ti ref3 metadata (`hcount=4, refsys=11, ascale=1, refcovcn=4.4385, refcn=6.6375,
   refq=1.4528`) all match dftd4 — only distinction from refs 0/1/2 is `refsys=11` (Na).
8. Native `d4_sscale_data` had no key 11 → `sscale=0` → subtraction skipped. dftd4 has
   `sscale(11)=0.5`.

## The fix

`d4_corrections_data.cpp`, in `initialize_d4_corrections()`:
```cpp
d4_sscale_data[10] = 0.50000000000000;   // Ne (refsys=10, unused but matches dftd4)
d4_sscale_data[11] = 0.50000000000000;   // Na (refsys=11) — the missing entry
```
Pure data addition; the correction/gradient machinery is unchanged. It only alters the
4th reference C6 of Sc/Ti/V/Zr/Nb/Hf/Ta (the only elements with `refsys=10/11`).

## Before / after (native GFN2 total energy vs tblite-gfn2, `release_tblite/curcuma`)

| Structure          | before (native) | after (native) | tblite         | after |Δ| |
|--------------------|-----------------|----------------|----------------|-----------|
| **PR40** (Ti)      | −25.07640018    | **−25.07442271** | −25.0744227052 | **4.8e-9** |
| PR40 Dispersion    | −0.01699139     | **−0.01501253**  | (xtb −0.0150126) | ~1e-8 |

## Controls — all still exact (native GFN2 vs tblite-gfn2), unchanged by the fix

| Structure       | native        | tblite          | |Δ|     |
|-----------------|---------------|-----------------|---------|
| ED40a (Ti)      | −11.90961224  | −11.9096122353  | 4.7e-9  |
| ED40b (Ti)      | −13.00089993  | −13.0008999271  | 2.9e-9  |
| ED02 (Fe)       | −27.85077571  | −27.8507757107  | 7.0e-10 |
| ED03 (Ni)       | −23.34102052  | −23.3410205149  | 5.1e-9  |
| PR22 (Ir)       | −73.05222509  | −73.0522250907  | 7.0e-10 |
| I2              | −7.65233253   | −7.6523325294   | 6.0e-10 |
| PR07            | −87.20777968  | −87.2077796562  | 2.4e-8  |
| PR09            | −107.07956359 | −107.0795635840 | 6.0e-9  |
| ED31            | −71.83128795  | −71.8312879433  | 6.7e-9  |
| H2O             | −5.07036982   | −5.0703698219   | 1.9e-9  |
| CH4             | −4.17507458   | −4.1750745757   | 4.3e-9  |
| C6H6            | −15.87853118  | −15.8785311796  | 4.0e-10 |

(All diffs are at the 8-dp print-precision floor; Fe/Ni/Ir are elements 26/28/77, none in
the affected `{Sc,Ti,V,Zr,Nb,Hf,Ta}` set, so their numbers are provably untouched.)

## Regression suite (release_tblite, `ctest`)

- `gfnff_val_*` — **18/18 pass** (incl. caffeine, complex, polymer, triose): the shared
  D4 generator is undisturbed for GFN-FF (its organic validation set has no
  Sc/Ti/V/Zr/Nb/Hf/Ta references). GFN-FF H2O = −0.32726266 Eh (unchanged).
- `d4_dedq` — pass (D4 q-response ∂E/∂q · ∂q/∂x chain rule intact).
- `xtb_gradient_{H2O,CH4,NH3}` — pass (native GFN2 analytic gradient).

## Caveats

- Machine-tested only (automated `ctest`). Not human production-tested.
- The fix is validated for Ti via PR40/ED40a/ED40b. The same missing entry also affected
  **Sc, V, Zr, Nb, Hf, Ta** (4th reference); those elements are now corrected toward the
  dftd4 reference but were not independently exercised by a test structure here.
- The `sscale[10]` (Ne) entry is dead in practice (no element uses `refsys=10`); it is
  added only for completeness with dftd4's table.
