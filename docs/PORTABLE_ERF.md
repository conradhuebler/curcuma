# Portable erf/acos/exp/log — Wine vs. native-Windows CRT divergence

> ⚠️ AI-implemented, machine-tested only (build succeeds, ctest suite unchanged
> with the flag OFF, vendored functions verified to ~1e-16 vs `std::` and
> against a second independent fdlibm source). **Not yet human production
> tested** — in particular, the original triggering report (a different bond
> term for the MOR-testset structure PR38 between a real Windows run and a
> Wine run of the same `curcuma.exe`) has not been re-confirmed fixed on a
> real Windows machine, because none is available to this project.

## Symptom

The same Windows binary (`curcuma.exe`, built once by GitHub Actions/MinGW)
produced a different GFN-FF bond term for the same input geometry when run
natively on Windows versus under Wine on Linux (used as a Windows-testing
proxy, since this project has no Windows machine).

## Root cause

`curcuma.exe` dynamically imports `erf`, `acos`, `exp`, `log`, `pow`, `sin`,
`cos`, `sqrt`, `round`, ... from `api-ms-win-crt-math-l1-1-0.dll` (the
Universal CRT math API set), confirmed via `objdump -p curcuma.exe`. Wine
cannot legally redistribute Microsoft's `ucrtbase.dll`, so it ships its own
reimplementation of these functions — which is not bit-identical to
Microsoft's for transcendental functions (a documented, long-standing Wine
characteristic; different rational/polynomial approximations round
differently in the last bit, same "table maker's dilemma" the C standard
never mandates a single answer for).

GFN-FF/EEQ build several **hard-coded classification thresholds** directly on
these functions' output:

- `erf()` — the CN (coordination number) sum is `0.5*(1+erf(kn*dr))` summed
  over neighbours; CN is compared against bond/hybridisation/pi-membership
  cutoffs throughout the topology builder.
- `acos()` — bond angles (`acos(cos θ)*180/π`) are compared against
  fixed-degree thresholds: 150° (carbene detection), 170° (linear
  N/dihedral), ~344-355° (planarity → sp2 vs sp3), to pick hybridisation,
  bond type, or HB donor/acceptor class.
- `exp()`/`log()` — CN's log-compression formula
  `log(1+e^cnmax) - log(1+e^(cnmax-cn_raw))` is later `std::round()`-ed to an
  integer that gates further hybridisation branches. `round()` itself is
  exact (IEEE 754 requires basic operations like round/floor/ceil to be
  correctly rounded for every finite input — no conformant implementation,
  Wine's included, can legally disagree there), so the divergence risk is
  entirely in `exp()`/`log()` feeding it, not in `round()` itself.

A single ULP difference in any of these can flip a discrete branch and
produce a genuinely different, qualitatively-different bond term for the same
molecule — not just a tiny energy difference.

## Fix

Four self-contained, statically-compiled functions ported from **fdlibm**
(the classic Sun Microsystems math library, permissively licensed —
`Copyright (C) 1993/2004 Sun Microsystems, "Permission to use, copy, modify,
and distribute this software is freely granted, provided that this notice is
preserved"`), via musl libc's `src/math/{erf,acos,exp,log}.c`, transcribed and
cross-verified against a second, independent source (FreeBSD
`lib/msun/src/e_{acos,exp,log}.c` / `s_erf.c`) — every numeric coefficient
matches bit-for-bit across both sources:

- `src/core/portable_erf.h` — `CurcumaMath::portable_erf()`
- `src/core/portable_acos.h` — `CurcumaMath::portable_acos()`
- `src/core/portable_exp.h` — `CurcumaMath::portable_exp()`
- `src/core/portable_log.h` — `CurcumaMath::portable_log()`

These compile straight into the binary — no DLL import, no dependency on the
OS/Wine's CRT, so they execute identical machine code regardless of platform.

`src/core/math_compat.h` is the single dispatch point: `curcuma_erf()`,
`curcuma_acos()`, `curcuma_exp()`, `curcuma_log()` call the vendored
implementation when `CURCUMA_PORTABLE_MATH` is defined, otherwise fall back to
`std::erf/acos/exp/log` (today's default behaviour, byte-for-byte).

### Build flag

```
cmake -DUSE_PORTABLE_MATH=ON ...
```

Default **OFF** — keeps existing validated golden-value ctests bit-for-bit
unchanged. `scripts/build_windows.bat` turns it **ON**, since the Windows
nightly is where the divergence was observed.

### Call sites converted

Every location identified (by a combination of an automated multi-agent audit
and manual grep+read follow-up, covering all files under
`src/core/energy_calculators/ff_methods/` and `.../dispersion/`) where one of
these functions' output feeds a discrete classification decision:

| Function | Files | Sites |
|---|---|---|
| `erf` | `gfnff_method.cpp`, `eeq_solver.cpp`, `ff_workspace_gfnff.cpp`, `forcefieldthread.cpp`, `cn_calculator.cpp`, `d4_charge_model.cpp`, `d4_ncoord.cpp` | 22 |
| `acos` | `gfnff_method.cpp` (5), `eeq_solver.cpp` (2, one is a shared helper used by 3 more hybridisation branches), `gfnff_torsions.cpp` (1) | 8 direct + 3 via shared helper |
| `exp`+`log` (CN log-compression only) | `cn_calculator.cpp`, `gfnff_method.cpp`, `d4_charge_model.cpp`, `xtb_vulkan_context.cpp` (Vulkan D4-EEQ host code, same formula) | 6 formula instances |

**Not converted** (audited, confirmed safe): the D3-style exponential CN
counting function (`1/(1+exp(arg))`, `cn_calculator.cpp`), all `exp`/`pow`
terms in continuous bond/angle/torsion/dispersion/HB/XB energy and gradient
formulas across `ff_workspace_gfnff.cpp`, `ff_workspace_uff.cpp`,
`forcefieldthread.cpp`, `gfnff_advanced.cpp`, `gfnff_inversions.cpp`,
`huckel_solver.cpp`, `d3param_generator.cpp` — a ULP difference there only
shifts an energy/gradient value by ~1e-16 relative, which is harmless. Also
not converted: two `log10()`-based neighbour-list cutoffs
(`gfnff_method.cpp` HB thresholds, `cn_calculator.cpp` CN cutoff) — these take
a **fixed accuracy constant** as input, not molecule geometry, so any
divergence is a deterministic-per-build broad-phase distance-filter shift far
below the noise floor of any realistic geometry, not a per-molecule risk.
`round`/`lround`/`llround` are not vendored anywhere (see Root cause above —
provably not a divergence risk given a deterministic input, which the erf/
acos/exp/log fix now provides).

## Verification

- Each vendored function checked against `std::` over a dense sweep plus
  explicit boundary/edge cases (domain limits, NaN, ±Inf, subnormals):
  `erf` max abs error 1.11e-16 (600k samples, [-30,30]); `acos` max abs error
  4.44e-16 ([-1,1]); `exp` max relative error 2.2e-16 ([-50,50]); `log` max
  relative error 2.22e-16 ((0,100]). The full GFN-FF CN log-compression
  formula, exercised end-to-end with the vendored `exp`+`log` together,
  differs from the `std::`-based formula by ≤8.9e-16.
- `ctest -R "gfnff|eeq" -E gpu` (60-73 tests, CUDA-linking failures on this
  machine excluded as pre-existing/unrelated): identical pass/fail result
  with the flag OFF before and after every change in this work.
- With the flag ON, the same suite passes; three representative molecules
  (`complex` — 231 atoms, the largest in the reference set; `caffeine`;
  `triose`) give bit-identical GFN-FF single-point energies at 8 decimal
  places between the OFF and ON builds.

## What this does *not* prove

- No real Windows machine was available to confirm the original PR38
  Wine-vs-native divergence is actually closed by this fix — only that the
  identified erf/acos/exp/log call sites are now architecturally immune to
  CRT-DLL-level divergence (deterministic machine code on any platform).
- The systematic audit for *other* CRT-math-fed threshold sites covered
  `src/core/energy_calculators/ff_methods/` and `.../dispersion/` exhaustively
  (all 23 CRT-imported symbol names grepped, every hit's context inspected);
  it did **not** extend to the native xTB (GFN1/GFN2) SCF code path or other
  parts of the codebase outside GFN-FF/EEQ/D4 — those were out of scope for
  the reported bug (a GFN-FF bond term) and have not been checked for the
  same failure pattern.
