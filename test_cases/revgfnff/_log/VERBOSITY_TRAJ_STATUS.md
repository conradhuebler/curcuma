# VERBOSITY_TRAJ_STATUS — why `-verbosity 4` changed a react trajectory (RESOLVED)

Closes the "Not determined" item in `BREAK_TAIL_STATUS.md`. Binary before **79bb76bb6720278c9da8b35007ed50ee**,
after **e7b6c62ace154d815a244645c29c3798** (`build_rev/curcuma`; `release/` untouched).
Harness `<scratch>/bt/run_one.sh`, `NOFILT=1`.

## 1. Reproduction — c2h6 / T=2000 K / frame 16, seed 42, 5000 fs

| run verbosity | calculator verbosity | rebuilds | max \|dE_jump\| |
|---|---|---:|---:|
| 1 / 2 / 3 | 0 / 1 / 2 | **90** | **+0.179425 Eh = +471.08 kJ/mol** |
| 4 | 3 | **136** | **+0.000212 Eh = +0.56 kJ/mol** |

Exactly the recorded numbers. ch4_H/T2000/f10 unaffected (18, +310.4 kJ/mol at both).
Count rebuilds as `max(REACT rebuild #N)`, not `grep -c`: from calculator verbosity 2 on,
GFN-FF prints the line itself **and** SimpleMD flushes it, so a naive count doubles.

## 2. The switch is the CALCULATOR level, 2 -> 3

A temporary cap on the thread verbosity `EnergyCalculator::CalculateEnergy` installs
(`CURCUMA_VCAP`, reverted): `-verbosity 4` with cap 0/1/2 -> 90, cap 3 -> 136. SimpleMD's own
verbosity-4 output is irrelevant; the FF at level 3 is not.

## 3. Divergence point

High-precision per-call trace (temporary `CURCUMA_EDUMP`, reverted): bit-identical for **4776
energy calls**, first difference at **call 4777 = step 4777, t = 1194.25 fs**, **1 ulp**
(`-1.00714392006672893` vs `...915` Eh), same `m_react_calls`, same bond count. A per-term dump
puts all of it in the **dispersion** term; a dispersion-parameter dump narrows it to the **C6
sum** (CN, EEQ charges, `r0_squared`, `zetac6`, `r4r2ij` all bit-identical). It first appears
just after the first react rebuild, where the C6 half-contraction table is rebuilt. Printed
6-decimal Epot rows then diverge at step 8767, the topology-event sequence at event ~10495.

## 4. Cause

`src/core/energy_calculators/dispersion/d4param_generator.cpp:1414`
(`D4ParameterGenerator::getChargeWeightedC6`):

```cpp
if (m_c6_half_valid && CurcumaLogger::get_verbosity() < 3
    && static_cast<int>(atom_i) < m_c6_half_natoms) {   // Lever-3 half-contraction fast path
```

Its own comment states the hazard — the fast path *"reassociates the FP sum vs the flat loop
below (~1e-16)"* — and it was skipped at verbosity >= 3 **so the `C6_DEBUG` print in the flat
path would still fire**. That makes a *numerical* path selected by the print level: the two
summation orders differ in the last ulp, so C6, the dispersion energy and the forces depend on
`-verbosity`. Category (b) of the brief, not a print side effect: disabling *every* `>= 3`
verbosity gate in `gfnff_method.cpp` (80), `eeq_solver.cpp` (39), `d4param_generator.cpp` (25)
and 9 further FF files at once did **not** remove the effect — this gate is an inverted `< 3`.
`grep -rn "get_verbosity() < \|m_verbosity < " src/` shows it is the only numerical gate of its
kind; every other hit is output-only.

## 5. Fix — `d4param_generator.cpp`, 18+/5-

Fast path taken unconditionally; the `C6_DEBUG` fall-through is opened by **`CURCUMA_C6DEBUG`**
instead (the `CURCUMA_BONDDUMP` / `CURCUMA_HUCKELDUMP` convention), read once into a
`static const bool` that also drives `log_c6`. Marked "Claude Generated (Sep 2026)" with the
measured numbers in the comment. Default verbosity (< 3) already took the fast path, so **no
default-verbosity result changes**.

## 6. Re-verification (e7b6c62a)

* c2h6/T2000/f16 at `-verbosity 1, 2, 3, 4`: **90 rebuilds, +0.179425 Eh (+471.08 kJ/mol)** at
  all four — default trajectory unchanged, and now verbosity-independent.
* ch4_H/T2000/f10 at `-verbosity 1` and `4`: 18, +310.4 kJ/mol (unchanged).
* `ctest -R "gfnff|sqm_val|react"` in `build_rev`: **95/98**. The 3 failures
  (`cli_simplemd_16_gfnff_rev_h2_smooth`, `_18_gfnff_rev_nve_vs_gfnff`,
  `_19_gfnff_rev_form_refuses_hbond`) are **pre-existing**: a rebuild of the unmodified source
  (md5 79bb76bb…) fails with byte-identical reasons (`T_max 7946.215041`, both slopes,
  `mean Epot -0.3995 kcal/mol`).

## 7. For future measurements

All c2h6 evidence in `BREAK_TAIL_STATUS.md` was taken at verbosity 3, where 471.1 reproduces —
that analysis stands. But any `-verbosity 4` replay of a react cell made **before** binary
e7b6c62a is on a different trajectory and must not be compared with a lower-verbosity run.
