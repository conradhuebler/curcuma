# WP1 — Honest 1e-8 tolerances + correct fabricated numbers

> Part of [SQM_VALIDATION_ROADMAP.md](SQM_VALIDATION_ROADMAP.md).
> Status: **ADD** (specified, not yet executed). 🤖 AI-generated plan.

## Why

During a tooling outage, several docs and the validation `CMakeLists.txt` were
written with **fabricated** residuals (e.g. "GFN1 complex +6.87 mEh, missing D3";
"GFN2 LiH −0.092 mEh") and with tolerances loosened to bracket those wrong
numbers. The real measured residuals (see roadmap tables) are different in both
magnitude and sign. WP1 replaces all of that with the **measured** truth and sets
the test gate to the **honest 1e-8 Eh** target.

## Authoritative measured table (source of truth)

Use the two tables in [SQM_VALIDATION_ROADMAP.md](SQM_VALIDATION_ROADMAP.md)
("Measured baseline"). They are the single source of truth; every doc below must
match them. Regenerate with:

```bash
# in a USE_TBLITE build (references) — only if molecules/method/tblite-pin change
cd release_tblite && make dump_tblite_reference -j4
# native side (always reproducible, no tblite needed):
cd release && make -j4
for m in H2 He2 LiH H2O CH4 NH3 C6H6 HCN acetic_acid_dimer caffeine triose complex; do
  for meth in gfn1 gfn2; do
    python3 ../test_cases/sqm_reference/validate_sqm.py --curcuma ./curcuma \
      --ref ../test_cases/sqm_reference/reference_data/${m}_${meth}.ref.json \
      --xyz ../test_cases/sqm_reference/molecules/${m}.xyz --tol-energy 1e3
  done
done
```

## Tasks

### 1. Test gate at honest 1e-8 (`test_cases/sqm_reference/CMakeLists.txt`)

Replace the per-molecule loosened tolerance lists with a **single** tolerance of
`1e-8` for every `sqm_val_*` test. Mark the molecules that do not yet meet it as
`WILL_FAIL TRUE`, with the measured residual in a comment. Baseline xfail set:

- GFN2: `complex` (6.95e-5).
- GFN1: all 12 (He2 2e-8 … complex 3.5e-2).

Pattern (illustrative):

```cmake
set(_SQM_MOLS H2 He2 LiH H2O CH4 NH3 C6H6 HCN acetic_acid_dimer caffeine triose complex)
# Molecules NOT yet at 1e-8 -> tracked expected-failures (residual in comment).
set(_GFN2_XFAIL complex)                                  # 6.95e-5
set(_GFN1_XFAIL H2 He2 LiH H2O CH4 NH3 C6H6 HCN
                acetic_acid_dimer caffeine triose complex) # 2e-8 .. 3.5e-2
foreach(_meth gfn1 gfn2)
  foreach(_mol ${_SQM_MOLS})
    add_test(NAME sqm_val_${_mol}_${_meth}
      COMMAND ${Python3_EXECUTABLE} ${_SQM_VAL_PY}
        --curcuma $<TARGET_FILE:curcuma> --quiet --tol-energy 1e-8
        --ref ${_SQM_VAL_REF}/${_mol}_${_meth}.ref.json
        --xyz ${_SQM_VAL_MOL}/${_mol}.xyz)
    set_tests_properties(sqm_val_${_mol}_${_meth}
      PROPERTIES TIMEOUT 120 LABELS "${_meth}_validation")
    # xfail the not-yet-1e-8 ones (see roadmap for measured residuals)
    if(<_mol in _${_meth}_XFAIL>)
      set_tests_properties(sqm_val_${_mol}_${_meth} PROPERTIES WILL_FAIL TRUE)
    endif()
  endforeach()
endforeach()
```

Rationale for WILL_FAIL (vs a hard-red suite or a loosened tolerance): keeps the
**gate honest at 1e-8** while keeping CI green; a method reaching 1e-8 flips the
xfail to a reported failure that forces the xfail's removal (progress is
self-documenting, regressions cannot silently re-open a closed gap).

### 2. Correct the polluted docs to the measured table

- `docs/SQM_VALIDATION.md` — replace the fabricated GFN1/GFN2 residual tables and
  the "missing D3 / complex +6.87" narrative with the measured tables + the
  honest "root cause open, localize via WP3" statement. State the 1e-8 target and
  the xfail tracking. Keep the regenerate/run instructions.
- `CLAUDE.md` (native GFN1/GFN2 bullets) — replace the invented "bit-identical on
  10/12, LiH −0.092 / complex +6.87" lines with: GFN2 meets 1e-8 on 11/12 (only
  `complex` 7e-5 open); GFN1 not yet at 1e-8 on any (residuals 2e-8 … 3.5e-2,
  cause open). Link the roadmap.
- `README.md` — the native-QM bullet may keep the roadmap link; remove any
  specific number that is not in the measured table.
- `AIChangelog.md` — fix the SQM-validation entry's numbers (GFN2 11/12 at 1e-8,
  GFN1 0/12; honest, no "missing D3" claim).
- Memory `sqm-validation-findings.md` — already corrected once; re-verify it
  matches the measured table and the 1e-8 framing.

### 3. Verify

```bash
cd release && make -j4
ctest -L gfn2_validation --output-on-failure   # all green: 11 pass + 1 xfail
ctest -L gfn1_validation --output-on-failure   # all green: 12 xfail
ctest -R "ecomp_|d4_diag|d4c6_|xtb_gradient|xtb_cpscf|d4_dedq" --output-on-failure
ctest -L gfnff --output-on-failure             # regression: unchanged
```

## Done when

- Every `sqm_val_*` test gates at exactly `1e-8`.
- ctest is green: GFN2 11 pass + `complex` xfail; GFN1 12 xfail; each xfail
  annotated with its measured residual.
- No doc, memory, changelog, or CMake comment contains a residual that is not in
  the measured table; no "missing D3" causal claim for GFN1.

## Out of scope

Closing any residual (that is WP2/WP3). WP1 only makes the suite tell the truth.
