# rev-PARAM delivery audit: does a plain `-sp` run honour the rev CLI flags? (2026-09-13, sp-audit)

Read-only audit, no src/ change, no build. Frozen binary: `build_rev/curcuma`, mtime
**2026-09-13 21:08:15**, md5 **b9c80176f8f411172ea46b6ddb04e049** (copy kept at
`.../scratchpad/audit/curcuma`) - this is the binary the audited measurements were taken on; the
tree moved during the audit and `build_rev/curcuma` was rebuilt at ~21:4x (md5
`bbefeece9fc85376c8135917ceb0cb40`, see §2).

## Verdict

**The premise is false, and the observed inertness was an argv-tokenisation artefact.**
`setupRevSettings()` is *not* ctor-only: it sits inside `loadParameterOverrides()`
(`gfnff_method.cpp:3315` at HEAD), and that function is called by **both** the constructor and
`setParameters()` (`:3277`), which flattens the CLI `gfnff` sub-config into `m_parameters` first.
So the whole `m_rev_settings` struct is re-read after the CLI merge. Measured: the same flag+value
passed as **two argv elements** works in `-sp` *and* in `-batch`; passed as **one argv element**
(`"-gfnff.rev_valence_share false"`) it is ignored in `-sp` *and* in `-batch`. The `-sp`/`-batch`
distinction in the original measurement was a coincidence of how the two loops were written.

| arm (frozen binary, `-sp … -gfnff.cache_topology false`) | E [Eh] |
|---|---|
| `-gfnff.rev_valence_share false`, 2 argv elements / 1 argv element (joined) | **0.82730126** / 1.17232892 |
| the same two arms under `-batch true` | 0.82730126 / 1.17232892 |
| replay of the orchestrator's own argv (joined) / split, their geometry | 1.17562701 / 0.82889302 |

The last row reproduces the two anchor numbers exactly, with the binary they were taken on.

## 1. Affected flags (enumeration from code, not from the task text)

41 registered `rev_*` PARAMs (`gfnff.h`). **40** are read (only) inside `setupRevSettings()` and
are therefore "affected" under the premise; **1** (`rev_over_p`) is also read lazily
(`revOverP():12033`, `fillRevPerAtom():12246`) and is not. Affected: `rev_enabled, rev_bond_weight,
rev_term_weights, rev_blend_repulsion, rev_over_coord, rev_valence_share, rev_h_not_sp, rev_blend,
rev_bo_* (center/width/form/break/2_*/3_*/4_*/5_*/13_form), rev_form_switch, rev_over_k,
rev_over_shift, rev_over_preset, rev_tr_* (begin/end/revert/prebreak), rev_demote_cooldown,
rev_max_transitions, rev_charge_model, rev_sqe_bmin, rev_sqe_kappa_{H,C,N,O,F,Cl}` (the six kappa
names are read through a name table, so a naive string scan misses them).
**Dead at HEAD** (registered, read nowhere): `rev_valence_share`, `rev_h_not_sp`. Dead in the
frozen 21:08 binary: `rev_h_not_sp` only. The working tree has since re-added both reads (`:12062`,
`:12066`; the other agent's comment at `:12063-12065` states the same finding). A dead flag is dead
in **-sp, -batch and MD alike** - it cannot produce a "-sp inert, -batch works" pattern.

## 2. Mechanism settled by measurement (frozen binary, split tokens)

- **`-gfnff.topology_mode react -verbosity 1` banner** (20 values, one run): with a 23-flag CLI set
  the banner reports `form switch weight`, `bo2 > 0.111`, `w > 0.071`/`< 0.031`, `bo3 1.61x k -8.1`,
  window `0.02..0.82`, `shift 0.510 k 11.0`, `valence share off` - defaults are `order`, `0.100`,
  `0.050/0.020`, `1.60x/-8.0`, `0.02..0.80`, `0.500/10.0`, `on`.
- **Energy moves in a plain `-sp`**: `rev_over_coord false` 1.17232892 -> 1.10354329;
  `rev_over_p 0.9` -> 1.30990018; `rev_over_shift 0.7` -> 1.12615438; also `rev_over_k`,
  `rev_bond_weight false`, `rev_bo_width`, `rev_bo_center`, `rev_bo_form`, `rev_bo4_width`;
  `rev_charge_model sqe + rev_sqe_kappa_H 0.1` -> 1.17688819. **23 of the 40 affected flags were
  positively observed reaching the model in `-sp`**; the remaining 17 (react/MD-only quantities:
  `bo5_*`, `bo13_form`, `demote_cooldown`, `max_transitions`, `blend*`, `term_weights`, `sqe_bmin`)
  are inert *for a static single point* by construction and share the one code path above.
- `rev_h_not_sp` shows `on` in the banner under `-gfnff.rev_h_not_sp false` **in the 21:08 binary**
  - consistent with the dead read, not with a path defect. **Direct confirmation after the rebuild**:
  the binary `build_rev/curcuma` rebuilt during this audit (md5 **bbefeece9fc85376c8135917ceb0cb40**,
  i.e. the other agent's missing-read fix) reports `H never sp off` for the same plain `-sp` call.
  The fix that was needed is the read, not the call order - and joined-token inertness is unchanged
  in the new binary (joined 1.17232892, split 0.82730126, as above).

## 3. Per-verification classification

| source | command path (quoted in the file) | verdict |
|---|---|---|
| `PRESET_STATUS.md` / `b4a3e0a4` | plain `-sp` single points in `release/`, `-dump_gradient`, `-gfnff.rev_enabled true`, fresh dir + `-gfnff.cache_topology false` | **b, re-measured, holds** |
| `REACT_JOIN_STATUS.md` / `4601be27` | MD: `-md`, water dimer 1 ps dt 0.25, H4 5 ps; + `ctest cli_simplemd_19` | **a (and positively confirmed)** |
| `NVE_TEST_STATUS.md` | MD / ctest (`cli_simplemd_18/19`); §4 `static` row = SP of one fixed geometry, no rev flag toggled | a |
| `WALL_ETOT_STATUS.md` | MD (NVE), 5 Etot assembly sites; no rev CLI toggle | a |
| `CLASSA_FROZENCN.md` | `-batch true` in all three modes (mode 3 also `-batch_reuse_topology false`) | a |
| `POLY_JUMP_BASELINE.md` | `-md … -gfnff.topology_mode react` grid, 5 ps (array-tokenised `run_one.sh`) | a |
| `R0_FIX_STATUS.md` / `e36d9925` | `-sp scan.xyz -batch true -batch_reuse_topology …`; `fdcheck.py` fresh dir + batch | a |
| `CIJ_STATUS.md` §2/§3 | arms built as **source variants** (`cij/build_hsp.log`, `build_nolam.log`, …) | a |
| `CIJ_STATUS.md` §4 (2x2) | `-batch true` CLI toggles (`toggle2.sh`); share axis effective, **H-not-sp axis was a dead flag** in this binary | **b/c, H-not-sp axis vacuous - re-derive** |
| `OUTLIER_STATUS.md` / `c66ffa13` | `-batch true -batch_reuse_topology …` | a |
| `EEQ_DRIFT_STATUS.md` | `-batch true -batch_reuse_topology true`, class-A protocol | a |
| `AGENT_STATE.md` "blockers (1)/(2)" | (1) `-batch true`, split tokens (`b-nh4-*`): share axis effective; (2) is this audit's premise | (1) holds, (2) **refuted** |

No verification in the list was a plain `-sp` run whose conclusion actually depended on an
undelivered rev flag. The two `-sp` entries (`PRESET_STATUS`) were re-measured and hold:

- `preset fit2026-09-12` == `-gfnff.param_file params/rev_over_fit_2026-09-12.json` (identical to
  8 dp on self-made H3 and CH5; the file's own H3/CH5 rows reproduce in structure exactly);
- an explicit `-gfnff.rev_over_shift 0.7` **takes effect in `-sp`** (CH5 -0.34146896 -> -0.34060998,
  i.e. 8.6e-4 Eh; banner shows `shift 0.700`), and `0.5` == fit default == the documented
  limitation; `stage1a` == no flag == default.
- `4601be27`'s toggles are safe for a stronger reason: an inert flag cannot change an event count.
  `rev_form_switch weight` reproduced the pre-change counts (water dimer 18/15/33 at 300 K,
  H4 94/94/240) against the new default's 0/0/0 and 17/17/50, and test 19's negative control
  flips 0/0/0 -> 6/3/9.

## 4. Answers to the two explicit questions

1. **Join work (`4601be27`, `REACT_JOIN_STATUS.md`)** - every run behind it is MD (`-md … -gfnff.
   topology_mode react`, water dimer 1 ps dt 0.25; H4 5 ps; the 6-formed/9-rebuild figure is
   `cli_simplemd_19`'s CSVR MD sub-run), plus ctest. No plain `-sp`/`-batch` single point carries
   any part of that conclusion, and the toggle demonstrably changed results (0/0 vs 6/9; 17/17/50 vs
   94/94/240) - an inert flag cannot. **Nothing to re-derive.**
2. **Preset work (`b4a3e0a4`, `PRESET_STATUS.md`)** - the runs were plain `-sp` single points, i.e.
   the nominally affected path, but they were **not** inert: re-measured here, the preset flag, the
   equivalent `param_file` and an explicit `rev_over_shift` all reach the model in `-sp` (numbers in
   §3). The 12-digit preset==param_file claim therefore stands; it is the *default/shift-equivalence*
   limitation, not flag delivery, that bounds it.

## 5. What must be re-derived

- `CIJ_STATUS.md` §4's 2x2: its H-not-sp axis was a dead flag in the 21:08 binary, so that axis was
  never exercised (a trivially-zero `dE` is expected either way). Re-run after the read at `:12066`
  is in a build. §2/§3's arm table is unaffected (source-variant builds).
- `AGENT_STATE.md` blocker (2) and the valfix task "defect 2": as written it is not a defect - the
  struct *is* re-read after the CLI merge. A "fix" there will show no before/after difference on any
  rev flag, which must not be read as "the fix did not take". The real defects were the two
  unassigned reads (`rev_valence_share` at HEAD, `rev_h_not_sp` at HEAD and in the 21:08 binary),
  already addressed in the working tree.
- Method note: the Bash tool's shell is zsh and does not word-split, so every `-flag value` must be
  its own argv element (array or bash script file) - the trap already recorded in
  `NVE_TEST_STATUS.md` §7, which the 21:09 `-sp` runs fell into.

## 6. Reproduce

All runs are in `.../scratchpad/audit/run/` (`t3.sh` joined-vs-split, `t5.sh` verbatim replay of the
orchestrator's argv, `banner.sh`, `sweep2.sh`, `ch5.sh`) plus the originals in `.../verify_cij/`.
