# RUNAWAY_STATUS — the default's ">= 50 kJ" events are explosions driven by the hydrogen valence budget (2026-09-15, orchestrator)

Binary: frozen copy of `curcuma-head/build/curcuma` at commit 30c186a4 (md5 `6febc06f…`; caffeine
revgfnff -4.67352165 / gfnff -4.67273707 Eh = the QP_STATUS records; the two cells reproduce
471.1 and 310.4 kJ/mol). Harness: the BREAK_TAIL copy, `CURCUMA_VERB=3`, `-md.print_frequency 1
-md.dump_frequency 1`, `CURCUMA_REACTSCAN=1`, then `CURCUMA_SHAREDUMP=1`. Scripts `rows.py`,
`terms2.py`, `markers.py`, `perbond.py` and the raw logs in
`/tmp/claude-1000/-home-conrad-src-curcuma-branches-curcuma/3443e495-…/scratchpad/trace{,1,2}/`.
Caveat: the scan trace itself changes the c2h6 trajectory count (135 rebuilds vs 90 without it) —
the same class as the `-verbosity 4` finding under investigation in `VERBOSITY_TRAJ_STATUS.md`; the
event studied here reproduces bit-identically (+471.1, bond +451.0) in every variant.

## 1. The scan is not the problem

`revgfnff` forces `react_check_every = 1` (method_factory.cpp:387): one scan per MD step. The
four cadence arms are bit-identical (`SCAN_CADENCE_STATUS.md`). The per-scan trace of pair 4-5 in
c2h6/T2000/f16: call 12836 r 1.7797 c 0.9053, call 12837 (0.25 fs later) r 3.0009 c 0.0000. A
displacement of 1.22 a0 in one step is not dynamics. ch4_H/T2000/f10, pair 4-6: 1.9195 -> 2.4246 a0
in one step.

## 2. The energy has already exploded before the "hard swap"

Per-step status rows (Time in ps, Eh):

| c2h6 t | Epot | Ekin | T / K | ch4_H t | Epot | Ekin | T / K |
|---|---:|---:|---:|---|---:|---:|---:|
| 3.20700 | -0.9611 | 0.0920 | 3229 | 1.42300 | -0.6386 | 0.0398 | 2096 |
| 3.20725 | **-1.3145** | 0.0949 | 3329 | 1.42325 | -0.6646 | 0.1535 | 8081 |
| 3.20750 | -1.2768 | 0.0711 | 2496 | 1.42350 | **-0.9147** | 0.2858 | 15040 |
| 3.20775 | -1.3075 | 0.0981 | 3441 | 1.42375 | -0.5818 | 0.0800 | 4210 |
| 3.20800 | -1.1561 | 0.2835 | 9945 | 1.42400 | -0.9918 | 0.4426 | 23292 |
| 3.20825 | **-0.1683** | 0.1395 | 4896 | 1.42425 | -0.5770 | 0.4131 | 21743 |
| 3.20850 | -0.9637 | **1.7699** | 62098 | 1.42575 | 11.87 | 744.6 | 3.9e7 |
| 3.20875 | hard swap +471.1 | | | | | | |

Per-step swings of 0.35-1.0 Eh in the Bond and bonded-Repulsion terms with **no** REACT event and
no topology re-perception in those steps (markers.py: the fingerprint-mismatch / rebuild lines sit
at 3.20675 and 3.20825+, not at 3.20725). The dE_jump statistic counts only topology events, so
none of this was ever visible in the 22-cell tables.

## 3. Cause: the valence budget grants a bridging hydrogen two full wells

`CURCUMA_SHAREDUMP=1`, c2h6, the same corner (8 bonds) in two adjacent steps:

| t / ps | r(H4-H5) | Val(H4) | c(C1-H4) | c(C1-H5) | c(H4-H5) | E(C1-H4)+E(C1-H5)+E(H4-H5) |
|---|---:|---:|---:|---:|---:|---:|
| 3.20700 | 1.7236 | 1.058 | 0.505 | 0.505 | 0.010 | -0.276 Eh |
| 3.20725 | 1.4731 | **1.928** | **0.992** | **0.993** | **0.985** | **-0.742 Eh** |
| 3.20800 | 1.6281 | 1.469 | 0.727 | 0.722 | 0.449 | -0.482 |
| 3.20825 | **0.5586** | 2.000 | 1.000 | 1.000 | 1.000 | -0.665 (brep +1.08) |
| 3.20850 | 1.7797 | 1.012 | 0.500 | 0.500 | 0.000 | -0.245 |

ch4_H, pair 4-6: Val(H4)/Val(H6) 1.151/1.160 at 1.42325 -> 1.996/1.999 at 1.42350, c(C1-H4)
0.529 -> 0.999, c(H4-H6) 0.065 -> 1.000, r(H4-H6) 1.694 -> 1.130 -> 0.781 a0.

Mechanism: `Val_i = Val_Z + softplus_50(sum_k shareClip(2 b_ik - 1) - Val_Z)` (Pass 4 of
`prepareValenceShare`). When the transient H-H tight bond order crosses the settled window, the
hydrogen's budget goes 1 -> 2 and BOTH its wells (C-H and H-H) go from half to full share: -0.47 Eh
over a 0.25 a0 change of one distance, a force of ~2 Eh/a0 pulling the H2 together. The system
runs into that artificial minimum (r(H-H) 0.56 a0, bonded repulsion +1.08 Eh), is thrown apart,
and the scan then finds the pair beyond its window: the hard swap and its +451 kJ neighbour-fc
re-derivation (BREAK_TAIL) are the aftermath. A hydrogen has one valence; a bridging H is a 3c-2e
bond whose two partial wells share it. The old fixed `Val = revValence(Z)` did exactly that, which
is why it measured 10.7 kJ / 0 events on the same grid — mechanism, not exposure (this corrects the
"exposure" reading in BREAK_TAIL and the vault note's Nachtrag 6).

## 4. Consequences

- The smoothness metric must include **max |dEpot| per MD step outside rebuild events**; the
  rebuild dE_jump missed the whole runaway.
- Fix under test (`HBUDGET_STATUS.md`, worktree `curcuma-head`): `rev_budget_fix_h` — hydrogen keeps
  Val = 1 exactly, derivative channel zero. Expected: runaway gone on both cells, hypervalent set /
  BF4- / equilibria bit-identical (no H there has two partners), rkt06 ~2.71.
- Noted, not addressed: the same rule gives the carbon of CH4 + H a budget of 4.95 (five full C-H
  wells in the H-H-unbonded corner) — probably the root of the class-B "H + CH4 TS 78 kcal/mol too
  low" question; a hypervalent budget should be limited to elements that can be hypervalent.
