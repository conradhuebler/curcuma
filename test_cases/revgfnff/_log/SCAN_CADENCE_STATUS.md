# SCAN_CADENCE_STATUS — does scanning more often remove the BREAK_TAIL hard swaps?

**Frozen binary** `curcuma_frozen` (worktree `curcuma-head/build`, commit 30c186a4), md5
**6febc06f4a72edd6db1968e2c0744e1d** on every one of the 88 `wall.txt` records below. Provenance
check before the sweep: `-sp caffeine.xyz` gives revgfnff **-4.67352165 Eh**, gfnff
**-4.67273707 Eh** — exact match to the QP_STATUS records. Harness `<scratch>/bt/` (BREAK_TAIL
copy, `W=`/joblist repointed into this session, `joblist22_fixed.txt` = 22 unique cells), runs in
`<scratch>/bt/runs_{A,B,C,D}/`. `CURCUMA_VERB=3`. 2/22 cells (`ch4_H/T1000_f16`, `ch4_H/T2000_f16`)
exit rc=1/0 rebuilds in every arm (pre-existing harness defect) — 20 live cells throughout.

**Provenance gate (arm A vs QP_STATUS A.2): PASSED**, bit-exact: 960 rebuilds, max 471.10,
5/478 hard swaps (s>=0.99), T_max 2036.57 K (=2037 rounded). Mean wall time 0.81 s/cell (small
molecules, GFN-FF, 20000 steps at dt 0.25 fs) — genuine full 5 ps runs, confirmed by `REACT
summary` lines reaching t=4998-4999 fs and rc=0 in every live cell.

## Four-arm table (20 live cells)

| arm | rebuilds | median | p99 | max &#124;dE&#124; | frac&lt;1kJ | n&ge;50kJ | T_max | hard/begin | wall (s) |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| A default | 960 | 0.0 | 1.8 | 471.1 | 0.9812 | 3 | 2036.6 | 5/478 | 16.21 |
| B `react_check_every 1` | 960 | 0.0 | 1.8 | 471.1 | 0.9812 | 3 | 2036.6 | 5/478 | 16.15 |
| C `react_check_every 2` | 1012 | 0.0 | 1.4 | 392.5 | 0.9852 | 3 | 2036.6 | 4/506 | 15.95 |
| D `every 1 + disp 0.05` | 960 | 0.0 | 1.8 | 471.1 | 0.9812 | 3 | 2036.6 | 5/478 | 16.13 |

**A, B, D bit-identical** on every stat and event (same rebuild indices, `s`, `w_ws`, `r_ws`, terms).
Only **C** (coarser than default) differs.
## Top-3 events

**A/B/D** (identical): 1) `c2h6/T2000_f16` #42 `begin_break` **+471.1**, s=1.00, bond +451.0
angle +28.9 brep -12.1 nbrep +3.7 (ij 3 4, H4-H5). 2) same pair #83 `complete` +471.1 (marker,
terms 0). 3) `ch4_H/T2000_f10` #4 `begin_break` **+310.4**, s=1.00, bond +276.4 angle +57.2
brep -37.9 nbrep +14.8 (ij 3 5).

**C**: 1) `ch3nh2/T2000_f16` #94 `begin_break` **+392.5**, s=1.00, bond +372.0 angle +41.3
brep -33.4 nbrep +12.5 (ij 4 6) — different pair/molecule than any A/B/D top event. 2) same pair
#187 `complete` +392.5. 3) `ch4_H/T2000_f10` #6 `begin_break` **+218.0**, s=1.00, bond +200.1
angle +18.0 (ij 3 4) — smaller than A/B/D's own ch4_H event, different rebuild.
## c2h6/T2000_f16 window (H4-H5)

**A/B/D**: bit-identical to the BREAK_TAIL reproduction — H4-H5 still hard-swaps at **s=1.00,
t=3208.8 fs**, broken in the *same* single scan step (blend 0.12->0.02 within one 2 fs interval;
begin and complete are the same rebuild, 0 intervening scans), largest single-step &#124;dE&#124;
= 471.1 kJ/mol.

**C**: diverges from A already at **rebuild #2** (dE -0.000013 vs -0.000017 Eh, first few fs) —
chaotic MD, no shared reference survives to t=3208.8 fs, so "same event" is not meaningful here
(BREAK_TAIL's own caveat). In C's own trajectory H4-H5 never reaches s=1.00 — no hard swap in this
cell, largest event **-0.7 kJ/mol** (#47). Divergence moved the trajectory away from the collision;
it did not resolve it — the 20-cell count still has 4 hard swaps under C, one (`ch3nh2/T2000_f16`)
newly at s=1.00, +392.5 kJ/mol.

## Wall-time cost

Negligible: 15.95-16.21 s total/arm over 20 cells (0.80-0.81 s/cell); C is marginally *cheaper*
than A/B/D, not more expensive.

## Verdict

**Cadence alone does not remove the hard swaps, and "scan more often" is not an available lever:**
the default already scans every single energy call for `revgfnff` (`react_check_every` is
force-set to 1 in `method_factory.cpp`'s preset, ahead of any CLI value), so `-gfnff.react_check_every
1` (B) and a tighter 0.05 Bohr displacement trigger (D) are no-ops — bit-identical to A on every
event. The only cadence value that changes anything is going *coarser* (`react_check_every 2`,
arm C), and it does not fix the problem: `begin_*` events rise (478->506), the hard-swap count
drops only 5->4, and a different pair (`ch3nh2` C-H) gets its own s=1.00 swap at +392.5 kJ/mol.
C2H6's own improvement is chaotic trajectory divergence moving that cell away from the collision,
not the coarser scan resolving the crossing — globally the tail is **moved, not removed**. Cost is
not the constraint (~16 s/arm either way). Points away from `react_check_every`/`react_check_disp_bohr`
and toward the transition-window mechanics itself (BREAK_TAIL: the window can be crossed faster
than any per-scan cadence can catch, at dt=0.25 fs).
