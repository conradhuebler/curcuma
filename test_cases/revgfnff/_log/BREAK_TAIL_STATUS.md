# BREAK_TAIL_STATUS — what the three `begin_break` hard swaps actually cost

Binary `build_rev/curcuma`, md5 **79bb76bb6720278c9da8b35007ed50ee** (unchanged; **no src edit** was
needed — `CURCUMA_SHAREDUMP`, `CURCUMA_BONDDUMP` and the verbosity-3 `BOND_FACTORS` line already
carry everything). Harness `<scratch>/bt/` (QP copy, `W=` repointed, `md5sum "$BIN" >> wall.txt`
added, `NOFILT=1` keeps the raw log); runs in `<scratch>/bt/{runs,dump,bp,bf,noshare}/`.
**Reproduction, exact**: c2h6/T2000/f16 `+471.1` (rebuild #42), ch4_H/T2000/f10 `+310.4` (#4) and
`+51.6` (#6) — identical in the plain, the `SHAREDUMP` and the `BONDDUMP` replays.

## Verdict

**None of (a), (b) or (c).** Not the swapped pair's own well, not the share/budget, not the dynamic
CN: it is the **static bond force constant `fc` of the OTHER bonds**, re-derived from the new
topology at the list edit. Per event, in kJ/mol:

| event | bond term | (a) own well | (b) share `c`/`Val` | (c) shape r0/alpha | **`fc` re-derivation** | resid |
|---|---:|---:|---:|---:|---:|---:|
| c2h6 f16 `+471.1` | +451.0 | 0.00 | 0.00 | −0.46 | **+451.41** | −0.05 |
| ch4_H f10 `+310.4` | +276.4 | 0.00 | 0.00 | +0.68 | **+275.72** | −0.00 |
| ch4_H f10 `+51.6`  | +56.2  | 0.00 | 0.00 | +0.00 | **+56.19**  | −0.01 |

(a) is 0 because the swapped pair's own well is already exactly 0 at the swap: c2h6 `w = 0.000000`
at r = 3.0009 a0, ch4_H(310) `c = 0.000000`, ch4_H(51.6) r = 24.86 a0. (b) is 0 because every
`share`/`shareD` column — `w, b, sig, g, sum, Val, u, f, c` — is **bit-identical before and after**
(c2h6 all `c = 1.000000`; ch4_H `0.823339 / 0.476335 / 0.000000` in both). (c) is 0 by construction:
`m_d3_cn` comes from `calculateGFNFFCNWithNeighbors`, a pure distance cutoff, and the geometry is
identical across the two evaluations — measured, `exp(-a dr^2)` agrees to 4 dp on every other bond.

## Per-bond table, c2h6/T2000/f16 (the +451.0)

| bond | r [a0] | E before | E after | dE [kJ] | fc before | fc after |
|---|---:|---:|---:|---:|---:|---:|
| C1-H4 | 2.4648 | −0.23485237 | −0.14327291 | **+240.44** | 0.272567 | 0.166691 |
| C1-H5 | 2.9351 | −0.18645593 | −0.11455791 | **+188.77** | 0.272567 | 0.166691 |
| C1-H3 | 2.1170 | −0.16844943 | −0.16049144 | +20.89 | 0.174956 | 0.166691 |
| C1-C2, C2-H6/7/8 | — | — | — | +0.85 total | ~0.1668 | ~0.1667 |
| **H4-H5 (swapped pair)** | 3.0009 | −0.00000000 | −0.00000000 | **0.00** | — | — |

ch4_H(+310.4), same shape: C1-H4 +111.5, C1-H6 +112.8 (the two H of the broken pair),
C1-H2/H3/H5 +18.1/+17.3/+16.7, swapped pair 0.00. ch4_H(+51.6): +19.4/+19.2/+17.6 on the three
surviving C-H, swapped pair (C1-H2, r = 24.9 a0) 0.00.

## Why `fc` moves — named factors

`BOND_FACTORS` (verbosity-3 line in `getGFNFFBondParameters`) on the ch4_H rebuild gives the whole
product. While the H-H bond is in the list: `bstr` **1.3234** instead of 1.0000 — the hydrogen has
2 neighbours, so it is hybridised **sp** and its C-H bond takes the `bsmat[sp][sp3]` entry;
`ringf` **1.1800** — C-H-H is perceived as a **3-membered ring** (both apply to the two ring C-H
only); and `fxh` **1.0500**, the 3-ring C-H rule (`gfnff_method.cpp:5003`), which fires on **every**
C-H of that carbon including the untouched ones — that is the whole "sibling" column above.

Product `1.3234 x 1.18 x 1.05` = 1.6397, times the small `fqq` shift = **1.6353** — measured
`0.272567/0.166691 = 1.635163` (c2h6), `0.272862/0.166720 = 1.636648` (ch4_H); sibling
`1.05 x fqq = 1.049594` vs measured `1.049583`. Five digits, both molecules. (c2h6's factors are
arithmetic, not printed: `-verbosity 4` changes that cell's trajectory, below.) So while the
transient H2 exists the surrounding C-H wells are **64 % deeper** than a real C-H bond.

## Pre-event REACT timeline

**c2h6/T2000/f16** — H4-H5 formed t=3169.8 fs (r_scan 2.0074), blend completed at r 1.1888; break
begun 3184.5 (r 2.0146) and **reverted** 3186.2; begun 3206.2 (r 2.0398) and **reverted** 3206.8
(r 1.8847); **3208.8 fs hard break at r_scan 3.0009, w_scan 0.0206, s = 1.00**. The pair crossed
1.88 -> 3.00 a0 inside one 2 fs scan interval, so the blend window was already behind it and the
swap completed in one step. Last (re-)formation 2.0 fs before the event, original formation 39 fs
before. 0.2 fs later C1-H5 breaks too (−0.1 kJ).

**ch4_H/T2000/f10** — H2-H4 formed 1375.5, reverted 1377.0 (+29.0 kJ); H3-H5 formed 1392.8,
reverted 1394.8; **H4-H6 formed 1419.8 (r 2.0245), kept, broken 1424.5 at r 2.4246 → +310.4** —
formed **4.7 fs** earlier, so yes, under 5 fs. That jump injects 0.118 Eh into 3 hydrogens and the
system flies apart: H2-H4 re-forms 1425.2, and the **+51.6** at 1426.0 is C1-H2 leaving the list at
**r = 24.9 a0** — a bookkeeping break in an already-exploded system, not a chemical event.

## `rev_valence_share false` (item 5)

There is **no switch for the softplus budget alone**; `Val_i` is used only by the share, so the only
lever is `-gfnff.rev_valence_share false`, which removes both. It does **not** remove the jump
class: c2h6 f16 goes to max **+383.3** / min **−641.3** kJ and ch4_H f10 to max **+453.6**
(n = 2 cells, one seed) — consistent with the share contributing exactly 0 above. That arm is *not*
the `10.7 kJ / 0 events` reference of VALFIX §6, which is share **on** with the old
`Val = revValence(Z)` (a source change, no runtime flag). So the budget did not create this
discontinuity; it changed the exposure to one that is present in every arm.

## Not determined

* Why the old `Val = revValence(Z)` avoided these configurations in all 22 cells — that is a
  trajectory/exposure question, not measurable from the event itself.
* **`-verbosity 4` changes the c2h6/T2000/f16 trajectory** (136 rebuilds, max +0.6 kJ, against 90
  and +471.1 at verbosity 1/2/3; deterministic on repeat). The only `>= 4` branch in the FF is a log
  line (`gfnff_method.cpp:1108`), so this is not root-caused. All c2h6 evidence above was taken at
  verbosity 3, where the 471.1 reproduces exactly. ch4_H is unaffected (its events reproduce at 4).
