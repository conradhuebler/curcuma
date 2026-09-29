# rev-gfnff stage 3a — the valence share, the budget rule and the bond-well form

🤖 AI-generated, machine-tested only. **Human production testing pending** (the operator removes
this line). Every number below was measured with one frozen binary per package; the run logs and
the per-cell tables are in `test_cases/revgfnff/_log/WORK_STATUS.md`.

Stage 3a modifies the GFN-FF **bond term** so that it survives a bond being made or broken. All
three parts are on by default since Sep 19, 2026:

| part | what it changes | flag | default |
|---|---|---|---|
| 3a(i) | the pair's own CN contribution is taken out of its own r0 | — | on |
| 3a(ii) | a valence SHARE factor multiplies the well | `-gfnff.rev_valence_share` | on |
| 3a(ii) | which share formula | `-gfnff.rev_share_form` | **`conserving`** (was `delivered`) |
| 3a(ii) | a dative donor may grow a budget | `-gfnff.rev_share_donor_rule` | **on** (new) |
| 3a(ii) | hydrogen keeps one valence | `-gfnff.rev_budget_fix_h` | on — but a **no-op** under `conserving` |
| 3a(iii) | the SHAPE of the well | `-gfnff.rev_well_form` | **`mg3`** (was `mg`, before that `gauss`) |

**Absolute energies moved with the well-form flips and relative ones did not.** An MG well carries
its own fitted depth `D = s |k_b|`, so `revgfnff` is further from plain `gfnff` in absolute terms
— caffeine **-4.789915106585** under the current `mg3` default, -4.546943898047 under `mg`,
against `gauss` -4.673521653477 and `gfnff` -4.672737068614. The pooled 167-reaction conformer/S66
guard moves only 1.0341 (`gauss`) -> 1.0439 (`mg`) -> 1.0547 (`mg3`) kcal/mol MAD. Do not compare
an absolute `revgfnff` energy with an absolute `gfnff` one.

---

## 1. The share: why a pairwise well needs one

A GFN-FF bond well is pairwise and knows nothing about the atom's other bonds. At an exchange
transition state X···H···Y both wells are evaluated at full depth, so the pair sum is two wells
where one bond's worth of valence exists. The share multiplies each well by a factor `c_ij` built
from how much of each end's valence its partners already claim.

### 1.1 `delivered` — the left-over rule (the default)

    S_i   = sum_k w_ik g_ik                     the claim of i's partners (w = the term weight)
    Val_i = Val_Z + G(settled_i - Val_Z)        G = softplus, settled = the tight bond order
    f_i   = clip((Val_i - sum_{k != j} w_ik) / w_ij)
    c_ij  = 1 - g_ij (1 - (f_i + f_j)/2)

`f_i` asks what the OTHER partners left over. A saturated equilibrium atom and a lone bond both
clip to exactly 1, so **the equilibrium term is bit-identical** to plain GFN-FF's.

### 1.2 The hydrogen budget (`rev_budget_fix_h`, DEFAULT ON since Sep 18, 2026)

`Val_i = Val_Z + G(settled_i - Val_Z)` is right for a hypervalent centre and wrong for hydrogen: a
bridging H — a just-formed H2 still bonded to its carbon — is a 3c-2e bond with ONE valence, but
the softplus lets its budget reach 2 as soon as the second partner's tight bond order crosses the
settled window, and BOTH of its partial wells go from half share to full share inside one 0.25 fs
step with **no topology event**. With the flag on, `Val_H = 1` exactly and its derivative channel
is zero. **Under the `conserving` share below the flag does nothing**: hydrogen's excess-budget
cap there is `X = 0` by element, so its budget cannot grow whatever the flag says. The numbers in
this section are therefore measured with `-gfnff.rev_share_form delivered -gfnff.rev_well_form
gauss`, which is also how `cli_simplemd_20` now pins its two arms.

Measured on the 20-cell react-MD grid (5 ps, dt 0.25 fs, CSVR, `-threads 1`), on vs off:

| metric | fix_h on (default) | off |
|---|---:|---:|
| max per-step abs(dEpot) outside rebuild intervals | **59.26** kJ/mol | 2593.60 |
| events >= 50 kJ/mol | 487 | 324 |
| max rebuild dE_jump / events >= 50 | 48.6 / **0** | 471.1 / 3 |
| hard swaps (`s >= 0.99` at a `begin_*`) | **0 of 591** | 5 of 478 |
| T_max | **8 306** K | 1.594e8 |
| min r(H-H) ever, c2h6/T2000_f16 | **1.770** a0 | 0.559 |

No falsifier moves: the six hypervalent ions, the equilibrium toggle set, rkt06 (rms 2.7140) and
the FD gradient at rkt06 point 10 (1.194e-08 Eh/A) are bit-identical in both arms, and `gfnff` is
untouched. Regression test: `cli_simplemd_20_gfnff_rev_h_budget` (two-armed — the `false` control
must violate both bounds).

### 1.3 `conserving` — the valence-conserving share (DEFAULT since Sep 19, 2026)

The left-over rule **forfeits** valence: with `w ~ 1` on every partner, one partner too many makes
every pair of that atom see "nothing left", so a five-coordinate carbon hands out 0 of its 4
valences and only the growing budget papers over it. The conserving rule spreads instead:

    f_i   = min(1, Val_i / S_i)                 PER ATOM, C1 smooth min (exact on both sides)
    c_ij  = f_i f_j                             the PRODUCT, not the mean
    Val_i = Val_Z + min(G(S_i - Val_Z), X_i)

so `sum_j f_i w_ij = min(Val_i, S_i)` exactly. `X_i`, the excess budget, is granted by **charge**,
not by element: 0 for H and F, 1 for group 13, `6 - Val_Z` for period >= 3 groups 15-17, and
`clip(Q_i)` otherwise, with `Q_i` the topological (phase-1 EEQ) charge of the atom plus that of its
H partners in this corner. What separates NH4+ (four real bonds) from NH3 + H (no bond) is the
electron count, and the charge is the only electron-count information a force field has.

Two implementation points carry the accuracy:
- the budget is built from the same **wide** `S_i` the share divides by, which is what makes NH4+
  come out at `Val = S` exactly (residual 0.0013 kcal/mol);
- the smooth min is exactly 1 for `Val >= S` (so an equilibrium atom keeps a literal 1.0 and the
  term stays bit-identical) and exactly `Val/S` below `1 - width` (so the conservation identity
  holds where the share bites).

**What it buys** — the class-C radical-adduct falsifier, model minus r2SCAN-3c, both curves
referred to their own d = 3.00 A point:

| scan | delivered, dev min / rms | conserving, dev min / rms |
|---|---|---|
| CH4 + H | -87.0 / 41.8 | **-1.5 / 12.1** |
| NH3 + H | -107.0 / 47.1 | **-1.4 / 5.1** |
| H2O + H | -54.6 / 24.3 | **-3.0 / 13.3** |
| N2H4 + H | -90.4 / 41.0 | **+0.0 / 7.0** |

and on the grid, 487 -> **72** per-step events >= 50 kJ/mol, with `ch4_H/T1000_f10` going
388.4 -> 33.5 kJ/mol and 40 -> 0 events (the carbon-budget snap of the first MD step disappears).
rkt06 rms 2.7140 -> 2.7614; equilibria, ClO4- and BF4- unmoved; FD gradients <= 1.2e-07 Eh/A.

**What it cost until Sep 19, 2026**: a dative or ylidic neutral has a donor with four partners at
a group charge well below 1, so `Val ~ 3.2` against `S ~ 4` and all four of its wells were scaled
by ~0.8 — where the delivered share is inert (0.002-0.017 kcal/mol). The missing physics is that a
dative bond puts a full valence into the acceptor's empty orbital, which the donor's EEQ charge
(~+0.2) cannot express.

### 1.4 The donor rule (`-gfnff.rev_share_donor_rule`, DEFAULT ON, `conserving` only)

An atom is granted `X_i >= 1` if, **in the corner being evaluated**, it has a partner that is
either

- a **group-13** element (B, Al, ...) — an empty p orbital, which no bond count can reveal; or
- an atom carrying **fewer partners than its own nominal sigma valence**, i.e. a free
  coordination site: the amine oxide's 1-coordinate O, the N-ylide's 3-coordinate C.

It is a `max()` against the charge rule, never a replacement, and it applies only where the charge
rule decides — group 13 and the period >= 3 octet expansion already carry a larger cap. Both tests
read the corner's own bond list, so the cap stays a per-corner constant: no new chain-rule term,
and a change is carried by the existing s-blend.

Measured (dev against the share-off arm, kcal/mol; geometries gfn2-optimised):

| system | delivered | conserving, no donor rule | conserving + donor rule |
|---|---:|---:|---:|
| H3N-BH3 | +0.02 | +94.4 | **+0.00** |
| H3N-O (amine oxide) | +0.00 | +73.4 | **+0.00** |
| H3N-CH2 (N-ylide) | +0.00 | +109.5 | **+0.00** |
| H5O2+ (Zundel) | +137.7 | +120.5 | +120.5 (unchanged) |
| N2H7+ | +174.1 | +157.8 | +157.8 (unchanged) |

All three dative/ylide neutrals are bit-identical to the share-off energy with the rule on. The
two proton-shared dimers are not touched by it and stay 16-24 kcal/mol *closer* to the pinned
`gfnff` value than the delivered share is — the bridging H there genuinely is a 3c-2e case, which
is what the share exists for.

**Nothing else moves.** With the rule on, every falsifier is bit-identical to conserving without
it: the four class-C adducts, rkt06 (rms 2.7614), the six hypervalent ions, the 20-cell grid
(902 rebuilds / 216.62 kJ per-step max / 72 events >= 50 / 0 of 449 hard swaps / T_max 11139 K),
the equilibrium toggle set and `gfnff` itself. FD gradient at the new rule's own geometry
(H3N-BH3) is 1.7e-08 Eh/A.

**Scope**: sulfoxides and phosphine oxides need no donor rule — DMSO and Me3P=O are inert in every
arm (+0.00 kcal/mol), because the period >= 3 octet expansion already caps their donor. The rule
is not exercised by any metal (no rev-gfnff reference set contains one) and the proton-shared
dimers remain an open residual: their charge is split over two groups, so neither the charge rule
nor the donor rule grants them a full budget.

---

## 2. The bond-well form (`-gfnff.rev_well_form`, DEFAULT `mg3` since Sep 22, 2026)

> The default was `gauss` until Sep 19, 2026, `mg` from Sep 19 to Sep 22, and is **`mg3`** since.
> This section describes the `gauss` -> `mg` step (it is the one that introduced the MG well);
> section 2.3 describes `mg2`/`mg3` and the flip to `mg3`.

The delivered well is `k_b exp(-alpha (r - r0)^2)` times the reactive term weight. Against the
class-A r2SCAN-3c bond scans its break-side RMS is **19.2 kcal/mol** (median over 32 bond types):
a Gaussian has no tail.

Both alternatives are `E = -D (2y - y^2)`:

    mg        y = exp(-(a x + beta x^2)),                    a = sqrt(K / (2D))       closed form
    erfmorse  y = erfc((x - u)/sigma) / erfc(-u/sigma),      u from 2 D h(u/sigma)^2/sigma^2 = K

with `x = r - r0` and `K = 2 alpha |k_b|`, the delivered Gaussian's own force constant. Both have
`E(r0) = -D` and `E'(r0) = 0` identically and reproduce r_min and the curvature **by
construction**, so only the depth scale `s = D/|k_b|` and the tail are fitted — per element pair,
by `scripts/revgfnff_wellfit.py`, on the break side against a **charge-frozen** rest (the Coulomb
drift's sign flips between the static and the react protocol, so a depth fitted on the raw rest
absorbs a charge-model error that is not the well's). The table is
`src/core/energy_calculators/ff_methods/rev_well_table.h`; an element pair with no class-A data
keeps the Gaussian and says so at verbosity 2.

The term weight is **not** applied to these wells — they decay by themselves (fitted `abs(E_pair)`
at the last grid point: median 0.000, max 0.69 kcal/mol), so multiplying by it would truncate the
tail the fit just put there. The inner side is capped at `y = 2` (C1, exact below y = 1.6), which
bounds the well in `[-D, 0]` exactly as the Gaussian is in `[k_b, 0]`; the repulsive wall stays
the repulsion term's job.

| row | gauss | mg | erfmorse |
|---|---:|---:|---:|
| class-A median break rms (offline fit, per system) | 19.24 | **2.07** | **2.13** |
| class-A median rms (harness, per-element-pair table) | 24.49 | **19.50** | **19.21** |
| class-A median dev D_e / r90 | -25.53 / -0.330 | **-12.74 / -0.058** | **-13.13 / -0.068** |
| guard, pooled MAD over 167 conformer/S66 reactions | **1.0341** | 1.0439 | 1.0427 |
| max equilibrium bond-length shift | — | 0.0063 A | 0.0067 A |
| class D dE_MAD / grad_RMS | **4.980** / 16.599 | 5.152 / **16.376** | 5.184 / 16.494 |
| 20-cell grid, events >= 50 kJ / max dE_jump | 487 / 48.6 | 458 / 49.7 | 567 / **58.3** |
| wall time over the grid (setup cached) | 17.39 s | 17.03 s | 17.33 s |

**The two forms are indistinguishable on the data** (the fitted rms differs by at most 0.85 over
32 curves) and **MG is the cheaper one**: its `a` is a closed form, while erf-Morse needs a
bisection for `u`. That bisection is a per-bond SETUP cost and must be cached — before the cache
it cost **1.42x** the whole react-MD wall time against MG's 1.02x. That is why `mg` is the
default and `erfmorse` is not.

The table above is measured with each piece ALONE, i.e. `mg` against the then-default `delivered`
share; it was re-measured on the Sep 19 binary and reproduces row for row.

---

## 2.1 The react-MD tail — the share, not the combination (root-caused Sep 19, 2026)

> **Every time in this section is in the OLD time scale.** It was measured before the MD
> time-step unit fix of Sep 2026, so multiply every "fs" and every duration by **1.9516144** to
> read what was actually simulated: the "dt 0.25 / 0.125 / 0.0625 fs" scan below was really
> 0.488 / 0.244 / 0.122 fs and the "5 ps" cells were 9.76 ps. The mechanism and every relative
> statement are unaffected; the true-fs re-measurement is in section 2.2.

Packages 3 and 4 measured the conserving share and the MG well independently, and on the 20-cell
react-MD grid the two defaults together looked worse than either alone (per-step max 800.12 kJ/mol
against 216.62 and 396.88, T_max 28 392 K). **Package 7 shows that reading was a 20-cell sampling
artefact.** Over 130 cells of the same three systems (every frame, both temperatures, one frozen
binary) the heavy tail belongs to `conserving` with EITHER well form:

| arm (130 cells) | sum reb | median step / kJ | p90 | max / kJ | cells > 100 kJ | T_max / K |
|---|---:|---:|---:|---:|---:|---:|
| gauss + delivered | 9473 | 45.5 | 88.4 | 222.5 | 9 | 8 306 |
| gauss + conserving | 9022 | 42.9 | 116.2 | **1071.3** | 16 | **46 601** |
| mg + delivered | 8636 | 46.8 | 80.6 | 362.6 | 8 | 11 035 |
| **mg + conserving (DEFAULT)** | 9708 | 49.4 | 213.5 | 800.1 | **27** | 28 392 |

`gauss + conserving` reaches 1071.3 kJ/mol on `c2h6/T2000_f7`, a cell the 20-cell grid does not
contain. MG adds **frequency**, not mechanism: 27 cells above 100 kJ against 16 (on c2h6 at 2000 K,
20/25 against 9/25), while its contribution to the mechanism is ~0.6 kJ/mol (below).

**Mechanism**, measured per bond at the worst cell's own event geometry. A transient geminal H2
forms inside the molecule (both hydrogens still on the same carbon). The delivered left-over rule
gives that pair `c = 0.000000` exactly; the conserving rule gives it `c = f_H f_H = 0.5031^2 =
0.2531` on a 0.204 Eh well, i.e. **-136 kJ/mol**. Both rules give the two C-H bonds the same
`c ~ 0.50`, so that one pair is the whole difference. Bond-term energy of the bridged corner
relative to the unbridged one, same geometry: **+124.4** (gauss+delivered) / **+128.6**
(mg+delivered) / **-10.6** (gauss+conserving) / **-11.2** (mg+conserving) kJ/mol. Under
`delivered` the artefact costs 124 kJ/mol and the dynamics is pushed out of it; under `conserving`
it is free, so the molecule visits it 1.5-2.8x as often (125 H-H formations against 52 on 25 c2h6
cells).

**The jump is a blend-window resolution failure, not a step in the potential.** `CURCUMA_BLENDDUMP=1`
shows the break transition of that H-H pair moving `s` from 0.000000 to 0.515777 — half its window —
in ONE 0.25 fs step, for a distance change of 0.09 a0, because the window is a fixed interval in the
bo3 ORDER and the bo3 switch is steepest for the smallest covalent sum: for H-H it is only ~0.175 a0
wide in distance, about two steps for a hydrogen at 2000 K. The corner gap it has to carry is a
median 153 and up to 390 kJ/mol. Halving the time step removes the event completely
(`c2h6/T2000_f0`: 800.1 -> 38.7 -> 13.5 kJ/mol and T_max 28 392 -> 5 579 -> 5 441 K at
dt = 0.25 / 0.125 / 0.0625 fs), and the static potential along the same path changes by only
+5.9 kJ/mol where the MD jumps by +393 — a discontinuity would survive a smaller step.

Context: `-gfnff.rev_valence_share false` on the same 130 cells is median **744.4** kJ/mol, max
7738.5, 87/130 cells above 100 — the share of either form is worth a factor ~15, so this is a choice
between two second-order failure modes. The two rules fail on disjoint motifs: every worst cell of
the delivered arms is `ch4_H` (the artificial radical adduct, 282 -> 20 H-H formations when the
share is flipped), every worst cell of the conserving arms is `c2h6` (the geminal H2). Full
measurement, the falsified hypotheses and a costed option list: `WORK_STATUS.md` package 7.

**Everything else is unaffected by the combination** and equals the MG-alone row: guard 1.0439,
class D 5.152 / 16.376, rkt06 2.72, the equilibrium toggle set 20/20 at dE = 0.000000000,
`gfnff` bit-identical, FD gradients 1.4e-08 to 1.24e-07 Eh/A over five geometries. The four
class-C adducts re-confirm package 4.3's strongest claim — **the conserving share fully
compensates the deeper MG well where the delivered share does not**: dev min / rms
-1.5 / 11.8, -1.4 / 5.1, -3.0 / 13.3, +0.0 / 2.9 under the default, against -89.4 / 43.1,
-110.1 / 48.7, -88.5 / 39.0, -93.6 / 42.8 for `mg + delivered`.

Class A moves **19.50 -> 20.52** median rms (dev D_e -12.74 -> -16.09, r90 -0.058 -> -0.120) under
the combination, and that is **two bond types out of 32**: `ncl3_N-Cl` 20.10 -> 23.57 (worse) and
`hocl_O-Cl` 35.86 -> 33.30 (better), the other 30 bit-identical; the mean rms moves 22.48 ->
22.51. It is the conserving share, not the well form and not the donor rule (identical with
`-gfnff.rev_share_donor_rule false`). The median is a fragile statistic under one large mover.

### 2.2 Recommended MD settings (operating recommendation, not a default change)

> **Every time in this section is a REAL femtosecond.** Sections 2.1 and 2.3 above, and everything
> packages 1 and 6-9 measured, were taken **before** the MD time-step unit fix of Sep 2026
> (`CurcumaUnit::Constants::MD_TIME_UNIT_FS`, `AIChangelog`): `SimpleMD` handed the requested step
> to the integrator without converting femtoseconds into the integrator's own time unit
> `sqrt(amu*A^2/Eh) = 1.9516144 fs`, so **every nominal "fs" in those sections is really
> 1.9516144 fs** and every duration is stretched by the same factor. Section 2.1's dt scan
> "0.25 / 0.125 / 0.0625" was really 0.488 / 0.244 / 0.122 fs, and its "5 ps" cells were 9.76 ps.
> The numbers below replace the recommendation in true femtoseconds.

**Use `-md.time_step 0.0625` for quantitative react-mode MD with the conserving share.** This is
still mitigation option 2 of `WORK_STATUS.md` 7.8, adopted by the operator on 2026-09-20 as a
**recommendation only**: no default value changed, and the transition-window code of section 2.1
was deliberately left alone (redefining the window in distance instead of bond order is a separate,
deferred decision). Only the number changed, and it was re-measured rather than rescaled.

Why 0.0625 fs: the jump of section 2.1 is a resolution failure of the integrator, not a step in the
potential, so a smaller step removes it. The old recommendation of "0.125" was calibrated on a
SINGLE cell and was really 0.244 fs. Re-measured on package 7's **full 130 cells** (c2h6 / ch3nh2 /
ch4_H, every frame, 1000 and 2000 K, 9.758 ps each = the same physical exposure as package 7,
shipped defaults `conserving` + `mg` + donor rule, `-threads 1`, one frozen binary; per-step
|dEpot| with rebuild-containing intervals excluded):

| true `-md.time_step` / fs | sum rebuilds | median / kJ | p90 | max | cells > 100 kJ | cells > 200 kJ | T_max / K |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 1.0 | 9209 | 97.8 | 347.1 | 1133.1 | 64 | 31 | 5.5e7 |
| 0.5 | 10475 | 51.8 | 218.4 | 453.1 | 25 | 16 | 18 928 |
| 0.25 (the effective default, see below) | 9677 | 27.5 | 64.8 | 483.0 | 11 | 9 | 11 948 |
| 0.125 | 13597 | 13.1 | 61.5 | 204.2 | 10 | 1 | 10 405 |
| **0.0625 (recommended)** | 20922 | **6.3** | **16.1** | **78.5** | **0** | **0** | 9 321 |

0.0625 fs is the first true step at which **no** cell exceeds 100 kJ/mol in a single step. It is
also the old threshold rescaled (0.125 / 1.9516 = 0.0640), so the arithmetic prior and the
measurement agree — but the measurement is what the recommendation rests on, because the old value
was a one-cell number. Cost: 4x wall time per picosecond relative to the 0.25 fs cap. No falsifier
moves (the equilibrium, class-A, class-C and gradient checks are time-step independent by
construction), and nothing about the model changes.

**The effective default is 0.25 fs, and that is now a genuine 0.25 fs** (before the unit fix the
same setting integrated 0.488 fs). `-md.time_step` itself defaults to 1.0 fs, but `-method
revgfnff` clamps it to `-md.rev_dt_cap` (default **0.25 fs**, stage 1,
`docs/REV_GFNFF_STAGE1.md`) with a warning. The clock fix therefore improved the default by about
a factor of two in every robust statistic for free — on the same 130 cells the effective default
went from median 51.8 / 25 cells above 100 kJ / T_max 18 928 K (what 0.25 used to integrate) to
27.5 / 11 / 11 948 — but 0.25 fs still sits **outside** the band where the sample is bounded, so
the warning still fires for a plain react run with no flags. `rev_dt_cap`'s default was **not**
changed here; whether to lower it is an operator decision.

`SimpleMD` prints a one-time startup warning (Claude Generated, Sep 2026,
`simplemd.cpp::Initialise`) when **all** of these hold: `-method revgfnff` (or `gfnff-rev`),
`-gfnff.topology_mode react`, `-gfnff.rev_share_form conserving` (the default), and an effective
time step above **0.0625 fs**. The gate is the share form alone — the well form is not the
mechanism (section 2.1) — so `-gfnff.rev_share_form delivered` does not warn with any well form,
and neither does a non-react run. Verified: fires at 0.25 and 0.125, silent at 0.0625 and 0.05,
silent for `delivered` and for plain `gfnff`. The warning is advisory: the spikes are smooth and
bounded, the thermostat recovers, hard swaps stay 0 and every rebuild reports `dE_jump 0.000000`.

**The worst cell of section 2.1 at the corrected clock.** `c2h6/T2000_f0` with the shipped
defaults, the SAME discrete dynamics as package 7 (CSVR per-step ratio 0.025, COM removal every
400 steps, 9.758 ps), varying only the physical step length:

| true dt / fs | rebuilds | per-step max / kJ/mol | steps >= 50 kJ/mol | T_max / K |
|---:|---:|---:|---:|---:|
| 0.4879 (what package 7 called "0.25") | 76 | 800.1 | 62 | 28 392 |
| 0.25 | 102 | 47.8 | 0 | 5 206 |
| 0.125 | 170 | 20.1 | 0 | 5 520 |
| 0.0625 | 82 | 13.2 | 0 | 5 602 |

The first row is reproduced **exactly** on the post-fix binary by scaling every dt-derived setting
by 1.9516144204 (`-md.rev_dt_cap 0 -md.time_step 0.4879036051 -maxtime 9758.072102 -md.coupling
19.516144204 -md.remove_com_motion 195.16144204`), which is the proof that the unit fix changed the
clock and nothing else.

---

## 2.3 Free curvature, r0 re-solve and the bond-order dimension (`mg3` is the DEFAULT since Sep 22, 2026)

Stage 3a(iii) step 2 and stage 3b. Delivered opt-in on Sep 20, 2026; **`mg3` became the default on
Sep 22, 2026** (operator decision, `test_cases/revgfnff/_log/WORK_STATUS.md` package 12).
Same MG well, but `x = r - (r0_model + dr0)` and four parameters per key instead of two:

    s    = D / |k_b|                   depth scale      (fitted)
    beta = the MG tail                                  (fitted)
    ca   = a / sqrt(alpha / s)         curvature scale, K_well = ca^2 * 2 alpha |k_b|   (SOLVED)
    dr0  = offset on the model's own dynamic r0 [A]                                     (SOLVED)

`ca = 1, dr0 = 0` is exactly `mg`, bit-for-bit — the verification the code path was landed with.
`mg2` is keyed on the element pair, `mg3` additionally on the CONTINUOUS bond order
`1 + pibo * (both ends sp ? 2 : 1)` (`Bond::rev_order`), interpolated linearly between the
single/double/triple fits of that pair. The order is a topology constant inside one energy call, so
it adds no geometry derivative and no switch; benzene's C-C reads 1.666 and gets a well between the
single and the double fit. Table: `rev_well_table_v2.h`, generated by the same fit script.

**Neither the curvature nor the offset is fitted, and that is a measured decision.** A free `ca` on
the break-side objective runs to zero (a flat, quartic well bottom) on 8 of 32 bonds; on the full
grid the five compressed points, where the repulsion carries the error, drag the well off. `ca` is
therefore solved in closed form from the measured curvatures, `ca^2 = 1 + (k_ref - k_model)/K_gauss`.
`dr0` solved against the class-A rigid-scan minimum moved the optimised water O-H to 0.9179 A
(reference 0.9644); it is solved against the **relaxed** bond length instead.

| row | gauss | mg (default) | mg2 | mg3 |
|---|---:|---:|---:|---:|
| class-A harness rms, median | 24.68 | 22.15 | **15.83** | **13.22** |
| class-A dev D_e / dev r90 / dev k, median | -25.79 / -0.341 / -169 | -14.32 / -0.066 / -158 | -16.17 / -0.101 / **-44** | **-6.86** / -0.131 / -82 |
| median \|b_model - b_r2SCAN-3c\| over 32 class-A bonds [A] | 0.0293 | 0.0246 | **0.0060** | **0.0044** |
| guard, pooled MAD (167 reactions) | **1.0341** | 1.0439 | 1.0543 | 1.0547 |
| max equilibrium bond-length shift vs gauss (4 molecules) | — | 0.0063 A | 0.0212 A | 0.0212 A |
| class D dE_MAD / grad_RMS (20 systems, mean) | 4.980 / 14.575 | 5.152 / 14.428 | **4.590 / 11.397** | **4.555 / 11.313** |
| rkt06 path rms | 2.3840 | 2.3788 | **2.2665** | **2.2665** |
| FD gradient, worst of 4 points | — | 1.02e-07 | 1.14e-07 | 1.14e-07 |
| 130-cell grid: median step / p90 / max / cells > 100 kJ | 44.1 / 154.3 / 417.2 / 19 | 50.3 / 215.7 / 446.8 / 27 | 51.2 / 214.9 / 432.6 / 21 | 47.4 / 254.2 / 441.3 / 24 |
| 130-cell grid: max dE_jump / n >= 50 | 259.4 / 35 | 272.3 / 14 | 359.8 / 26 | **456.3 / 89** |
| hard swaps (20 cells, CURCUMA_VERB=3) | — | 0 / 623 | 0 / 497 | 0 / 469 |

The free curvature buys the INNER branch (dev k -158 -> -44) and the bond-order split buys the
depth (dev D_e -14.32 -> -6.86). The equilibrium bond lengths move more than `mg`'s **towards the
reference**: the median error against r2SCAN-3c falls by a factor 5-7.

> **The two 130-cell rows above are in the OLD time scale** (see the box in section 2.2): their
> "dt 0.25 fs" was really 0.488 fs. Re-measured in TRUE femtoseconds (package 10, same 130 cells,
> 9.758 ps each, one frozen binary), **the smoothness ordering does not survive** — `mg3` is no
> longer the worst arm, and the absolute severity collapses for all three:
>
> | arm | rebuilds | max dE_jump / kJ | n(jump) >= 50 | rate per 1000 rebuilds | step median | step max | cells > 100 kJ |
> |---|---:|---:|---:|---:|---:|---:|---:|
> | **old clock (real 0.488 fs)** ||||||||
> | mg | 14551 | 304.4 | 24 | 1.65 | 49.4 | 800.1 | 27 |
> | mg2 | 12882 | 352.8 | 37 | 2.87 | 48.3 | 460.9 | 27 |
> | mg3 | 12459 | 440.5 | 72 | **5.78** | 47.3 | 700.1 | 22 |
> | **true 0.25 fs (the effective default)** ||||||||
> | mg | 14506 | 260.5 | 9 | 0.62 | 27.5 | 483.0 | 11 |
> | mg2 | 11801 | 300.6 | 21 | **1.78** | 25.2 | 447.8 | 14 |
> | mg3 | 11780 | 283.5 | 16 | 1.36 | 25.7 | **316.5** | **9** |
> | **true 0.0625 fs (the recommended point)** ||||||||
> | mg | 31371 | 317.2 | 11 | 0.35 | 6.3 | 78.5 | 0 |
> | mg2 | 23695 | 254.4 | 2 | **0.08** | 6.7 | 81.6 | 0 |
> | mg3 | 26712 | 354.0 | 6 | 0.22 | 6.9 | 79.5 | 0 |
>
> At the corrected clock `mg3`'s jump-event rate falls from 3.5x `mg`'s to 2.2x and it swaps places
> with `mg2`; on the per-step statistics `mg3` is the **best** of the three at true 0.25 fs
> (max 316.5 vs 483.0/447.8, 9 cells above 100 kJ vs 11/14). At the recommended 0.0625 fs all three
> have **zero** cells above 100 kJ/mol per step and the ordering is `mg2 < mg3 < mg`, i.e. noise.
> **The "mg3 has a clearly worse tail" argument was an artefact of the too-coarse clock and should
> not weigh in the adoption decision.** The class-A / guard / class-D rows above are time-step
> independent and are unaffected.

**The two costs, stated plainly** (in the old time scale, see the box): the conformer/S66 guard
opens from 1.0341 to 1.0543/1.0547,
twice the +0.0098 the `mg` flip cost; and `mg3`'s rebuild `dE_jump` tail grows to 89 events above
50 kJ/mol against `mg`'s 14. The tail is **not** the new dimension — over 259 rebuilds in five
independent cells, 0 of the 10 rebuilds above 50 kJ/mol involves any change of `rev_order`, while
the 4 rebuilds that do move an order carry at most 0.1 kJ/mol. Forcing the order with the new
diagnostic `-gfnff.rev_well_order_override` and scanning 1.00 -> 3.00 gives a strictly linear
response (largest step over d(order) = 0.05 is 0.7930 kcal/mol against a linear prediction of
0.7929, `|dE/d(order)| <= 15.9` kcal/mol per unit order).

Every number: `test_cases/revgfnff/_log/WORK_STATUS.md` packages 9a / 9b, including the rejected
r0-solve variant and the `dr0 = 0` variant (guard 1.0508/1.0490, but the relaxed bond lengths stay
0.020 A from the reference).

### 2.3.1 The default flip to `mg3` (Sep 22, 2026) and what the alternatives are for

The operator made `mg3` the default on the package-9 accuracy rows, which are time-step
independent, once package 11 had removed the smoothness argument from both sides: over **11 700
paired-replicate react-MD trajectories (780 per arm and per true time step)** neither `mg2` nor
`mg3` is distinguishable from `mg` on the per-step |dEpot| tail at any time step, and the single
trajectory that had been read as "mg3's tail is clearly worse" is a singleton (same cell and arm,
six replicates: 256.8 / 43.4 / 35.8 / 49.9 / 44.7 / 38.5 kJ). What remains is the package-9
accuracy/guard trade, restated for the flip that was actually made, `mg` -> `mg3`:

| row | `mg` (old default) | `mg3` (new default) |
|---|---:|---:|
| class-A harness rms, median | 22.15 | **13.22** |
| class-A dev D_e, median [kcal/mol] | -14.32 | **-6.86** |
| median \|b_model - b_r2SCAN-3c\| over 32 class-A bonds [A] | 0.0246 | **0.0044** |
| class D dE_MAD / grad_RMS (20 systems, mean) | 5.152 / 14.428 | **4.555 / 11.313** |
| rkt06 path rms (conserving / delivered) | 2.3788 / 2.3447 | **2.2665 / 2.2599** |
| guard, pooled MAD over 167 conformer/S66 reactions | **1.0439** | 1.0547 |
| max equilibrium bond-length shift vs `gauss` (4 molecules) | **0.0063 A** | 0.0212 A |

The guard is the one row that gets worse, by +0.011 kcal/mol on a 1.04 kcal/mol MAD; the
equilibrium shift is larger but moves **towards** the reference (the row above it). Every
alternative stays available:

- **`mg`** — the Sep 19 - Sep 22 default, the two-parameter MG well with the curvature and the
  minimum pinned to the delivered Gaussian's. Keep it to reproduce anything measured in that
  window, or when the +0.011 kcal/mol guard cost matters more than the class-A/class-D gain.
- **`mg2`** — the same four-parameter family as `mg3` but keyed on the element pair only, without
  the bond-order dimension. It costs the same guard (1.0543) and is equal or worse than `mg3` on
  every package-9 row, so it is not preferred; it exists to separate "free curvature + r0" from
  "bond order" when attributing a change.
- **`erfmorse`** — the erf-Morse alternative of the same curvature-pinned family as `mg`;
  indistinguishable from `mg` on the data (fitted rms differs by at most 0.85 over 32 curves) but
  more expensive (a per-bond bisection for `u` instead of a closed form). Kept as the independent
  check that the MG functional form is not itself doing the work.
- **`gauss`** — the delivered GFN-FF Gaussian, i.e. the state before stage 3a(iii). It has the
  best guard (1.0341) and no tail at all on the break side (class-A median rms 24.68). Use it to
  reproduce pre-Sep-19 behaviour or to attribute anything to the well form as a whole.

**Open against this default, DIAGNOSED (package 13) — not a defect, the test's statistic is
misleading on this bath.** `cli_simplemd_18_gfnff_rev_nve_vs_gfnff` fails on its dt = 0.125 arm
at the `mg3` default (|slope| 2.95e-3 against a 2.5e-3 floor); the measurement reproduces (also
independently by the orchestrator on two fresh replicates). But `Etot(t)` on this 12-H2/8000K NVE
bath is NOT a drift ramp — it is a step function: one H2 dissociates within the first ~1 ps,
injecting an amount of energy set by its bond well's depth, and the trajectory is then FLAT
(late-window slope consistent with 0) for the remaining 9 ps. The test's fixed-window OLS slope
is therefore a one-time plateau HEIGHT divided by the window, not an ongoing rate, and `mg3`'s
deeper H-H well (`D_e` larger by 0.0141 Eh, 98 % of the observed amplitude difference) makes that
single release bigger — exactly reproducing the "1.6-1.7x higher slope" without any higher
dissipation. Measured separately (dt-scaling, n = 6 replicates/cell), the actual energy-injection
rate scales as dt^2 and is 19-34 % LOWER for `mg3` than `mg` at every dt below 0.25, never higher.
The 2.5e-3 floor also fails the DELIVERED `gauss` well at small dt (1.21e-2 at 0.03125 fs, the
worst of all four arms) and is non-monotone in dt for every arm — it is not diagnosing `mg3`
specifically, or well-posed for any of the new forms. No source change; the test still carries the
old measurement and awaits an operator decision on a better statistic, not a defect fix.
Full detail: `test_cases/revgfnff/_log/MG3_DISSIPATION_STATUS.md`.

---

## 3. What was tested, what was not, what is not implemented

**Tested** (all on the reference sets under `test_cases/revgfnff/ref/`, the GMTKN55 conformer/S66
guard and a 20-cell react-MD grid of c2h6 / ch3nh2 / ch4+H at 1000 and 2000 K): single-point
energies to 12 digits, analytic gradients against central finite differences at four geometries
(equilibrium, an exchange TS, a react-MD runaway frame and a radical-approach point), NVE-relevant
per-step energy continuity, the class-A bond scans, class-D MD frames against r2SCAN-3c, and the
class-C radical approaches.

**NOT tested**: anything containing a metal (the conserving share's element rules do not cover the
d block — a transition metal deliberately keeps the delivered growth, and that carve-out is not
measured; the donor rule is likewise unexercised by any metal); periodic systems; charged species
beyond the six hypervalent ions and the two proton-shared dimers; any well form on a
hydrogen-bonded X-H (the HB alpha modulation is not applied in the new forms, and the class-A set
contains no hydrogen bond); long-time MD stability beyond 5 ps per cell; every combination of the
flags except the ones tabulated above.

**OPEN, measured**:
- the react-MD tail of section 2.1 — **root-caused** (the conserving share makes a transient
  geminal H2 free instead of +124 kJ/mol unfavourable, and the stage-1b break window for an H-H
  pair is only ~2 time steps wide at dt = 0.25 fs), **not fixed**: every remedy is a redesign of
  either the share rule or the transition window. `dt = 0.125 fs` removes it and is now the
  documented recommendation (section 2.2) plus a startup warning; the window redesign is deferred.
  Five costed options in `WORK_STATUS.md` 7.8;
- `ncl3_N-Cl`, the one class-A bond type the conserving share makes worse (rms 20.10 -> 23.57);
- the proton-shared dimers H5O2+ and N2H7+, where the charge is split over two groups so that
  neither the charge rule nor the donor rule grants a full budget — both modes stay 120-158
  kcal/mol from the share-off arm, and `conserving` is only the smaller of the two errors;
- the compressed-BF4- probe, where the share costs +33.4 kcal/mol under `mg` against +18.6 under
  `gauss`. That geometry is documented as a perception question, not a share question
  (`FABLE_BOND_STATE.md`), and is not a regression gate. **Superseded (Sep 29, 2026)**: the
  "+33.4/+18.6" figures were stale (`FABLE_BOND_STATE_2.md` section 0 re-measured the current
  shipped default at +343.6 kcal/mol vs share-off, +895.7 vs the corner's own 4-bond evaluation) -
  and the perception question itself is now RESOLVED, not merely documented: see below.

### 2.4 The pair-validity gate (Sep 29, 2026) resolves the compressed-BF4- perception question

`FABLE_BOND_STATE_2.md` section 2.1-rev's VALID(i,j,b) rule - a per-corner, geometry-free test of
whether a perceived bonded pair is a genuine bond or a closed-shell repulsion (BF4-'s six spurious
F...F contacts; a hot react-MD trajectory's geminal H...H "bond") - is implemented as
`-gfnff.rev_pair_validity` (default off). An invalid pair's corner is regenerated from the reduced
bond list (hybridisation/angles/torsions/pi-systems/EEQ all freshly derived), reusing the existing
`m_forced_bonds` corner-forcing seam every react-mode corner already goes through - no new
architecture needed. With the gate on, the compressed-BF4- probe reproduces its own 4-bond
topology's energy **bit-for-bit** (`-1.25025415` Eh), and over the full reference set (GMTKN55
2462 + MOR41 285 + S30L-CI 90 structures) exactly the 5 structures the offline sweep predicted
move (`PX13/hf_2_ts` + 4 named `MB16-43` clusters), nothing else - bit-identical with the flag
off. Full acceptance table, the react-mode MD smoothness measurement (max per-step |dEpot| 70.5 ->
14.3 kJ/mol on the extracted geminal-H...H frame) and a bug found+fixed on a live react-mode
corner: `test_cases/revgfnff/_log/PAIR_VALIDITY_IMPL_STATUS.md`. Implementation:
`src/core/energy_calculators/ff_methods/gfnff_pair_validity.cpp`.

**Implemented since Sep 20, 2026, OPT-IN (section 2.3)**: the second step of stage 3a(iii)
(`rev_well_form mg2`, free curvature + r0 re-solve) and the **bond-order-resolved well table**
(`mg3`). The shipped default `mg` is still keyed on the element pair, so under it C-C, C=C and C#C
share one median (per-system `s` 1.18 / 1.01 / 0.91) — that remains true until the operator flips a
default.

**NOT implemented / not covered by mg2/mg3**: the HB alpha modulation in any new form; a second
Newton iteration of the r0 re-solve (residual median 0.0060 A); anything the class-A set does not
contain (metals, hydrogen bonds, charged species, periodic systems); and the C-H / H-O keys still
average two systems each, which is a HYBRIDISATION difference (sp3 vs sp C-H, water vs methanol
O-H) that a bond-ORDER dimension cannot separate.

### 2.5 Q5, the hydrogen-perception rule set (Sep 29, 2026, opt-in)

`FABLE_BOND_STATE_2.md` section 7/7.6 ("Q5"): a bridging hydrogen (an X-H-Y 3c-4e/3c-2e bond, a
migrating H) is given `hyb=1` ("sp") by the topology perception, which lets it trigger four rules
meant for a genuine sp centre and not for an exchange intermediate - `-gfnff.rev_h_scope`
(default off) removes them, one PARAM per rule (`rev_h_scope_h1/h2/r1`, all default true, read
only when the master is on):
- **H1**: every atom with `Z==1` reads the terminal-H bond-strength column `bsmat[hyb_X][0]`
  instead of the "sp" column, whatever its partner count (tests the element, never the periodic
  group, so Li/Na/K are untouched). Implemented narrowly, at the bond-strength lookup and the
  `is_bridge` 0.30-scaling gate only (`getGFNFFBondParameters`) - NOT as a blanket override of the
  hybridisation array, which was tried first and reverted: it also changes the r0 "shift"
  correction (an unrelated hybridisation-keyed term), moving the exactly-collinear rkt06 H+H2 path
  by ~0.013 Eh where the falsifier requires 0.00. Subsumes the narrower, pre-existing
  `-gfnff.rev_h_not_sp` for real hydrogen (kept, unchanged, as the mechanism reachable when
  `rev_h_scope_h1` is off).
- **H2**: no GFN-FF angle term is ever centred on a `Z==1` atom (`generateAnglesNative`). Needed
  together with H1, which was tried alone once before and rejected for exactly this reason (a
  spurious tetrahedral `theta0`; see the note at `determineHybridizationFortran`'s return).
- **R1**: ring enumeration excludes every `Z==1` atom entirely, on an H-free copy of the metal-free
  adjacency list fed to `findSmallestRings` - clears `ringf` AND the 3-ring `fxh` correction on
  every bond of a bridged carbon at once (bridging and sibling C-H alike), since `fxh` is keyed on
  the carbon's own ring membership.
- **P1**: a `Z==1` atom never counts as an sp/sp2 "picon" neighbour (`detectPiSystems`); a
  structural consequence of H1 that still needs its own explicit skip since H1 does not touch the
  hybridisation array (no separate PARAM - there is no independent code path to gate).

Bit-identical when off (GMTKN55 2462 + MOR41 285 + S30L-CI 90 structures, 0 moved); on, exactly
the 80 structures with a genuinely 2+-coordinate hydrogen in the EVALUATED (pass-2, never pass-1)
topology move (79 + `MB16-43/34`, a mu3-bridging hydride where the rule is a near-total no-op).
FHF- De -120.6 -> -77.6 kcal/mol (predicted -76.9); the exactly-collinear rkt06 path is 0.00 at
every point measured, both the symmetric TS and an asymmetric point (H-H bond strength is a pure
`Z==1&&Z==1` check, independent of hybridisation). One of the seven pass-1-only structures the
design's own offline verification named as required-bit-identical genuinely is NOT
(`WATER27/OHmH2O`) - root-caused to the unrelated `-gfnff.frag_charge_model ensemble` default
(Sep 24, 2026), which spawns independent sub-`GFNFF` topology evaluations that Q5 correctly and
consistently applies to as well; with `-gfnff.frag_charge_model reference` that structure is
bit-identical too. Full acceptance table, the CH5+/B2H6 probe numbers and the sweep methodology:
`test_cases/revgfnff/_log/H_SCOPE_IMPL_STATUS.md`. Implementation: three call sites in
`gfnff_method.cpp` (`determineHybridizationFortran`, `generateAnglesNative`,
`calculateTopologyInfoOnce`'s ring-building block, `detectPiSystems`); PARAMs and `RevSettings`
fields in `gfnff.h`/`ff_workspace.h`.

---

## 4. Files

| what | where |
|---|---|
| the share and the well kernels | `src/core/energy_calculators/ff_methods/ff_workspace_gfnff.cpp` (`calcBonds`, `prepareValenceShare`, `prepareConservingShare`, `prepareWellForms`) |
| the settings struct and the smooth helpers | `src/core/energy_calculators/ff_methods/ff_workspace.h` (`RevSettings`, `shareClip`, `shareExcess`, `shareMinOne`) |
| the flags | `src/core/energy_calculators/ff_methods/gfnff.h` (PARAM block), read in `gfnff_method.cpp::setupRevSettings` |
| the AI-fitted well table | `src/core/energy_calculators/ff_methods/rev_well_table.h` (generated, `mg`/`erfmorse`) |
| the free-curvature / bond-order table | `src/core/energy_calculators/ff_methods/rev_well_table_v2.h` (generated, `mg2`/`mg3`) |
| the fit | `scripts/revgfnff_wellfit.py` |
| the class-A harness | `scripts/revgfnff_classa.py` (`--mode kept --extra "-gfnff.topology_mode react"`) |
| diagnostics | `CURCUMA_SHAREDUMP=1` (per-pair `share`/`shareD` rows and, in `conserving`, a per-atom `shareA` row whose `cap` column shows the donor grant), `CURCUMA_BLENDDUMP=1` (per energy call a `blendD` row per stage-1b transition: pair, forming/tight, window `[w_a, w_b]`, r, coordinate c, corner weight s), `CURCUMA_REVDUMP=1` (the resolved settings, incl. `share_form` / `share_donor_rule` / `well_form`), `CURCUMA_BONDDUMP=1` (per perceived bond the dynamic r0/fc/alpha/fqq/CN **and its continuous bond order**), `-gfnff.rev_well_order_override` (forces that order so `dE/d(order)` is measurable; `mg3` only, off by default) |
| measurements | `test_cases/revgfnff/_log/WORK_STATUS.md`, `HBUDGET_STATUS.md`, `RUNAWAY_STATUS.md`, `FABLE_REVIEW_2.md` |
| regression tests | `cli_simplemd_20_gfnff_rev_h_budget`, `cli_gfnff_03_rev_adduct_falsifier`, `cli_gfnff_04_rev_well_form`, `cli_gfnff_07_pair_validity_gate`, `cli_gfnff_08_h_scope` |
