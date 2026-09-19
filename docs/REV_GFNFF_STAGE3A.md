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
| 3a(iii) | the SHAPE of the well | `-gfnff.rev_well_form` | **`mg`** (was `gauss`) |

**Absolute energies moved with the `mg` flip and relative ones did not.** An MG well carries its
own fitted depth `D = s |k_b|`, so `revgfnff` is now further from plain `gfnff` in absolute terms
— caffeine -4.546943898047 against `gauss` -4.673521653477 and `gfnff` -4.672737068614, and
-133.6 kcal/mol on the acetic-acid dimer of `cli_simplemd_19` against -0.56 before. The pooled
167-reaction conformer/S66 guard moves only 1.0341 -> 1.0439 kcal/mol MAD. Do not compare an
absolute `revgfnff` energy with an absolute `gfnff` one.

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

## 2. The bond-well form (`-gfnff.rev_well_form`, DEFAULT `mg` since Sep 19, 2026)

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

**Use `-md.time_step 0.125` for quantitative react-mode MD with the conserving share.** This is the
mitigation option 2 of `WORK_STATUS.md` 7.8, adopted by the operator on 2026-09-20 as a
**recommendation only**: no default value changed, and the transition-window code of section 2.1
was deliberately left alone (redefining the window in distance instead of bond order is a separate,
deferred decision).

Why 0.125 fs: the jump of section 2.1 is a resolution failure of the integrator, not a step in the
potential, so halving the step removes it. Measured on the worst known cell, `c2h6/T2000_f0` with
the shipped defaults (`conserving` + `mg`), 5 ps, `-threads 1`:

| `-md.time_step` / fs | rebuilds | per-step max / kJ/mol | steps >= 50 kJ/mol | T_max / K |
|---:|---:|---:|---:|---:|
| 0.25 (the effective default, see below) | 76 | **800.1** | 62 | **28 392** |
| **0.125 (recommended)** | 102 | **38.7** | **0** | 5 579 |
| 0.0625 | 44 | 13.5 | 0 | 5 441 |

Cost: 2x wall time per picosecond. No falsifier moves (the equilibrium, class-A, class-C and
gradient checks are time-step independent by construction), and nothing about the model changes.

**The effective default is 0.25 fs, so this applies to a plain react run with no flags.**
`-md.time_step` itself defaults to 1.0 fs, but `-method revgfnff` clamps it to `-md.rev_dt_cap`
(default **0.25 fs**, stage 1, `docs/REV_GFNFF_STAGE1.md`) with a warning. 0.25 fs is exactly the
step at which the tail of section 2.1 lives.

`SimpleMD` therefore prints a one-time startup warning (Claude Generated, Sep 2026,
`simplemd.cpp::Initialise`) when **all** of these hold: `-method revgfnff` (or `gfnff-rev`),
`-gfnff.topology_mode react`, `-gfnff.rev_share_form conserving` (the default), and an effective
time step above 0.125 fs. The gate is the share form alone — the well form is not the mechanism
(section 2.1) — so `-gfnff.rev_share_form delivered` does not warn with any well form, and neither
does a non-react run. The warning is advisory: the spikes are smooth and bounded, the thermostat
recovers, hard swaps stay 0 and every rebuild reports `dE_jump 0.000000`.

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
  (`FABLE_BOND_STATE.md`), and is not a regression gate.

**NOT implemented** relative to the plan: the second step of stage 3a(iii) (freeing the curvature
with an r0 re-solve); a **bond-order-resolved well table** — the current one is keyed on the
element pair, so C-C, C=C and C#C share one median (their per-system fits differ by `s`
1.18 / 1.01 / 0.91), and that is exactly why the class-A harness median lands at 19.50 rather than
at the 2.07 a per-system fit reaches. That stage-3b work is the remaining half of the MG flip's
benefit.

---

## 4. Files

| what | where |
|---|---|
| the share and the well kernels | `src/core/energy_calculators/ff_methods/ff_workspace_gfnff.cpp` (`calcBonds`, `prepareValenceShare`, `prepareConservingShare`, `prepareWellForms`) |
| the settings struct and the smooth helpers | `src/core/energy_calculators/ff_methods/ff_workspace.h` (`RevSettings`, `shareClip`, `shareExcess`, `shareMinOne`) |
| the flags | `src/core/energy_calculators/ff_methods/gfnff.h` (PARAM block), read in `gfnff_method.cpp::setupRevSettings` |
| the AI-fitted well table | `src/core/energy_calculators/ff_methods/rev_well_table.h` (generated) |
| the fit | `scripts/revgfnff_wellfit.py` |
| the class-A harness | `scripts/revgfnff_classa.py` (`--mode kept --extra "-gfnff.topology_mode react"`) |
| diagnostics | `CURCUMA_SHAREDUMP=1` (per-pair `share`/`shareD` rows and, in `conserving`, a per-atom `shareA` row whose `cap` column shows the donor grant), `CURCUMA_BLENDDUMP=1` (per energy call a `blendD` row per stage-1b transition: pair, forming/tight, window `[w_a, w_b]`, r, coordinate c, corner weight s), `CURCUMA_REVDUMP=1` (the resolved settings, incl. `share_form` / `share_donor_rule` / `well_form`) |
| measurements | `test_cases/revgfnff/_log/WORK_STATUS.md`, `HBUDGET_STATUS.md`, `RUNAWAY_STATUS.md`, `FABLE_REVIEW_2.md` |
| regression tests | `cli_simplemd_20_gfnff_rev_h_budget`, `cli_gfnff_03_rev_adduct_falsifier`, `cli_gfnff_04_rev_well_form` |
