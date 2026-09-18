# rev-gfnff stage 3a — the valence share, the budget rule and the bond-well form

🤖 AI-generated, machine-tested only. **Human production testing pending** (the operator removes
this line). Every number below was measured with one frozen binary per package; the run logs and
the per-cell tables are in `test_cases/revgfnff/_log/WORK_STATUS.md`.

Stage 3a modifies the GFN-FF **bond term** so that it survives a bond being made or broken. It has
three parts, and only the first is on by default:

| part | what it changes | flag | default |
|---|---|---|---|
| 3a(i) | the pair's own CN contribution is taken out of its own r0 | — | on |
| 3a(ii) | a valence SHARE factor multiplies the well; hydrogen keeps one valence | `-gfnff.rev_valence_share`, `-gfnff.rev_budget_fix_h`, `-gfnff.rev_share_form` | on, on, `delivered` |
| 3a(iii) | the SHAPE of the well | `-gfnff.rev_well_form` | `gauss` |

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
is zero.

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

### 1.3 `conserving` — the valence-conserving share (OPT-IN, `-gfnff.rev_share_form conserving`)

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

**What it costs, and why it is not the default**: a dative or ylidic neutral has a donor with four
partners at a group charge well below 1, so `Val ~ 3.2` against `S ~ 4` and all four of its wells
are scaled by ~0.8 — where the delivered share is inert (0.002-0.017 kcal/mol):

| system | conserving minus share-off |
|---|---:|
| H3N-BH3 | **+94.4** kcal/mol |
| H3N-O (amine oxide) | **+73.4** |
| H3N-CH2 (N-ylide) | **+109.5** |
| H5O2+ (Zundel) | -17.2 (closer to the pinned-topology `gfnff` value) |
| N2H7+ | -16.3 (likewise) |

The missing physics is that a dative bond puts a full valence into the acceptor's empty orbital,
which the donor's EEQ charge (~+0.2) does not express. Closing that is the open item.

---

## 2. The bond-well form (`-gfnff.rev_well_form`, DEFAULT `gauss`)

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
it cost **1.42x** the whole react-MD wall time against MG's 1.02x.

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
measured); periodic systems; charged species beyond the six hypervalent ions and the two
proton-shared dimers; any well form on a hydrogen-bonded X-H (the HB alpha modulation is not
applied in the new forms, and the class-A set contains no hydrogen bond); long-time MD stability
beyond 5 ps per cell; every combination of the three flags except the ones tabulated above.

**NOT implemented** relative to the plan: the second step of stage 3a(iii) (freeing the curvature
with an r0 re-solve); a bond-order-resolved well table (the current one is keyed on the element
pair, so C-C, C=C and C#C share one median — their per-system fits differ by `s` 1.18 / 1.01 /
0.91); a donor rule for the conserving budget, which is what the dative/ylide regression needs.

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
| diagnostics | `CURCUMA_SHAREDUMP=1` (per-pair `share`/`shareD` rows and, in `conserving`, a per-atom `shareA` row), `CURCUMA_REVDUMP=1` (the resolved settings) |
| measurements | `test_cases/revgfnff/_log/WORK_STATUS.md`, `HBUDGET_STATUS.md`, `RUNAWAY_STATUS.md`, `FABLE_REVIEW_2.md` |
| regression tests | `cli_simplemd_20_gfnff_rev_h_budget`, `cli_gfnff_03_rev_adduct_falsifier`, `cli_gfnff_04_rev_well_form` |
