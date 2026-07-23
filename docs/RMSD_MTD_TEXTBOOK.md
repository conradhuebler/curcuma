# RMSD metadynamics in SimpleMD: textbook (well-tempered) MTD vs. the counter scheme

Status: design note / specification. 🤖 AI-generated analysis, no code changes yet. Not a
feature description — it derives the reference method, analyses the current implementation against
it, and specifies the intended change plus its validation. All code references are to the
`confsearch` branch (`src/capabilities/simplemd.cpp`, `simplemd.h`).

## 0. TL;DR

SimpleMD's `-rmsd_mtd` bias is a memory-efficient, region-quantised approximation of
metadynamics, not textbook metadynamics. It differs in three independent places:

1. **Deposition cadence** — the bias (force + height growth) is evaluated **every MD step**
   (`rmsd_mtd_pace` default 1), not once per deposition interval.
2. **Occupancy gate** — a hill's height counter is incremented for **any** hill within
   `~1.1-1.4 Å` RMSD of the walker (`rmsd_econv = 1e8`), far outside the width where that hill
   contributes to the force. Growth is decoupled from real occupancy.
3. **Height damping** — hill height `W_i = k · counter_i` grows **linearly and unbounded**; the
   well-tempered damping (`rmsd_mtd_dt`) is applied only to a reported energy, never to the force.

Together these are the analytic root of the "temperature runaway" that the `freeze_inherited`,
`max_height`, `temp_abort` and `max_gaussians` safeguards were added to contain. This note specifies
the replacement scheme and its two-phase roadmap (Section 0.1). The Jul-2026 Gaussian-cutoff **screen**
(`-rmsd_mtd_screen`, see `CONFSEARCH_MTD_SCREEN.md`) is orthogonal — it is a physics-preserving
speed optimisation and none of the changes here touch it.

## 0.1 Roadmap — two phases

The work splits into two phases with a shared engine. **Phase 1 is implemented first and stands on
its own**; Phase 2 builds on its machinery.

- **Phase 1 — robust + fast ConfSearch** (Sections 4-7). Replace the counter/`econv` scheme with the
  *strided scheme*: decoupled force/deposition cadences (held + smoothed force), an interpretable
  deposition spacing `r_dep`, a soft residence-weighted counter (deletes the runaway-causing gate),
  a displacement gap-guard, and provenance diagnostics. Goal: exploration that is faster (the
  expensive Kabsch fleet runs ~10-20x less often) and more robust (no runaway, no coverage gaps, no
  impulsive kicks). Well-tempering is built in only as a *hook* (the soft increment can be damped) —
  not required for exploration.

- **Phase 2 — RMSD-family path CV for dissociation thermodynamics** (Section 8). Turn the WT hook on
  and add a path collective variable (Branduardi s, z) on an RMSD-family metric (RMSD vs. MSD is an
  open Phase 2 decision, 8.1) to compute dissociation free energies of organic complexes (e.g. a
  sugar unbinding from a receptor). Target: relative `ΔΔG` between similar guests first (much
  cancels, robust), absolute `ΔG` later. Reuses Phase 1's exploration to *discover* the unbinding
  path, then runs WT-MTD along it. Validated against PLUMED WT-MTD/PATHMSD.

## 1. Textbook metadynamics (reference)

### 1.1 General collective variable

Pick a collective variable (CV) `s(x)` — a function of the nuclear coordinates (a dihedral, a
distance, a coordination number, an RMSD, ...); it may be a vector `s = (s_1, ..., s_d)`. The
history-dependent bias deposits one Gaussian of **fixed** height every `tau_G` steps at the
current CV value:

```
              t' = tau_G, 2 tau_G, ... ; t' < t
V(s, t)  =   SUM   W * exp( - SUM_d (s_d - s_d(x(t')))^2 / (2 sigma_d^2) )
```

- `W`     fixed hill height
- `sigma_d` hill width per CV
- `tau_G`  deposition interval ("pace")
- `s(x(t'))` CV value at the deposition time `t'`

The sum runs over **deposition times**, not regions. Every hill persists; the hill count grows
linearly with simulation time; overlapping hills in an oft-visited basin stack up on their own.

### 1.2 Force

The bias acts on the atoms through the CV via the chain rule:

```
F_k = - dV/dx_k = - (dV/ds) . (ds/dx_k)
```

The first factor is analytic from the Gaussian sum; the second is the CV gradient (for RMSD:
`dRMSD/dx`).

### 1.3 Convergence to the free energy

In the long-time limit the bias fills the wells until the biased CV distribution is flat, and

```
lim (t -> inf)  V(s, t)  =  - F(s) + const
```

with `F(s)` the free energy (PMF) along `s`. Recovering `F(s)` is the point of metadynamics; it
requires the fixed-height accumulation.

### 1.4 Well-tempered MTD

Fixed height never fully converges (V oscillates around `-F`). Well-tempered MTD damps each new
hill by the bias already present at its deposition point:

```
W(t')  =  W_0 * exp( - V(s(x(t')), t') / (k_B * Delta_T) )

  =>  lim (t -> inf) V(s, t)  =  - [ Delta_T / (T + Delta_T) ] * F(s) + const
```

Where bias already sits, new hills shrink; the bias grows bounded and converges smoothly. This is
why well-tempered MTD does not heat the dynamics.

### 1.5 RMSD as the CV

**(a) Single fixed reference, 1-D CV.** `s(x) = RMSD(x, x_ref)` (scalar, after optimal
superposition). Hills are laid along the RMSD axis at the visited `s(t')`.

**(b) Configuration-space flooding (multi-reference).** Every visited snapshot `x_i` becomes its
own centre; the CV is the distance to it:

```
V(x) = SUM_i W_i * exp( - alpha * RMSD(x, x_i)^2 )
```

The textbook form sets `W_i = W` (fixed) and adds a **new** centre `i` every `tau_G` steps.
Curcuma's scheme is variant (b) with a different height rule — see below.

## 2. The counter scheme as implemented (as-is)

Local (per-walker) path: `BiasThread::execute()`. Shared-pool (ConfSearch) path:
`SimpleMD::ApplyRMSDMTD()`. Both compute

```
V(x) = SUM_i W_i * SUM_p exp( - alpha * RMSD(x, x_{i,p})^2 ),   W_i = k * counter_i
```

(`p` over identity + symmetry images), and the force is its exact negative gradient.
Height `W_i = k * counter_i`: `simplemd.cpp:181` (local), `:3296` (shared).

Instead of depositing a fresh fixed-height Gaussian per pace, the scheme keeps **one** reference
per resolved region and raises its counter on every visit. Defaults (`simplemd.h:650-661`):
`k = 0.01 Eh`, `alpha = 10`, `rmsd_econv = 1e8` (legacy, not in PARAM block), `pace = 1`,
`Delta_T = 2000 K`, `wtmtd = false`.

### Deviation 1 — cadence (pace)

`ApplyRMSDMTD()` is called under `if (m_step % m_mtd_steps == 0)` (`:2747` Verlet, `:2993`
RATTLE). `m_mtd_steps` is loaded from `rmsd_mtd_pace` (`:322`), whose registry default is **1** —
so the bias is evaluated **every step**. (The header initialiser `m_mtd_steps = 10`,
`simplemd.h:477`, is dead: it is overwritten at construction. The `rmsd_mtd_pace` help text calls
the parameter "unused"; that is true only for *new-structure deposition*, which is gated by the
bias level — the value still controls how often force + counter update run.)

### Deviation 2 — occupancy gate

A hill's counter is bumped only if it is in the `visited` set (`:214-215` local, `:3381` shared),
gated by

```
expr * rmsd_econv > N          (expr = exp(-alpha * RMSD^2), N = pool size)
```

With `alpha = 10`, `rmsd_econv = 1e8`, this fails (hill *not* counted) only when
`RMSD^2 > (8 ln10 - ln N) / alpha`, i.e. an effective cutoff of

| pool size N | counter stops growing above RMSD |
|-------------|----------------------------------|
| 1           | 1.36 Å |
| 100         | 1.18 Å |
| 1000        | 1.07 Å |

But the Gaussian itself has meaningful weight only well below that:

| RMSD | `exp(-10 RMSD^2)` |
|------|-------------------|
| 0.3 Å | 0.41 |
| 0.5 Å | 0.082 |
| 0.7 Å | 0.0074 |
| 1.0 Å | 4.5e-5 |

So the counter grows out to `~1.2 Å` while the hill contributes force only below `~0.5 Å`. The
code comment reads "Visited if the walker sits **inside** this Gaussian" (`:204`), but the
threshold sits far outside it. In a compact basin many references lie within `~1 Å` of each other,
so a single step bumps the counters of the **whole local cluster** at once, regardless of which
hill the walker is actually in. (Numbers above are analytic from the defaults, not measured.)

### Deviation 3 — height damping

`W_i = k * counter_i` grows linearly and without bound. The well-tempered weight `factor`
(`:216-217`, `:3312`) uses the correct `exp(-V/(k_B Delta_T))` damping but only feeds a separate
*reported* energy `current_bias_wt` — never the force. Confirmed by the `rmsd_mtd_dt` help text:
"only used ... for the reported well-tempered energy -- it never affects the force or the
exploration."

### Consequence and existing safeguards

Combine the three: every step (D1) × whole nearby cluster (D2) × undamped (D3). Over an
`M`-step run a hill can reach, as an analytic upper bound under the loose gate,

```
W_i = k * counter_i  <=  k * M  =  0.01 * M  Eh     (M = 10000 -> up to ~100 Eh)
```

which drives the walker's kinetic energy up run-over-run. The current mitigations —
`rmsd_mtd_max_height` (cap `counter_i` in the force), `rmsd_mtd_freeze_inherited` (freeze
inherited heights), `temp_abort` (`simplemd.h:679-681`), and the `max_gaussians` pool cap — bound
the symptom after the fact rather than removing the cause.

## 3. The three levers

| Lever | Textbook (WT-)MTD | Current counter scheme | Knob to change |
|-------|-------------------|------------------------|----------------|
| Cadence | one hill per `tau_G` | every step (`pace = 1`) | `rmsd_mtd_pace` |
| Where height grows | at the deposited CV point | every hill within `~1.2 Å` | the `visited` gate `expr * econv > N` |
| Height rule | fixed `W`, or WT-damped `W * exp(-V/k dT)` | linear `k * counter_i`, undamped | `W_i` formula in force |

## 4. Design options

### Option A — opt-in well-tempered RMSD-MTD mode

A rigorous well-tempered variant selected by `-rmsd_mtd_wt`, layered on the strided scheme (it only
changes the height increment — see 5.3):

- Deposit one hill per `pace` steps at the current configuration (multi-reference flooding,
  Section 1.5b), width `sigma` (reuse `alpha` via `alpha = 1/(2 sigma^2)`).
- Height in the **force**: `W_i = W_0 * exp( - V(x_i, t_dep,i) / (k_B Delta_T) )` evaluated at
  deposition time and stored per hill (constant thereafter), or the exact per-step WT form.
- `V -> -[Delta_T/(T+Delta_T)] F` recoverable; no runaway by construction; `temp_abort` etc.
  become unnecessary in this mode.
- Cost: a second, parallel bias code path and its bookkeeping (deposition, per-hill stored height).

### Option B — in-place fixes to the counter scheme

Keep one code path, remove the runaway cause:

- **Gate**: replace the binary `expr * econv > N` with either a value gate `expr > c`
  (e.g. `c = 0.5` -> RMSD < 0.26 Å at `alpha = 10`) or a soft counter `counter += expr` — so only
  hills the walker is actually in grow, which is also closer to continuous deposition.
- **Damping**: use `W_i = k * counter_i * exp(-V/(k_B Delta_T))` in the force (not only in the
  reported energy), bounding the height.
- Minimal change, one path; departs from strict textbook MTD but fixes the physics that matters.

### Recommendation — mapped onto the two phases

**Phase 1** is **B** built out as the concrete *strided scheme* (Section 5), with the WT damping of
**A** wired in only as a dormant hook on the height increment. It removes the runaway at its root,
replaces the opaque `econv` with interpretable spacing/cadence parameters, and makes ConfSearch
faster and more robust — without needing WT. **Phase 2** then activates **A** (WT) and adds the path
CV (Branduardi s, z; RMSD/MSD metric an open Phase 2 decision, 8.1) for thermodynamics (Section 8).
One engine, no second code path.

## 5. Phase 1 — the strided scheme (robust + fast ConfSearch)

This scheme **replaces** the current counter/`econv` behaviour of `-rmsd_mtd` — it is the new
default, not an opt-in mode. Its purpose is exploration that is **faster** (the expensive Kabsch
fleet runs once per `deposit_stride` instead of every step — ~10-20x fewer fits, on top of the
existing screen) and **more robust** (soft counter + deleted gate remove the runaway; the gap-guard
removes coverage holes; held + smoothed force removes impulsive kicks). WT is present only as the
dormant hook of 5.3, used in Phase 2. The old path is kept transiently as a deprecated
`rmsd_mtd_scheme = legacy` switch so validation can A/B against it, removed after operator sign-off
(Section 6). Four decoupled parts. Defaults below: `k = 0.01 Eh`, `alpha = 10` (`sigma = 1/sqrt(2
alpha) = 0.22 Å`, `FWHM = 2.355 sigma = 0.53 Å`).

### 5.1 Two cadences (force vs. deposition)

The current code fuses force evaluation and deposition into one call gated by `pace`
(`simplemd.cpp:2746-2749`, `:2992-2995`), so raising `pace` turns the bias into impulsive kicks
(the force is *overwritten* fresh each step at `:4091`, never held). Split them:

- **Expensive evaluation** (the Kabsch fleet, deposition test, counter growth) runs every
  `deposit_stride` (default **10 fs**; converted to steps via `dt`), or earlier if the displacement
  trigger fires (5.4).
- **Bias force** acts **every** step over the *frozen* hill set, held and smoothly ramped between
  evaluations: store `F_old`, `F_target` (3N vectors); each step add
  `F_applied = smoothstep(F_old, F_target; lambda)` to `m_eigen_gradient`; at each evaluation set
  `F_old <- F_target`, recompute `F_target`, reset `lambda`. `smoothstep(a,b;l) = a + (b-a)(3l^2 -
  2l^3)` gives C1-continuous force. `transition_fraction f` (default **1.0**) sets the ramp length
  as a fraction of the stride; `f = 1` glides continuously (a ~one-stride causal lag, negligible at
  small `deposit_stride`); `f < 1` ramps fast then holds (fresher magnitude, small stale plateau).
- **Honesty note**: during a ramp the applied force is not exactly `-grad V(x)`. This is acceptable
  — MTD is already non-conservative and thermostatted — and, crucially, the **recorded** potential
  `V(x) = SUM_i k*counter_i*exp(...)` stays exact, so FES reconstruction (A) is unaffected; only the
  sampling force is smoothed.

### 5.2 Deposition threshold (replaces `econv`, deposition half)

Deposit a new hill at `x` when

```
V(x) < V_min = k * exp(-alpha * r_dep^2)
```

`r_dep` (Å) is the interpretable hill spacing in RMSD space; default ties to the width,
`r_dep = FWHM(alpha) = 2.355 / sqrt(2 alpha)` (**0.5 Å** at `alpha = 10`) so neighbouring hills
overlap. `V_min` is constant in pool size `N` (unlike the old `current_bias * econv < N`, whose
effective spacing was ~1 Å and drifted with `N`). Coarser conformer granularity should be obtained
by lowering `alpha` (widening hills) *together with* raising `r_dep`, not by widening `r_dep` alone.

### 5.3 Counter growth (replaces `econv`, gate half — the gate is deleted)

Per evaluation, grow every hill by its own Gaussian weight:

```
counter_i += expr_i          (expr_i = exp(-alpha * RMSD(x, x_i)^2), already computed for the force)
```

`counter` becomes `double`. No `visited` threshold survives (`:205`/`:3381` removed): growth is
`~1` on a hill, `~exp(-alpha r_dep^2)` at spacing, `~0` far away — smooth, self-limiting
(`SUM expr_i` per step ≈ number of overlapping hills), and Voronoi-flip-free. This is the fixed-grid
continuum form of MTD deposition, and it is the WT hook: for **A**, multiply the increment by
`exp(-V(x)/(k_B Delta_T))`. `rmsd_mtd_max_height` remains an optional backstop on the double.

### 5.4 Displacement trigger + the `deposit_stride`/`r_dep` coupling

`pace = 1` accidentally guaranteed the walker never moved far between deposition checks. With
`deposit_stride = 10-20 fs` that protection is gone: if the walker's RMSD moves more than `r_dep`
in one stride it can skip a region, leaving a coverage gap. Safe regime:

```
deposit_stride * v_rmsd  <~  r_dep
```

Rather than tune `deposit_stride` against an unknown `v_rmsd`, force an early evaluation by
displacement: fire the deposition test when the walker has moved `> r_dep` from the last deposited
hill since the last evaluation. Implement as **one** Kabsch/step to the last-deposited hill
(`O(1)`, hard guarantee); optionally use the screen's rotation-invariant sigma lower bound
(`CONFSEARCH_MTD_SCREEN.md`) as a free additional trigger (upper-bounds `V(x)`, so a provable
`V < V_min` can advance the deposit). Exposed as `rmsd_mtd_gap_guard` (Bool, default **true**).

### 5.5 Parameters (register in both `simplemd` and `confsearch`, snake_case)

| flag | type | default | meaning |
|------|------|---------|---------|
| `rmsd_mtd_scheme` | String | `strided` | new default; `legacy` kept transiently for A/B validation, then removed |
| `rmsd_mtd_deposit_stride` | Double (fs) | `10` | expensive-eval / deposition cadence |
| `rmsd_mtd_transition_fraction` | Double | `1.0` | force-ramp length as fraction of stride, in `(0,1]` |
| `rmsd_mtd_r_dep` | Double (Å) | `-1` = auto `FWHM(alpha)` | deposition spacing / `V_min` |
| `rmsd_mtd_gap_guard` | Bool | `true` | displacement trigger against coverage gaps |
| `rmsd_mtd_diag` | Bool | `true` | write provenance CSV + gnuplot scripts (BMT-routed); `false` = COLVAR only (see Section 7) |

`econv` and `rmsd_mtd_pace` are retired — accepted as deprecated no-op aliases that warn (so
existing command lines still parse), removed together with the `legacy` path. Existing `alpha`,
`k`, `max_height`, `wtmtd`, `dt` (Delta_T) are reused unchanged.

### 5.6 Code touch points (confsearch line numbers)

- Cadence split: `:2746-2749` / `:2992-2995` (call site) and the force overwrite at `:4091`; new
  `F_old`/`F_target`/`lambda` state on `SimpleMD`, blended into `m_eigen_gradient` every step.
- Local path `BiasThread::execute()`: height `:181`, gate `:205-206` (delete), counter `:214-217`
  (-> `+= expr`).
- Shared path `SimpleMD::ApplyRMSDMTD()`: height `:3296`, gate `:3381` (delete), `registerVisits`,
  deposition `:3428` / `:3542` (-> `V < V_min`).
- Screen: `eps` at `:141` re-derived from `V_min` (or `cutoff_tol` alone) instead of `N/econv`.

### 5.7 Migration (replacement, not opt-in)

The `strided` scheme is the default; there is no bit-identical fall-back in normal use. Migration
points:
- **Deprecated flags**: `econv`, `rmsd_mtd_pace` parse but warn and no-op; the warning names the
  replacements (`r_dep`/`V_min`, `deposit_stride`).
- **Restart format**: the checkpoint changes (`counter` int -> double; `deposit_stride`/`r_dep`
  stored; `econv`/`mtd_steps` dropped). Bump the checkpoint version and refuse an old counter-scheme
  restart with a clear message rather than misinterpret it.
- **ConfSearch safeguards**: `temp_abort`, `rmsd_mtd_freeze_inherited`, `rmsd_mtd_max_height` were
  defaults-ON to contain the old runaway. Once check 5 confirms the runaway is gone by construction,
  revisit those defaults (expected: OFF) and update the ConfSearch docs.
- **`legacy` switch**: retained only for the Section 6 A/B comparison, removed after sign-off.

## 6. Phase 1 validation plan

Each claim below must be pinned by a reproducible check before the mode is called correct
(machine-tested only until the operator runs it on real problems). Phase 2 (thermodynamics) has its
own validation in Section 8.

1. **Analytic gradient** — finite-difference `F_target` (the per-evaluation bias force) with the
   soft-counter height rule on a small molecule (butane), max component error `< 1e-6` relative.
2. **Characterise vs. the old scheme** — since we replace (not preserve) the old behaviour, run
   butane `-md -rmsd_mtd` under `strided` vs. `legacy` (same seed, screen off) and report how they
   differ (hill count, nearest-neighbour spacing, `<T>`), showing the new one covers RMSD space more
   evenly and does not run away. The `legacy` switch enables this comparison and is removed
   afterwards. (`legacy` itself stays bit-identical to the pre-change binary, per
   `CONFSEARCH_MTD_SCREEN.md` "Verification", so it is a trustworthy reference.)
3. **Held-force accuracy** — `strided` vs. a force-every-step reference driven by the *same*
   deposition sequence: sweep `deposit_stride` in {5, 10, 20} fs, report the trajectory RMSD
   divergence and confirm the same basins are escaped. Establishes how large `deposit_stride` may
   grow before the held+ramped force degrades sampling; quote the worst case.
4. **Coverage / no gaps** — histogram of nearest-neighbour RMSD among deposited hills should center
   near `r_dep` with no cluster above it; on a fast-escape case, verify `rmsd_mtd_gap_guard` keeps
   the max RMSD jump between consecutive deposits `<= r_dep + tol` (and that disabling it reproduces
   the gap, confirming the guard is what closes it).
5. **Runaway removed** — the `-startT 500` ConfSearch case that historically ran away: with
   `strided` (soft counter + `deposit_stride`) the running-mean `<T>` stays bounded **without**
   `temp_abort`/`freeze_inherited`/`max_height`. Report the max `<T>` observed vs. the legacy path.
6. **Reproducibility** — every quoted number carries its command line and seed; put the FD,
   legacy-identity, and coverage checks into `test_cases/cli/simplemd/` so they stay true.

(FES recovery is a Phase 2 check — Section 8.5.)

## 7. Diagnostics: structure provenance (CSV + gnuplot)

Make the origin and fate of every deposited bias structure inspectable graphically. All files route
through `outputPath()` (BMT-mandatory) and follow the existing analysis naming schema
(`basename.type.csv` + auto-generated `basename.type.gnu` -> PNG — the same pattern as the scattering
analysis, reusing `TrajectoryWriter`/the analysis output handlers where practical). Console output
stays ASCII-only; unicode is fine inside the files. Controlled by `rmsd_mtd_diag` (default on at
verbosity >= 1).

### 7.1 `basename.mtd_hills.csv` — one row per deposited hill (the provenance table)

```
# index  step  time_fs  energy_Eh  rmsd_ref  trigger  counter_final  temperature_K  cycle  persistent
```

- `trigger` — why this hill was born: `initial` | `bias_below_vmin` | `displacement` (the two
  deposition paths of 5.2 / 5.4). This is the core "where did it come from" column.
- `rmsd_ref` — position in RMSD space (RMSD to hill 0).
- `energy_Eh` — physical potential energy of the structure at deposition (tracks the search reaching
  lower minima over time).
- `counter_final` — accumulated soft height at run end (how strongly the region was revisited).
- `cycle` / `persistent` — shared pool / ConfSearch only: which temperature cycle deposited it and
  whether it is a fed-back optimised minimum. Single-walker MD writes `0` / `false`.

### 7.2 Height evolution — `basename.mtd_counter.csv` (long format)

```
# step  time_fs  hill_index  rmsd  expr  counter  bias_contribution_Eh
```

One blank-line-separated block per hill (so gnuplot `plot ... index` works). Supersedes the current
per-index `COLVAR_<index>` files (same data, one tidy file). Shows which hills dominate the bias and
when they grew.

### 7.3 Coverage — `basename.mtd_coverage.csv` (doubles as validation check 4)

```
# index  nn_rmsd
```

Nearest-neighbour RMSD among deposited hills; the distribution should sit near `r_dep`. A
`basename.mtd_coverage_statistics.csv` companion gives min / mean / max and the count above `r_dep`
(i.e. gaps).

### 7.4 gnuplot scripts (auto-generated, like the scattering `.gnu`)

- `basename.mtd_hills.gnu` -> **deposition map**: `rmsd_ref` (y) vs `step` (x), point colour =
  `counter_final` (or `energy_Eh`), point type by `trigger`. The single plot that shows where
  structures come from and how tall they grow.
- `basename.mtd_counter.gnu` -> `counter` vs `step`, one line per hill.
- `basename.mtd_coverage.gnu` -> nearest-neighbour RMSD histogram with an `r_dep` marker line.

Each emits a PNG when gnuplot is present (mirroring `*.scattering_plot.png`); otherwise the `.gnu` +
`.csv` are left for the user to run. All quoted diagnostics are reproducible from these files.

## 8. Phase 2 — RMSD-family path CV for dissociation thermodynamics

Builds on Phase 1. Goal: dissociation free energies of organic complexes (sugar-receptor and the
like), where COM distance (curved exit, pose vs. separation, multiple channels) and coordination
number (degenerate, plateaued) both fail, but the RMSD-family metric keeps the pose/contact pattern.

### 8.1 The metric and the CV (precise)

Not a single RMSD. Use the **MSD metric** `D(x, x_i)` = optimally-aligned mean-square displacement
(`= RMSD^2`, length^2) to an ordered set of `P` reference frames, and the Branduardi path variables:

```
s(x) = SUM_i i*exp(-lambda*D_i) / SUM_i exp(-lambda*D_i)     progress, in [1, P]
z(x) = -(1/lambda)*ln SUM_i exp(-lambda*D_i)                 distance from path, length^2
lambda ~ 2.3 / <D between adjacent frames>
```

RMSD-*family* (same optimal-alignment displacement RMSD measures) so it keeps pose sensitivity, but
bounded with both endpoints defined. We already compute the aligned displacement for Phase 1 RMSD;
Phase 2 adds the metric-to-ordered-frames, `lambda`, and the `s, z` values + analytic gradients
(chain rule through the exp/ln).

**Open decision (deferred to Phase 2): RMSD metric vs. MSD metric for `D_i`.** The `s, z` formulas
work with either. MSD (`RMSD^2`, PLUMED `PATHMSD`) has a non-singular gradient (`dRMSD/dx ~ 1/RMSD`
diverges at coincidence) and is the canonical choice; the RMSD metric (PLUMED generic `PATH`) is the
quantity we already compute and is more intuitive (Å), at the cost of regularising the gradient near
`RMSD=0`. Not decided here — needs more reading/experience. Phase 1 is unaffected (it uses RMSD as
today); this only fixes the Phase 2 CV metric.

### 8.2 CV definition for unbinding

- **References**: ordered frames bound -> released along the exit path.
- **Alignment**: align on the host (receptor), displace on the guest (separate ALIGN/DISPLACE atom
  sets) — the binding displacement, not whole-system rigid motion.
- Bias `s` with WT-MTD; add a wall on `z` to keep sampling near the path.

### 8.3 Path from exploration (reuse Phase 1)

Run Phase 1 multi-reference flooding to discover the exit; extract an ordered, roughly MSD-equidistant
subset of the deposited hills as the path frames; set `lambda` from the adjacent-frame MSD. Pipeline:
`explore (flooding) -> extract ordered frames -> PathMSD s -> WT-MTD -> F(s), Delta_G`.

### 8.4 WT core and free energy

Activate the Phase 1 hook: increment `*= exp(-V(s)/(k_B*Delta_T))`; `Delta_T` a parameter,
`gamma = (T+Delta_T)/T`. At convergence `V(s) = -[Delta_T/(T+Delta_T)] F(s) + C`, so

```
F(s) = -[(T+Delta_T)/Delta_T] * V(s) + C
Delta_G = F(released) - F(bound)
```

- **Relative `Delta-Delta-G` first** (similar guests): standard state and common path segments
  cancel -> robust; usually the chemically interesting quantity.
- **Absolute `Delta_G` later**: needs the standard-state / volume correction for the released state
  (define it at a finite `s`/`z`, add the analytic correction). Flagged as the harder, error-prone
  step.

### 8.5 Phase 2 validation (reference-match, no novelty in the WT/FES math)

1. **PLUMED cross-check** — reproduce the paper's own alanine-dipeptide path-CV FES against PLUMED
   PATHMSD + WT-MTD: same frames, same `lambda`, same `Delta_T`; agreement within a stated tolerance.
   This is the hard reference bar — the WT/FES machinery must match, not innovate.
2. **`s, z` gradient** — finite-difference check, `< 1e-6` relative.
3. **Convergence** — bias deposition rate -> small; FES-vs-time overlays; block-average error bar on
   `Delta_G`. A number without a converged FES and an error bar is not reported.
4. **Relative `Delta-Delta-G`** — a pair of sugar analogues: report value + error; the ordering must
   match the known chemistry.

### 8.6 Diagnostics (extends Section 7)

FES table `basename.fes.csv` (`s  V_Eh  F_Eh`) + gnuplot; `s(t)`/`z(t)` in COLVAR; per-frame
occupancy. All BMT-routed, same pattern as Section 7.

## References

- Laio & Parrinello, PNAS 99 (2002) 12562 — metadynamics.
- Barducci, Bussi & Parrinello, PRL 100 (2008) 020603 — well-tempered metadynamics.
- Branduardi, Gervasio & Parrinello, J. Chem. Phys. 126 (2007) 054103 — path collective variables
  (s, z) on the MSD metric (PLUMED PATHMSD).
- `docs/CONFSEARCH_MTD_SCREEN.md` — the orthogonal Gaussian-cutoff screen and pool cap.
