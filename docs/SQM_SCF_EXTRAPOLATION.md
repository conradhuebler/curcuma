# Multi-step SCC extrapolation (native GFN1/GFN2)

> Status: 🤖 AI-generated, ⚙️ machine-tested only. Opt-in; the default path is
> unchanged. Human production testing pending.

## What it is

The native GFN SCF reuses the converged charge state of the previous geometry
step as the initial guess for the next one (the **1-step warm-start**,
`-warm_start true`, on by default). For a smooth trajectory the density is
*moving* in a predictable direction, so a **multi-step extrapolation** of the
self-consistent-charge (SCC) vector from several past steps lands closer to the
new fixpoint than the previous step alone — fewer SCF iterations, most
noticeably in MD (uniform timestep).

The extrapolated quantity is the full SCC mixing vector:
- **GFN1**: the shell charges `q_sh`.
- **GFN2**: `[q_sh; dp_at; qp_at]` (shell charges + atomic dipoles + quadrupoles).

The prediction is `P_pred = Σ_j w_j · P(n-1-j)` over the history of converged
steps, with **charge-conserving** weights (`Σ w_j = 1`).

This is **opt-in** and disabled by default (`scf_extrapolation = none` reproduces
the existing 1-step warm-start byte-for-byte).

## Flags

| Flag | Default | Values | Meaning |
|------|---------|--------|---------|
| `-scf_extrapolation` | `none` | `none` \| `aspc` \| `gauss` | Predictor scheme |
| `-scf_extrapolation_order` | `3` | integer ≥ 0 | ASPC order *k* (history *k+2*); Gauss polynomial degree |
| `-scf_extrapolation_apply` | `guess` | `guess` \| `xlbomd` | How the prediction couples to the SCF |
| `-scf_xlbomd_correctors` | `1` | integer ≥ 1 | Corrector SCF maps per step (`xlbomd` only) |

All four route to the `xtb` config scope (flat CLI flags auto-route). They apply
to `-method gfn1` and `-method gfn2`, in `-sp` (no effect: no history), `-opt`,
and `-md`.

### Predictor schemes

- **`aspc`** — Always Stable Predictor-Corrector (Kolafa, *J. Comput. Chem.*
  2004). Fixed coefficients `w_j = (-1)^j C(k+2, j+1)` over the last `k+2`
  converged steps (k=0 → `[2, -1]`, the 2-point linear predictor; k=3 →
  `[5, -10, 10, -5, 1]`). Time-reversible for a fixed MD timestep — the
  recommended scheme for MD.
- **`gauss`** — least-squares polynomial of degree `order` fit to the SCC-vector
  history and extrapolated to the next step. More robust than high-order ASPC for
  the *irregular* step lengths of a geometry optimisation, where ASPC can
  overshoot.

### SCF coupling

- **`guess`** (default, **safe**) — the prediction only seeds the SCF initial
  guess; the SCF still iterates to `scf_threshold`, so the converged energy and
  density are unchanged within the threshold (same guarantee as the 1-step
  warm-start). A poor prediction can only cost a few extra iterations, never
  accuracy.
- **`xlbomd`** (**experimental**) — extended-Lagrangian Born-Oppenheimer MD: the
  SCC density is a **time-reversibly propagated auxiliary variable** (Verlet +
  Niklasson dissipation, `J. Chem. Phys.` **130**, 214109 (2009)) that seeds the
  SCF each step. The SCF then **converges normally** (the auxiliary is a good
  guess, so this is few iterations), so the energy is **exact**, and the
  time-reversibility of the auxiliary lowers long-MD energy drift. The dissipation
  order is `scf_extrapolation_order` clamped to 3..7 (default 3). The first K+1
  steps bootstrap with a normal SCF seeding the auxiliary trajectory.

  Note: a naive "few bare SCF maps, no convergence" corrector is **not
  contractive** for tight-binding SCF and diverges, so the implemented corrector
  converges — the iteration savings come from the good guess (like `guess` mode),
  while XL-BOMD's distinct value is the time-reversible auxiliary for MD energy
  conservation. `scf_xlbomd_correctors` is a minimum-corrector-iteration floor
  (default 1 = no effect; raise for tighter MD forces). MD only; use with
  `scf_extrapolation=aspc` and a fixed timestep. **Energy drift is unvalidated.**

## Usage

```bash
# MD with ASPC charge extrapolation (largest win at a fixed timestep)
curcuma -md mol.xyz -method gfn2 -scf_extrapolation aspc -scf_extrapolation_order 3

# Optimisation — prefer gauss (or a lower aspc order) for irregular LBFGS steps
curcuma -opt mol.xyz -method gfn1 -scf_extrapolation gauss -scf_extrapolation_order 2

# XL-BOMD MD (experimental — time-reversible auxiliary density, converged corrector)
curcuma -md mol.xyz -method gfn2 -scf_extrapolation aspc -scf_extrapolation_apply xlbomd
```

At `-verbosity 2` (`-verbosity 3` in opt/MD iterative mode) the per-step log
shows `SCF initial guess: aspc extrapolation over N steps`.

## Correctness and caveats

- **`guess` mode converges to the same fixpoint** as the 1-step warm-start within
  `scf_threshold`. Validated: caffeine GFN2 `-opt` reaches the same minimum
  (−42.1538368 Eh) with `none`/`aspc`/`gauss`. The
  `scf_extrapolation_caffeine` ctest (`ctest -L scf_extrapolation`) drives a
  smooth GFN1+GFN2 trajectory and asserts equal final energy **and** that
  aspc/gauss never cost more iterations than `none`.
- **Same loose-SCF gradient caveat as the warm-start**: at a loose
  `scf_threshold` (default 1e-5) the converged point sits anywhere in the
  convergence ball, so the gradient can shift at the ~1e-5 level. Tighten
  `-scf_threshold` for high-precision forces — independent of extrapolation.
- **MD vs opt**: the win is largest in MD (uniform timestep → ASPC is exact to
  its order). In `-opt` the LBFGS steps are irregular and high-order ASPC can
  overshoot; prefer `gauss` or a lower `scf_extrapolation_order`, or leave it off.
- **`xlbomd` is experimental** — energy conservation / drift is not yet
  validated. Treat results as exploratory.

## What was tested / not tested

- **Tested**: caffeine GFN1+GFN2 along a smooth interpolated trajectory.
  `guess` (`aspc`/`gauss`): same converged energy as `none`, iteration count
  ≤ `none`. `xlbomd`: energy bit-identical to `none` (converged corrector) with
  fewer iterations; a short NVE MD stays bounded. Default-path regression (`none`)
  byte-unchanged across the `gfn1_validation`/`gfn2_validation`/`xtb_gradient`/
  `native_xtb`/`xtb_cpscf` suites. All in `ctest -L scf_extrapolation`.
- **Not tested**: long production MD energy drift / conservation (the distinct
  XL-BOMD benefit — only smoke-tested for boundedness); systems with large
  per-step geometry changes / bond breaking; charged / open-shell systems; the GPU
  resident path (it consumes the same `m_wfn` guess, so it benefits transparently,
  but this was not separately benchmarked).

## Citations

When a scheme is actually applied, curcuma registers its reference in the
citation system (printed in the end-of-run summary and written to
`*_citations.bib`):

- **`aspc`** → Kolafa, J. *J. Comput. Chem.* **25**, 335–342 (2004),
  doi:10.1002/jcc.10385 (citation key `aspc` / `kolafa2004aspc`).
- **`gauss`** → Pulay, P.; Fogarasi, G. *Chem. Phys. Lett.* **386**, 272–278
  (2004), doi:10.1016/j.cplett.2004.01.069 (key `density_extrapolation` /
  `pulay2004fockdynamics`).
- **`xlbomd`** → Niklasson, A. M. N. *Phys. Rev. Lett.* **100**, 123004 (2008),
  doi:10.1103/PhysRevLett.100.123004 (key `xlbomd` / `niklasson2008xlbomd`).
