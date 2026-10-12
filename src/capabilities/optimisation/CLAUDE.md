# CLAUDE.md - Optimisation Algorithms

Native algorithms behind `-optimizer lbfgs|diis|rfo` and the RF solver shared with ANCOpt; loop, convergence and output
belong to `OptimizerDriver` ([../CLAUDE.md](../CLAUDE.md)). Educational-first: keep the mathematics readable here.

## Files

- `lbfgs.*`: class `LBFGS` (all three methods), driven by `NativeLBFGS/DIIS/RFOAdapter` in `../native_optimizer_adapters.*`
- `rf_solver.*`: `RFSolver::calculateRFStep()` (Lanczos if nvar+1 >= 50, else dense), `lanczosLowestEigenpair()`; used by `RFOStep()` and ANCOpt
- `modern_optimizer_simple.*` was removed on 2026-10-01 (not compiled, included nowhere); its `modern_optimizer` PARAM module (11 parameters) left the registry with it. `main.cpp` still lists `modern_optimizer` among the scope modules, which has no effect

## Algorithms (`lbfgs.cpp`)

- **L-BFGS** `LBFGSStep()`: two-loop recursion (Nocedal & Wright Alg. 7.4), H0 = I, pair stored only if y.s > 1e-10;
  line search backtracking (Armijo c1 = 1e-4, halving, max 30) or `lineSearchStrongWolfe()` (N&W Alg. 3.5/3.6, unreachable via `-opt`)
- **DIIS** `diisExtrapolation()`: Pulay, CPL 73, 393 (1980); gradients as error vectors, B_ij = e_i.e_j + Lagrange row, QR solve,
  SVD only for the condition check (> 1e12 drops the oldest vector); L-BFGS steps until `diis_start`
- **RFO** `RFOStep()`: Banerjee et al., JPC 89, 52 (1985); mass-weighted lowest eigenpair of [[H, g], [g^T, 0]]; trust-radius clamp,
  step 1 / 0.5 / 0.1 or reject (radius x0.5), radius x1.2 after a decrease
- **SR1** `sr1Update()` after every RFO step: H + (y-Hs)(y-Hs)^T / ((y-Hs).s), skipped below 1e-12 (Broyden 1970, Dennis & Moré 1977)
- No caller: `lineSearchRFO()`, `line_search_backtracking()`, `updateHessian()`
- `EIGEN_USE_LAPACKE` is defined per file (`lbfgs.cpp`, `rf_solver.cpp`) only with `HAVE_LAPACKE`; no variable may be named `I` there

## Status

> **All native optimizers (L-BFGS, DIIS, RFO, SR1) are 🤖 AI-generated, not human-tested or approved.
> Human production testing pending - do not use for production calculations.** No ctest passes `-optimizer`.

| Method | Status | Tested on | Not tested |
|--------|--------|-----------|------------|
| Native L-BFGS | 🤖 AI-generated, ⚙️ compiles | water, ethane (UFF) | large systems, QM methods, constrained opt |
| Native DIIS | 🤖 AI-generated, ⚙️ compiles | small UFF molecules | convergence stability, near-degenerate cases |
| Native RFO | 🤖 AI-generated, ⚙️ compiles | small UFF molecules | saddle point searches, transition states |
| SR1 update | 🤖 AI-generated, ⚙️ compiles | indirectly via DIIS | standalone correctness vs. reference |

## Open items

- L-BFGS stalls once backtracking is exhausted (takes the ~1e-9 step, keeps its history); driver stall detection ends the run ([TODO.md](../../../TODO.md))
- The adapter never calls `setConfig`/`setLineSearch`/`setRFOSolver`: `native_lbfgs.lbfgs_line_search` and `rfo_solver` have no effect, verbosity is forced to 0
- Atom constraints: `initialize()` resets `m_constraints` after the adapter set them (`native_optimizer_adapters.cpp:53`, `:73`) and
  `getEnergyGradient()` indexes the 3N vector by atom (`lbfgs.cpp:674-676`); code reading, not run
- `lbfgs_line_search` is an Int (1-4) in module `opt` (LBFGSpp) but a String in `native_lbfgs`
- No numerical gradient check, step control not compared to xtb/Gaussian, no DIIS convergence guarantee, RFO mode following untested

## Recorded design decisions (Claude Generated, kept verbatim)

1. **`-optimizer` Parameter**: Bestimmt Optimierungsalgorithmus (nicht `-method`)
2. **Native First**: Curcuma's eigene Implementierungen bevorzugt
3. **Educational Output**: Literaturzitate und wissenschaftliche Parameter sichtbar
4. **Automatic Fallback**: Legacy system für unbekannte Optimizer
5. **Erweiterbar**: System kann einfach um weitere Optimizer erweitert werden

Code today: `auto` selects LBFGSpp and the `-opt` help recommends lbfgspp/ancopt (2); unknown names fall back to LBFGSpp (4).

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future tasks and visions defined by operator/programmer*

---

Previous version (removed 2026-10-01): [docs/archive/OPTIMISATION_NOTES_2026-10.md](../../../docs/archive/OPTIMISATION_NOTES_2026-10.md)
