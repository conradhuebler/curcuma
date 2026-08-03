# Convergence presets for geometry optimisation (`-convergence_preset`)

> 🤖 AI-generated, machine-tested (Aug 2026). Human production testing pending.

`-convergence_preset loose|normal|tight|verytight` sets all four convergence parameters at once.
The values live in **one** place, `src/capabilities/optimisation/convergence_presets.h`, which the
three call sites (`OptimizationContext::fromJson`, `buildPresetConfig`, `CurcumaOpt::LoadControlJson`)
read — they used to be three hand-copied blocks, identical by accident rather than by construction.

| preset | ΔE [kJ/mol] | ΔRMSD [Å] | \|g\| [Eh/Bohr] | max steps (N = atoms) |
|---|---|---|---|---|
| `loose` | 1 | 0.1 | 5e-3 | max(50, 1·N) |
| `normal` (default) | 0.1 | 0.01 | 5e-4 | max(500, 10·N) |
| `tight` | 0.01 | 1e-3 | 5e-5 | max(1000, 25·N) |
| `verytight` | 0.001 | 1e-4 | 5e-6 | max(5000, 100·N) |

One decade per step in every criterion, so the name states the accuracy to expect. All three
thresholds must hold in the same step (`-convergence_count 7`, a bit field: 1 = energy, 2 = RMSD,
4 = gradient). An explicit `-max_iterations` / `-energy_threshold` / `-rmsd_threshold` /
`-gradient_threshold` overrides the preset; a named step cap also switches the size scaling off and
is used verbatim.

## Why the step cap is the interesting number

Measured on seven MD snapshots of a 107-atom peptide (GFN-FF, identical starting geometries), the
**thresholds alone** cost almost nothing: 656 steps at `loose` versus 750 at `verytight` — three
decades of gradient threshold for 13 % more work. A cartesian L-BFGS spends its steps getting *near*
the minimum, not on the last digits.

The consequence for the old presets was that they did not distinguish anything. Eight GFN2
optimisations of that peptide, started from GFN-FF minima:

| old preset | wall time | result |
|---|---|---|
| `loose` | 131 s | within 0.06–0.17 kJ/mol of the others, identical ranking |
| `normal` | 137 s | " |
| `tight` | 204 s | " |

Three names, one optimisation. The step cap is therefore what separates the presets today.

## Why the cap scales with the system

The number of steps needed to reach a minimum grows with the degrees of freedom. Measured (GFN-FF,
steps to full convergence):

| system | atoms | steps |
|---|---|---|
| butane | 14 | 134 |
| caffeine | 24 | 232 |
| helicene | 42 | 149 |
| peptide (from MD snapshots) | 107 | 748 |

Roughly ten steps per atom — which is where the `normal` factor comes from. With a **fixed** cap of
500, `normal` truncated the peptide: 6 of 7 structures stopped at the cap, 16.9 kJ/mol above the
converged energy, and the ranking degraded (ρ = 0.86 instead of 1.00). With the size term the same
run reaches **0.096 kJ/mol with nothing truncated**, while butane, caffeine and helicene are
bit-identical — their caps do not move.

**What the size term does not fix is stiffness.** A water octamer needs ~1100 steps at 24 atoms
because it is floppy and hydrogen-bonded; no function of the atom count predicts that, so `tight`
still truncates 2 of 5 starts there and lands 5.2 kJ/mol short. Any truncation is reported
(`not converged after N`), never silent.

## Accuracy actually delivered

GFN-FF, five perturbed starting geometries per system, deviation from the `verytight` result in
kJ/mol (mean over the set):

| system | `loose` | `normal` | `tight` |
|---|---|---|---|
| butane, 14 | 1.03 | 0.013 | 0.011 |
| caffeine, 24 | 2.14 | 0.93 | 0.001 |
| helicene, 42 | 2.31 | 0.003 | 0.000 |
| water octamer, 24 | 252 | 208 | 5.2 |
| peptide, 107 (MD snapshots) | 114 | 0.096 | 0.000 |

`loose` is a *pre-optimisation*, not a result: on compact molecules it delivers roughly the 1 kJ/mol
its name claims, on floppy or large systems it stops far short by construction. That is exactly what
the crude stage of the two-stage refinement wants — eight GFN2 optimisations in 15 s instead of
120 s, 29 kJ/mol from converged, but ranking at ρ = 0.88 against the accurate result, where the
GFN-FF surface it replaces ranks at ρ = −0.13 (see
[CONFSEARCH_DUAL_METHOD.md](CONFSEARCH_DUAL_METHOD.md)).

## Known weakness

The gradient criterion is the **norm of the full 3N vector**, so it is system-size dependent:
5e-4 means 7.7e-5 per component for butane but 2.8e-5 for a 107-atom molecule — a large molecule is
converged more strictly than a small one at the same preset. ORCA and Gaussian use size-independent
max- and RMS-component criteria instead. Kept as is for now because every stored convergence
behaviour and every golden test value rests on it.
