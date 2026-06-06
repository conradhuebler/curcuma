# Native GFN large-system modes

> Status: 🤖 AI-generated, ⚙️ machine-tested (June 2026). **Approximate, opt-in**
> research paths — NOT human production tested. The dense default is unchanged.
> Validate against the dense reference before research use.

The native GFN1/GFN2 SCF (`curcuma::xtb::XTB`) is dense: the wall above ~1000
atoms is the **O(N³) eigensolve**. The three `-large_system_mode` modes trade a
controlled, *measurable* accuracy loss for the ability to run large systems by
exploiting locality. They are opt-in via `-large_system_mode`; the default
`none` is the exact dense path and is byte-for-byte unchanged.

```
-large_system_mode none        # default: exact dense O(N^3) eigensolve
-large_system_mode fragments   # disconnected-fragment SCF
-large_system_mode dc          # divide-and-conquer overlapping fragments
-large_system_mode sparse      # sparse + non-orthogonal density purification (0 K, gapped)
```

All three **reuse the validated dense `XTB` engine** (SCF, gradient, dispersion,
energy evaluation); only the partition/assembly or the eigensolve replacement is
new. Energies are evaluated through the existing, tblite-validated
`evaluateComponentsAtFixedDensity`, so a large-system mode's correctness reduces
to producing a density that converges to the dense one — which the accuracy
tables below confirm.

---

## Combination with `-eigensolver`

The three modes propagate `-eigensolver` differently. `-eigensolver` itself
(`mkl` / `native` / `purify` / `lobpcg`, see [SQM_EIGENSOLVE_GPU.md](SQM_EIGENSOLVE_GPU.md))
governs **how the eigenvalue problem is solved**; `-large_system_mode` governs
**how the system is partitioned**. They compose like this:

| `-large_system_mode` | `-eigensolver` propagation | notes |
|----------------------|----------------------------|-------|
| `none` (default) | applies to the single dense matrix | unchanged from the default path |
| `fragments` | applies **per fragment** (each fragment is a small dense SCF) | `-eigensolver=purify` on every fragment forces T=0 |
| `dc` | `purify` and `lobpcg` **fall back to dense GES** per sub-block; `native` applies directly | purify produces projector eigenvalues (0/1), not orbital energies; lobpcg gives only a partial spectrum; DC needs the full spectrum for Fermi occupation / chemical-potential bisection |
| `sparse` | **ignored** — sparse IS the eigensolver replacement | `-eigensolver` is consulted for nothing; the driver always runs purification |

**Hard rule:** `-large_system_mode=fragments|dc + -eigensolver=purify` requires
`-electronic_temperature 0`. Purification is 0-K by construction (integer
occupation, no Fermi smearing). The wrapper enforces this as a **hard error in
`setMolecule`** — not a warning, not a silent fallback. Example:

```bash
# OK: T=0, purify runs per fragment (fragments mode — each fragment is a full SCF)
./curcuma -sp cluster.xyz -method gfn2 \
          -large_system_mode fragments -eigensolver purify -electronic_temperature 0

# SOFT FALLBACK: purify in DC falls back to dense GES per sub-block
# (DC needs full spectrum for Fermi occupation; purify can't provide it)
./curcuma -sp complex.xyz -method gfn2 \
          -large_system_mode dc -eigensolver purify -electronic_temperature 0

# HARD ERROR: T=300 incompatible with purify
./curcuma -sp complex.xyz -method gfn2 \
          -large_system_mode dc -eigensolver purify -electronic_temperature 300
# NativeXtbMethod: -eigensolver=purify is 0-K only; set -electronic_temperature 0
```

---

## `fragments` — disconnected-fragment SCF

Partition by **bond connectivity** (`Molecule::GetFragments`), run one dense SCF
per fragment, **sum** energies and **scatter** each fragment's charges/gradient
onto its atoms (block-diagonal). Exact for non-interacting fragments; **neglects
all inter-fragment coupling**. A connected molecule is a single fragment →
transparent exact dense fallback. **Energy + gradient** (the only mode with a
gradient — geometry optimisation works).

Use for: non-bonded molecular clusters / many-molecule assemblies. Charged
multi-fragment systems are not supported (each fragment treated as neutral; warns).

**Validation** — 8 well-separated waters (2×2×2 grid), error vs dense = the
neglected inter-water interaction, which vanishes with separation:

| spacing | dense / Eh | fragments / Eh | error |
|--------:|-----------:|---------------:|------:|
| 5 Å | -40.56363345 | -40.56296612 | 0.667 mEh |
| 8 Å | -40.56328980 | -40.56296612 | 0.324 mEh |
| 12 Å | -40.56309210 | -40.56296612 | 0.126 mEh |
| 20 Å | -40.56299735 | -40.56296612 | 0.031 mEh |

Connected molecules (caffeine, triose) reproduce the dense energy **exactly**
(single-fragment fallback). Optimisation converges (grad norm → 1e-5).

---

## `dc` — divide-and-conquer

The general method for **connected** large systems (Yang & Lee, *J. Chem. Phys.*
**103** (1995) 5674). Spatial cubic-cell core partition
(`-large_system_cell_bohr`) + a buffer shell (`-large_system_buffer_bohr`); each
core+buffer is a sub-system. Each outer SCF iteration builds the **global** Fock
from the current density (so severed covalent bonds still see their real
neighbours), diagonalises only the **sub-blocks** (using the user-selected
`-eigensolver` on each), finds a shared chemical potential by global
electron-count bisection over all sub-spectra, and assembles a global density by
Yang's core-projection. Broyden charge mixing (the same mixer the dense SCF
uses for stiff systems). **Energy-only** first cut (no analytic gradient).

**Accuracy knob = `-large_system_buffer_bohr`**: larger buffer → converges
monotonically to the dense energy.

**Validation** — `complex` (231 atoms), cell 10 Bohr, vs dense -329.52707823 Eh:

| buffer | max sub | dc / Eh | error | mEh/atom |
|-------:|--------:|--------:|------:|---------:|
| 6 Bohr | 62 | -329.23565150 | 291.4 mEh | 1.26 |
| 8 Bohr | 101 | -329.49909519 | 27.98 mEh | 0.121 |
| 10 Bohr | 145 | -329.51312658 | 13.95 mEh | 0.060 |
| 12 Bohr | 178 | -329.52373560 | 3.343 mEh | 0.014 |
| 14 Bohr | 200 | -329.52650459 | 0.574 mEh | 0.0025 |
| 16 Bohr | 225 | -329.52697970 | 0.099 mEh | 0.0004 |

**Validation** — `polymer` (1410 atoms), cell 10 Bohr, vs dense -2088.25340678 Eh
(dense single point: 70.3 s):

| buffer | max sub | dc / Eh | error | mEh/atom | time |
|-------:|--------:|--------:|------:|---------:|-----:|
| 8 Bohr | 145 | -2088.16258674 | 90.8 mEh | 0.064 | 35.8 s |
| 12 Bohr | 328 | -2088.25105909 | 2.35 mEh | 0.0017 | 60.1 s |
| 16 Bohr | 576 | -2088.25330745 | 0.099 mEh | 0.0001 | 175.8 s |

At buffer 8 Bohr DC is **2× faster than the dense SCF** at 0.064 mEh/atom (well
under a 1 mEh/atom target). The speedup is modest for this *condensed* polymer
(a small buffer already pulls in ~145 atoms); DC pays off more for elongated /
sparse systems and for much larger N where the dense O(N³) becomes prohibitive.
A huge cell+buffer (whole molecule per sub-system) reproduces dense to ~3e-7 Eh
(machinery check).

---

## `sparse` — sparse + non-orthogonal purification

The O(N) path for **gapped** systems at **0 K** (`-electronic_temperature 0`).
Replaces the eigensolve with non-orthogonal **canonical density-matrix
purification** (Palser & Manolopoulos, *Phys. Rev. B* **58** (1998) 12704,
generalised to the S metric: products P² → PSP, traces Tr(P) → Tr(PS)). It
purifies the projector `M = P·S` (plain matrix powers; the metric enters only in
the M₀ initial guess and the AO density recovery `D = 2·M·S⁻¹`). The converged
projector is pruned at `-large_system_sparse_threshold` and the pruned density
feeds the SCF, so the charges feel the sparsity. **Gapped only** — purification
needs a HOMO–LUMO gap; on non-convergence it warns and **falls back to the
eigensolver** (energy stays correct, sparsity unavailable). **Energy-only.**

> First cut: **dense Eigen storage + thresholding** — it measures the energy
> error and the achievable **nnz fraction** vs the threshold (the validation
> contract). True sparse storage / S⁻¹-free O(N) work is **deferred**, so on
> very large N this is currently slow (dense O(N³) per purification step).

**Accuracy knob = `-large_system_sparse_threshold`** (energy → dense as threshold
→ 0; nnz decreases as the threshold grows).

**Validation** — 0 K, vs dense-0K. Energy error / projector nnz fraction:

`complex` (231 atoms, dense-0K -329.52707823 Eh):

| threshold | error | nnz(M) |
|----------:|------:|-------:|
| 0 | 0.057 mEh | 1.000 |
| 1e-6 | 0.057 mEh | 0.842 |
| 1e-5 | 0.057 mEh | 0.604 |
| 1e-4 | 0.074 mEh | 0.330 |
| 1e-3 | 0.474 mEh | 0.136 |

`caffeine` / `triose` converge to the dense-0K energy to ~1e-4 mEh at threshold 0
and sparsify with the threshold (triose: nnz 1.00 → 0.38 over 0 → 1e-3).
Larger / more insulating systems sparsify more (density-matrix nearsightedness,
Kohn 1996) — at 1e-5 `complex`'s projector is already 60 % sparse for a 0.057 mEh
error. The 0.057 mEh floor at threshold 0 is the purification's idempotency
tolerance vs the eigensolver (≈0.25 µEh/atom).

**Scope of the validation.** Purification + sparsity is demonstrated on
**gapped systems up to `complex` (231 atoms)**. On the 1410-atom `polymer` the
purification does **not** converge and `sparse` **falls back to the eigensolver**
(the energy is still exact, -2088.25340678 Eh): the frontier-orbital gap is too
small for the first-cut canonical-PM watchdog at that scale, and the dense-storage
purification is in any case O(N³) per step (~3 min) until true sparse storage
lands. So `sparse` today is a **validated-but-small-scale** prototype; `dc` is
the mode to use for the largest connected systems.

---

## Out of scope (first-cut prototypes)

- Analytic DC / sparse gradients (`dc` / `sparse` are energy-only).
- True sparse storage / cell lists / FMM-Ewald long-range Coulomb (the Coulomb
  `gamma·q` stays an O(N²) matvec).
- Metallic / gapless systems (`dc` smearing helps; `sparse` requires a gap and
  falls back to the eigensolver).
- Periodic boundary conditions; inter-fragment electronic polarisation in
  `fragments` (that is what `dc`'s buffer adds).

## Reproduce

```bash
cd release
# dense reference
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfn2
# dc sweep
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfn2 \
          -large_system_mode dc -large_system_cell_bohr 10 -large_system_buffer_bohr 12 -verbosity 1
# sparse sweep (0 K, gapped)
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfn2 \
          -large_system_mode sparse -large_system_sparse_threshold 1e-5 -electronic_temperature 0 -verbosity 1
# fragments + native eigensolver: many-fragment cluster, MKL-free pipeline
./curcuma -sp cluster.xyz -method gfn2 \
          -large_system_mode fragments -eigensolver native
# fast ctest
ctest -L large_system --output-on-failure
```
