# Native GFN large-system modes (C1)

> Status: 🤖 AI-generated, ⚙️ machine-tested (June 2026). **Approximate, opt-in**
> research paths — NOT human production tested. The dense default is unchanged.
> Validate against the dense reference before research use.

The native GFN1/GFN2 SCF (`curcuma::xtb::XTB`) is dense: the wall above ~1000
atoms is the **O(N³) eigensolve**. The three **C1** modes trade a controlled,
*measurable* accuracy loss for the ability to run large systems by exploiting
locality. They are opt-in via `-c1_mode`; the default `none` is the exact dense
path and is byte-for-byte unchanged.

```
-c1_mode none        # default: exact dense O(N^3) eigensolve
-c1_mode fragments   # C1c: disconnected-fragment SCF
-c1_mode dc          # C1a: divide-and-conquer overlapping fragments
-c1_mode sparse      # C1b: sparse + non-orthogonal density purification (0 K, gapped)
```

All three **reuse the validated dense `XTB` engine** (SCF, gradient, dispersion,
energy evaluation); only the partition/assembly or the eigensolve replacement is
new. Energies are evaluated through the existing, tblite-validated
`evaluateComponentsAtFixedDensity`, so a C1 mode's correctness reduces to
producing a density that converges to the dense one — which the accuracy tables
below confirm.

---

## C1c — disconnected-fragment SCF (`-c1_mode fragments`)

Partition by **bond connectivity** (`Molecule::GetFragments`), run one dense SCF
per fragment, **sum** energies and **scatter** each fragment's charges/gradient
onto its atoms (block-diagonal). Exact for non-interacting fragments; **neglects
all inter-fragment coupling**. A connected molecule is a single fragment →
transparent exact dense fallback. **Energy + gradient** (the only C1 mode with a
gradient — geometry optimisation works).

Use for: non-bonded molecular clusters / many-molecule assemblies. Charged
multi-fragment systems are not supported (each fragment treated as neutral; warns).

**Validation** — 8 well-separated waters (2×2×2 grid), error vs dense = the
neglected inter-water interaction, which vanishes with separation:

| spacing | dense / Eh | C1c / Eh | error |
|--------:|-----------:|---------:|------:|
| 5 Å | -40.56363345 | -40.56296612 | 0.667 mEh |
| 8 Å | -40.56328980 | -40.56296612 | 0.324 mEh |
| 12 Å | -40.56309210 | -40.56296612 | 0.126 mEh |
| 20 Å | -40.56299735 | -40.56296612 | 0.031 mEh |

Connected molecules (caffeine, triose) reproduce the dense energy **exactly**
(single-fragment fallback). Optimisation converges (grad norm → 1e-5).

---

## C1a — divide-and-conquer (`-c1_mode dc`)

The general method for **connected** large systems (Yang & Lee, *J. Chem. Phys.*
**103** (1995) 5674). Spatial cubic-cell core partition (`-c1_cell_bohr`) + a
buffer shell (`-c1_buffer_bohr`); each core+buffer is a sub-system. Each outer
SCF iteration builds the **global** Fock from the current density (so severed
covalent bonds still see their real neighbours), diagonalises only the
**sub-blocks**, finds a shared chemical potential by global electron-count
bisection over all sub-spectra, and assembles a global density by Yang's
core-projection. Broyden charge mixing (the same mixer the dense SCF uses for
stiff systems). **Energy-only** first cut (no analytic gradient).

**Accuracy knob = `-c1_buffer_bohr`**: larger buffer → converges monotonically
to the dense energy.

**Validation** — `complex` (231 atoms), cell 10 Bohr, vs dense -329.52707823 Eh:

| buffer | max sub | C1a / Eh | error | mEh/atom |
|-------:|--------:|---------:|------:|---------:|
| 6 Bohr | 62 | -329.23565150 | 291.4 mEh | 1.26 |
| 8 Bohr | 101 | -329.49909519 | 27.98 mEh | 0.121 |
| 10 Bohr | 145 | -329.51312658 | 13.95 mEh | 0.060 |
| 12 Bohr | 178 | -329.52373560 | 3.343 mEh | 0.014 |
| 14 Bohr | 200 | -329.52650459 | 0.574 mEh | 0.0025 |
| 16 Bohr | 225 | -329.52697970 | 0.099 mEh | 0.0004 |

**Validation** — `polymer` (1410 atoms), cell 10 Bohr, vs dense -2088.25340678 Eh
(dense single point: 70.3 s):

| buffer | max sub | C1a / Eh | error | mEh/atom | time |
|-------:|--------:|---------:|------:|---------:|-----:|
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

## C1b — sparse + non-orthogonal purification (`-c1_mode sparse`)

The O(N) path for **gapped** systems at **0 K** (`-electronic_temperature 0`).
Replaces the eigensolve with non-orthogonal **canonical density-matrix
purification** (Palser & Manolopoulos, *Phys. Rev. B* **58** (1998) 12704,
generalised to the S metric: products P² → PSP, traces Tr(P) → Tr(PS)). It
purifies the projector `M = P·S` (plain matrix powers; the metric enters only in
the M₀ initial guess and the AO density recovery `D = 2·M·S⁻¹`). The converged
projector is pruned at `-c1_sparse_threshold` and the pruned density feeds the
SCF, so the charges feel the sparsity. **Gapped only** — purification needs a
HOMO–LUMO gap; on non-convergence it warns and **falls back to the eigensolver**
(energy stays correct, sparsity unavailable). **Energy-only.**

> First cut: **dense Eigen storage + thresholding** — it measures the energy
> error and the achievable **nnz fraction** vs the threshold (the validation
> contract). True sparse storage / S⁻¹-free O(N) work is **deferred**, so C1b on
> very large N is currently slow (dense O(N³) per purification step).

**Accuracy knob = `-c1_sparse_threshold`** (energy → dense as threshold → 0; nnz
decreases as the threshold grows).

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

**Scope of the C1b validation.** Purification + sparsity is demonstrated on
**gapped systems up to `complex` (231 atoms)**. On the 1410-atom `polymer` the
purification does **not** converge and C1b **falls back to the eigensolver** (the
energy is still exact, -2088.25340678 Eh): the frontier-orbital gap is too small
for the first-cut canonical-PM watchdog at that scale, and the dense-storage
purification is in any case O(N³) per step (~3 min) until true sparse storage
lands. So C1b today is a **validated-but-small-scale** prototype; C1a is the
mode to use for the largest connected systems.

---

## Out of scope (first-cut prototypes)

- Analytic DC / sparse gradients (C1a/C1b are energy-only).
- True sparse storage / cell lists / FMM-Ewald long-range Coulomb (the Coulomb
  `gamma·q` stays an O(N²) matvec).
- Metallic / gapless systems (C1a smearing helps; C1b requires a gap and falls
  back to the eigensolver).
- Periodic boundary conditions; inter-fragment electronic polarisation in C1c
  (that is what C1a's buffer adds).

## Reproduce

```bash
cd release
# dense reference
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfn2
# C1a sweep
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfn2 \
          -c1_mode dc -c1_cell_bohr 10 -c1_buffer_bohr 12 -verbosity 1
# C1b sweep (0 K, gapped)
./curcuma -sp ../test_cases/molecules/larger/complex.xyz -method gfn2 \
          -c1_mode sparse -c1_sparse_threshold 1e-5 -electronic_temperature 0 -verbosity 1
# fast ctest
ctest -L c1 --output-on-failure
```
