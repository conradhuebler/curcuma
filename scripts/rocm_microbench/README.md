# ROCm micro benchmarks (Oct 13, 2026)

Standalone HIP programs behind [docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md](../../docs/ROCM_PERFORMANCE_ANALYSIS_2026-10-13.md). Not part of the build.

```
hipcc -O2 --offload-arch=gfx1200 eig_solver_steps.cpp   -o eig1 -L/opt/rocm/lib -lrocsolver -lrocblas
hipcc -O2 --offload-arch=gfx1200 eig_reduction_paths.cpp -o eig2 -L/opt/rocm/lib -lrocsolver -lrocblas
hipcc -O3 --offload-arch=gfx1200 coulomb_farfield.cpp   -o coul   && ./coul test_cases/structures/bulk/peo201x2-water1500.gfnff-min.xyz
hipcc -O2 --offload-arch=gfx1200 erf_saturation.cpp     -o erf    && ./erf
```
Matrices are synthetic (random symmetric F, near-identity S), so timings show solver cost, not SCF behaviour. `coulomb_farfield.cpp` uses
random charges and alpha values, no periodic boundary; it is a copy of `k_coulomb_dense` plus the far-field variant.
