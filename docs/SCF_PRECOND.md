# Fragment SCF preconditioner for the native GFN2 and GFN1 Broyden mixer (`-scf_precond frag`)

Status: AI-generated, machine-tested only, human production testing pending. Opt-in; the default
(`-scf_precond none`) leaves the SCF path bit-identical to the state before this feature.
Code: `src/core/energy_calculators/qm_methods/xtb_precond.{h,cpp}` (Oct 2026).

## What it does

The native SCF iterates x -> G(x) on the SCC vector x = [q_sh | dp | qp] and mixes with modified
Broyden (`-scf_mode broyden`, the default). Broyden starts every inverse-Jacobian estimate from the
scalar step alpha*F (alpha = `scf_damping` = 0.4). `-scf_precond frag` replaces that seed by M*F,
where M approximates -J^-1 molecule by molecule. The convergence test still acts on F = G(x) - x,
so the fixed point and the energy do not change; only the number of SCF iterations does.

## Mathematics

- Residual F(x) = G(x) - x, Jacobian J = G'(x) - I; Newton step dx = -J^-1 F.
- Modified Broyden (Johnson, PRB 38, 12807 (1988); Eyert, J. Comput. Phys. 124, 271 (1996)) updates an
  inverse-Jacobian estimate H_n from the history, with the seed H_0 = alpha*I. Here H_0 = M.
- Fragments: connected components of the bond graph r_ij < 1.25 (r_i + r_j), Pyykko-Atsumi covalent radii.
- For each fragment f: own neutral closed-shell GFN1/GFN2 SCF (method of the parent) of the isolated fragment (FP64), then
  J_f by central differences (h = 1e-4, 2 dim_f SCC evaluations) at its converged state; B_f = -J_f^-1.
- M = sum_f R_f^T D_f B_f R_f + C (additive Schwarz form; R_f gathers the fragment's shells, dipoles
  and quadrupoles from x, D_f = 1 because the blocks do not overlap, C = alpha on every index without a block).
- M*v costs O(sum_f dim_f^2); no dense dim x dim matrix is formed. Overlapping blocks: only the restricted form of
  the segmentation below (zero rows outside the own atoms), no weighted partition of unity.
- M is built at the first SCF of an XTB object that uses it and reused by every later SCF of that object
  (MD, optimisation steps); rebuilt for a new molecule or when the RMS displacement from the build geometry
  exceeds `-scf_precond_rebuild_rms`. Different XTB objects do not share M, except the worker calculators of one Hessian
  (section Hessian and MD); ConfSearch or parallel optimisations each build their own.
- Guard (default `g2:4:best`): with f_n = ||F_n||_2 of the full residual, if f_4 / f_0 >= 0.15, M is dropped
  for the rest of that SCF and the scalar step restarts from the stored input x_k with the smallest f_k
  (x_k + alpha F_k, no new map evaluation). `g1:<tau>` instead fires at n >= 2 if f_n > tau min_{k<n} f_k;
  `keep` keeps the Broyden history (u_k rewritten to the scalar form u_k - M dF_k + alpha dF_k), otherwise it is reset.
  The 0.15 is the worst scalar (no M) reduction over 24 prototype SCFs (6 systems x 4 geometries), rounded up.

## Options

| Parameter | Default | Meaning |
|---|---|---|
| `-scf_precond` | `none` | `none` or `frag`; other values: warning, `none` |
| `-scf_precond_guard` | `g2:4:best` | `g2:<n0>[:<ratio>][:keep][:best]`, `g1:<tau>[:keep][:best]` or `off`; only with `frag` |
| `-scf_precond_rebuild_rms` | `0` | rebuild M when the RMS displacement from its build geometry exceeds this (Angstrom); 0 = never |
| `-scf_precond_max_fragment_atoms` | `100` | molecules above this keep the scalar seed (one warning) or are segmented; 0 = no limit, negative: warning, 100. A block costs 2*dim SCC-map evaluations (see below) |
| `-scf_precond_segment_atoms` | `0` | > 0: segment molecules above `-scf_precond_max_fragment_atoms`, at least this many heavy atoms per segment, whole rings and multiple bonds (section Segmentation); 0 = off |
| `-scf_precond_segment_radius` | `2.0` | domain radius of a segment (Angstrom), default from the radius scan (section Segmentation); <= 0: warning, value kept |
| `-scf_precond_jacobian` | `fd` | cluster Jacobian: `fd` central differences, `response` linear response (section Response Jacobian); other values: warning, `fd` |

Developer environment variables (not for production runs):

| Variable | Effect |
|---|---|
| `CURCUMA_SCF_PRECOND_FILE=<file>` | dense M (dim x dim float64, row-major) from a file; wins over `-scf_precond` |
| `CURCUMA_SCF_PRECOND_GUARD=<spec>` | guard spec, wins over `-scf_precond_guard` (also used with the file) |
| `CURCUMA_SCF_JAC=<dir>` | after convergence, write J (central and forward differences), x*, G(x*), gamma, eta, multipole kernels, orbitals to `<dir>`; `CURCUMA_SCF_JAC_H`, `_LINTEST`, `_LINSCALE` tune it |
| `CURCUMA_SCF_JAC_COMPARE=1` | every cluster Jacobian of `-scf_precond frag` by both methods; stderr line `SCFJACCMP` (max abs and relative Frobenius deviation of J and of J^-1, times) and `SCFJACRESP` (phase times of the response) |
| `CURCUMA_SCF_TRACE=1` | per CPU SCF iteration: residual norms of the q and multipole blocks and phase times on stderr; `SCFTRACE_GUARD` when the guard fires |

`-verbosity 2` prints `SCF precond: frag, N fragments, dim D, built in T ms` and the block count, the
atoms left on the scalar seed and the memory of M; with segmentation also the number of segments, the mean cluster size
(atoms including caps) and the skipped segments.

## When it is ignored (one warning, scalar seed)

GFN1 on GPU paths other than the CUDA and ROCm resident loops; GPU SCF paths other than the fully device-resident CUDA and ROCm loops (host-driven device loops: Vulkan,
CUDA without the resident loop, ROCm with `-xtb.rocm_resident_loop false` or `-xtb.rocm_resident_precond false`); a multi-GPU SCF split (eigensolver or density over several devices at this basis
size); `-large_system_mode` other than `none`; `-scf_mode` other than `broyden`; total
charge != 0; open shell or odd electron count; no fragment qualifies for a block. Per fragment, single
atoms, fragments with an odd number of valence electrons and fragments with |sum of EEQ charges| > 0.5 e
(ions) keep the scalar seed. Fragments with elements other than H, C, N, O are used, with one warning
per run, because they are unmeasured. An unusable `CURCUMA_SCF_PRECOND_FILE` is a hard error; on a GPU path the
file is ignored with a warning. An ignored switch on a GPU run is not a counted GPU fallback, so `-gpu_strict true`
does not stop on it (nothing moves to the CPU).

## Measurements (all n=1 per entry)

Conditions: AMD Ryzen 9 9950X3D, `release/` build of this branch, `-sp -method gfn2 -threads 8`, CPU only,
default `scf_threshold` 1e-5. Iteration counts depend on the thread count. Structures: cuts of the
library structure `urea400-water1000`; `_relax` = GFN2 relaxed, `_ffrelax` = relaxed with GFN-FF.

| System (atoms) | none | frag | file M (prototype) | energy | build of M | memory of M (dense) |
|---|---|---|---|---|---|---|
| urea3-water8.gfn2-relaxed (48) | 11 | 6 | 6 | same 8 decimals | 15 ms | 0.23 MB (2.0 MB) |
| ureawater197_relax (197) | 12 | 8 | 8 | same 8 decimals | 19 ms | 1.05 MB (33.9 MB) |
| ureawater197, unrelaxed (197) | 24 | 33 (guard fired at 4) | - | same 8 decimals | 18 ms | 1.05 MB |
| ureawater599_ffrelax (599) | 12 | 8 | 8 | same 8 decimals | 49 ms | 2.9 MB (312 MB) |
| ureawater599, unrelaxed (599) | 58 | 62 (guard fired at 4) | 61 (guard) | same 8 decimals | 48 ms | 2.9 MB |

- Unrelaxed 599 atoms over thread counts 4/6/8/12/16: none 63/138/58/58/58, frag 61/62/62/55/61. The
  count at a single thread count is not a stable measure on this structure; without the guard the
  prototype needed 114 (none 58).
- Guard variants on the unrelaxed 599 (`-threads 8`): g2:4:keep 58, g1:3:best 58, g2:4:best 62,
  g2:4:keep:best 69, g2:4 109, g1:3:keep:best not converged within the iteration limit (hard error, no
  wrong result).
- NVE MD (`-time_step 1 -temperature 300 -thermostat none -seed 42 -threads 8`, guard default, mean iterations of
  SCFs 2..n, warm start): urea3-water8 50 fs 8.55 -> 5.00 (M built once; `-scf_precond_rebuild_rms 0.02`: 27 builds,
  mean 5.00); ureawater197 unrelaxed 50 fs 10.76 -> 7.34; ureawater197_relax 100 fs 9.04 -> 5.91; ureawater599
  unrelaxed 20 fs 13.35 -> 11.50. First-step energies identical to 8 decimals; later energies differ by up to
  4.3e-6 Eh from diverging trajectories, not from the fixed point.
- Hessian (`-hessian`, ureawater48_relax, `-threads 1`, `-scf_threshold 1e-7`, 288 SCFs, follow-up SCF mean / median /
  max): none 10.05/10/12, frag 5.80/6/7; 144 frequencies, max deviation from none 0.1 cm-1 (printed resolution).
- Default path: iteration tables and energies of water51, ureawater48, nacl64 (`-verbosity 3`) identical to
  the binary before the feature.
- Prototype (dense M from external scripts, same algorithm): MD 9.1 -> 6.0 iterations; Hessian warm starts
  6-7 -> 3-4; M from the start geometry over 500 MD steps without clear ageing (RMS displacement up to
  0.69 A); valid up to 0.10 A random displacement.

## GPU (CUDA resident loop)

Status: AI-generated, machine-tested only (Oct 11, 2026), human production testing pending.
On `-gpu cuda` the fully device-resident GFN2 loop mixes on the device (`XtbGpuContext::runBroydenUpdate`). With
`-scf_precond frag` the host builds M exactly as above (same fragment code), flattens the blocks
(`SchwarzPreconditioner::flatten`), remaps the indices to the device packing of the SCC vector (moments atom-major)
and uploads them once per build (`broydenSetPreconditioner`, reused by tag in later SCFs of the object, e.g. MD).
The device replaces alpha*F by M*F in the linear step and the seed and alpha*dF by M*dF in u (kernels `k_pc_diag`,
`k_pc_blocks`: one CUDA block per fragment block, block vector in shared memory, one warp per row; blocks up to 6144
entries). The guard keeps the CPU rule (`ScfPrecondGuard::firesNorm` as a host callback) on ||F||_2 of each update
(one `cublasDnrm2` scalar), keeps x_best/F_best as device copies and drops M like `dropStepHook`; no extra map
evaluation. Without the switch the device Broyden runs its unchanged kernels.

Measurements, RTX 5080, `-sp -method gfn2 -gpu cuda -threads 8 -verbosity 3`, default mixed precision, n=3 per cell
(iterations identical in all three runs; times median; energies identical to 8 decimals within each row and between
none and frag):

| System | none: iter / SCF ms / total s | frag: iter / SCF ms / total s | frag, guard off | build of M |
|---|---|---|---|---|
| ureawater599_ffrelax | 13 / 426 / 1.09 | 9 / 330 / 1.04 | 9 / 331 / 1.04 | 53-56 ms |
| ureawater1499_ffrelax | 13 / 1940 / 3.96 | 9 / 1504 / 3.72 | 9 / 1504 / 3.69 | 162-164 ms |
| ureawater599, unrelaxed | 60 / 4487 / 5.11 | 47 / 3480 / 4.18 (guard fired at 4) | not converged in 151 (hard error) | 52-56 ms |

- Upload of M: 0.7-3.9 ms (2.9 MB for 599 atoms). A repeat of the table on a busier machine (load 4-6) gave the same
  iterations and energies, SCF times up to 20 % higher.
- Default path unchanged: iteration tables (energy to 10 decimals, max|dq|) of the relaxed structures identical to
  the binary before the change in all runs (n=3 before, n=6 after); the unrelaxed 599 table varies from run to run
  already before the change (GPU not run-to-run reproducible), with 60 iterations and the same final energy every time.
- With mixed precision off (`-scf_mixed_precision false`, n=1): ureawater599_ffrelax 12 -> 8 and ureawater1499_ffrelax
  12 -> 8, i.e. the CPU counts (12 -> 8, 12 -> 7); unrelaxed 599 with guard 61 (CPU 62). The remaining iteration on the
  default GPU path comes from the FP32 phase.
- NVE MD 5 fs on ureawater599_ffrelax (`-max_time 5 -time_step 1 -temperature 300 -thermostat none -seed 42`, n=3):
  iterations per SCF none 13, 11, 11, 10-11, 10-11, 11; frag 9, 7, 7, 7, 7, 7 (M built once); first energy identical
  (-1037.75117514 Eh), later energies differ by up to 7e-7 Eh; total 4.45-4.65 s vs 3.98-4.08 s.
- Component test `sqm_cuda_broyden`: device Broyden with M vs host `BroydenMixer` with the step hook, with guard
  variants reset/keep/best, agree to 3e-15.
- ROCm (Oct 24, 2026; AI-generated, machine-tested on one gfx1200): the resident GFN1 and GFN2 loop applies the same M and guard (`-xtb.rocm_resident_precond`, default true; indices in the host `packSCC` order); iterations equal to the CPU loop in 15 of 16 test pairs (PR07 GFN2: 8 against 10), guard hits at the same iteration; numbers in [ROCM_RESIDENT_LOOP_WP.md](ROCM_RESIDENT_LOOP_WP.md) (WP5).
- CUDA GFN1 (Oct 25, 2026; AI-generated, machine-tested on one RTX A4500): the resident GFN1 loop takes M and the guard like GFN2 (`supportsResidentPreconditionerGfn1` in `xtb_gpu_method.cpp`; the device Broyden of the GFN1 loop has dimension nsh, M is layout-neutral). Closed-shell neutral clusters, `-gpu cuda -gpu_devices 0`, energies equal to the run without the switch and to the CPU run at 8 digits in 5 of 5 systems; SCF iterations plain / frag on the GPU: urea3-water8 (48 atoms) 11 / 8, water8 6 / 4, acetic acid dimer 12 / 7, urea-water 197 atoms 14 / 10, peo201-water500 (2910 atoms) 11 / 11 (CPU frag 7, 3, 6, 8, 10). Test `cli_sqm_13` checks GFN1 (device M, energy, 8 < 11 iterations). Not run: GFN1 with the guard firing, systems with charged fragments, multi-GPU SCF split.

## Tested and rejected (Oct 2026, all n=1)

On the unrelaxed ureawater197 `frag` does not help (24 -> 33 iterations, guard). Two extensions were measured and dropped:
- Fragment embedding with a frozen environment instead of isolated fragments: 33 -> 32 iterations.
- Two-level variant with a coarse level of one charge-transfer vector per molecule (Jc = Phi^T J Phi): with the exact J
  of the same geometry 13 instead of 31 iterations. Built by finite differences after the first SCF it ages within
  1 fs (MD ureawater197 unrelaxed 8.98, ureawater599 17.30 against frag 7.34 / 11.50 iterations per SCF). Built from the
  orbitals (sum over states) it is accurate with converged orbitals (1.7 % relative error of the coarse matrix) but
  useless with the orbitals of early iterations (single point ureawater197 unrelaxed 27 against 24 without M; MD 7.20 /
  11.60). Cause according to the measurements: HOMO-LUMO gap 0.23 eV with charge transfer between molecules, so the
  coarse matrix depends on the states at the Fermi level and changes with them. Code removed Oct 11, 2026; lab journal
  in the vault, 'curcuma GFN2-SCF-Vorkonditionierung', steps 7 to 10.

## Segmentation of large molecules (`-scf_precond_segment_atoms`, Oct 12, 2026)

Status: AI-generated, machine-tested only, human production testing pending. Code: `segmentFragment` and
`XTB::buildFragmentPreconditioner` in `xtb_precond.cpp`.

- Restricted additive Schwarz (Cai and Sarkis, SIAM J. Sci. Comput. 21, 792 (1999); Smith/Bjorstad/Gropp 1996): M = sum_i R_i^T D_i B_i R_i
  with overlapping domains, D_i = 1 on the components of the segment's own atoms and 0 elsewhere; own atoms are disjoint, so no
  component is counted twice.
- Only single bonds outside rings are cut (Oct 12, 2026, second version): a heavy-heavy bond is cuttable if it is a bridge of the heavy-atom
  bond graph (not in a ring) and r_ij / (r_cov,i + r_cov,j) >= 0.95 (measured ratios: single ~1.0, amide C-N 0.94, aromatic 0.92-0.93,
  C=O 0.88). Heavy atoms joined by uncuttable bonds form units (rings, C=O, amide groups) that are never split.
- Segments: units in breadth-first order (start at the smallest index), collected until a segment holds >= `scf_precond_segment_atoms`
  heavy atoms; each H joins its nearest heavy atom. Domain: every unit with a heavy atom within `scf_precond_segment_radius` of an own
  heavy atom, with the H. Cut bonds get an H cap on the bond axis at r_cov(atom) + r_cov(H).
- Each domain is a neutral cluster (real atoms, then caps) built like a fragment job (own XTB, central-difference Jacobian, all jobs in
  one thread pool, one result slot each); B_loc = -J_loc^-1 is restricted to the real atoms, rows of non-own atoms are set to zero
  (stored, so the CUDA copy is unchanged). Clusters with an odd electron count (valence electrons + one per cap) are skipped; their own
  atoms keep alpha (counted at verbosity 2; 0 in all runs of the second version). A domain that covers the whole molecule gives the exact block rows: with R = 100 A the
  SCF trace of PEO-58 is identical to the whole-molecule block (4 iterations).
- Measurements (n=1 each, 9950X3D, `-sp -method gfn2`, `-scf_precond_segment_atoms 3`, R = 3.0, `-scf_precond_max_fragment_atoms 20`):

| System | none | segmented | segments, mean cluster | build of M | energy |
|---|---|---|---|---|---|
| PEO 58 atoms (`peo8-chain.h-capped`), `-threads 2` / `1` | 11 | 5 / 5 | 8, 18.4 atoms | 0.18 s / 0.34 s | same 8 decimals |
| PEO 212 atoms, `-threads 8` / `1` | 12 | 6 / 6 | 30, 19.6 atoms | 0.20 s / 1.39 s | -311.85122954 Eh in all |
| PEO 212 atoms, `-gpu cuda` (RTX 5080), mixed precision / FP64 | 12 / 11 | 8 / 7 | as above | 0.20 s | same 8 decimals |

  PEO chain 1410 atoms (`peo201-chain`), `-threads 12`, n=3 alternating: none 12 iterations, SCF 10.2 to 10.5 s, wall 11.9 to 12.4 s; segmented 7
  iterations, SCF 5.7 to 5.9 s, M built in 0.95 s (202 segments), wall 8.4 to 8.5 s (-30 %), energy -2088.25340678 Eh in both. Short MD (5 fs, 12 threads,
  n=2): wall 72.8 / 71.1 s without, 50.4 / 49.2 s with (-30 %). Also at a single point of 212 atoms the build (0.2 s) costs more than the SCF saves (68 ms).
  The whole 212-atom block gave 4 iterations after 129.5 s of build (section Not measured).
- First version (cut any bond): polystyrene 20 units (325 atoms) and a 20-residue peptide (309 atoms) gave no gain (12 -> 13, 15 -> 14
  iterations), most segments skipped for an odd electron count (scalar seed on 196 / 156 atoms); cause: cuts through aromatic rings,
  C=O and amide C-N bonds (counted on the segments of k = 3, R = 3.0). Second version, `-threads 12`, n=1 each, machine under load
  (load average 28-37, build times only indicative); energies identical to the run without M to 8 decimals in every row; no segment
  skipped, no cluster SCF failed:

| System (atoms), none | k / R (A) | iterations | segments | mean cluster (atoms incl. caps) | build of M (s) |
|---|---|---|---|---|---|
| PEO 212, none 12 | 3 / 3.0, 4 / 3.0, 4 / 4.0 | 6, 6, 6 | 30, 23, 23 | 19.6, 20.8, 44.0 | 0.37, 0.44, 7.5 |
| polystyrene 325, none 12 | 3 / 3.0, 4 / 3.0, 4 / 4.0 | 5, 5, 5 | 21, 20, 20 | 42.0, 43.0, 69.7 | 6.3, 4.1, 26.9 |
| peptide 309, none 15 | 3 / 3.0, 4 / 3.0, 4 / 4.0 | 7, 7, 6 | 39, 32, 32 | 37.6, 41.0, 65.9 | 7.1, 15.5, 92.7 |
| peptide 309 | 8 / 3.0, 8 / 4.0 | 7, 6 | 18, 18 | 51.0, 80.0 | 16.1, 118.9 |
| peptide 309 | 12 / 3.0, 12 / 4.0 | 7, 6 | 13, 13 | 58.7, 90.3 | 30.8, 164.7 |

  Larger segments (k = 8, 12, about one residue) do not lower the iteration count further; R = 4.0 saves one iteration on the peptide at
  5 to 10 times the build cost. Rings and amide groups make the clusters about twice as large as on PEO, so the build costs seconds.
  Total cost (library structures `polystyrene20.gfnff-min`, `pep20-peptide.gfnff-min`, k = 3, R = 3.0, `-threads 12`, machine load about 29 from other jobs, n=3
  single points alternating, n=1 MD): single point polystyrene wall 1.2 s -> 3.2 to 3.9 s (SCF 1.0 -> 0.4 s, M 2.5 to 3.3 s), peptide 1.4 to 1.6 s -> 4.8 to 5.9 s
  (SCF 1.15 to 1.33 -> 0.54 to 0.62 s, M 4.0 to 5.1 s): a single point does not pay, break-even about 6 to 7 SCFs. MD (1 fs steps): 10 fs wall 11.8, 10.8 s ->
  10.3, 9.9 s (polystyrene, n=2) and 12.2 s -> 12.7, 13.0 s (peptide, n=2); 40 fs polystyrene 41.0 -> 26.6 s (-35 %), peptide 44.6 -> 32.3 s (-28 %). Energies equal.
  Radius scan (Oct 12, 2026, k = 2 or 3, R = 2.0 to 3.5, n=1, `-threads 12`, load 16 to 29; all energies equal): the build cost falls steeply with R, the iterations
  rise by at most 1. Polystyrene k = 3: R 3.0 -> 2.0 gives 5 -> 6 iterations, build 2.7 -> 0.32 s; peptide k = 3: R 3.0 -> 2.0 gives 7 -> 8, build 5.2 -> 0.87 s; peo212 R 2.0: 7, 0.09 s.
  Default radius changed from 3.0 to 2.0 A. Confirmed with R = 2.0, k = 3 (n=3 alternating single points, wall none -> segmented): polystyrene 1.35, 1.42, 1.35 s -> 1.20, 1.24,
  1.10 s (a single point already pays); peptide 1.59, 1.60, 1.43 s -> 1.95, 1.89, 1.83 s (build 0.87 to 0.92 s, break-even about 1.5 SCFs). MD 40 fs (n=1): peptide 44.6 -> 25.6 s,
  polystyrene 41.0 -> 22.5 s (-43 / -45 %). PEO 1410 atoms, R = 2.0 (n=2, loaded machine): 12 -> 8 iterations, build 0.54 s.
  Build cost split of a 58-atom cluster (1 thread, per SCC-map evaluation): eigensolver `dsyevd` 43 %, its reduction and back-transformation 15 %, potential 13 %, Fock 10 %, rest 11 %.
  Not measured: molecules other than PEO, this polystyrene and this peptide (both only as GFN-FF minima, not GFN2 relaxed); Hessian with a segmented M; charge-separated molecules.

## Response Jacobian (`-scf_precond_jacobian response`, Oct 12, 2026)

Status: AI-generated, machine-tested only, human production testing pending. Code: `XTB::computeScfJacobianResponse` in
`xtb_response.cpp` (independent of the preconditioner; GFN2 and GFN1, closed shell). `CURCUMA_SCF_JAC=<dir>` also writes `J_response.bin`.

- One eigensolve at v(x*), then per column k: dv = dv/dx_k by central differences of the potential alone (no eigensolve; exact for the
  linear terms and the quadratic third-order term, O(h^2) for D4), dF' = C^T dF C, dP'_pq = w_pq dF'_pq with
  w_pq = (f_p - f_q)/(eps_p - eps_q), for |eps_p - eps_q| < 1e-8 the mean of f'_p = -f_p (1 - f_p/2)/kT (Fermi-Dirac, f in [0,2]), on
  the diagonal dP'_pp = f'_p (dF'_pp - dmu) with dmu = sum f'_p dF'_pp / sum f'_p (fixed electron count), dP = C dP' C^T, Mulliken and
  multipole projection. Columns are batched into two GEMMs (C^T M and dP' C^T). Coupled-perturbed SCF form: Pople, Krishnan, Schlegel,
  Binkley, Int. J. Quantum Chem. Symp. 13, 225 (1979).
- Accuracy against central differences (h = 1e-4), n=1 per system: PEO 58 atoms (one block) max |dJ| 2.8e-9, relative Frobenius 5.7e-10,
  of J^-1 2.2e-10; urea3-water8 (48 atoms, whole system, dump) 6.7e-9 / 1.4e-9 / 6.0e-10; all segment clusters (k = 3, R = 2.0) of
  polystyrene 325 atoms (21 clusters, 15-24 atoms), peptide 309 atoms (39 clusters, 9-37 atoms; the 20-atom cluster 8.6e-9 / 1.6e-9 / 5.7e-10)
  and PEO 212 (30 clusters): worst 2.4e-8 / 4.1e-9 / 1.6e-9. The deviation is at the level of the finite-difference error of fd itself.
- SCF (`-threads 12`, k = 3, R = 2.0, n=1, machine under load): polystyrene fd 6 / response 6 iterations, build 296 / 107 ms; peptide 8 / 8,
  789 / 278 ms; PEO 212 7 / 7, 78 / 38 ms; energies identical to 8 decimals.
- Build time of one whole-molecule block (`-threads 1`, n=3 alternating, load average 18-20, so only indicative): PEO 58 atoms fd 2.3-5.0 s,
  response 1.2-1.4 s; PEO 100 atoms (cut of `peo201-chain`, not in the library) fd 27.8-47.5 s, response 5.3-5.6 s. Phase times of the
  response for PEO 58 (nao 130, dim 604): GEMMs 416 ms, potential differences 189 ms, rest 170 ms.

## GFN1 (Oct 12, 2026)

Status: AI-generated, machine-tested only, human production testing pending. GFN1 mixes only the shell charges (x = q_sh, dim = nsh;
`packSCC`), so a block holds only the shell indices of its atoms; fragments, segmentation and the response Jacobian are otherwise the same.
The GFN1 potential in x is Coulomb (gamma) plus the atom-resolved third order (v_at); the halogen-bond and D3 terms do not depend on the
charges and are not part of the map. CPU loop only: on a GPU path the switch is ignored with one warning.

- Response against central differences (`CURCUMA_SCF_JAC_COMPARE=1`, n=1, relative Frobenius of J / of J^-1): PEO 58 (one block)
  3.1e-10 / 2.2e-10; urea3-water8 whole system (dump, dim 96) 1.2e-9 / 7.2e-10, its 11 fragments <= 1.4e-9 / 9.6e-10; segment clusters
  (k = 3, R = 2.0) of polystyrene 325 (21) <= 1.2e-9 / 5.1e-10 and of the peptide 309 (39) <= 2.0e-9 / 9.0e-10.
- SCF (`-threads 12`, n=1, machine under load, energies identical to 8 decimals in every row), iterations none / frag fd / frag response,
  build fd / response: polystyrene 10 / 6 / 7, 248 / 147 ms; peptide 11 / 8 / 8, 489 / 247 ms; PEO 212 9 / 8 / 8, 163 / 165 ms (segments
  k = 3 with `-scf_precond_max_fragment_atoms 20`); urea3-water8 (11 fragments) 10 / 7 / 7, 94 / 120 ms.

## Census: frag against none on GMTKN55 and MOR41 (Oct 12, 2026)

Same binary, `-threads 1`, n=1 per structure, 2555 structures (open shells with `-spin`), 14 processes in parallel; single point
of ONE molecule each (a molecule is one fragment, its block is the Jacobian of the whole system). 624 structures ignore the switch
with a warning (odd electrons or open shell 233, charge 1 186, charge -1 128, charge 2 10, no fragment qualifies 67); 3 do not
converge in either mode (G21IP b+, be+, c+, single-atom cations). Energies: |E_frag - E_none| median 0, maximum 1.0e-8 Eh, none above 1e-6.

| atoms | n | SCF iterations median none / frag | wall sum none / frag (s) | build median (ms) | break-even (SCFs, median build / median saved per SCF) |
|---|---:|---|---|---:|---:|
| < 10 | 783 | 10 / 4 | 20.2 / 25.3 (1.25x) | 14 | ~2 |
| 10-30 | 958 | 11 / 4 | 28.0 / 74.7 (2.7x) | 37 | ~8 |
| 30-100 | 252 | 13 / 4 | 12.5 / 280.8 (22.5x) | 409 | ~62 |

The iteration count drops by more than 2 in 1888 of 2552 structures and rises by more than 2 in 3. A single point of one molecule is
slower with frag in almost all cases (faster in 138 of 1928 structures with a block); the build of one block of ~60 atoms takes 7 s
(C60ISO, PR08). MOR41: 18 -> 6 iterations, wall 5.6 -> 105.8 s. frag therefore pays off only for repeated SCFs (MD, optimisation,
Hessian, ConfSearch) of the same molecule, not for single points; it is not the default. Data: vault `_artefacts/2026-10-12_curcuma_scf-census2_*`.

## Not measured or not implemented

- Not measured: other elements than H, C, N, O; solvation (fragments are computed in gas phase); systems
  above 599 atoms; Hessian with the built-in M (each displaced-geometry object builds its own M; only the prototype with a file matrix was measured, up to 197 atoms); MD longer than
  500 fs; RMS drift above 1 A; GPU other than one RTX 5080 (section GPU); multi-GPU nodes.
- Fragment size (`-scf_precond_max_fragment_atoms`, default 100, 0 = no limit): a fragment Jacobian costs 2*dim evaluations of its SCC map
  (dim = nsh + 9*atoms for GFN2, nsh for GFN1). Measured on one thread, PEO chains as ONE molecule (`CURCUMA_SCF_JAC`, same algorithm, n=1): 58 atoms 1.4 s,
  212 atoms 120 s (about N^3.4; 100 atoms about 10 s and 300 atoms several minutes are extrapolated, not measured). A large molecule pays
  off only in a long MD or Hessian run. Measured once (n=1, `-threads 8`, PEO 212 atoms as ONE fragment, `-scf_precond_max_fragment_atoms 300`):
  M built in 129.5 s (one fragment is built serially), SCF 12 -> 4 iterations (164 -> 55 ms), i.e. break-even after about 1200 SCFs. The
  block of a single-fragment system is the exact Jacobian of the whole system. The iteration gain for ordinary fragments was measured up to ~30 atoms.
- Not implemented: GFN1 on CUDA; Vulkan and the host-driven CUDA and ROCm loops; multi-GPU SCF split; open shells; charged systems and ionic fragments (NaCl needs
  overlapping neighbourhoods of about R = 5 A for its soft dipole mode); overlapping Schwarz blocks between molecules (segmentation
  overlaps only inside one molecule).
- The guard default is a rule from the prototype test set (24 scalar SCFs); on the unrelaxed 599-atom
  cluster it costs 4 iterations against `none` at 8 threads.

## Hessian and MD with segmentation (Oct 12, 2026)

Campaign by an agent on a frozen binary (commit `773f8e17`, GFN2, `-threads 12`, machine load 17 to 28, wall times indicative only, n=2 per variant, one structure per system, seed 7);
iteration counts re-derived from the raw SCFTRACE files by the operator. Segmentation k = 3, R = 2.0.

- MD NVE 200 fs (`-thermostat none -print 1`): peptide 309 atoms 10.99 -> 6.00 SCF iterations per step, wall 722 / 466 s -> 286 / 302 s; polystyrene 325 atoms 9.97 -> 5.00, wall 406 / 420 s -> 264 / 259 s.
  No ageing of M over 200 fs (blocks of 40 steps constant at 6.0 and 5.0). Energies: seg against none differ by at most 1.1e-5 Eh in Etot (peptide) and 2e-6 Eh (polystyrene,
  output resolution 1e-6); the 4 to 9 mEh fluctuation of Etot is the same in both variants (dt 1 fs, start at a GFN-FF minimum).
- Hessian (semi-numerical, 12 threads): iterations per displaced SCF 6.6 to 7.2 -> 3.2 to 3.8; wall time urea3-water8 -6 to -8 %, 82-atom polystyrene -15 to -23 %, 67-atom peptide no gain
  (seg 8.9 / 9.1 s against 7.2 / 8.5 s). Frequencies differ from none by 0.3 to 0.8 cm^-1 (SCF threshold effect; none against none is exactly 0.0).
- Finding: the SCF count of the Hessian runs shows one build of M per worker thread (segmented PS5: 60 extra SCFs = 12 threads x 5 segments; PEP5: 96 = 12 x 8), not one build per Hessian.
  The M of a worker is built from the same geometry, so sharing M across the workers would remove 11 of 12 builds.
- Fixed (Oct 12, 2026, AI-generated, machine-tested): the Hessian builds M once at the reference geometry before the workers start (one extra
  SCF, not on GPU runs) and keeps that calculator alive; the workers' XTB objects take M over from a process-wide cache
  (`XTB::adoptSharedPreconditioner`, used only by objects the Hessian marks with the internal key `xtb._scf_precond_share`, only for
  `gfn1`/`gfn2`, so other parallel users of XTB keep their own M and do not depend on thread timing; key = method, charge, atoms and every build parameter; geometry within `scf_precond_rebuild_rms`, or
  0.01 A when that is 0; weak references, so an entry lives only as long as an object holds M). Measured (`-threads 12`, n=3 alternating, old
  binary `773f8e17`-based frozen copy): PS5 60 -> 5 cluster SCFs, wall 4.04 / 4.22 / 4.34 -> 3.80 / 4.05 / 3.80 s; PEP5 96 -> 8, wall 3.18 / 3.19 ->
  2.59 / 2.59 s (third pair under a load spike: 6.0 / 19.4 s); iterations per displaced SCF unchanged (3.18 / 3.73). Frequencies old against new
  0.8 (PS5) and 0.6 (PEP5) cm^-1 at the default threshold (M now from the reference geometry instead of each worker's first displacement; the
  same size as none against seg); with `-scf_threshold 1e-7` PEP5 old against new 0.000 cm^-1. Test `cli_sqm_12` (l).
- Not measured: other conformers; GPU Hessian (there each worker still builds its own M).

## GPU with segmentation (Oct 12, 2026)

Campaign by an agent (frozen binary and CUDA plugin of `c19047e3`, RTX 5080, resident loop, `-gpu_strict true`, n=3 alternating, machine load 26 to 164 from other jobs, so wall times are noisy by up to a factor 3).
Iterations (none -> seg): urea3-water8 12 -> 8 (frag), peo212 12 -> 8, polystyrene 325 atoms 14 -> 8, peptide 309 atoms 16 -> 9, PEO 1410 atoms 12 -> 10 (seg and seg+response). All energies equal to 1e-8 Eh, no fallback,
no "ignored" message. M on the device 0.22 to 20.4 MB (chain), device memory unchanged within the nvidia-smi resolution (2730 against 2750 MiB). M is built on the CPU, 0.3 to 7 s depending on load: a single point is slower in wall time
(peptide 2.4 -> 5.3 s, chain 15.4 -> 22.5 s under load). In the one low-load pair (chain, n=2) the SCF time fell by 11 % (1621 / 1592 ms -> 1413 / 1420 ms seg, 1428 / 1421 ms seg+response) and the build took 750 to 804 ms (fd) or 288 to 297 ms
(response). MD 40 fs: chain 501 s -> 527 s (n=1 each, load 35), peptide no difference (n=5): 40 fs is too short to amortise the build. Not measured: a quiet machine, MD longer than 40 fs, building M on the device, Hessian on the GPU
(the Hessian builds M once per worker there, see above), GFN1 on the GPU (not supported).

Repeat on a quiet machine (Oct 13, 2026, current HEAD, RTX 5080, resident loop, `-gpu_strict true`, load 0.2 to 3 for the single points, n=3 alternating, all energies equal to 1e-8 Eh, no fallback):
iterations none -> seg: peo212 12 -> 8, polystyrene 14 -> 8, peptide 16 -> 9, chain 12 -> 10; SCF time 170 -> 146, 266 -> 203, 281 -> 206, 1559 -> 1373 ms; build of M 54 / 23 (peo212), 199 to 221 / 75 to 78 (polystyrene),
561 to 568 / 153 to 167 (peptide), 330 / 138 ms (chain) with fd / response. Wall time of a single point with seg+response against none: peo212 0.49 / 0.50 to 0.54 s, polystyrene 0.65 / 0.66 s, peptide 0.73 to 0.76 / 0.65 to 0.67 s (+15 %),
chain 3.05 to 3.16 / 3.14 to 3.30 s (up to -7 %). MD 100 fs (seg+response, n=2, load 3 to 14): chain 285.7, 271.9 s -> 259.6, 242.0 s (-9 / -11 %), peptide 25.9, 26.0 s -> 20.4, 20.4 s (-21 %); end energies equal within 1e-4 Eh (chain, GPU scatter)
and 1e-6 Eh (peptide). The earlier statement "no gain in wall time" refers to the loaded campaign above and is superseded for MD.
