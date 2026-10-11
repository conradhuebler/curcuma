# Fragment SCF preconditioner for the native GFN2 Broyden mixer (`-scf_precond frag`)

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
- For each fragment f: own neutral closed-shell GFN2 SCF of the isolated fragment (FP64), then
  J_f by central differences (h = 1e-4, 2 dim_f SCC evaluations) at its converged state; B_f = -J_f^-1.
- M = sum_f R_f^T D_f B_f R_f + C (additive Schwarz form; R_f gathers the fragment's shells, dipoles
  and quadrupoles from x, D_f = 1 because the blocks do not overlap, C = alpha on every index without a block).
- M*v costs O(sum_f dim_f^2); no dense dim x dim matrix is formed. Overlapping blocks (with weights D_f)
  fit the same storage, not implemented.
- M is built at the first SCF of an XTB object that uses it and reused by every later SCF of that object
  (MD, optimisation steps); rebuilt for a new molecule or when the RMS displacement from the build geometry
  exceeds `-scf_precond_rebuild_rms`.
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

Developer environment variables (not for production runs):

| Variable | Effect |
|---|---|
| `CURCUMA_SCF_PRECOND_FILE=<file>` | dense M (dim x dim float64, row-major) from a file; wins over `-scf_precond` |
| `CURCUMA_SCF_PRECOND_GUARD=<spec>` | guard spec, wins over `-scf_precond_guard` (also used with the file) |
| `CURCUMA_SCF_JAC=<dir>` | after convergence, write J (central and forward differences), x*, G(x*), gamma, eta, multipole kernels, orbitals to `<dir>`; `CURCUMA_SCF_JAC_H`, `_LINTEST`, `_LINSCALE` tune it |
| `CURCUMA_SCF_TRACE=1` | per CPU SCF iteration: residual norms of the q and multipole blocks and phase times on stderr; `SCFTRACE_GUARD` when the guard fires |

`-verbosity 2` prints `SCF precond: frag, N fragments, dim D, built in T ms` and the block count, the
atoms left on the scalar seed and the memory of M.

## When it is ignored (one warning, scalar seed)

GFN1; GPU-resident SCF; `-large_system_mode` other than `none`; `-scf_mode` other than `broyden`; total
charge != 0; open shell or odd electron count; no fragment qualifies for a block. Per fragment, single
atoms, fragments with an odd number of valence electrons and fragments with |sum of EEQ charges| > 0.5 e
(ions) keep the scalar seed. Fragments with elements other than H, C, N, O are used, with one warning
per run, because they are unmeasured. An unusable `CURCUMA_SCF_PRECOND_FILE` is a hard error.

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

## Not measured or not implemented

- Not measured: other elements than H, C, N, O; solvation (fragments are computed in gas phase); systems
  above 599 atoms; Hessian above 48 atoms (each displaced-geometry object builds its own M); MD longer than
  500 fs; RMS drift above 1 A; GPU.
- Not implemented: GFN1; GPU/ROCm paths; open shells; charged systems and ionic fragments (NaCl needs
  overlapping neighbourhoods of about R = 5 A for its soft dipole mode); overlapping Schwarz blocks.
- The guard default is a rule from the prototype test set (24 scalar SCFs); on the unrelaxed 599-atom
  cluster it costs 4 iterations against `none` at 8 threads.
