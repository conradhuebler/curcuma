# Umbrella restraints in SimpleMD (`-simplemd.umbrella_q`, `-simplemd.umbrella_cv`)

AI-generated (Oct 2026), machine-tested only (unit test + ethane GFN-FF runs); human production testing pending.

## `umbrella_q`: q = d(i,j) - d(i,k)

- **Purpose**: harmonic umbrella on q = d(i,j) - d(i,k), e.g. the antisymmetric stretch of an identity transfer (i = C, j/k = leaving/entering O); an independent reference for NEB-MD path-CV free energies ([NEB_MD.md](NEB_MD.md)).
- **Format**: `-simplemd.umbrella_q "i,j,k,q0,kappa"`; atoms 1-based and distinct, q0 in Angstrom, kappa > 0 in kcal/mol/A^2 (converted to Eh/A^2 with `CurcumaUnit::Energy`). Bad input stops `Initialise()` with an error.
- **Potential**: W = kappa/2 (q - q0)^2, analytic gradient (`SimpleMD::umbrellaQEnergy`), installed through `SimpleMD::setExtraPotential`, so it acts in both velocity half-kicks and is part of `trueGradient()`.
- **Output**: `<basename>.umbrella_q.csv` (BMT dir or CWD with `-no_bmt`), columns `step,time_fs,q_A,d_ij_A,d_ik_A,W_Eh,Epot_Eh,Ekin_Eh`, every `umbrella_q_frequency` steps (default 20) plus step 0. Epot excludes W; Epot + Ekin + W is the NVE conserved quantity.
- **WHAM**: run one MD per window q0 (same kappa), take the `q_A` column of each window and unbias with w_i(q) = kappa/2 (q - q0_i)^2 (same 1/2 convention as Grossfield's `wham`); discard equilibration lines.
- **Limits**: plain `-md` only; NEB-MD's own extra potential (z wall) and umbrella_q exclude each other (error). One restraint per run. The csv is overwritten on restart. Not tested: GPU path, restart files, thermostats other than CSVR/none, RATTLE on the restrained atoms.
- **Test**: `ctest -R simplemd_umbrella_q` (gradient vs central FD, W = grad = 0 at q = q0, parser rejections).

## `umbrella_cv`: sum of two distances, signed dihedral (Oct 20, 2026)

- **Formats** (atoms 1-based, W = kappa/2 (q - q0)^2, kappa > 0):
  - `sum,i,j,k,l,q0,kappa`: q = d(i,j) + d(k,l) (i != j, k != l; a shared atom is allowed), q0 in Angstrom, kappa in kcal/mol/A^2; e.g. the synchronous coordinate d(C1-C6) + d(C4-C5) of a concerted retro-Diels-Alder.
  - `dihedral,i,j,k,l,phi0,kappa`: signed dihedral phi(i,j,k,l), phi0 in degrees, kappa in kcal/mol/rad^2; dphi = phi - phi0 is wrapped to [-180, 180] deg (periodic). Atoms distinct.
  - `diff,i,j,k,q0,kappa`: synonym of `umbrella_q` (same code, same numbers).
- **Dihedral convention**: IUPAC/Klyne-Prelog, as `gfnff_geometry.h` `calculateDihedralAngle`: looking along j -> k, phi > 0 if bond j-i turns clockwise onto k-l; cis 0, gauche +60/-60, trans 180. Distances are mirror invariant, the sign of phi is not: only the dihedral separates enantiomers. `Molecule::CalculateDihedral` uses a different convention.
- **Gradient**: analytic (`SimpleMD::umbrellaCVEnergy`, `SimpleMD::dihedralAngle`, Blondel-Karplus form); a dihedral with a bond angle i-j-k or j-k-l whose |sin| < 1e-6 stops the run with an error.
- **Output**: `<basename>.umbrella_cv.csv` via `outputPath()`, every `umbrella_cv_frequency` steps (default 20) plus step 0; columns `step,time_fs,cv_A,d_ij_A,d_kl_A` (sum; `d_ik_A` for diff) or `step,time_fs,cv_deg,dphi_deg` (dihedral; cv in [-180, 180]), then `W_Eh,Epot_Eh,Ekin_Eh`.
- **Collisions**: `umbrella_cv` and `umbrella_q` exclude each other (error); both exclude NEB-MD's extra potential (z wall), `setExtraPotential` throws. One restraint per run.
- **Use as 1-D reference**: one MD per window q0 (same kappa), WHAM on the `cv` column with w_i = kappa/2 (q - q0_i)^2; for a dihedral, convert to radians and
  take q - q0_i wrapped (periodic WHAM). The resulting 1-D free energy is an independent check of an NEB-MD path-CV profile, not a replacement.
- **Limits**: a dihedral umbrella needs windows around the whole circle (or the range of interest plus overlap at the ends). For a helicene enantiomerisation the dihedral must run through the inner rim; which four atoms is the user's choice, the code does not check that the coordinate separates the two enantiomers.
- **Test**: `ctest -R simplemd_umbrella_q` also covers `umbrella_cv` (FD gradient on random 6-atom structures, W = grad = 0 at the minimum, wrap at 175/-175 and +-180, sign convention, collinear guard, parser rejections). Not tested: GPU path, restart, thermostats other than CSVR/none (BMT output only smoke-tested, 100 fs).
