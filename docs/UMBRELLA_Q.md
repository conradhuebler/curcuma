# Umbrella restraint on a distance difference (`-simplemd.umbrella_q`)

AI-generated (Oct 2026), machine-tested only (unit test + ethane GFN-FF runs); human production testing pending.

- **Purpose**: harmonic umbrella on q = d(i,j) - d(i,k), e.g. the antisymmetric stretch of an identity transfer (i = C, j/k = leaving/entering O); an independent reference for NEB-MD path-CV free energies ([NEB_MD.md](NEB_MD.md)).
- **Format**: `-simplemd.umbrella_q "i,j,k,q0,kappa"`; atoms 1-based and distinct, q0 in Angstrom, kappa > 0 in kcal/mol/A^2 (converted to Eh/A^2 with `CurcumaUnit::Energy`). Bad input stops `Initialise()` with an error.
- **Potential**: W = kappa/2 (q - q0)^2, analytic gradient (`SimpleMD::umbrellaQEnergy`), installed through `SimpleMD::setExtraPotential`, so it acts in both velocity half-kicks and is part of `trueGradient()`.
- **Output**: `<basename>.umbrella_q.csv` (BMT dir or CWD with `-no_bmt`), columns `step,time_fs,q_A,d_ij_A,d_ik_A,W_Eh,Epot_Eh,Ekin_Eh`, every `umbrella_q_frequency` steps (default 20) plus step 0. Epot excludes W; Epot + Ekin + W is the NVE conserved quantity.
- **WHAM**: run one MD per window q0 (same kappa), take the `q_A` column of each window and unbias with w_i(q) = kappa/2 (q - q0_i)^2 (same 1/2 convention as Grossfield's `wham`); discard equilibration lines.
- **Limits**: plain `-md` only; NEB-MD's own extra potential (z wall) and umbrella_q exclude each other (error). One restraint per run. The csv is overwritten on restart. Not tested: GPU path, restart files, thermostats other than CSVR/none, RATTLE on the restrained atoms.
- **Test**: `ctest -R simplemd_umbrella_q` (gradient vs central FD, W = grad = 0 at q = q0, parser rejections).
