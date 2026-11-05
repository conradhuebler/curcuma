# Native QM SCF: level shift and wavefunction stability analysis (Nov 4, 2026)

> 🤖 AI-generated, ⚙️ machine-tested. One machine (4 vCPU container, no GPU), n = 1 per number, not reviewed by the
> operator. Human production testing is pending.

Native ab-initio SCF (`QMEngine`: `hf`, `hf-3c`, `lda`/`pbe`/`b3lyp`/`tpss`/`scan`, closed shell `qm_scf.cpp`, open
shell `qm_uscf.cpp`). Both options are opt-in; with them off the SCF is bitwise the one before (checked below).

## Level shift (`-qm.scf_level_shift b`, `-qm.scf_level_shift_off g`, `-qm.scf_mode level-shift`)

- Diagonalised matrix F + b (S - S P_s S) (closed shell S - 1/2 S P S; Saunders/Hillier 1973), applied after the
  DIIS/ADIIS extrapolation (as PySCF). DIIS error, energy and the three convergence tests use the unshifted F.
- When a shifted iteration meets every test, the unshifted (extrapolated) F of that iteration is diagonalised once more
  and its orbitals, orbital energies and density are returned (the rule of the run without shift). A first version
  continued the loop unshifted instead: undamped Roothaan steps (ozone) left the converged point again.
- `scf_level_shift_off g > 0` removes the shift once ||g||_F < g. `scf_mode level-shift` = diis + b = 0.2 Eh.
- Host loops only. With `-gpu` the device-resident RHF loop is skipped and counted (`qm_scf`); not run (no GPU here).
- Code: `ScfLevelShift` in `qm_scf_internal.h`, used by `runSCF` and `runUSCF`.

Result does not depend on b (dump_qm_energy, scf_threshold 1e-10, def2-SVP, 2 threads):

| case | b (Eh) | iterations | E - E(b = 0) (Eh) | max dg (Eh/Bohr) |
|---|---|---|---|---|
| water HF, gradient | 0 / 0.2 / 0.5 / 1.0 | 15 / 20 / 26 / 40 | 0 / 0 / 0 / -4e-14 | - / 1.2e-12 / 1.3e-12 / 1.3e-12 |
| water B3LYP sg1, gradient | 0 / 0.5 | 12 / 35 | 4e-14 | 1.4e-13 |
| CH3 UHF, gradient | 0 / 0.5 | 18 / 27 | -3e-14 | 1.2e-15 |
| ozone HF | 0 / 0.3 / 1.0 | 24 / 30 / 53 | 0 / -5e-13 / 0 |  |

Where it helps and where not (iterations to scf_threshold 1e-10, "fail" = not converged in the limit):

| case | without | b = 0.3 | b = 0.5 | b = 1.0 |
|---|---|---|---|---|
| ozone HF, `scf_mode plain`, `scf_damping 1` (undamped Roothaan), 150 it. | fail | 76 | - | 127 |
| N2 at 2.0 A HF, plain, damping 1, 150 it. | fail | fail | - | 51 |
| water, O-H 1.8 A HF, plain, damping 1, 200 it. | fail | - | fail | fail |
| ozone UHF (`-qm.reference uhf -qm.uhf_guess_mix true`), DIIS, 200 it. | 60 | 99 | 56 / 73 / 70 (off 1e-1 / 1e-2 / 1e-3) | fail |
| ozone RHF, DIIS | 24 | 30 | 24 / 24 / 26 (off 1e-1 / 1e-2 / 1e-3) | 53 |

With the default DIIS no tested case failed without the shift, and the shift always cost iterations unless it was
removed early (`scf_level_shift_off`). It rescues oscillating undamped Roothaan iterations. No DIIS failure case was found
in this test set, so the benefit for DIIS runs is not demonstrated.
