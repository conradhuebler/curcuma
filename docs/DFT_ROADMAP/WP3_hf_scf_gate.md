# WP3 — HF-SCF (rung 666) — hartes ERI-Gate

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP2
- **Validierung:** ORCA HF + xcDFT rung=666 (zwei Referenzen)

## Ziel

Vollständiger geschlossenschaliger SCF mit Löwdin `X=S^-1/2`, verallg.
Eigenproblem (Cholesky wie `xtb_scf.cpp:199`), Dichte, Fock = Hc+J+K, [F,P]_S-
Konvergenz, DIIS (`diis_accelerator.h`), Energiekomponenten (ET/EV/EJ/Ex/EKS+ENuc).
Validiert den ganzen ERI/J/K-Pfad gegen ORCA HF — entscheidendes Gate.

## Deliverables

- [ ] `dft_scf.cpp`: SCF-Schleife (xcDFT `RKS.f90` als Struktur-Vorbild, zitiert);
      eigensolver-Dispatch via `curcuma::eigsolver` (mkl/native/purify/lobpcg)
      analog `xtb_scf.cpp:229/302/327/351`; DIIS mit konfigurierbarer History
      (xcDFT hat nur n_diis=1 → erweitern); Energiekomponenten + Orbitalenergien.
- [ ] `DFT::Calculation(gradient=false)`: HF-Pfad; `getEnergyDecomposition`,
      `getOrbitalEnergies/Occupations`.

## Erfolgskriterien

- [ ] ORCA `! HF def2-SVP TightSCF` auf He/Be/Ne/H2O/CH4: E_KS+ENuc ≤1e-6 Eh,
      jede Komponente (ET/EV/EJ/Ex) ≤1e-6; HOMO/LUMO ≤1e-5 Eh.
- [ ] xcDFT rung=666 auf He VDZ: ≤1e-8 (gespiegelte Basis).
- [ ] SCF konvergiert <50 Iter auf DIIS; Plain-SCF (kein DIIS) ≤1e-7 als Fallback.

## Quellen & Lektüre

- curcuma SCF-Vorbild (Hauptlektüre): `xtb_scf.cpp:1` — Löwdin/Cholesky `:199`,
  `solveEigen` `:75` (LAPACK `dsyevd_` `:35-37` vs `Eigen::SelfAdjointEigenSolver`
  `:99`), eigensolver-Dispatch `:229` (purify), `:287` (`m_external_eigensolver`
  GPU-Hook), `:302-327` (native/lobpcg), `:351` (`solveSymmetric`);
  MKL-Scope `xtb_native.h` (`MklSerialScope`/`MklThreadScope`).
- Wiederverwendbare Eigenlöser: `native_eigensolver.h:38` — `solveSymmetric:56`,
  `purifyDensity:81`, `lobpcgLowest:118` (generisch, direkt nutzbar).
- DIIS/Mixer: `diis_accelerator.h`, `broyden_mixer.h` (statt xcDFT
  `DIIS_extrapolation.f90` n_diis=1).
- ComputationalMethod-Hooks: `computational_method.h:70` (`calculateEnergy`),
  `:211/217/223` (OrbitalEnergies/Occupations/NumElectrons), `:241`
  (`getEnergyDecomposition`), `:282-285` (Warm-Start/IterativeMode).
- xcDFT-SCF-Vorbild (Struktur + Komponenten, zitiert): `RKS.f90:1` —
  Per-Iter-Flow `:111-200`, Dichte `:131` (`P=2·c·cᵀ`), Fock-Build `:154`
  (`F=Hc+J+Fx+Fc`), Konvergenz `[F,P]_S` `:158-159`, Energiekomponenten
  `:169-189` (ET/EV/EJ/Ex/Ec/EKS), `:208-218` Abbruch, `:222` `print_RKS`;
  `orthogonalization_matrix.f90:1` (Löwdin `X=S^-1/2`); `wrap_lapack.f90:1`
  (`dsyev`), `utils.f90` (`trace_matrix`, `AtDA/ADAt`).

## Validierungsergebnisse (nach Ausführung eintragen)

```
ORCA HF def2-SVP:  E-Komponent (He/Be/Ne/H2O/CH4):
  ΔE_total:        [ ]
  ΔET/EV/EJ/Ex:    [ ]
  ΔHOMO/LUMO:      [ ]
xcDFT rung=666 He VDZ:  Δ: [ ]
SCF Iter (DIIS):   [ ]  Plain-SCF fallback Δ: [ ]
```