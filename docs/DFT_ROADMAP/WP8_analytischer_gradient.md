# WP8 — Analytischer Gradient

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP7 (alle Functionals) + WP2/WP4
- **Validierung:** ORCA EnGrad + FD

## Ziel

Analytischer Kerngradient (Eh/Bohr) für `-opt`/`-md`; Pulay 1e + 2e-ERI-
Ableitung + Grid-XC-Ableitung + Becke-Gewichtsableitung.

## Deliverables

- [ ] `dft_gradient.cpp`: dS/dR, dT/dR, dV/dR (1e-Pulay, Obara-Saika-Derivat),
      d(ERI)/dR (MD raised-Hermite / Obara-Saika-Derivat → dJ/dR, dK/dR für Hybrid),
      dV_xC/dR (AO-Gradienten aus WP4 + Funktionalableitung + Becke-Gewichtsableitung,
      Becke 1988), Kernabstoßung. `hasGradient()=true`; `copyGradientTo`.
- [ ] FD-Validierung gegen `ComputationalMethod::NumGrad`.

## Erfolgskriterien

- [ ] FD-Check (NumGrad dx=1e-4): max|g_analytisch − g_FD| ≤1e-6 Eh/Bohr für
      LDA/PBE/B3LYP auf H2O/CH4.
- [ ] ORCA `! <functional> def2-SVP EnGrad` Gradientkomponenten ≤1e-5 Eh/Bohr.
- [ ] `-opt`-Smoke: H2O optimiert, konvergiert (E sinkt monoton), Gradientennorm <Schwelle.

## Quellen & Lektüre

- curcuma analytischer Gradient (Hauptvorbild): `xtb_gradient.cpp` — Aufbau
  (Repulsion + H0/Pulay + Coulomb + CN + Multipole), `calculateGradient()`
  füllt `m_gradient` in Eh/Bohr (`xtb_native.h:985`); Wrapper-Hooks
  `computational_method.h:80` (`getGradient`), `:121` (`copyGradientTo` —
  überschreiben), `:142` (`NumGrad` FD-Fallback), `:129` (`hasGradient`).
- 1e-Pulay: Ableitung von S/T/V nach Kernposition — Obara-Saika-Derivat-Rekursion
  (gleiche Quelle wie WP1); Theorie: Yamaguchi et al. *A New Dimension to QC*,
  Pulay, Mol. Phys. 17, 197 (1969) (Pulay-Kräfte).
- 2e-Ableitung: ERI-Ableitung via MD raised-Hermite-Indizes oder Obara-Saika
  Transfer-Derivat (Helgaker/Jørgensen/Olsen Kap. 9.4); J/K-Gradient wie HF-Pulay.
- Grid-XC-Ableitung: AO-Gradienten aus WP4 (`AO_values_grid.f90:78-90` portiert);
  Becke-Gewichtsableitung (Becke 1988, d w_i/d R_A); Funktionalableitung nach ρ/∇ρ
  (Standard KS-DFT-Gradient, Johnson/Gill-Pople, JCP 98, 5612 (1993)).
- xcDFT hat **keine** Gradienten (Doku vermerken als curcuma-native Erweiterung);
  einzige Ableitungs-Spuren: kartesische AO-Gradienten `AO_values_grid.f90:78-90`.
- `-opt`-Integration: `main.cpp:1709` `EnergyCalculator`, `:1724` Warm-Start
  (DFT hinzufügen), Optimizer-Konsum von `getGradient`.

## Validierungsergebnisse (nach Ausführung eintragen)

```
FD vs analytisch (LDA/PBE/B3LYP, H2O/CH4) max Δg: [ ]
ORCA EnGrad Komponenten max Δ: [ ]
-opt Smoke H2O konvergiert: [ ]
```