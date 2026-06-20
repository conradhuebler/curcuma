# WP5 — LDA (Slater-Dirac + VWN5)

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP4
- **Validierung:** ORCA LSDA + xcDFT rung=1 (zwei Referenzen) + FD

## Ziel

LDA-Exchange + LDA-Korrelation (VWN5, **xcDFTs fehlerhafte Formel korrigiert**);
V_XC-Matrix auf Grid; LDA-SCF.

## Deliverables

- [ ] `dft_xc.hpp/.cpp`: Slater-Dirac-Exchange (xcDFT `lda_exchange_*.f90`-Port,
      zitiert); VWN5-Korrelation (Vosko-Wilk-Nusair 1980 — ersetzt xcDFTs
      fehlerhaftes `lda_correlation_*.f90`, dokumentiert warum); Ex/Ec-Energie +
      Vxc-Potentialmatrix (ρ^(1/3), ρ^(4/3)-Skalare, Summe über Grid).
      Functional-Dispatch `select_rung`-Modell.
- [ ] SCF LDA-Pfad in `dft_scf.cpp` (Fock = Hc + J + Vxc). LDA = Method-Name
      `-method lda` (kein `-dft.functional`).

## Erfolgskriterien

- [ ] ORCA `! LSDA VWN def2-SVP` auf He/Be/Ne/H2O/CH4: E ≤1e-6 Eh, Komponenten
      (Ex/Ec) ≤1e-6.
- [ ] xcDFT rung=1 He VDZ: ≤1e-8 (Komponenten) bei gespiegelter Basis + SG-3.
- [ ] `dE/dρ`-Konsistenz: numerische Ableitung Ex vs analytische Vxc (FD) ≤1e-6.

## Quellen & Lektüre

- xcDFT LDA-Exchange (portieren): `lda_exchange_energy.f90:12-18` (Dirac-Koeff
  `C=-0.5^(1/3)·(3/2)·(3/(4π))^(1/3)`, `Ex=ΣC·ρ^(4/3)·w`),
  `lda_exchange_potential.f90:21-31` (`Fx_μν=Σ AO_μ·C·(4/3)ρ^(1/3)·w·AO_ν`).
- xcDFT LDA-Korrelation (**fehlerhaft — ersetzen, dokumentieren**):
  `lda_correlation_energy.f90:12-19` (`Ec=-a·Σρw/(1+bρ^(-1/3)w)` — `w` im
  Nenner ist falsch), `lda_correlation_potential.f90:22-44` (hand-derive,
  keine saubere Ableitung der Energie). Ersatz: VWN5
  (Vosko-Wilk-Nusair, Can. J. Phys. 1980) parametrisierte RPA-Korrelation —
  in Doku als Fix vermerken + warum.
- Dispatch-Vorbild (portieren, um GGA/Hybrid zu erweitern): `select_rung.f90:1`
  (rung-Labels), `exchange_potential.f90:37-79` (Fx-Dispatch, rung 2/4 Stubs),
  `correlation_potential.f90:35-71`, `exchange_energy.f90:34-76`,
  `correlation_energy.f90:30-66`.
- SCF-Einbindung: `RKS.f90:154` (`F=Hc+J+Fx+Fc`), `:138` (`gradient_density`
  in xcDFT auskommentiert — hier aktivieren für WP6).
- Theorie-Quellen (Doku): Dirac 1930 (Slater-Austausch), Vosko-Wilk-Nusair 1980
  (VWN5), Slater 1951; XC-Matrix-Build auf Grid = Standard KS-DFT (z.B. Koch &
  Holthausen *A Chemist's Guide to DFT*).

## Validierungsergebnisse (nach Ausführung eintragen)

```
ORCA LSDA VWN def2-SVP: ΔE/ΔEx/ΔEc (He/Be/Ne/H2O/CH4): [ ]
xcDFT rung=1 He VDZ:    ΔE/ΔEx/ΔEc: [ ]
FD dE/dρ vs Vxc:        [ ]
```