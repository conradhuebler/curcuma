# WP6 — PBE-GGA (+ ∇ρ, rung 2)

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP5
- **Validierung:** ORCA PBE + FD der ∇ρ-Terme

## Ziel

GGA-Functionals via ∇ρ auf Grid; PBE-Exchange + PBE-Korrelation.

## Deliverables

- [ ] `dft_xc.cpp` erweitert: `gradient_density` (xcDFT-Port, in WP4 vorhanden)
      aktiviert; PBE-Exchange (Perdew-Burke-Ernzerhof 1996), PBE-Korrelation;
      reduced-density-gradient s = |∇ρ|/(2·(6π²)^(1/3)·ρ^(4/3)); V_xC mit ρ- und
      ∇ρ-Termen. Functional `pbe` → Method-Name `-method pbe`.

## Erfolgskriterien

- [ ] ORCA `! PBE def2-SVP` auf H2O/CH4/CH3OH: E ≤1e-6 Eh; Ex/Ec ≤1e-6 (Grid medium).
- [ ] Gitterkonvergenz PBE: SG-1→SG-2 Δ<1e-4 Eh.
- [ ] `dE/d(∇ρ)`-Konsistenz (FD) ≤1e-5.

## Quellen & Lektüre

- xcDFT-GGA-Stub (zeigt, wo angesetzt wäre — war auskommentiert):
  `exchange_potential.f90:54`, `correlation_energy.f90:47` (kommentierte
  GGA-Aufrufe), `RKS.f90:138` (`gradient_density`-Aufruf auskommentiert),
  `gradient_density.f90:1` (in WP4 bereits aktiviert).
- PBE primär aus Literatur (xcDFT hat keine GGA-Implementierung):
  Perdew-Burke-Ernzerhof, PRL 77, 3865 (1996); PBE-Revision (revPBE/RPBE) nicht
  nötig. Standard-Implementierung: Koch & Holthausen Kap. 6; reduced-gradient s
  und Enhancement-Faktor F_x(s)=1+κ−κ/(1+μs²/κ).
- curcuma-Einbindung: `dft_xc.cpp` (aus WP5) Dispatch `exchange_potential`/
  `correlation_potential` um `case(2)` erweitern; ∇ρ aus `dft_grid.cpp`
  `gradient_density` (WP4).
- Theorie-Quellen (Doku): PBE 1996; Becke 1988 (B88, Basis für B3LYP in WP7);
  LYP 1988; GGAs auf Jacob's-Ladder-Rung 2.

## Validierungsergebnisse (nach Ausführung eintragen)

```
ORCA PBE def2-SVP: ΔE/ΔEx/ΔEc (H2O/CH4/CH3OH): [ ]
Gitterkonvergenz PBE SG-1→SG-2: [ ]
FD dE/d(∇ρ):      [ ]
```