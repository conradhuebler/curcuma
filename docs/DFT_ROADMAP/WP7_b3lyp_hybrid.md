# WP7 — B3LYP-Hybrid (+ exakter Austausch, rung 4)

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP6 + WP2 (K_HF)
- **Validierung:** ORCA B3LYP

## Ziel

Hybrid-Functional B3LYP: 0.20·HF + 0.80·B88 + 0.72·LYP + 0.19·VWN5_1
(Stephens/Becke/Frisch 1994); K-Build-Kopplung im SCF.

## Deliverables

- [ ] `dft_xc.cpp` erweitert: B88-Exchange (Becke 1988), LYP-Korrelation
      (Lee-Yang-Parr 1988), B3LYP-Mischungsform. Functional `b3lyp` → Method-Name
      `-method b3lyp`.
- [ ] SCF-Hybrid-Pfad: Fock = Hc + (1-a)·V_xGGA + a·K_HF + V_c; K_HF aus WP2-ERI
      (kein CP-KSCF nötig für geschlossene Schale, direkter K-Build wie HF).

## Erfolgskriterien

- [ ] ORCA `! B3LYP def2-SVP` auf H2O/CH4/CH3OH/C6H6: E ≤1e-6 Eh, Komponenten ≤1e-6.
- [ ] HF-Anteil (0.20·Ex_HF) gegen WP3-Ex reproduzierbar ≤1e-6.
- [ ] Konvergenz <80 Iter (Hybrid langsamer); DIIS stabil.

## Quellen & Lektüre

- xcDFT-Hybrid-Stub (war unvollständig): `exchange_energy.f90:58-59`
  (`cX=0.20, aX=0.72`), `exchange_potential.f90:61-62`,
  `correlation_energy.f90` (`aC=0.81`); alle Sub-Calls auskommentiert → 0. Nur als
  Hinweis auf Mischungs-Koeffizienten-Ort, nicht als Implementierung.
- B3LYP-Form aus Literatur: Stephens, Devlin, Chabalowski, Frisch, JPC 98,
  11623 (1994); B88 Becke PRA 38, 3098 (1988); LYP Lee-Yang-Parr PhysRevB 37, 785
  (1988); VWN5_1-Form aus WP5. Standard-Mischung
  `E_xc[B3LYP] = 0.20·Ex[HF] + 0.80·Ex[B88] + 0.72·Ec[LYP] + 0.19·Ec[VWN5]_1 + 0.81·Ec[LYP]`
  (Lypsche α=0.81-Formulierung) — in Doku exakt festhalten, Toleranz gegen
  ORCA-Variante prüfen.
- K_HF: Wiederverwendung aus WP2 (`fock_exchange_potential.f90:22-30` Schema,
  in curcuma `qm_integrals.cpp`); HF-SCF-Pfad aus WP3 zeigt den K-Build-Loop.
- Hybrid-SCF: `xtb_scf.cpp` als Gerüst (geschlossene Schale = direkter K, kein
  CP-KSCF); `xtb_response.cpp` nur Referenz, falls später offene Schale.
- Theorie-Quellen (Doku): Stephens/Becke/Frisch 1994; Becke 1993 (ursprüngliche
  3-Parameter-Hybride); Hybrid = Rung 4.

## Validierungsergebnisse (nach Ausführung eintragen)

```
ORCA B3LYP def2-SVP: ΔE/ΔKomponenten (H2O/CH4/CH3OH/C6H6): [ ]
0.20·Ex_HF vs WP3-Ex: [ ]
SCF-Iterationen:      [ ]
```