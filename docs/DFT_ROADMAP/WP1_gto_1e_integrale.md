# WP1 — GTO-1e-Integrale (Overlap, kinetisch, Kernanziehung)

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP0
- **Validierung:** ORCA 1e-Komponenten + Symmetrie-Checks

## Ziel

Native 1e-Integrale für kontrahierte kartesische Gauß-Basis; Basis-Laden via
`BasisSetParser`; `MakeOverlap`/`MakeH` implementiert.

## Deliverables

- [ ] `dft_integrals.hpp/.cpp`: kontrahierte Gauß-Overlap (`GTOIntegrals.hpp`
      erweitern), kinetische Energie (Kartesisch-Gauß-Rekursion), Kernanziehung
      (Obara-Saika für V), je mit Literaturzitaten (Obara-Saika 1986).
- [ ] `def2-SVP.dat` für H/He/Li/Be/B/C/N/O/F/Ne (gespiegelt zu ORCA `def2-SVP`).
- [ ] `DFT::InitialiseMolecule` baut Basis + 1e-Integrale; `MakeOverlap`,
      `MakeH` liefern S, Hc = T+V; Kernabstoßung.

## Erfolgskriterien

- [ ] S, T, V gegen ORCA-1e-Ausgabe (`! HF def2-SVP` + Hcore/Komponenten) bzw.
      Python-Kontrollskript auf ≤1e-10.
- [ ] `1==Tr(P·S)` für besetzte Dummy-Dichte (Symmetrie/Orthogonalität).
- [ ] Alle drei Matrizen symmetrisch, Dimension = #Basisfunktionen.

## Quellen & Lektüre

- curcuma-Basis: `basissetparser.hpp:15` (`#include GTOIntegrals`), `:27`
  `BasisSetParser`, `:47-57` `ElementBasis`/`BasisSetMap`; gelieferte Basis
  `qm_methods/def2-SV(P).dat` als `.dat`-Format-Referenz.
- 1e-Primitive: `GTOIntegrals.hpp` (Overlap-Primitive; kinetisch/V nach Obara-Saika
  als Erweiterung prüfen — ggf. `STOIntegrals.hpp`/`integrals/MNDOIntegrals.hpp`
  für Rekursions-Stil, nicht für Physik).
- QMDriver-Hooks: `qm_driver.h:59-60` (`MakeOverlap(Basisset&)`, `MakeH(...)`),
  Speicher `m_H,m_S,m_mo,m_energies` `:63-69`.
- xcDFT-1e-Lesen (nur Format/Erwartung, Integrale werden **selbst** berechnet, nicht
  gelesen): `read_integrals.f90:1` (liest `int/{Ov,Kin,Nuc}.dat`, baut
  `Hc = T+V` `:59`); `read_basis.f90:1` (Basis-Format, Shell-Buchstaben `:57-77`,
  Kontraktion `:40-89`); `read_geometry.f90:31-38` (Kernabstoßung);
  `NormCoeff.f90:27-29` (Gauß-Normalisierung, zu portieren).

## Validierungsergebnisse (nach Ausführung eintragen)

```
S vs ORCA max abs diff:    [ ]
T vs ORCA max abs diff:    [ ]
V vs ORCA max abs diff:    [ ]
Tr(P·S):                   [ ]
Symmetrie-Check:           [ ]
```