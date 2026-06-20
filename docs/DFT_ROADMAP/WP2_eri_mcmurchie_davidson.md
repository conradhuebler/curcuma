# WP2 — 4-Zentren-ERI-Engine (McMurchie-Davidson)

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP1
- **Validierung:** ORCA 2e / unabhängiges ERI-Skript / xcDFT-`examples/*.dat`

## Ziel

Vollständige 4-Zentren-ERI für kontrahierte kartesische Gauß-Basis; 8-fache
Permutationssymmetrie; performant genug für ≤~30 Basisfunktionen.

## Deliverables

- [ ] `dft_integrals.cpp` erweitert: Hermite-Gauß-Koeffizienten-Rekursion (MD),
      elektrostatischer Hermite-Moment-Build, Kontraktion, 8-fache
      Symmetrie-Speicherung (chemists'-Notation (μν|λσ)). Didaktisch kommentiert
      (McMurchie-Davidson 1978).
- [ ] Coulomb-Build `J_μν = Σ P_λσ (μλ|νσ)` und HF-Exchange `K` aus ERI (xcDFT
      `hartree_coulomb.f90` / `fock_exchange_potential.f90` als Vorbild, zitiert).

## Erfolgskriterien

- [ ] ERI(μν|λσ) gegen ORCA-2e-Ausgabe / unabhängiges Python-Skript
      (pylibcint/hand-rolled) auf ≤1e-10 (He VDZ, H2O def2-SVP).
- [ ] 8-fache Symmetrie: max|ERI(a,b,c,d)−ERI(perm)| ≤1e-12.
- [ ] J symmetrisch; `Tr(P·J) == 2·Tr(P·K)` für geschlossenschalige Dummy-P.

## Quellen & Lektüre

- curcuma-ERI-Bauplatz: `integrals/MNDOIntegrals.hpp` (diatomarer Rahmen —
  Stil-Referenz, nicht Algorithmus; GTO-ERI ist **neu**); `integrals.h`.
- xcDFT ERI-Lesen + Symmetrisierung (Format/Notation, **nicht** die Berechnung):
  `read_integrals.f90:28` (liest `int/ERI.dat`), `:66-81` (8-fache
  Permutationssymmetrie in chemists'-Notation — direkt als Speicher-Schema
  portieren).
- J/K-Build-Vorbild: `hartree_coulomb.f90:26` (`J_μν = Σ P_λσ ERI(μ,λ,ν,σ)`),
  `fock_exchange_potential.f90:22-30` (`K`, -½-Faktor + Indexreihenfolge),
  `fock_exchange_energy.f90:23` (`Ex = ½ Tr(P·K)`).
- Theorie (externe Referenz, zu zitieren): McMurchie-Davidson, JCP 1978;
  Obara-Saika, PRA 1986; Helgaker/Jørgensen/Olsen *Molecular Electronic-Structure
  Theory* Kap. 9 — als didaktische Primärquelle in der Doku nennen.
- Validierungs-Referenz: `examples/{Ov,Kin,Nuc,ERI}.He.VDZ.dat` (vorberechnete
  He-VDZ-ERI zum Abgleich auf ≤1e-10); `examples/*.Be.VDZ.dat`, `*.Ne.VDZ.dat`.

## Validierungsergebnisse (nach Ausführung eintragen)

```
ERI vs Referenz max abs diff:   [ ]
8-fold symmetry max abs diff:   [ ]
Tr(P·J) vs 2·Tr(P·K):           [ ]
```