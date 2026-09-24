# WP1 — GTO-1e-Integrale (Overlap, kinetisch, Kernanziehung)

- **Status:** ⚠️ AI-generated / ⚙️ machine-tested only (kein ✅ — Human-Production-Test ausstehend)
- **Abhängigkeit:** WP0
- **Validierung:** unabhängiger Python-Obara-Saika-Zeuge (Kernel, 1e-10) + ORCA kleinster S-Eigenwert (sphärisch, 1e-4) + interne Konsistenz

## Ziel

Native 1e-Integrale für kontrahierte kartesische Gauß-Basis; Basis-Laden via
`BasisSetParser`; `MakeOverlap`/`MakeH` implementiert.

## Deliverables

- [x] `qm_integrals.hpp/.cpp`: Obara-Saika 1986 Rekurrenz (J. Chem. Phys. 84,
      3963) für kontrahierte kartesische Gauß-Overlap, kinetische Energie
      (Gradienten-Identität aus verschobenen Overlap-Primitiven) und
      Kernanziehung (McMurchie-Davidson Hermite + Boys-Funktion), je mit
      Literaturzitaten. Getrennt von `GTOIntegrals.hpp` (ODR-sicher, WP1-isoliert).
- [x] `def2-SVP.dat` für H/He/Li/Be/B/C/N/O/F/Ne (gespiegelt zu ORCA `def2-SVP`,
      BSE Turbomole-Format, pre-normalisierte Primitive + Schalen-Renormierung).
- [x] `DFT::InitialiseMolecule` baut Basis + 1e-Integrale; `MakeOverlap`,
      `MakeH` liefern S, Hc = T+V (Kernabstoßung separat in
      `calculateCoreRepulsionEnergy`). Sphärisch-5d-Transform via `cartesian_d=false`.

## Erfolgskriterien

- [x] S, T, V gegen unabhängigen Python-Obara-Saika-Zeuge (selbe `def2-SVP.dat`,
      gleiche AO-Ordnung) auf ≤1e-10 (Kernel-Gate, kartesisch).
- [x] `1==Tr(P·S)` für besetzte Dummy-Dichte aus `S^{-1/2}`-Spalte (≤1e-9).
- [x] Alle drei Matrizen symmetrisch (`max|M-Mᵀ|≤1e-12`), Dimension = #Basisfunktionen.
- [x] ORCA-Kreuzcheck: kleinster **sphärischer** S-Eigenwert vs ORCA `def2-SVP`
      ≤1e-4 (AO-invariant, basis-sensitiv — bestätigt `def2-SVP.dat` == ORCA-intern
      UND korrekte kartesisch→sphärisch-d-Transform).

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

## Validierungsergebnisse

`ctest -L qm_1e` — 10/10 PASS (H2, He, LiH, BeH2, BH, CH4, NH3, H2O, HF, Ne).
Drei unabhängige Gates pro Molekül:

**(a) Kernel-Gate (kartesisch curcuma vs Python-Zeuge, tol 1e-10) — alle 10 ≤1e-14:**
```
              S          T          V          H
H2      2.8e-17    0.0e+00    2.2e-16    2.2e-16
H2O     1.1e-16    0.0e+00    7.1e-15    7.1e-15   (repräsentativ, nbf=25)
```
Der Zeuge teilt die exakte AO-Ordnung und liest dieselbe `def2-SVP.dat`, also ist
dies ein direkter, integralkernel-weiser Beweis — unabhängig von der d-Konvention.

**(b) Interne Konsistenz (kartesischer Dump) — alle 10:**
- `max|M-Mᵀ| ≤ 1e-12` für S, T, V, H (Symmetrie)
- `max|H-(T+V)| = 0.0e+00` (Hcore-Zerlegung exakt)
- `max|S_ii-1| ≤ 4.4e-16` (renormalisierte kontrahierte AOs)
- `Tr(P·S) = 1.000000000000` (Dummy-Dichte aus `S^{-1/2}`-Spalte, ≤1e-9)
- nbf stimmt mit Zege überein

**(c) ORCA-Kreuzcheck (sphärischer Dump vs ORCA `def2-SVP`, tol 1e-4) — alle 10:**
```
       curc(sph)   ORCA       diff
Ne    1.9108e-01  1.911e-01  2.3e-05
HF    4.8614e-02  4.861e-02  3.7e-06
H2O   3.5828e-02  3.583e-02  2.0e-06
BH    3.1426e-02  3.143e-02  4.5e-06
CH4   1.3844e-02  1.384e-02  3.5e-06
```
ORCA druckt den kleinsten S-Eigenwert nur mit ~4 signifikanten Stellen, daher
die 1e-4-Toleranz. **Wichtig:** ORCA rechnet def2-SVP sphärisch (5d); der
Vergleich gebruikt den sphärischen Dump (ohne `--cartesian_d`). Der kartesische
6d-Satz trägt eine zusätzliche s-artige (dxx+dyy+dzz)-Kombination, die mit dem
s-Block nahezu linear abhängig ist — sein kleinster Eigenwert ist systematisch
kleiner (z.B. Ne kartesisch 7.9e-2 vs sphärisch 1.9e-1) und wäre ein falscher
Vergleich. Die 6d→5d-Transformation ist eine 5×6-Projektion (nicht quadratisch),
also unterscheiden sich die Spektren tatsächlich.

### Analytische Stichproben (Primitive, 1e-10)
- `T_ss = 3α/2` (exakt, s-Shell auf dem Zentrum)
- `V = 2π/γ·F_0(T)` (on- und off-center, exakt) — bestätigt gegen Handrechnung
- l-Erhaltung: `<p|T|d> = 0` für Einzelatom (kinetisch ist Skalar) — war der
  Schlüssel zum Auffinden des la/lb-Kopierfehlers im Zeugen (curcuma korrekt =0)

### Was NICHT validiert wurde (konservative Selbsteinschätzung)
- Atome > Ne (`def2-SVP.dat` trägt H–Rn, kernel validiert nur H–Ne)
- f-Schalen (def2-SVP H–Ne hat keine)
- general contraction (`numContractions>1`; Code korrekt, aber ungetestet)
- spin-polarisiert/offen/Anionen/Kationen (kein SCF in WP1, `m_num_electrons` ungenutzt)
- ORCA MO-Spektrum vs curcuma Hcore-Spektrum: **ungültig in WP1** — ORCA
  `OrbitalEnergy` = SCF-konvergierte Fock-Eigenwerte (Hc+2J-K), curcuma hat bei
  WP1 nur das 1e-Hcore-Spektrum (keine 2e/SCF). Vollständiger MO-Vergleich = WP3+.
- SCF/XC/Gradient (WP3+)
- Numerische Stabilität: OS-Rekurrenz validiert nur l≤2, Exponenten ≤ ~1e3

### Referenz-Regeneration (ORCA installiert)
```bash
python3 scripts/qm_1e_reference.py test_cases/qm_1e/H2O.xyz   # -> H2O.orca_ref.json
```