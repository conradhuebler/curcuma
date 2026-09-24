# WP2 — 4-Zentren-ERI-Engine (McMurchie-Davidson)

- **Status:** ⚙️ machine-tested (`ctest -L qm_2e` 10/10, Juni 2026) — AI-generated, noch nicht human production tested
- **Abhängigkeit:** WP1
- **Validierung:** unabhängiges Pure-stdlib-Python-MD-Witness (`scripts/qm_2e_python_ints.py`) auf 10 def2-SVP-Molekülen; (ss|ss)-Stichprobe gegen geschlossene Form; optionaler xcDFT-He-VDZ-Abgleich (vorbereitet, Referenz nicht im Repo).

## Ziel

Vollständige 4-Zentren-ERI für kontrahierte kartesische Gauß-Basis; 8-fache
Permutationssymmetrie; performant genug für ≤~30 Basisfunktionen.

## Deliverables

- [x] `qm_integrals.{hpp,cpp}` erweitert: `ERITensor` (flaches n⁴, chemists'
      (μν|λσ)), `buildERI` (MD: `hermiteCoeffs` + neuer ERI-R-Block mit
      `R^n_{000}=(-2ρ)^n F_n(T)` und (P−Q)-Verschiebung, Hermite-Hermite-Kontraktion
      mit (−1)^(τ+υ+ω), kanonische Quartet-Schleife mit 8-facher Füllung via
      `set8`). Didaktisch kommentiert (McMurchie-Davidson 1978, Helgaker Kap. 9.9).
- [x] Coulomb `J_μν = Σ_λσ P_λσ (μν|λσ)` und Exchange `K_μν = Σ_λσ P_λσ (μλ|νσ)`
      (`buildCoulomb`/`buildExchange`) aus ERI; 4-Index-Sphärisch-Transform
      `applySphericalTransformERI` für WP3. xcDFT-Vorbild zitiert (siehe unten);
      curcuma speichert chemists' direkt (xcDFT speichert physicists' ⟨ij|kl⟩).
- [x] Lazy `DFT::cartesianERI()` (gebaut on-demand, NICHT im Scaffold-`-sp`-Pfad);
      `ctest -L qm_2e` 10 Moleküle via `dump_qm_2e` + `diff_qm_2e.py`.

## Erfolgskriterien

- [x] ERI(μν|λσ) vs unabhängiges Python-MD-Witness (Pure-stdlib, gleiche AO-Reihenfolge)
      auf 10 def2-SVP-Molekülen: max|curc−witness| = **2.1e-14** (worst, CH4/NH3/H2O;
      He/LiH exakt 0), Ziel ≤1e-10. (ss|ss)-Stichprobe gegen geschlossene Form
      `2π^{5/2}/(pq√(p+q))·K_AB·K_CD·F_0(T)` = 2.2e-16 (inkl. T=0-Gleichzentrum).
- [x] 8-fache Symmetrie: max|ERI(a,b,c,d)−ERI(perm)| = **0.0** (exakt, Ziel ≤1e-12).
- [x] J und K symmetrisch (≤1e-12); `Tr(P·J) == Tr(P·K)` für Dummy-P (Rang-1,
      P=2cc^T) = **1.8e-13** (worst Ne; Identität gilt für jedes 4-Tensor per
      Dummy-Index-Umordnung b↔c). Die Roadmap-Form `Tr(P·J)==2·Tr(P·K)` ist die
      Energie-Relation E_Coulomb=½Tr(PJ)=2·E_exchange=2·¼Tr(PK) für ein Orbital
      (Tr(PJ)==Tr(PK) auf Matrix-Spur-Niveau, da J_cc==K_cc für ein Orbital).

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
ERI vs Referenz max abs diff:   2.1e-14  (Python MD witness, 10 def2-SVP mols; worst CH4/NH3/H2O)
8-fold symmetry max abs diff:   0.0      (exakt)
Tr(P·J) vs Tr(P·K):             1.8e-13  (Rang-1 Dummy-P, worst Ne; == Energie-Relation E_C=2·E_x)
(ss|ss) closed form:            2.2e-16  (inkl. T=0-Gleichzentrum)
ctest -L qm_2e:                10/10    (H2 He LiH BeH2 BH CH4 NH3 H2O HF Ne)
ctest -L qm_1e (Regression):   10/10
```

Anmerkung: der optionale xcDFT-He-VDZ-Abgleich (`--xcDFT-ref ERI.He.VDZ.dat`)
ist im `diff_qm_2e.py`-Gate (c) vorbereitet, wird aber automatisch übersprungen,
wenn die Referenzdatei nicht im Repo liegt (Checkout ohne Referenz läuft (a)+(b)).
Human-Produktionstest steht aus.