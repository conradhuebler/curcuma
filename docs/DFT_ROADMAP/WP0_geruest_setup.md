# WP0 — Gerüst, Quellenangabe, Setup

- **Status:** ⚙️ machine-tested (Gerüst läuft; ✅ TESTED/APPROVED nur durch Mensch)
- **Abhängigkeit:** keine
- **Validierung:** manueller Smoke-Test + ORCA-Aufrufbarkeit

## Ziel

Leeres Engine/Wrapper-Gerüst im MethodFactory registriert; Doku-Gerüst mit
Quellenangabe; ORCA-Referenzumgebung nutzbar.

## Deliverables

- [ ] `dft.h`/`dft.cpp` — Engine `DFT : public QMDriver`, `DFTFunctional`-Enum
      (HF/LDA/PBE/B3LYP), `Calculation` liefert erst Kernabstoßung + 0.
- [ ] `dft_method.h/.cpp` — Wrapper, delegiert, `hasGradient()=false` vorerst;
      ctor nimmt `DFTFunctional`.
- [ ] `BEGIN_PARAMETER_DEFINITION(dft)` mit basis/grid/threads/scf_* (functional
      ist **kein** Parameter — wird durch Method-Namen festgelegt).
- [ ] `method_factory.cpp`: **kein** `dft`-Sammel-Method; stattdessen pro Funktional
      ein eigener Method-Name, jeweils `DFTMethod(DFTFunctional::..., config)`:
      `hf`→HF, `lda`→LDA, `pbe`→PBE, `b3lyp`→B3LYP. Eintrag in `getAvailableMethods`;
      `CMakeLists.txt`: neue `.cpp`; `make GenerateParams` sauber.
- [ ] `docs/NATIVE_QM_IMPLEMENTATION.md`-Gerüst mit Herkunfts-Sektion.
- [ ] CLAUDE.md-Status-Zeile (qm_methods + Haupt), README-Link.
- [ ] ORCA-Check (`ORCA_PATH=/opt/orca_6_1/orca --version`); Test-Geometrien
      festgelegt (He/Be/Ne-XYZ erzeugen, H2O/CH4 vorhanden).

## Erfolgskriterien

- [ ] `curcuma -sp <He.xyz> -method hf` (bzw. `-method lda`) läuft (Energie =
      Kernabstoßung, Meldung „native DFT — nur Gerüst").
- [ ] `make -j4` ohne Warning; `make GenerateParams` ohne Validierungswarnung.
- [ ] `curcuma -methods` listet `hf`/`lda`/`pbe`/`b3lyp` unter Quantum Methods.
- [ ] Doku-Gerüst existiert mit Quellen-Sektion (xcDFT-TCCM-2019).

## Quellen & Lektüre

- Template Wrapper: `nddo_method.h:21`, `nddo_method.cpp:42-88` (delegieren +
  getDefaultConfig); Engine-Gerüst: `nddo.h:55` (`NDDO : public QMDriver`),
  `qm_driver.h:39` (`QMDriver`, `MakeOverlap`/`MakeH` `:59-60`),
  `abstract_interface.h:24` (`QMInterface::InitialiseMolecule` `:30-40`).
- Factory-Registrierung: `method_factory.cpp:352` (`create()`), `:366` native
  Block, `:509` `getAvailableMethods`, `:555` `getMethodInfo`,
  `:633` `printAvailableMethods`; ADR `:324-351`.
- Parameter: `parameter_macros.h:38-43`; PARAM-Vorbild `xtbinterface.h:33-87`;
  Flow `qm_methods/QM_ARCHITECTURE.md:306-363` (Skeleton `:326-343`,
  Parameter-Kette `:382-403`).
- Doku-Vorbild: `docs/NATIVE_QM_METHODS_IMPLEMENTATION.md` (Attribution-Template
  `:27-49`, Theorie-Zugänglichkeitsregeln `:57-80`).
- ORCA-Anbindung (kein Code-Change): `orcainterface.cpp:201` `findOrcaExecutable`
  (Priorität `orca_executable`→`ORCA_PATH`→`ORCA_ROOT`→`which orca`);
  `orca_method.cpp:230` Methodenliste.
- CMake: `CMakeLists.txt:428-464` (`curcuma_core_SRC`); CLI-Gruppierung
  `main.cpp:2309`, Hilfe `:1636`, Warm-Start `:1724`.
- xcDFT (nur zum Verstehen des Gesamtflusses): `xcDFT.f90:1` (Driver),
  `RKS.f90:1` (SCF-Struktur), `read_options.f90:14` (rung/SGn),
  `print_RKS.f90:23-57` (Komponenten-Ausdruck, zu replizieren).

## Validierungsergebnisse (nach Ausführung eintragen)

```
make -j4:           [x]  (rc=0, no DFT-related warnings; pre-existing -Winline Mol::~Mol in test_xtb_cpscf unrelated)
make GenerateParams:[x]  (clean; no validation warning mentioning the dft module)
curcuma -methods (hf/lda/pbe/b3lyp): [x]  (listed under Quantum Methods)
ORCA --version:     [x]  (ORCA 6.1.0 at /opt/orca_6_1/orca)
```

Smoke-test (2026-06-20, release/ build):
- `curcuma -sp test_cases/he.xyz -method {hf,lda,pbe,b3lyp}` -> prints
  "native DFT -- nur Geruest", E = 0.00000000 Eh (single atom, E_nn = 0).
- `curcuma -sp test_cases/water.xyz -method pbe` -> E_nn = 9.64357925 Eh
  (physically sensible nuclear repulsion for H2O).
- Full `ctest`: pre-existing failures (d4_diag_*, ecomp_* / gfn{1,2}_align,
  sqm_scf_*_gfn2) are NOT regressions — reproduced identically on a clean
  baseline (stash + rebuild). None touch DFT; GFN2/D4/SCF paths are outside
  the additive WP0 scaffold.