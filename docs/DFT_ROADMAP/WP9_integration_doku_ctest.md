# WP9 — Integration, CLI, Parameter, Doku, CTest (Finalisierung)

- **Status:** ⚠️ AI-generated / offen — ⚙️ machine-tested (nach Erfüllung)
- **Abhängigkeit:** WP1–WP8
- **Validierung:** gesamte CTest-Suite + Round-Trip

## Ziel

Vollständige CLI/Factory/Parameter-Integration, Doku final, CTest-Suite,
Method-Funktional-Aliase, Aufräumen.

## Deliverables

- [ ] `method_factory.cpp`: die pro-Funktional-Method-Namen `hf`/`lda`/`pbe`/`b3lyp`
      (bereits in WP0 registriert) in `getMethodInfo`/`printAvailableMethods`
      aufnehmen; keine Sammel-`dft`.
- [ ] `main.cpp`: Warm-Start-Bedingung um `hf`/`lda`/`pbe`/`b3lyp` erweitern (`:1724`),
      Hilfe/Beispiele (`:1636`/`:1678`), `-methods`-Gruppierung (`:2309`).
- [ ] Parameter final: `make GenerateParams` sauber; `-export_run`/`-import_config`
      Round-Trip für `-method hf`/`lda`/`pbe`/`b3lyp` geprüft.
- [ ] `docs/NATIVE_DFT_IMPLEMENTATION.md` vollständig (Theorie, Gleichungen,
      Literatur, Quellen-Map, Test-Sektion: getestet/nicht getestet/nicht implementiert).
- [ ] CLAUDE.md-Status (⚠️/⚙️, kein ✅), Haupt-CLAUDE.md Kap. 1, README, AIChangelog.
- [ ] CTest `cli_dft_*`: sp LDA/HF/PBE/B3LYP auf H2O/CH4 gegen gespeicherte
      Referenz-Energien; `cli_dft_opt` Smoke.

## Erfolgskriterien

- [ ] Alle WPs: `make -j4` ohne Warning; `ctest` komplett grün (inkl. `cli_dft_*`).
- [ ] `curcuma -sp H2O.xyz -method b3lyp -dft.basis def2-svp -dft.grid medium`
      liefert Referenz-tolerante Energie; `-export_run run.json` + Replay identisch.
      (Functional über Method-Namen `b3lyp`, nicht über `-dft.functional`.)
- [ ] Doku enthält vollständige Quellenangabe xcDFT und Selbstbewertungs-Caveats.

## Quellen & Lektüre

- Factory/Aliase: `method_factory.cpp:352/509/555/633` (wie WP0); CLI-Gruppierung
  `main.cpp:2309`, Methode-Lesen `:1181/1606/1705`, Warm-Start `:1724`, Hilfe
  `:1636`, Beispiele `:1678`.
- Parameter-Round-Trip: `docs/CLI_ROUND_TRIP.md` (`-export_run`/`-import_config`,
  `_command`/`_input`); Flat-Flag-Routing `main.cpp` (CLAUDE.md JSON-Controller-Sektion).
- Doku-Vorbild: `docs/NATIVE_QM_METHODS_IMPLEMENTATION.md`,
  `docs/SQM_VALIDATION.md` (Toleranz-Tabellen),
  `docs/SQM_WP3_component_validation.md` (Komponenten-Validierung).
- CTest-Vorbild: `test_cases/test_orca_interface.cpp`; CLI-Tests `cli_rmsd_*`/
  `cli_curcumaopt_*` Muster in `test_cases/CLAUDE.md`.
- Status-Label-Regeln: CLAUDE.md „AI-Generated Content and Validation Policy" +
  „Conservative Self-Assessment Rules" (kein selbst vergebenes ✅ TESTED/APPROVED).

## Validierungsergebnisse (nach Ausführung eintragen)

```
make -j4 / make GenerateParams: [ ]
ctest (cli_dft_*):              [ ]
-export_run/-import_config:      [ ]
Doku Quellen-Sektion:           [ ]
```