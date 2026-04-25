# AP 5 — Konkreter Implementierungsplan: Validierung und Regressionstests

**Status:** Bereit zur Umsetzung (nach AP 4)
**Erstellt:** 2026-04-25
**Vorgängerdokument:** [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md), Arbeitspaket 5
**Vorbedingung:** AP 1–4 abgeschlossen
**Implementation:** AI-generiert; Validierungsergebnisse erfordern Operator-Sichtung (`✅ TESTED` darf nur der Operator setzen)

---

## Ziel

Nach AP 1–4 ist die neue Implementation **funktional** — sie kompiliert, läuft und liefert Energien plus Gradienten. AP 5 prüft, **ob die Werte richtig sind**, und macht das Ergebnis durch CTests regressions­fest. Der AP wird damit der erste, in dem der Operator `✅ TESTED` für die nativen `gfn1`/`gfn2` setzen kann (oder eben dokumentiert, dass das nicht passiert ist und warum).

Drei Validierungs­achsen:

1. **Energie** — Native gegen TBLite-Referenz
2. **Gradient** — Analytisch gegen Finite-Differenz UND analytisch gegen TBLite-Referenz
3. **End-to-end** — Optimierung, Single-Point, Trajektorien — Verhalten gegenüber bestehender TBLite-Pipeline

---

## Vorbedingungen

```bash
cd release && make -j4 curcuma 2>&1 | tail -3
ctest --output-on-failure 2>&1 | tail -5
./curcuma -sp ../test_cases/sqm_reference/molecules/H2O.xyz -method gfn2 2>&1 | grep -i energy
./curcuma -opt ../test_cases/sqm_reference/molecules/H2O.xyz -method gfn2 2>&1 | grep -iE "converg|cycles|grad"
```

Erwartet:
- build grün
- CTests passen (3 als DISABLED markierte gfn1/gfn2-Goldens aus AP 3)
- gfn2-Single-Point liefert endliche Energie
- gfn2-Optimierung konvergiert (AP 4 abgeschlossen)

USE_TBLITE muss aktiv sein, damit die Referenzwerte für AP 5 erzeugt werden können:

```bash
grep "USE_TBLITE" release/CMakeCache.txt
```

Erwartet: `USE_TBLITE:BOOL=ON`. Falls nicht: rekonfigurieren.

---

## Test-Molekül-Suite

| Molekül | Datei | Zweck | Element-Coverage |
|---|---|---|---|
| H₂ | `test_cases/sqm_reference/molecules/H2.xyz` | minimaler Edge-Case | H |
| He₂ | `He2.xyz` | closed-shell ohne kovalente Bindung | He |
| LiH | `LiH.xyz` | gemischtes Element, polar | Li, H |
| H₂O | `H2O.xyz` | kanonisches Test­molekül | O, H |
| CH₄ | `CH4.xyz` | rein kovalent, sp³ | C, H |
| NH₃ | `NH3.xyz` | lone pair | N, H |
| C₆H₆ | `C6H6.xyz` | π-System, Ringe | C, H |

Alle existieren bereits unter `test_cases/sqm_reference/molecules/`. Optional erweiterbar nach AP 5 mit Cu/Zn-haltigen Molekülen für Element-Reichweite.

---

## Aufgabenliste

### Aufgabe 5.1 — Energie-Referenztabelle aufbauen

**Skript:** Neu unter `scripts/native_xtb_energy_validation.sh`. Variante:

```bash
#!/usr/bin/env bash
set -e
MOLDIR=test_cases/sqm_reference/molecules
RESULTS=test_cases/reference_data/sqm/native_xtb_energy_validation.json
mkdir -p $(dirname $RESULTS)

declare -a MOLS=(H2 He2 LiH H2O CH4 NH3 C6H6)
echo "[" > $RESULTS
first=1
for mol in "${MOLS[@]}"; do
    for method in gfn1 gfn2 ipea1; do
        # ipea1 dient als „TBLite-aktiv"-Beweis; nicht der Vergleichsmaßstab.
        e_native=$(./release/curcuma -sp $MOLDIR/${mol}.xyz -method ${method} 2>/dev/null \
            | grep -oP 'Total Energy:\s*\K[-+0-9.eE]+')
        # TBLite-Referenz über explizit als ipea1/xtb-gfn2 verfügbar
        :
    done
done
```

Die exakte Form ist eine Skript-Frage, das **Wesentliche** ist die Tabelle als CSV/JSON für AP-5-Reporting:

| Molekül | Methode | E_native (Eh) | E_TBLite (Eh) | ΔE (Eh) | ΔE (kJ/mol) |
|---|---|---|---|---|---|
| H₂O | GFN2 | … | … | … | … |
| ... | ... | ... | ... | ... | ... |

**Datenpfad für TBLite-Referenz:** zwei Optionen

(a) **Direkt aus laufendem Curcuma**, indem temporär `gfn2-tblite` als zusätzliche, alternative Methode in MethodFactory verfügbar gemacht wird. Aufwand niedrig — eine neue Methode `tblite-gfn2`/`tblite-gfn1`, die direkt `TBLiteMethod` instanziiert. **Empfehlung.**

(b) **Externes tblite-CLI**, gestartet aus `external/tblite/`. Reproduzierbar, aber mehr Setup.

**Empfehlung:** Option (a). Anpassung in `method_factory.cpp` (in AP 5, nicht AP 3 vorgezogen):

```cpp
if (method == "tblite-gfn2") return std::make_unique<TBLiteMethod>("gfn2", config);
if (method == "tblite-gfn1") return std::make_unique<TBLiteMethod>("gfn1", config);
```

Damit ist der Vergleich `gfn2` (native) vs. `tblite-gfn2` (TBLite) im selben Prozessbaum möglich.

**Akzeptanzschwellen:**

| Bereich | ΔE pro Molekül |
|---|---|
| 1–4 Atome (H₂, He₂, LiH, H₂O) | < 1e-4 Eh |
| 5–8 Atome (CH₄, NH₃) | < 5e-4 Eh |
| > 8 Atome (C₆H₆) | < 1e-3 Eh |

Falls einer dieser Werte überschritten wird: **Eintrag in „Schwierigkeiten / Blocker"**, plus klare Dokumentation der wahrscheinlichen Ursache (z. B. „GFN2 Multipole-Implementierung simplified, AP 6").

---

### Aufgabe 5.2 — Analytischer vs. numerischer Gradient

**Skript:** Erweitert `scripts/native_xtb_energy_validation.sh` um Gradient-Sektion oder als separates `scripts/native_xtb_gradient_validation.sh`.

```bash
for mol in "${MOLS[@]}"; do
    for method in gfn1 gfn2; do
        # CLI-Flag aus AP 4 (-grad-numerical-check) — vergleicht analytisch mit FD intern,
        # gibt max_abs_error pro Molekül zurück.
        ./release/curcuma -sp $MOLDIR/${mol}.xyz -method ${method} -grad-numerical-check 1 \
            >> grad_validation.log
    done
done
```

Falls dieser CLI-Flag nicht in AP 4 angelegt wurde: in AP 5 nachholen (in `src/main.cpp` als Single-Point-Variante mit Gradient-Vergleich), oder als Standalone-Test im sqm_reference-Setup belassen.

**Akzeptanzschwellen** (passend zu AP 4):

| Methode | Bereich | max ‖∇E_an − ∇E_num‖ |
|---|---|---|
| GFN1 (kein Multipole) | alle Moleküle | < 1e-4 Eh/Bohr |
| GFN2 (mit Multipole) | alle Moleküle | < 1e-3 Eh/Bohr |

---

### Aufgabe 5.3 — Analytischer Gradient vs. TBLite-Referenz

Schwieriger als 5.2: TBLite muss seinen Gradient ausgeben. `tblite-gfn2`-Methode aus 5.1 nutzen, plus:

```cpp
// Im TBLiteMethod-Wrapper getGradient() bereits implementiert.
auto m_native = MethodFactory::create("gfn2", {});
auto m_tblite = MethodFactory::create("tblite-gfn2", {});
m_native->setMolecule(mol);
m_tblite->setMolecule(mol);
m_native->calculateEnergy(true);
m_tblite->calculateEnergy(true);
double max_diff = (m_native->getGradient() - m_tblite->getGradient()).cwiseAbs().maxCoeff();
```

**Akzeptanzschwellen** dieselben wie 5.2 — wir messen praktisch dieselbe Größe (TBLite ist die Referenz, FD ist Selbstvalidierung).

---

### Aufgabe 5.4 — Optimierungs-End-to-End

```bash
for mol in H2 H2O CH4 NH3 C6H6; do
    for method in gfn1 gfn2; do
        # Beide Pfade: native vs. TBLite
        cp $MOLDIR/${mol}.xyz tmp_${mol}_native.xyz
        cp $MOLDIR/${mol}.xyz tmp_${mol}_tblite.xyz
        ./release/curcuma -opt tmp_${mol}_native.xyz -method ${method} > log_native_${mol}_${method}.txt
        ./release/curcuma -opt tmp_${mol}_tblite.xyz -method tblite-${method} > log_tblite_${mol}_${method}.txt

        # Final structure RMSD
        ./release/curcuma -rmsd tmp_${mol}_native.opt.xyz tmp_${mol}_tblite.opt.xyz
    done
done
```

**Akzeptanz:** finale RMSD zwischen native und TBLite < 1e-3 Å pro Atom (oder dokumentiert, warum nicht). Konvergenz beider Pfade in vergleichbarer Iterations­zahl (Faktor ≤ 2x).

---

### Aufgabe 5.5 — CTests neu eichen und reaktivieren

**5.5a — gfn1/gfn2-Goldens** aus AP 3 (DISABLED-Marker aufheben):

```bash
# In test_cases/cli/curcumaopt/CMakeLists.txt und test_cases/cli/sqm/CMakeLists.txt:
# - Die set_tests_properties(... DISABLED TRUE) aus AP 3 entfernen
# - Golden-Dateien (Energie-Werte in test_cases/cli/.../golden/) gegen die neuen
#   native-Energien aktualisieren
```

Die Golden-Werte stammen aus 5.1 (Native-Spalten). Toleranzen in den Test-Skripten anpassen, falls nötig (Tolerance ≤ 1e-4 Eh).

**5.5b — Neue End-to-end-Tests**:

```cmake
# test_cases/cli/sqm/CMakeLists.txt
add_cli_test(sqm 09_native_gfn2_h2o_singlepoint)
add_cli_test(sqm 10_native_gfn2_h2o_optimization)
add_cli_test(sqm 11_native_gfn1_h2o_singlepoint)
add_cli_test(sqm 12_native_gfn2_grad_consistency_check)  # FD vs analytical
```

Skripte unter `test_cases/cli/sqm/09_*.sh` … erstellen, jeweils mit Golden-Energie/-Gradient.

**5.5c — Cross-Provider-Konsistenztest**:

```cmake
add_cli_test(sqm 13_native_vs_tblite_h2o_energy)   # |E_native - E_tblite| < 1e-4 Eh
add_cli_test(sqm 14_native_vs_tblite_h2o_grad)     # |∇E_native - ∇E_tblite|_max < 1e-4 Eh/Bohr
```

Diese Tests laufen nur, wenn USE_TBLITE aktiv ist — entsprechend gating mit `if(USE_TBLITE) add_cli_test(...) endif()`.

---

### Aufgabe 5.6 — Performance-Vergleich

| Molekül | TBLite (s) | Native (s) | Verhältnis | Vertretbar? |
|---|---|---|---|---|
| H₂O | … | … | … | … |
| C₆H₆ | … | … | … | … |
| Aspirin (21 Atome) | … | … | … | … |

```bash
time ./release/curcuma -sp benzene.xyz -method tblite-gfn2 -loop 100
time ./release/curcuma -sp benzene.xyz -method gfn2 -loop 100
```

`-loop 100` ist hypothetisch — alternativ Skript-Schleife.

**Akzeptanzschwelle:** Native bis zu 5x langsamer als TBLite ist OK (TBLite ist Fortran + LAPACK + Optimierungen). Faktoren > 10x → in „Schwierigkeiten" eintragen, Profile Pflicht.

---

### Aufgabe 5.7 — Status-Doku schreiben

**`docs/NATIVE_XTB_STATUS.md`** komplett überarbeiten (existiert bereits, von AP 1–4 inkrementell aktualisiert):

```markdown
# Native xTB Implementation Status

**Last Updated**: <Datum>
**Implementation**: AI-generiert + machine-tested gegen TBLite (AP 5)
**Human production testing**: <ausstehend / Operator-Status>

## Test Results (AP 5)

### Energy Validation (nat-vs-TBLite)

| Molekül | GFN2 ΔE (μEh) | GFN1 ΔE (μEh) |
|---|---|---|
| H₂  | … | … |
| He₂ | … | … |
...

### Gradient Validation

| Molekül | GFN2 max ‖∇an−∇TBLite‖ (Eh/Bohr) | GFN1 |
|---|---|---|
...

### Optimization Convergence

| Molekül | RMSD nat vs TBLite (Å) | Iters nat | Iters TBLite |
|---|---|---|---|

### Performance

(Tabelle aus 5.6)

## Known Issues / Limitations

- Multipole-Gradient: simplified radial form (AP 4e); accuracy 1e-3 Eh/Bohr statt 1e-4
- d-Funktionen: nicht in Multipol-Integralen abgedeckt (Zn, Cu, Fe-haltige Moleküle untergenau)
- DIIS: nicht implementiert; Konvergenz langsam für stark polarisierte Systeme
```

**`CLAUDE.md`** Status-Tabelle aktualisieren:

```
| Native GFN2-xTB | X/Y vs TBLite | <Status nach 5.5c> |
| Native GFN1-xTB | X/Y vs TBLite | <Status nach 5.5c> |
```

`CLAUDE.md` „Known Issues" um die in 5.7 dokumentierten Limitationen ergänzen.

**`AIChangelog.md`**: eine Zeile.

```
2026-XX-XX: AP 5 abgeschlossen — native gfn1/gfn2 validiert gegen TBLite (X/Y Moleküle, ΔE < N μEh, ∇ < M Eh/Bohr). Status-Tabelle aktualisiert in NATIVE_XTB_STATUS.md.
```

---

### Aufgabe 5.8 — Operator-Übergabe

Nach Abschluss von 5.1–5.7 ist die Implementation **machine-tested**, aber noch **nicht human-tested**. Die endgültige `✅ TESTED`-Markierung darf nur der Operator setzen. AP 5 endet daher mit einem expliziten Übergabe-Punkt:

```markdown
## Status: Bereit für Operator-Review

[Datum]: AP 5 abgeschlossen — siehe `docs/NATIVE_XTB_STATUS.md` für vollständige
Validierungs­ergebnisse. Code ist machine-tested gegen TBLite-Referenz auf den
in der Suite enthaltenen Molekülen. **Operator wird gebeten, vor Setzen von
✅ TESTED:**
1. Validierungs­tabellen in NATIVE_XTB_STATUS.md sichten
2. Repräsentatives Produktions­molekül durchrechnen
3. Falls OK: ⚠️-Markierung in CLAUDE.md durch ✅ ersetzen
```

---

## Akzeptanzkriterien

### Numerisch

- [ ] **Energie:** Alle 7 Moleküle × 2 Methoden (gfn1, gfn2) — ΔE_native − E_TBLite| in Toleranz (siehe 5.1)
- [ ] **Gradient (FD):** Alle Moleküle — max-Abweichung in Toleranz (siehe 5.2)
- [ ] **Gradient (TBLite):** Alle Moleküle — max-Abweichung in Toleranz (siehe 5.3)
- [ ] **Optimierung:** Native gfn1/gfn2 konvergiert für alle Moleküle, finale Struktur ähnlich TBLite (5.4)

### CTest

- [ ] Alle bisherigen 26 CTests passieren weiter
- [ ] AP-3-DISABLED-Tests reaktiviert mit neuen Goldens (5.5a)
- [ ] Mindestens 4 neue native-spezifische CTests (5.5b)
- [ ] Mindestens 2 Cross-Provider-Tests (5.5c, gating mit USE_TBLITE)

### Doku

- [ ] `docs/NATIVE_XTB_STATUS.md` enthält vollständige Validierungs-Tabelle
- [ ] `CLAUDE.md` Status-Tabelle aktuell
- [ ] `AIChangelog.md` enthält Eintrag
- [ ] Bekannte Limitationen explizit dokumentiert

### Nicht-Ziele

- ❌ Operator-Tests an Produktions-Molekülen (=Aufgabe des Operators)
- ❌ Performance-Optimierung (separater AP)
- ❌ Erweiterung der Element-Coverage über die Test-Suite hinaus (AP 6)

---

## Risiken und Stolperfallen

| Risiko | Eintrittswahrscheinlichkeit | Gegenmaßnahme |
|---|---|---|
| Energieabweichungen > 1e-3 Eh ⇒ AP 4 hat versteckte Fehler | mittel | pro Molekül und Anteil isolieren — `getEnergyDecomposition()` aus AP 1 zeigt, welcher Term differiert |
| Gradient-Fehler > 1e-3 Eh/Bohr in einer Sub-AP | mittel | FD-Sweep über h ∈ {1e-3, 5e-4, 1e-4} — falls inkonsistent, ist das Problem im analytischen Gradient |
| TBLite-Wrapper-Pfad selbst war nie sauber validiert (Curcuma-Bridge) | niedrig | bei Diskrepanz: `./release/curcuma -sp X -method tblite-gfn2` mit external `tblite` CLI cross-checken |
| C₆H₆-Optimierung divergiert wegen fehlendem DIIS | mittel | Linear damping mit 0.4 zu konservativ, `damping=0.7` testen; ggf. als Limitation dokumentieren |
| Performance > 10x TBLite | mittel | Profile mit `perf` oder `gprof`, häufigste Hotspots: GammaMatrix-Aufbau (O(nsh²)) und Multipol-Setup |
| Operator akzeptiert nicht, weil Goldens-Differenzen vorhanden | hoch | klare Dokumentation der **erwarteten** Differenz (ja, das ist eine Re-Implementation, ja, sie ist nicht bit-genau) — TBLite selbst ist Referenz, aber „1:1" ist unrealistisch |

---

## Test-Plan

### Vollständiger Validierungs-Lauf

```bash
cd release && make -j4 curcuma
bash ../scripts/native_xtb_energy_validation.sh > validation_energies.txt
bash ../scripts/native_xtb_gradient_validation.sh > validation_gradients.txt
ctest --output-on-failure 2>&1 | tee ctest_ap5.log
```

### Reporting

Nach jedem Lauf wird `docs/NATIVE_XTB_STATUS.md` mit den Tabellen­werten aktualisiert. Reporting-Skript: optional ein `scripts/native_xtb_report.sh`, das aus den JSON-Logs die Markdown-Tabellen erzeugt.

---

## Fortschritt

| Datum | Aufgabe | Status | Notizen |
|-------|---------|--------|---------|
| 2026-04-25 | 5.1 Energie | Plan | — |
| | 5.2 Grad-FD | | |
| | 5.3 Grad-TBLite | | |
| | 5.4 Opt-E2E | | |
| | 5.5 CTests | | |
| | 5.6 Perf | | |
| | 5.7 Doku | | |
| | 5.8 Übergabe | | |

## Schwierigkeiten / Blocker

- *Noch keine dokumentiert*

---

## Referenzen

- **Vorgänger:** [`NATIVE_XTB_AP4_PLAN.md`](NATIVE_XTB_AP4_PLAN.md)
- **Roadmap-Übersicht:** [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md)
- **Status-Doku (Ziel der Aktualisierung):** [`NATIVE_XTB_STATUS.md`](NATIVE_XTB_STATUS.md)
- **Test-Infrastruktur:**
  - `test_cases/sqm_reference/` — Standalone-Tests (kernel-Niveau)
  - `test_cases/cli/sqm/` — CLI-CTests (end-to-end)
  - `test_cases/cli/curcumaopt/` — Optimierungs-CTests
  - `test_cases/reference_data/sqm/` — Referenzwerte
- **TBLite-Referenz:** über neue `tblite-gfn1`/`tblite-gfn2`-Methoden in MethodFactory (5.1)

---

## Änderungshistorie

| Datum | Autor | Änderung |
|-------|-------|----------|
| 2026-04-25 | Claude | Erstdokument |
