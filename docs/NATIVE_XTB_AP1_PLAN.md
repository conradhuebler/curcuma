# AP 1 — Konkreter Implementierungsplan: XTB-Klasse erweitern und konsolidieren

**Status:** Bereit zur Umsetzung
**Erstellt:** 2026-04-25
**Vorgängerdokument:** [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md), Arbeitspaket 1
**Implementation:** AI-generiert, nicht-getestet — `🤖 AI-generated`-Status nach Abschluss

---

## Ziel

Die Klasse `curcuma::xtb::XTB` (`xtb_native.h/cpp`) erhält die fehlende Funktionalität, die der `GFN2Method`/`GFN1Method`-Wrapper in AP 2 erwartet. **Energie- und Gradientenberechnung bleiben in AP 1 unverändert** — es geht ausschließlich um API-Vervollständigung, Cache-Hygiene und Persistierung von Berechnungsartefakten (Coordination Numbers, Energiekomponenten).

Nach AP 1 sind folgende Aufrufe auf einem `XTB`-Objekt zulässig und liefern korrekte Werte:

| Aufruf | Quelle | Liefert |
|---|---|---|
| `xtb.UpdateMolecule(geometry)` | bereits vererbt, hier nur Cache-Invalidierung | `bool` |
| `xtb.getCharges()` | bereits vorhanden | `Vector` (q_at) |
| `xtb.getCoordinationNumbers()` | **neu** | `Vector` (CN nach letztem `Calculation()`) |
| `xtb.getOrbitalEnergies()` | bereits vorhanden | `Vector` (eps) |
| `xtb.getMOCoefficients()` | bereits vorhanden | `Matrix` (C) |
| `xtb.getNumElectrons()` | **neu** | `int` (ganzzahlige nocc·2) |
| `xtb.getEnergyDecomposition()` | **neu** | `json` mit fixiertem Schema |
| `xtb.getHOMOEnergy/getLUMOEnergy/getHOMOLUMOGap()` | bereits vorhanden | `double` |
| `xtb.Gradient()` | von `QMInterface` vererbt | `Matrix::Zero` (Platzhalter; AP 4) |

---

## Vorbedingungen / Status-Check

Vor Beginn:

```bash
cd release && make -j4 curcuma 2>&1 | tail -20
ctest --output-on-failure 2>&1 | tail -20
```

Erwartet: build grün; CTests passieren wie zuvor (26/26). Falls nicht — **erst Baseline herstellen, dann AP 1**.

Die `XTB`-Klasse ist aktuell **nicht** in `MethodFactory` verdrahtet — siehe `NATIVE_XTB_STATUS.md`. AP 1 ändert das nicht.

---

## Aufgabenliste (in Bearbeitungs­reihenfolge)

### Aufgabe 1.1 — `UpdateMolecule()` überschreiben (no-arg-Variante)

**Ziel:** Cache-Invalidierung bei jedem Geometrie-Update, unabhängig vom verwendeten Setter-Overload.

**Begründung:** Heute baut `Calculation()` `m_S`, `m_H0`, `m_gamma` und `setupMultipole` *immer* neu, daher ist die Invalidierung **noch keine Korrekturmaßnahme**, sondern Design-Hygiene und Vorbereitung für AP 4 (Gradienten benötigen geometrieabhängige Zwischengrößen, die effizient gecached werden sollen).

**Rationale für no-arg-Override (Abweichung vom Roadmap):** `QMInterface` bietet fünf `UpdateMolecule()`-Overloads (Mol, Mol*, Matrix, double*, Vector). Alle delegieren am Ende an `UpdateMolecule()` (no-arg). Wenn wir nur `UpdateMolecule(const Matrix&)` überschreiben, verlieren wir die Invalidierung bei Mol-/Pointer-Updates. Die saubere Lösung ist die Geometrie-Setter unverändert zu lassen und nur die finale parameterlose Methode zu überschreiben.

**Datei:** `xtb_native.h` — Methodendeklaration (private Sektion, vor `MakeOverlap`):

```cpp
// QMInterface hook called after every geometry-setter overload.
// Invalidates all geometry-dependent cached state.
bool UpdateMolecule() override;
```

**Datei:** `xtb_native.cpp` — Implementierung (am Ende der Datei, vor schließendem Namespace):

```cpp
bool XTB::UpdateMolecule()
{
    // Invalidate geometry-dependent caches. Calculation() will rebuild them.
    m_S.resize(0, 0);
    m_H0.resize(0, 0);
    m_gamma.resize(0, 0);
    m_mp_initialized = false;
    // m_basis, m_h0 stay valid: they depend on element list, not geometry.
    // m_wfn populations are reset implicitly when SCF runs again.
    return true;
}
```

**Verifizierung:** kompilierbar; nach `xtb.UpdateMolecule(new_geometry)` ist `m_mp_initialized == false`.

---

### Aufgabe 1.2 — `m_h0.rad` befüllen

**Ziel:** Konsistenz mit TBLite-Datenstruktur und Vorbereitung für AP 4 (gradientenseitige Distanz-Polynom-Ableitung greift auf pro-Atom-Radien zu).

**Hinweis:** `getHamiltonianH0()` (xtb_h0.cpp:166) ruft heute direkt `atomic_rad_au(z)` auf. AP 1 ändert diesen Konsumenten **nicht** — Performance-/Refactoring-Schritt ist Teil von AP 4. Wir füllen nur die Datenstruktur für künftige Konsumenten.

**Datei:** `xtb_native.cpp`, Funktion `buildH0Data()`, direkt nach der `m_h0.rad.assign(...)`-Zeile:

```cpp
// Fill per-atom radii (consumed by future gradient kernels in AP 4).
// Keep in sync with atomic_rad_au() in parameters/xtb_params_extra.hpp.
for (int iat = 0; iat < m_basis.nat; ++iat) {
    m_h0.rad[iat] = atomic_rad_au(m_basis.z[iat]);
}
```

**Include hinzufügen** falls noch nicht vorhanden — ein kurzer Check zeigt, dass `parameters/xtb_params_extra.hpp` aktuell **nicht** in `xtb_native.cpp` inkludiert ist (nur in `xtb_h0.cpp`, `xtb_coulomb.cpp`, `xtb_multipole.cpp`). Daher:

```cpp
// xtb_native.cpp — bei den anderen parameter-Includes
#include "parameters/xtb_params_extra.hpp"
```

**Verifizierung:** Nach `InitialiseMolecule()` für H₂O gilt `m_h0.rad[0] ≈ atomic_rad_au(8)` (Sauerstoff) und `m_h0.rad[1] = m_h0.rad[2] ≈ atomic_rad_au(1)` (Wasserstoff).

---

### Aufgabe 1.3 — Veralteten `hscale`-Kommentar korrigieren

**Datei:** `xtb_native.cpp`, Zeilen 320–322 (im Anschluss an die `selfenergy`/`kcn`/`shpoly`-Schleife in `buildH0Data()`).

**Aktuell:**

```cpp
// hscale is Phase 3: its per-method construction rule is in
// tblite/xtb/{gfn1,gfn2}.f90 (get_hscale).
m_h0.hscale = Matrix::Zero(m_basis.nsh, m_basis.nsh);
```

**Ersetzen durch:**

```cpp
// hscale is computed on-the-fly in getHamiltonianH0() (xtb_h0.cpp), where
// per-pair (kpair * kshell * enscale * shpoly) factors are assembled directly
// from the parameter tables. m_h0.hscale stays zero-sized; the field is kept
// for future use (e.g. element-resolved caching when revisiting AP 4).
m_h0.hscale.resize(0, 0);
```

Der Wechsel von `Matrix::Zero(nsh, nsh)` auf `resize(0, 0)` macht klar, dass das Feld unbenutzt ist, und spart `nsh²` Speicher.

---

### Aufgabe 1.4 — Coordination Numbers persistieren

**Ziel:** Wrapper braucht `getCoordinationNumbers()`. Aktuell wird der CN-Vektor in `Calculation()` lokal berechnet und verworfen.

**Datei:** `xtb_native.h` — neuer privater Member-Vektor (im Bereich der gecachten Matrizen):

```cpp
Vector m_cn;   // coordination numbers from last Calculation()
```

**Datei:** `xtb_native.h` — neuer öffentlicher Getter (bei den anderen Property-Acces­soren):

```cpp
Vector getCoordinationNumbers() const { return m_cn; }
```

**Datei:** `xtb_native.cpp`, `Calculation()` Zeile 69:

```cpp
// 1. Coordination numbers
Vector cn = computeCoordinationNumbers();
m_cn = cn;   // persist for getCoordinationNumbers()
```

**Datei:** `xtb_native.cpp`, in `UpdateMolecule()` (aus 1.1):

```cpp
m_cn.resize(0);   // invalidate alongside other geometry caches
```

**Verifizierung:** Nach `Calculation()` für H₂O: `xtb.getCoordinationNumbers().size() == 3`, Werte ≈ `[1.99, 0.99, 0.99]` (CN für O ≈ 2, je H ≈ 1).

---

### Aufgabe 1.5 — `getNumElectrons()` hinzufügen

**Datei:** `xtb_native.h` — bei den Property-Accessoren:

```cpp
int getNumElectrons() const { return static_cast<int>(std::lround(m_wfn.nocc)); }
```

`m_wfn.nocc` ist `double` (in `buildReferenceOccupations()` als `total = -charge + Σ refocc` gefüllt). `std::lround` → exakte Ganzzahl (Closed-shell-Erwartung).

**Include** in `xtb_native.h`: `<cmath>` (typischerweise schon transitiv vorhanden, sicherheitshalber prüfen — `xtb_native.h` includet aktuell `<Eigen/Dense>`, `<array>`, `<memory>`, `<string>`, `<vector>`. `<cmath>` fehlt → hinzufügen).

**Verifizierung:** Für H₂O: `xtb.getNumElectrons() == 8` (O: 2s²2p⁴ = 6 Valenz; 2× H: 1s¹ = 2 Valenz; total 8).

---

### Aufgabe 1.6 — `getEnergyDecomposition()` mit fixiertem JSON-Schema

**Ziel:** Wrapper nutzt heute Keys `"electronic"`, `"repulsion"`, `"coulomb"`, `"dispersion"` (gfn2_method.cpp:140-143). Wir liefern dieselben Keys plus zusätzliche Komponenten, die die neue Implementierung sauber trennt.

**Schema (fest):**

```json
{
    "total":          double,   // m_E_total
    "electronic":     double,   // m_E_electronic (Tr(P·H0) + coulomb + third-order + multipole)
    "repulsion":      double,   // m_E_repulsion
    "coulomb":        double,   // m_E_coulomb_shell
    "third_order":    double,   // m_E_third_order
    "multipole":      double,   // m_E_multipole (0 für GFN1)
    "halogen_bond":   double,   // m_E_halogen_bond (0 für GFN2)
    "dispersion":     double    // m_E_dispersion (Stub: 0.0 — AP 4)
}
```

**Begründung der Schemawahl:**
- `"electronic"`, `"repulsion"`, `"coulomb"`, `"dispersion"` bleiben rückwärtskompatibel zum bestehenden Wrapper.
- `"third_order"`, `"multipole"`, `"halogen_bond"` sind Neuentwicklungen der modularen Architektur und erweitern die Decomposition ohne den Wrapper zu brechen.
- `"total"` als Top-Level-Sicherung ist Standard in vergleichbaren TBLite-Outputs.

**Datei:** `xtb_native.h` — bei den Property-Accessoren (Include `<nlohmann/json.hpp>` über `src/core/global.h` bereits vorhanden):

```cpp
nlohmann::json getEnergyDecomposition() const;
```

**Datei:** `xtb_native.cpp` — neue Methode (z. B. nach `getHOMOLUMOGap`):

```cpp
nlohmann::json XTB::getEnergyDecomposition() const
{
    return {
        {"total",        m_E_total},
        {"electronic",   m_E_electronic},
        {"repulsion",    m_E_repulsion},
        {"coulomb",      m_E_coulomb_shell},
        {"third_order",  m_E_third_order},
        {"multipole",    m_E_multipole},
        {"halogen_bond", m_E_halogen_bond},
        {"dispersion",   m_E_dispersion}
    };
}
```

**Verifizierung:** Nach `Calculation()` ist `decomp["total"]` ≈ Summe der übrigen Beiträge (bis auf Vorzeichenkonvention; siehe `Calculation()` Zeile 187-188: `total = electronic + repulsion + halogen_bond + dispersion`, wobei `electronic` bereits Coulomb/3rd-order/Multipole enthält).

**Hinweis zur Konsistenz:** `m_E_electronic` wird in `Calculation()` zweimal gesetzt — einmal pro Iteration mit `(P⊙H0).sum()` (Zeile 144), und einmal final mit `coulomb + third_order + multipole + Tr(P·H0)` (Zeile 180-182). Nach Abschluss von `Calculation()` enthält `m_E_electronic` die Summe der elektronischen Beiträge. Das Schema ist damit selbstkonsistent.

---

### Aufgabe 1.7 — Build- und Smoke-Test

**Build:**

```bash
cd release && make -j4 curcuma 2>&1 | tee build_ap1.log | tail -10
```

Erwartung: keine neuen Warnings, keine Errors.

**Kernel-Tests:**

```bash
ctest -R "test_xtb_overlap|test_xtb_h0|test_xtb_coulomb|test_xtb_scf_snapshot" --output-on-failure
```

Erwartung: Status wie vor AP 1 (siehe `NATIVE_XTB_STATUS.md` Zeilen 76-82). AP 1 ändert nichts an den numerischen Pfaden.

**End-to-end Smoke (nicht durch CTest abgedeckt — manuell):** Da `XTB` noch nicht in `MethodFactory` ist, wird ein kleines Standalone-Snippet im `sqm_reference`-Testsuite-Stil verwendet. Falls kein passender Test existiert: durch AP 2 abgedeckt.

---

## Akzeptanzkriterien

Markiere abgehakt, wenn lokal verifiziert:

- [ ] `release/` build kompiliert ohne neue Warnings (`build_ap1.log` sauber).
- [ ] `XTB::UpdateMolecule()` (no-arg) ist überschrieben; setzt `m_S`, `m_H0`, `m_gamma`, `m_cn` auf Größe 0 und `m_mp_initialized = false`.
- [ ] `m_h0.rad[i] == atomic_rad_au(m_basis.z[i])` nach `InitialiseMolecule()` für ein Test­molekül.
- [ ] Kommentar zu `m_h0.hscale` reflektiert tatsächliche On-the-fly-Berechnung in `xtb_h0.cpp`.
- [ ] `getCoordinationNumbers()` liefert nicht-leeren Vektor nach `Calculation()`.
- [ ] `getNumElectrons()` liefert für H₂O = 8.
- [ ] `getEnergyDecomposition()` liefert JSON mit den 8 spezifizierten Keys; `total ≈ electronic + repulsion + halogen_bond + dispersion` (Toleranz 1e-12).
- [ ] Bestehende Kernel-CTests (`test_xtb_overlap`, `test_xtb_h0`, `test_xtb_coulomb`) passieren weiter.

**Nicht-Ziele für AP 1 (gehören zu späteren Paketen):**
- ❌ MethodFactory-Wiring (AP 3)
- ❌ Wrapper-Umstellung (AP 2)
- ❌ Gradienten (AP 4)
- ❌ Refactoring von `getHamiltonianH0`, sodass `m_h0.rad` statt `atomic_rad_au` verwendet wird (Performance-Optimierung in AP 4 sinnvoller, wenn der Konsument feststeht)

---

## Risiken und Stolperfallen

| Risiko | Eintrittswahrscheinlichkeit | Gegenmaßnahme |
|---|---|---|
| `<cmath>` ist transitiv vorhanden, fällt bei einer Compiler-Version aber weg → `std::lround`-Fehler | niedrig | explizit `#include <cmath>` in `xtb_native.h` |
| `nlohmann::json` zirkulärer Include | niedrig | `global.h` bringt `json` bereits ein; siehe `xtb_native.h:30` |
| Falsche Override-Wahl bei `UpdateMolecule` lässt Mol-/Pointer-Updates ohne Cache-Invalidierung | mittel | no-arg-Override gewählt — alle Setter delegieren dorthin |
| `m_E_electronic` wird durch `Calculation()` zweimal überschrieben — Doppelzählung von Coulomb-Beiträgen in `total` | niedrig | aktuell ok, weil `total = electronic + repulsion + halogen + dispersion`; Coulomb steckt nur in `electronic` |
| `getCoordinationNumbers()` gibt vor `Calculation()` einen leeren Vektor zurück → Wrapper-Aufrufer muss damit leben | niedrig | Verhalten wie alte `GFN2`-Klasse (dort initial leer); im Wrapper wird das in AP 2 abgefangen |

---

## Test-Plan

### Standalone-Smoke (manuell ausgeführt nach AP 1)

H₂O als Referenzmolekül (`test_cases/water.xyz` ist im Repo vorhanden). Da AP 1 noch keine Factory-Integration bringt, läuft der Test nur über die kernel-tests.

**Trockenlauf:** Falls die Lust besteht, kann ein winziges main()-Programm in `test_cases/sqm_reference/` direkt eine `curcuma::xtb::XTB`-Instanz bauen und alle neuen Getter sequenziell aufrufen. Pseudocode:

```cpp
auto mol = readXYZ("water.xyz");
curcuma::xtb::XTB x(curcuma::xtb::MethodType::GFN2);
x.QMInterface::InitialiseMolecule(mol);
x.Calculation(false);
assert(x.getCoordinationNumbers().size() == 3);
assert(x.getNumElectrons() == 8);
auto d = x.getEnergyDecomposition();
assert(d.contains("electronic") && d.contains("dispersion"));
x.UpdateMolecule(mol.m_geometry);   // forces invalidation path
assert(x.getCoordinationNumbers().size() == 0);
```

Das Snippet ist **kein** Bestandteil von AP 1 als Lieferung, sondern ein optionaler Verifizierungs­schritt für den Operator.

### Regressionssperre

Nach AP 1 dürfen folgende CTests **nicht** schlechter abschneiden als vorher (Baseline siehe `NATIVE_XTB_STATUS.md`):

```bash
ctest -R "test_xtb_" --output-on-failure
ctest -R "cli_" --output-on-failure
```

---

## Fortschritt

| Datum | Status | Notizen |
|-------|--------|---------|
| 2026-04-25 | Plan erstellt | — |
| | | |

## Schwierigkeiten / Blocker

- *Noch keine dokumentiert*

---

## Referenzen

- **Roadmap-Übersicht:** [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md)
- **Status-Bestandsaufnahme:** [`NATIVE_XTB_STATUS.md`](NATIVE_XTB_STATUS.md)
- **Wrapper-Code (Konsument der hier hinzugefügten API):**
  - `src/core/energy_calculators/qm_methods/gfn2_method.cpp` — Zeilen 132-178 zeigen aktuell verwendete Getter
- **Basisklasse:** `src/core/energy_calculators/qm_methods/interface/abstract_interface.h` — `UpdateMolecule()`-Overload-Hierarchie
- **TBLite-Pendant:** `external/tblite/src/tblite/xtb/h0.f90` (`tb_hamiltonian%rad`)

---

## Änderungshistorie

| Datum | Autor | Änderung |
|-------|-------|----------|
| 2026-04-25 | Claude | Erstdokument |
