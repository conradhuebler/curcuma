# AP 2 — Konkreter Implementierungsplan: ComputationalMethod-Wrapper auf neue XTB-Klasse umstellen

**Status:** Bereit zur Umsetzung (nach AP 1)
**Erstellt:** 2026-04-25
**Vorgängerdokument:** [`NATIVE_XTB_ROADMAP.md`](NATIVE_XTB_ROADMAP.md), Arbeitspaket 2
**Vorbedingung:** AP 1 abgeschlossen (`getCoordinationNumbers`, `getNumElectrons`, `getEnergyDecomposition`, `UpdateMolecule()`-Override müssen in `curcuma::xtb::XTB` vorhanden sein)
**Implementation:** AI-generiert, nicht-getestet

---

## Ziel

Die Wrapper `GFN2Method` und `GFN1Method` (Konsumenten von `ComputationalMethod`) wechseln intern von der alten monolithischen `GFN2`/`GFN1`-Klasse auf die neue modulare `curcuma::xtb::XTB`-Klasse. Das öffentliche Interface der Wrapper bleibt unverändert — `EnergyCalculator`, `MethodFactory` und alle Capabilities sind nicht betroffen. **Energiewerte werden hier noch nicht numerisch validiert** — das ist Aufgabe von AP 5. AP 2 ist rein „API-Plumbing".

Nach AP 2:
- `GFN2Method` / `GFN1Method` halten `std::unique_ptr<curcuma::xtb::XTB>` statt `GFN2`/`GFN1`
- Alle bestehenden Wrapper-Methoden funktionieren mit dem neuen Backend
- `Gradient()` liefert weiterhin korrekt (Platzhalter `Matrix::Zero` aus AP 1 — echter Gradient kommt in AP 4)
- Die alten Klassen `GFN2`/`GFN1` sind **noch nicht gelöscht** (gehört zu AP 3 oder einem separaten Cleanup-AP); sie werden nur nicht mehr von den Wrappern verwendet.

---

## Vorbedingungen

```bash
cd release && make -j4 curcuma 2>&1 | tail -5
ctest --output-on-failure 2>&1 | tail -5
```

Erwartet: build grün, alle 26 CLI-Tests passen. Falls nicht — Baseline herstellen.

Sicherstellen, dass AP 1 abgeschlossen ist:

```bash
grep -n "getCoordinationNumbers\|getNumElectrons\|getEnergyDecomposition" \
    src/core/energy_calculators/qm_methods/xtb_native.h
```

Erwartung: alle drei Methoden sind dort deklariert.

---

## API-Mapping (alte GFN2/GFN1 → neue XTB)

Vor der Code-Änderung gilt es, die Namensänderungen zu fixieren. Tabelle als Referenz für die Aufgaben unten:

| Alt (`GFN2`/`GFN1`) | Neu (`curcuma::xtb::XTB`) | Bemerkung |
|---|---|---|
| `Calculation(bool)` | `Calculation(bool)` | identisch |
| `UpdateMolecule(geometry)` | `UpdateMolecule(geometry)` | von QMInterface vererbt; AP 1 hat no-arg-Override eingebaut |
| `Gradient()` | `Gradient()` | von QMInterface vererbt → `m_gradient` (Zero bis AP 4) |
| `getPartialCharges()` | `getCharges()` | Umbenennung |
| `Energies()` | `getOrbitalEnergies()` | Umbenennung |
| `MolecularOrbitals()` | `getMOCoefficients()` | Umbenennung |
| `NumElectrons()` | `getNumElectrons()` | AP 1 |
| `getCoordinationNumbers()` | `getCoordinationNumbers()` | AP 1, identischer Name |
| `getHOMOEnergy/LUMOEnergy/HOMOLUMOGap()` | identisch | identisch |
| `getEnergyDecomposition()` | `getEnergyDecomposition()` | **anderes Schema** — siehe unten |

**JSON-Schema-Diff bei `getEnergyDecomposition()`:**

| Alt-GFN2 (`gfn2.cpp:984`) | Neu (AP 1 fixiert) |
|---|---|
| `electronic` | `electronic` ✅ |
| `repulsion` | `repulsion` ✅ |
| `coulomb` | `coulomb` ✅ |
| `aes2` | **fehlt** — neu: `multipole`, `third_order` |
| `total` | `total` ✅ |
| `scf_converged` | **fehlt im neuen Schema** |
| — | `halogen_bond`, `dispersion` (neu) |

⇒ Wrapper-Code, der `aes2` oder `scf_converged` liest, **bricht**. Das betrifft aktuell nichts (gfn2_method.cpp:140-143 liest nur `electronic/repulsion/coulomb/dispersion`), muss aber in 2.5 verifiziert werden.

---

## Aufgabenliste

### Aufgabe 2.1 — Header `gfn2_method.h` umstellen

**Datei:** `src/core/energy_calculators/qm_methods/gfn2_method.h`

**Diff:**

```cpp
// Alt
#include "gfn2.h"
// ...
std::unique_ptr<GFN2> m_gfn2;
```

```cpp
// Neu
#include "xtb_native.h"
// ...
std::unique_ptr<curcuma::xtb::XTB> m_xtb;
```

Außerdem im Header: alle Doxygen/Hinweise auf "wrapped GFN2" zu "wrapped curcuma::xtb::XTB" anpassen. Konstruktor- und Destruktor-Signaturen bleiben.

**Nicht mehr benötigte private Helfer löschen:** `initializeGFN2()`, `updateGFN2Parameters()`, `handleGFN2Error()` waren auf die alte Klasse zugeschnitten. Falls Wrapper-Logik weiter benötigt wird (z. B. Element-Validation), ersetzen wir sie durch knappe Inline-Calls in `setMolecule()`/`calculateEnergy()`.

---

### Aufgabe 2.2 — Implementation `gfn2_method.cpp` umstellen

**Datei:** `src/core/energy_calculators/qm_methods/gfn2_method.cpp`

**2.2a — Konstruktor:**

```cpp
GFN2Method::GFN2Method(const json& config)
    : m_xtb(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    m_parameters = MergeJson(getDefaultConfig(), config);
    m_xtb = std::make_unique<curcuma::xtb::XTB>(curcuma::xtb::MethodType::GFN2);
}
```

**2.2b — `setMolecule`:**

```cpp
bool GFN2Method::setMolecule(const Mol& mol)
{
    try {
        if (!GFN2MethodUtils::isMoleculeSupported(mol)) {
            m_has_error = true;
            m_error_message = "Molecule contains unsupported elements for GFN2";
            return false;
        }
        m_molecule = mol;
        m_calculation_done = false;
        m_initialized = true;
        clearError();

        if (!m_xtb->QMInterface::InitialiseMolecule(mol)) {
            m_has_error = true;
            m_error_message = "GFN2 (native xTB) initialization failed";
            return false;
        }
        return true;
    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("GFN2 setMolecule failed: {}", e.what());
        return false;
    }
}
```

**2.2c — `updateGeometry`:**

```cpp
bool GFN2Method::updateGeometry(const Matrix& geometry)
{
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return false;
    }
    m_molecule.m_geometry = geometry;
    m_calculation_done = false;
    clearError();
    return m_xtb->UpdateMolecule(geometry);
}
```

**2.2d — `calculateEnergy`:**

```cpp
double GFN2Method::calculateEnergy(bool gradient)
{
    if (!m_initialized) {
        m_has_error = true;
        m_error_message = "Method not initialized - call setMolecule first";
        return 0.0;
    }
    try {
        clearError();
        m_last_energy = m_xtb->Calculation(gradient);
        m_calculation_done = true;

        if (m_parameters.value("print_orbitals", false)
            && CurcumaLogger::get_verbosity() >= 2) {
            json d = m_xtb->getEnergyDecomposition();
            CurcumaLogger::info("GFN2 (native xTB) Energy Decomposition:");
            CurcumaLogger::param("electronic",   fmt::format("{:.6f} Eh", d.value("electronic",  0.0)));
            CurcumaLogger::param("repulsion",    fmt::format("{:.6f} Eh", d.value("repulsion",   0.0)));
            CurcumaLogger::param("coulomb",      fmt::format("{:.6f} Eh", d.value("coulomb",     0.0)));
            CurcumaLogger::param("third_order",  fmt::format("{:.6f} Eh", d.value("third_order", 0.0)));
            CurcumaLogger::param("multipole",    fmt::format("{:.6f} Eh", d.value("multipole",   0.0)));
            CurcumaLogger::param("dispersion",   fmt::format("{:.6f} Eh", d.value("dispersion",  0.0)));
        }
        return m_last_energy;
    } catch (const std::exception& e) {
        m_has_error = true;
        m_error_message = fmt::format("GFN2 (native xTB) energy calculation failed: {}", e.what());
        return 0.0;
    }
}
```

`d.value(key, default)` schützt gegen fehlende Keys — wichtig, falls `getEnergyDecomposition()`-Schema in Zukunft erweitert wird.

**2.2e — Property-Getter:**

```cpp
Matrix GFN2Method::getGradient() const {
    if (!m_calculation_done || !m_xtb) return Matrix::Zero(m_molecule.m_number_atoms, 3);
    Matrix grad = m_xtb->Gradient();
    if (CurcumaLogger::get_verbosity() >= 3) {
        CurcumaLogger::param("GFN2Method_gradient_norm", fmt::format("{:.6e} Eh/Å", grad.norm()));
    }
    return grad;
}

Vector GFN2Method::getCharges() const {
    if (!m_calculation_done || !m_xtb) return Vector::Zero(m_molecule.m_number_atoms);
    return m_xtb->getCharges();
}

Vector GFN2Method::getOrbitalEnergies() const {
    if (!m_calculation_done || !m_xtb) return Vector{};
    return m_xtb->getOrbitalEnergies();
}

Matrix GFN2Method::getMolecularOrbitals() const {
    if (!m_calculation_done || !m_xtb) return Matrix{};
    return m_xtb->getMOCoefficients();
}

int GFN2Method::getNumElectrons() const {
    if (!m_calculation_done || !m_xtb) return 0;
    return m_xtb->getNumElectrons();
}

double GFN2Method::getHOMOEnergy()  const { return m_xtb && m_calculation_done ? m_xtb->getHOMOEnergy()  : 0.0; }
double GFN2Method::getLUMOEnergy()  const { return m_xtb && m_calculation_done ? m_xtb->getLUMOEnergy()  : 0.0; }
double GFN2Method::getHOMOLUMOGap() const { return m_xtb && m_calculation_done ? m_xtb->getHOMOLUMOGap() : 0.0; }

Vector GFN2Method::getCoordinationNumbers() const {
    if (!m_calculation_done || !m_xtb) return Vector{};
    return m_xtb->getCoordinationNumbers();
}

json GFN2Method::getEnergyDecomposition() const {
    if (!m_calculation_done || !m_xtb) return json{};
    return m_xtb->getEnergyDecomposition();
}
```

**2.2f — `try/catch` durchforsten:** Wegen `XTB::Calculation()` keine Exceptions mehr wirft (es liefert `0.0` und einen Logger-Warn bei Konvergenz-Fehler) sind die `try/catch` bei den Getter-Methoden unnötig. Klassische Ausnahme: Initialisierung. Daher: Wrapper-Code-Reduktion, lediglich `setMolecule` und `calculateEnergy` behalten ihre `try/catch`-Blöcke.

**2.2g — `setParameters`:** `updateGFN2Parameters()` löschen. `setParameters()` setzt nur `m_parameters` — die ConfigManager-Verteilung erfolgt heute via `EnergyCalculator`-Konstruktor, nicht im Wrapper-Setter.

```cpp
void GFN2Method::setParameters(const json& params) {
    m_parameters = MergeJson(m_parameters, params);
    // No live propagation: XTB-Parameter (scf_max_iter, threshold, damping) sind
    // im Konstruktor gesetzt. Live-Updates wären Zukunftsfeature in AP 4/5.
}
```

---

### Aufgabe 2.3 — `gfn1_method.h/.cpp` analog umstellen

**Header `gfn1_method.h`:**

```cpp
#include "xtb_native.h"
// ...
std::unique_ptr<curcuma::xtb::XTB> m_xtb;
```

**Implementation `gfn1_method.cpp`:**

Konstruktor:

```cpp
GFN1Method::GFN1Method(const json& config)
    : m_xtb(nullptr)
    , m_calculation_done(false)
    , m_last_energy(0.0)
{
    json full_config = getDefaultConfig();
    if (!config.empty()) full_config.merge_patch(config);
    m_xtb = std::make_unique<curcuma::xtb::XTB>(curcuma::xtb::MethodType::GFN1);

    if (CurcumaLogger::get_verbosity() >= 2) {
        CurcumaLogger::info("GFN1Method initialized with native xTB implementation");
    }
}
```

Übrige Methoden 1:1 nach dem Muster aus 2.2e (replace `m_gfn1` durch `m_xtb`, Methodennamen aus dem Mapping). `getCoordinationNumbers`, `getEnergyDecomposition` bleiben.

`saveToFile()` ist im aktuellen `GFN1Method` nur ein `return false`-Stub — nichts zu ändern.

---

### Aufgabe 2.4 — Build und Wrapper-Smoke-Test

**Build:**

```bash
cd release && make -j4 curcuma 2>&1 | tee build_ap2.log | tail -10
```

Erwartet: keine neuen Compiler-Fehler. Likely warnings: ungenutzte alte `GFN2`/`GFN1`-Includes via Headerketten — sind harmlos, werden in AP 6 (Cleanup) entfernt.

**Smoke-Test über `ngfn2`/`ngfn1`** (diese gehen via Factory direkt an `GFN2Method`/`GFN1Method`, ohne Priority-Chain — perfekter Test von AP 2 isoliert von AP 3):

```bash
cd release
./curcuma -sp ../test_cases/sqm_reference/molecules/H2O.xyz -method ngfn2 2>&1 | tail -20
./curcuma -sp ../test_cases/sqm_reference/molecules/H2O.xyz -method ngfn1 2>&1 | tail -20
```

Erwartet: Curcuma terminiert ohne Crash, gibt eine Energie aus. Energie-Wert ist **noch nicht validiert** — Validierung erfolgt in AP 5. Ziel hier ist: kein Segfault, keine Exception, sinnvoller Output (nicht NaN, nicht 0.0).

**Regression CTest:**

```bash
ctest -R "cli_curcumaopt|cli_sqm" --output-on-failure
```

⚠️ `cli_curcumaopt_02_gfn2_single_point` und `cli_sqm_05_gfn2_singlepoint` laufen heute über `gfn2` mit TBLite-Priority. Nach AP 2 läuft `gfn2` **immer noch über TBLite** (AP 3 stellt erst um). Daher dürfen diese Tests sich **nicht** ändern. `ngfn2`-/`ngfn1`-Pfade haben aktuell keine CTests — die kommen in AP 5.

---

### Aufgabe 2.5 — Konsumenten-Audit

Verifizieren, dass kein Aufrufer `getEnergyDecomposition()` mit dem alten Schlüssel `aes2` oder `scf_converged` liest:

```bash
grep -rn 'energy_decomposition\|getEnergyDecomposition' src/ test_cases/ | \
    grep -v 'qm_methods/gfn[12]' | \
    head -20
```

Erwartung: nur `gfn2_method.cpp` (von AP 2 selbst geändert) und `gfn1_method.cpp`. Falls weitere Konsumenten existieren, deren Output mit `d.value(key, 0.0)` defensiv lesen.

---

## Akzeptanzkriterien

- [ ] `release/` build kompiliert ohne neue Errors. Warnings nur „unused include"-Klasse.
- [ ] `gfn2_method.h/.cpp` enthält keinen `GFN2`-Typ mehr, ausschließlich `curcuma::xtb::XTB`.
- [ ] `gfn1_method.h/.cpp` analog.
- [ ] `./curcuma -sp H2O.xyz -method ngfn2` terminiert mit endlicher Energie (kein NaN/Inf, kein Crash).
- [ ] `./curcuma -sp H2O.xyz -method ngfn1` analog.
- [ ] Bestehende `cli_curcumaopt_*` und `cli_sqm_*` CTests passieren weiter (kein Regression auf `gfn2`-Pfad, weil AP 3 noch nicht aktiv).
- [ ] `getEnergyDecomposition()`-Konsumenten greifen nur via `value(key, default)` zu — kein direktes `[]`-Indexing für Keys, die im neuen Schema fehlen.

**Nicht-Ziele:**
- ❌ Numerische Übereinstimmung mit TBLite (AP 5)
- ❌ Funktionierende Optimierung mit `ngfn2` (Gradient = Zero ⇒ konvergiert nicht; AP 4)
- ❌ Default-Methode `gfn2` benutzt neue Implementation (AP 3)
- ❌ Löschung der alten `gfn2.cpp`/`gfn1.cpp` (Cleanup in AP 6)

---

## Risiken und Stolperfallen

| Risiko | Eintrittswahrscheinlichkeit | Gegenmaßnahme |
|---|---|---|
| `XTB`-Konstruktor nimmt enum, alte `GFN2()` nimmt nichts → Kompilierfehler beim Konstruktor-Init | hoch (offensichtlich) | Konstruktor explizit umbauen wie in 2.2a |
| `m_gfn2->Gradient()` liefert ≠ Zero in der Alt-Implementation, `m_xtb->Gradient()` liefert Zero → Optimierungstests scheitern | mittel | Tests, die `gfn2`-Optimierung erwarten, laufen weiter über TBLite-Priority (AP 3 noch nicht aktiv). `ngfn2`-Optimierungen scheitern bewusst — nicht Teil von AP 2 |
| `getEnergyDecomposition()` mit fehlendem Key wirft `json::type_error` | mittel | Defensive `d.value(key, default)` statt `d[key]` |
| `print_orbitals=true` mit Verbosity ≥ 2 enthält jetzt mehr Felder → User-Output ändert sich | niedrig (Output, kein API-Bruch) | dokumentieren in `AIChangelog.md` als Verbesserung |
| Compiler warnt vor ungenutztem `gfn2.h`-Include in der Headerkette | niedrig | belassen — wird in AP 6 (Cleanup) entfernt |

---

## Test-Plan

### Manueller Smoke-Lauf (für Operator)

```bash
cd release
for mol in H2 H2O CH4 NH3 LiH; do
    echo "=== $mol GFN2 (ngfn2 = neu) ==="
    ./curcuma -sp ../test_cases/sqm_reference/molecules/${mol}.xyz -method ngfn2 -verbosity 1 2>&1 | grep -i "energy\|error"
    echo "=== $mol GFN1 (ngfn1 = neu) ==="
    ./curcuma -sp ../test_cases/sqm_reference/molecules/${mol}.xyz -method ngfn1 -verbosity 1 2>&1 | grep -i "energy\|error"
done
```

Erwartung: für alle fünf Moleküle eine endliche Energie. Werte werden in AP 5 numerisch geprüft.

### Regressions-Sperre

```bash
ctest --output-on-failure
```

Erwartet: 26/26 wie vor AP 2.

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

- **Vorgänger:** [`NATIVE_XTB_AP1_PLAN.md`](NATIVE_XTB_AP1_PLAN.md)
- **Nachfolger:** [`NATIVE_XTB_AP3_PLAN.md`](NATIVE_XTB_AP3_PLAN.md)
- **Wrapper-Quelldateien:**
  - `src/core/energy_calculators/qm_methods/gfn2_method.h/.cpp`
  - `src/core/energy_calculators/qm_methods/gfn1_method.h/.cpp`
- **Neues Backend:** `src/core/energy_calculators/qm_methods/xtb_native.h`
- **Altes Backend (zur Kontrolle, noch nicht löschen):** `gfn2.h/.cpp`, `gfn1.h/.cpp`

---

## Änderungshistorie

| Datum | Autor | Änderung |
|-------|-------|----------|
| 2026-04-25 | Claude | Erstdokument |
