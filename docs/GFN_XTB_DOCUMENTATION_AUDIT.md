# GFN/xTB Dokumenten-Analyse: Redundanzen und Obsolete Teile

**Erstellt**: 2026-04-25  
**Umfang**: Alle GFN- und xTB-bezogenen Dokumente in `docs/`, `src/helpers/`, `external/docs/`, Root-Level

---

## Zusammenfassung

| Kategorie | Anzahl | Status |
|---|---|---|
| Aktuell und relevant | 4 | Beibehalten |
| Veraltet, aber archivwürdig | 5 | Nach `docs/archive/` verschieben |
| Direkt obsolet / redundant | 4 | Löschen oder nach `archive/` |
| Inkonsistent zwischen Dokumenten | 2 | Konsolidieren |

---

## 1. Native QM-Methoden (GFN1/GFN2) — Dokumente

### 1.1 `docs/NATIVE_QM_IMPLEMENTATION_STATUS.md` — VERALTET

**Status**: 🔴 Veraltet (Nov 2025), irreführend für aktuelle Entwicklung

**Probleme**:
- Beschreibt die **alte** `gfn2.cpp` (~1236 Zeilen) mit Status "0/7 vs TBLite"
- Die neue modulare `xtb_native.cpp` + Kernel-Module existiert seit Apr 2026 und ist nicht dokumentiert
- "Phase 4: Parameter-Integration" als TODO gelistet — wurde bereits durch `scripts/extract_xtb_params.py` gelöst (auto-generierte `parameters/gfn2_params.hpp`)
- "Phase 5: Analytische Gradienten" als TODO — die alte `gfn2.cpp` HAT tatsächlich Gradienten, die neue `xtb_native.cpp` hat sie noch nicht
- Tabelle "Implementierte Features" zeigt AES2 als "vereinfacht" — in der neuen Implementierung ist AES2 vollständig
- Keine Erwähnung der neuen `curcuma::xtb::XTB`-Klasse oder der TBLite-Architektur-Portierung

**Empfehlung**: Nach `docs/archive/` verschieben. Ein neues Status-Dokument für die neue Implementierung erstellen (oder `NATIVE_XTB_ROADMAP.md` erweitern).

---

### 1.2 `docs/NATIVE_QM_METHODS_IMPLEMENTATION.md` — TEILWEISE VERALTET

**Status**: 🟡 Teilweise veraltet (Nov 2025), 1600+ Zeilen

**Probleme**:
- Enthält extrem viel **Pseudocode** (~800 Zeilen), der eine hypothetische Klasse beschreibt, die nicht der aktuellen Implementierung entspricht
- Der Pseudocode in Sektion "Core Algorithm Structure" beschreibt ein monolithisches Design mit `class GFN2 : public QMDriver`, das der alten `gfn2.cpp` ähnelt
- Die tatsächliche neue Implementierung (`xtb_native.cpp`) ist modular und spiegelt TBLites Architektur wider (separate Dateien für H0, Coulomb, Multipole, SCF)
- Parameter-Struktur `GFN2Parameters` im Pseudocode ist vereinfacht — die tatsächliche Implementierung nutzt auto-generierte `gfn2_params.hpp` mit echten TBLite-Parametern
- Pseudocode für `calculateEffectiveCoulomb()` ist stark vereinfacht (keine shell-resolved Gamma-Matrix, kein Klopman-Ohno)
- "Implementation Roadmap" Phase 1-5 ist veraltet — die neue Architektur ist bereits implementiert
- Abschnitt "Phase 5: MethodFactory Integration" zeigt Lambda-Registry, die in `method_factory.cpp` nicht verwendet wird (GCC-15-Kompatibilität)

**Was behalten werden sollte**:
- Referenzabschnitt (Primär-Literatur, DOIs, Quellen)
- Vergleichs-Tabelle GFN1 vs GFN2
- Theoretischer Hintergrund (Tight-Binding, NDDO)

**Empfehlung**: Aufteilen — theoretische Referenzen behalten, Pseudocode und Roadmap entfernen oder nach `archive/` verschieben.

---

### 1.3 `docs/NATIVE_XTB_ROADMAP.md` — AKTUELL

**Status**: 🟢 Aktuell (Apr 2026)

**Bemerkung**: Gerade erstellt. Dies ist die einzige aktuelle Dokumentation der neuen `curcuma::xtb::XTB`-Implementierung.

---

### 1.4 `docs/GFN2_MULTIPOLE_IMPLEMENTATION.md` — AKTUELL

**Status**: 🟢 Aktuell und relevant

**Bemerkung**: Beschreibt korrekt die Multipol-Implementierung in `xtb_multipole.cpp`. Verweist auf TBLite-Module. Korrekt und nützlich.

---

### 1.5 `GFN1_FIXES_SUMMARY.md` (Root-Level) — ARCHIVWÜRDIG

**Status**: 🟡 Archivwürdig (historisch)

**Probleme**:
- Dokumentiert Fixes an der **alten** `gfn1.cpp` (self-energy sign, Hamiltonian construction)
- Die alte `gfn1.cpp` wird durch `xtb_native.cpp` mit `MethodType::GFN1` ersetzt
- Diese Fixes sind im neuen Code bereits korrekt implementiert (die `getHamiltonianH0()` in `xtb_h0.cpp` nutzt das korrekte Vorzeichen)

**Empfehlung**: Nach `docs/archive/` verschieben. Historischer Wert für das Verständnis der Entwicklung.

---

### 1.6 `external/docs/GFN2_PARAMETER_EXTRACTION.md` — VERALTET

**Status**: 🔴 Veraltet (Jan 2025)

**Probleme**:
- Beschreibt die **manuelle** Parameter-Extraktion aus TBLite Fortran-Source
- "Result: ~420 lines vs ~4000+ for manual element-by-element approach"
- Die aktuelle Implementierung nutzt `scripts/extract_xtb_params.py`, das direkt aus `external/tblite/src/tblite/xtb/gfn2.f90` extrahiert und `parameters/gfn2_params.hpp` generiert
- Der Ansatz "compact, data-driven" wurde durch den Auto-Generator überholt

**Empfehlung**: Nach `docs/archive/` verschieben oder löschen. Die Pipeline ist jetzt vollständig automatisiert.

---

## 2. GFN-FF Dokumente

### 2.1 `docs/GFNFF_STATUS.md` — AKTUELL

**Status**: 🟢 Aktuell (Apr 2026)

**Bemerkung**: Das zentrale Status-Dokument für GFN-FF. Wird regelmäßig aktualisiert. Enthält die aktuellsten Fehleranalysen und Fixes.

---

### 2.2 `docs/GFNFF_IMPLEMENTATION_HUB.md` — AKTUELL

**Status**: 🟢 Aktuell (Feb 2026)

**Bemerkung**: Gute Einstiegsseite mit Verweisen auf andere Dokumente. Struktur ist sauber.

---

### 2.3 `docs/GFNFF_GRADIENTS.md` — INKONSISTENT

**Status**: 🟡 Inkosistent mit `GFNFF_STATUS.md`

**Probleme**:
- Sagt "BATM grad TODO" in der Gradient-Tabelle
- `GFNFF_STATUS.md` (Mar 2026) sagt "ATM ✅ 100% — Separated to own GradientATM() component" und "Gradients ✅ 85%"
- "Validation Results" sind von Feb 2026 und wurden in `GFNFF_STATUS.md` übertroffen
- Einheitliche Dokumentation: "ForceField internally calculates gradients in Hartree/Bohr" — dies ist eigentlich ein generelles Konzept, nicht GFN-FF-spezifisch

**Empfehlung**: Mit `GFNFF_STATUS.md` konsolidieren oder nach `archive/` verschieben, wenn `GFNFF_STATUS.md` alle relevanten Gradienten-Informationen enthält.

---

### 2.4 `docs/GFNFF_PHYSICS_AUDIT.md` — ARCHIVWÜRDIG

**Status**: 🟡 Archivwürdig (Jan 2026)

**Probleme**:
- Dokumentiert historische Untersuchungen (Coulomb alpha-Parameter, C-H R0-Abweichung, EEQ Self-Energy)
- Alle drei Issues sind als "✅" markiert (behoben)
- Enthält keine aktuell offenen Probleme
- Die Informationen sind auch in `GFNFF_STATUS.md` (historische Abschnitte) enthalten

**Empfehlung**: Nach `docs/archive/` verschieben.

---

### 2.5 `docs/GFNFF_DISPERSION_FIX.md` — ARCHIVWÜRDIG

**Status**: 🟡 Archivwürdig (Jan 2026)

**Probleme**:
- Dokumentiert den D4-CN-only-Weighting-Fix und den D4-as-Default-Fix
- Beide Fixes sind in `GFNFF_STATUS.md` (Abschnitt "Previous: Dispersion WEIGHT_THRESHOLD Fix") ebenfalls dokumentiert
- Keine neuen Informationen gegenüber `GFNFF_STATUS.md`

**Empfehlung**: Nach `docs/archive/` verschieben.

---

### 2.6 `docs/GPU_GFNNF_DISCREPANCIES.md` — AKTUELL (spezifisch)

**Status**: 🟢 Aktuell (März 2026), aber GPU-spezifisch

**Bemerkung**: Enthält GPU-spezifische Untersuchungen, die nicht in `GFNFF_STATUS.md` abgedeckt sind. Relevant für CUDA-Entwicklung. Beibehalten.

---

## 3. Code-Level Artefakte (Nicht-Dokumente)

### 3.1 `gfn2_params_full.hpp` (Root-Level) — OBSOLET

**Status**: 🔴 Obsolet

**Probleme**:
- Alte Parameter-Header, manuell erstellt
- Die neue Implementierung nutzt `src/core/energy_calculators/qm_methods/parameters/gfn2_params.hpp` (auto-generiert)
- Root-Level Header wird nicht mehr eingebunden

**Empfehlung**: Löschen oder nach `archive/` verschieben.

---

### 3.2 `src/helpers/gfnff_helper.cpp` — UNKLAR

**Status**: 🟡 Status unklar

**Bemerkung**: Hilfsprogramm für GFN-FF. Muss geprüft werden, ob es noch genutzt wird.

---

### 3.3 `src/helpers/xtb_helper.cpp` — UNKLAR

**Status**: 🟡 Status unklar

**Bemerkung**: Hilfsprogramm für xTB. Muss geprüft werden, ob es noch genutzt wird.

---

## 4. Inkonsistenzen zwischen Dokumenten

### 4.1 Gradienten-Status GFN-FF

| Dokument | BATM Gradient | ATM Gradient |
|---|---|---|
| `GFNFF_GRADIENTS.md` | "📋 Planned" | Nicht erwähnt |
| `GFNFF_STATUS.md` | "~0.000 mEh ✅ FIXED" (Mar 6) | "✅ 100%" (Mar 12) |

**Empfehlung**: `GFNFF_GRADIENTS.md` aktualisieren oder nach `archive/` verschieben.

---

### 4.2 Native GFN2-Status

| Dokument | Implementierungsstatus | Teststatus |
|---|---|---|
| `NATIVE_QM_IMPLEMENTATION_STATUS.md` | "0/7 vs TBLite, major errors" | 0/7 |
| `src/core/energy_calculators/qm_methods/CLAUDE.md` | "0/7 vs TBLite, major errors" | 0/7 |
| `NATIVE_XTB_ROADMAP.md` | Neue modulare Implementierung in Arbeit | Noch nicht validiert |

**Empfehlung**: `NATIVE_QM_IMPLEMENTATION_STATUS.md` muss entweder aktualisiert oder archiviert werden. `CLAUDE.md` sollte auf die neue Implementierung verweisen.

---

### 4.3 Dispersion-Status

| Dokument | D4-Integration |
|---|---|
| `NATIVE_QM_IMPLEMENTATION_STATUS.md` | "⏸️ Stub" |
| `NATIVE_QM_METHODS_IMPLEMENTATION.md` | "Stub functions initially" |
| `gfn2.cpp` (alt) | `calculateDispersionEnergy()` gibt 0.0 zurück |
| `xtb_native.cpp` (neu) | `m_E_dispersion = 0.0; // D3/D4 — to be added` |
| `GFNFF_STATUS.md` | D4 vollständig implementiert und validiert |

**Bemerkung**: Die GFN-FF-D4-Implementierung ist vollständig, aber die GFN2-xTB-D4-Integration ist noch offen. Das ist korrekt, sollte aber klarer dokumentiert sein.

---

## 5. Redundante Inhalte über Dokumente hinweg

### 5.1 Theoretische Grundlagen

- `NATIVE_QM_METHODS_IMPLEMENTATION.md` enthält ~200 Zeilen Theorie (Tight-Binding, NDDO)
- `docs/theory/GFNFF_COMPLETE_GUIDE.md` enthält ebenfalls Grundlagen zur xTB-Theorie
- `docs/theory/GFNFF_EEQ_THEORY.md` spezialisiert auf EEQ

**Empfehlung**: Theoretische Grundlagen sollten in einem zentralen Dokument konsolidiert werden, nicht über 3+ Dateien verteilt.

---

### 5.2 TBLite-Referenzen

Fast jedes Dokument enthält die gleichen TBLite-Referenzen:
- `NATIVE_QM_IMPLEMENTATION_STATUS.md`: "TBLite: https://github.com/tblite/tblite"
- `NATIVE_QM_METHODS_IMPLEMENTATION.md`: "TBLite Library: https://github.com/tblite/tblite"
- `GFN2_MULTIPOLE_IMPLEMENTATION.md`: Verweise auf `tblite/xtb/...`
- `NATIVE_XTB_ROADMAP.md`: Verweise auf `external/tblite/src/tblite/...`

**Empfehlung**: Zentrale Referenzdatei für alle TBLite-Referenzen.

---

## 6. Empfohlene Aktionen

### Sofort (keine Risiken)

| Aktion | Datei(en) | Grund |
|---|---|---|
| **Archivieren** | `docs/NATIVE_QM_IMPLEMENTATION_STATUS.md` | Veraltet, irreführend |
| **Archivieren** | `GFN1_FIXES_SUMMARY.md` | Historisch, alte Implementierung |
| **Archivieren** | `external/docs/GFN2_PARAMETER_EXTRACTION.md` | Auto-Generator existiert |
| **Archivieren** | `docs/GFNFF_PHYSICS_AUDIT.md` | Alle Issues behoben |
| **Archivieren** | `docs/GFNFF_DISPERSION_FIX.md` | In `GFNFF_STATUS.md` enthalten |
| **Löschen** | `gfn2_params_full.hpp` (Root) | Nicht genutzt |

### Konsolidierung (mittelfristig)

| Aktion | Datei(en) | Grund |
|---|---|---|
| **Aktualisieren** | `docs/GFNFF_GRADIENTS.md` oder archivieren | Inkosistent mit `GFNFF_STATUS.md` |
| **Aktualisieren** | `docs/NATIVE_QM_METHODS_IMPLEMENTATION.md` | Pseudocode entfernen, aktuelle Architektur beschreiben |
| **Erweitern** | `docs/NATIVE_XTB_ROADMAP.md` | Fortschrittstabelle füllen |

### Überprüfung nötig

| Datei | Was prüfen |
|---|---|
| `src/helpers/gfnff_helper.cpp` | Wird es noch gebaut/genutzt? |
| `src/helpers/xtb_helper.cpp` | Wird es noch gebaut/genutzt? |
| `src/helpers/gfnff_term_validator.cpp` | Wird es noch gebaut/genutzt? |
| `test_cases/GFNFF_TESTING.md` | Ist es aktueller als `GFNFF_STATUS.md`? |
| `test_cases/GFNFF_TEST_RESULTS_2026-01-12.txt` | Hat es Ersatz durch CTest? |

---

## 7. Vorgeschlagene neue Dokumentenstruktur

```
docs/
├── NATIVE_XTB_ROADMAP.md              ← Aktuell (neue Implementierung)
├── NATIVE_XTB_STATUS.md               ← NEU: Status der neuen Implementierung
├── GFNFF_STATUS.md                    ← Aktuell (zentrales GFN-FF-Dokument)
├── GFNFF_IMPLEMENTATION_HUB.md        ← Aktuell (Einstiegspunkt)
├── GFN2_MULTIPOLE_IMPLEMENTATION.md   ← Aktuell (Multipol-Details)
├── GPU_GFNNF_DISCREPANCIES.md         ← Aktuell (GPU-spezifisch)
├── theory/
│   ├── GFNFF_COMPLETE_GUIDE.md
│   └── GFNFF_EEQ_THEORY.md
└── archive/
    ├── NATIVE_QM_IMPLEMENTATION_STATUS.md    ← Verschoben
    ├── NATIVE_QM_METHODS_IMPLEMENTATION.md   ← Verschoben (oder gekürzt)
    ├── GFN1_FIXES_SUMMARY.md                  ← Verschoben
    ├── GFNFF_GRADIENTS.md                     ← Verschoben (oder aktualisiert)
    ├── GFNFF_PHYSICS_AUDIT.md                 ← Verschoben
    ├── GFNFF_DISPERSION_FIX.md                ← Verschoben
    └── GFN2_PARAMETER_EXTRACTION.md           ← Verschoben
```

---

## Anhang: Dateien, die nicht analysiert wurden

| Datei | Grund |
|---|---|
| `test_cases/GFNFF_TESTING.md` | Konnte nicht gelesen werden |
| `test_cases/GFNFF_TEST_RESULTS_2026-01-12.txt` | Konnte nicht gelesen werden |
| `docs/theory/PHASE1.3_FORMULA_FIXES.md` | Vermutlich historisch, nicht gelesen |
| `docs/theory/PHASE2_TOPOLOGY_DETECTION.md` | Vermutlich historisch, nicht gelesen |
| `docs/theory/PHASE3_EEQ_CHARGES.md` | Vermutlich historisch, nicht gelesen |
| `docs/archive/*` | Bereits archiviert |
