# GFN-FF Native Implementation Status & Roadmap

**Goal**: Replace external Fortran GFN-FF library with native C++ implementation (`cgfnff`)

**NEW**: Comprehensive roadmap created → See `docs/GFNFF_NATIVE_ROADMAP.md`

---

## **WICHTIG: Neue Strategie (2025-11-10)** 🎯

### Externe Library als Referenz nutzen
Die externe Fortran-Implementierung ist jetzt verfügbar und dient als **Validierungs-Referenz**:

```bash
# Initialisiert (bereits erledigt in dieser Session):
git submodule update --init external/gfnff

# Externe Library befindet sich in:
external/gfnff/src/  # 42.849 Zeilen Fortran, 368 Funktionen
```

**Alle Rückfragen vom alten TODO sind beantwortet**:
1. ✅ **Referenz**: `external/gfnff/` ist die offizielle Spicher/Grimme Implementation
2. ✅ **Parameter**: Alle in `gfnff_param.f90` vollständig Z=1-86
3. ✅ **Formeln**: In `gfnff_engrad.F90` komplett dokumentiert
4. ✅ **Integration**: Bestehende D4/H4 Infrastruktur kann genutzt werden

---

## **Aktueller Status: ~50% Komplett** ⚠️

### Was funktioniert ✅
- ✅ **Bond Stretching**: Einfache anharmonische Implementation
- ✅ **Angle Bending**: Fourier-Expansion
- ✅ **Parameter-Arrays**: Alle Z=1-86 aus Fortran extrahiert (chi, gam, cnf, alp, rad, bond, angle)
- ✅ **Basis-Topologie**: Einfache CN, Hybridization (neighbor-based)
- ✅ **ForceField Backend**: Integration funktioniert
- ✅ **ConfigManager**: Parameter-System integriert
- ✅ **CurcumaLogger**: Verbosity-Control implementiert

### Was fehlt ❌ (Kritisch)
- ❌ **Torsions/Dihedrals**: Arrays leer, keine Energy/Gradient Implementation
- ❌ **Inversions/Out-of-Plane**: Arrays leer, keine Implementation
- ❌ **EEQ Charges**: Nur simple Placeholder (H=+0.1, C=-0.1, etc.) statt Matrix-Löser
- ❌ **Non-bonded/vdW**: D4 existiert aber nicht gekoppelt
- ❌ **Ring Detection**: Gibt nur Nullen zurück
- ❌ **Pi-System Detection**: Leer
- ❌ **CN Derivatives**: Nur Null-Matrizen (benötigt für Gradienten)
- ❌ **H-Bond Detection**: Gibt leeres Array

### Problematische Design-Entscheidungen ⚠️
```cpp
// PROBLEM 1: Equilibrium aus aktueller Geometrie statt Parametertabelle
params.equilibrium_distance = distance;  // gfnff.cpp:492

// PROBLEM 2: Willkürliche Scaling Factors
params.force_constant = angle_param * 0.001;  // gfnff.cpp:533

// PROBLEM 3: Hardcoded Werte
params.anharmonic_factor = -0.1;  // gfnff.cpp:495
```

---

## **Neue Dokumentation** 📚

Alle Details sind jetzt strukturiert dokumentiert:

1. **`docs/GFNFF_NATIVE_ROADMAP.md`** (NEU) ⭐
   - 8-Phasen Implementierungsplan
   - Detaillierte Aufgabenlisten pro Phase
   - Zeitschätzungen (~6-8 Wochen full-time)
   - Validierungsstrategien
   - Erfolgskriterien

2. **`docs/GFNFF_FORTRAN_FUNCTIONS.md`** (NEU) ⭐
   - Mapping Fortran → C++ Funktionen
   - Code-Porting-Patterns
   - Test-Strategien
   - Quick Reference für Entwicklung

3. **`scripts/validate_gfnff_native.py`** (NEU) ⭐
   - Automatische Validierung native vs. extern
   - Energie- und Gradienten-Vergleich
   - Report-Generierung

---

## **Roadmap Kurzfassung** (Details in ROADMAP.md)

### **Phase 1: Kritische Energie-Terme** (1-2 Wochen) ⚡
- [ ] Torsion Energy/Gradient (`egtors` portieren)
- [ ] Inversion Energy/Gradient (Out-of-Plane)
- [ ] ForceField Integration vervollständigen
- **Deliverable**: Bonds+Angles+Torsions+Inversions funktional

### **Phase 2: Topologie-Algorithmen** (1.5-2 Wochen) 🔍
- [ ] Ring Detection (DFS/BFS, 3-20 Ringe)
- [ ] Pi-System Detection (Konjugation)
- [ ] Erweiterte Hybridization Detection (Geometrie-basiert)
- **Deliverable**: Topologie-Informationen für Parameter-Korrekturen

### **Phase 3: EEQ Charge Calculation** (2-3 Wochen) ⚡ KOMPLEX
- [ ] Matrix-Setup (A·q = b)
- [ ] Linear Solver (Eigen LU)
- [ ] CN Derivatives (3D-Tensor)
- [ ] Fragment-Constraints
- **Deliverable**: Echte EEQ-Ladungen, korrekte Elektrostatik

### **Phase 4: Non-bonded Interaktionen** (1-2 Wochen) 🌐
- [ ] D4 Dispersion Integration (bestehende D4Interface nutzen)
- [ ] Repulsion Term (short-range)
- [ ] H-Bond Detection & Energy
- **Deliverable**: Vollständige Nicht-kovalente Wechselwirkungen

### **Phase 5: Parameter-Korrekturen** (1 Woche) 🔧
- [ ] Equilibrium-Werte aus Tabellen (nicht aktuelle Geometrie!)
- [ ] Scaling Factor Validierung
- [ ] Topologie-abhängige Parameter (Ring-Strain, Konjugation)
- **Deliverable**: Wissenschaftlich korrekte Parameter

### **Phase 6: Validierung & Testing** (1-2 Wochen) ✅
- [ ] 20+ Test-Moleküle mit Referenz-Daten
- [ ] Automatische Validierung (`validate_gfnff_native.py`)
- [ ] Accuracy Benchmarks (±0.5 kcal/mol Ziel)
- **Deliverable**: 95% Agreement mit Fortran-Referenz

### **Phase 7: Performance-Optimierung** (1 Woche) ⚡
- [ ] Profiling (gprof, perf)
- [ ] Neighbor Lists (O(N) statt O(N²))
- [ ] Parallelisierung (Bond/Angle Berechnungen)
- **Deliverable**: ≤2x langsamer als Fortran

### **Phase 8: Default Integration** (1 Woche) 📚
- [ ] MethodFactory Priority Update (Native > External > XTB)
- [ ] Dokumentation (README, CLAUDE.md)
- [ ] CMake: Native immer verfügbar, External optional
- **Deliverable**: cgfnff als Production-Default

**Total**: ~14 Wochen konservativ, 6-8 Wochen fokussiert

---

## **Nächste Schritte für Entwickler** 🚀

### Sofort starten:
```bash
# 1. Fortran-Code studieren
less external/gfnff/src/gfnff_engrad.F90  # Zeile 1041: Torsions

# 2. Test-Molekül vorbereiten
mkdir -p test_cases/gfnff_validation/hydrocarbons
# Butane.xyz erstellen für Torsion-Test

# 3. Referenz-Daten generieren
./build/curcuma -sp test_cases/.../butane.xyz -method gfnff > butane_ref.out

# 4. Torsion implementieren
vim src/core/energy_calculators/qm_methods/gfnff.cpp
# → calculateTorsionEnergy() hinzufügen

# 5. Testen
./build/curcuma -sp test_cases/.../butane.xyz -method cgfnff > butane_nat.out
python scripts/validate_gfnff_native.py
```

### Entwicklungs-Workflow:
1. **Feature aus Roadmap wählen** (z.B. Torsions)
2. **Fortran-Code analysieren** (`gfnff_engrad.F90`)
3. **Test schreiben** (externe Referenz generieren)
4. **C++ implementieren** (Energy → Gradient)
5. **Validieren** (validate_gfnff_native.py)
6. **Dokumentieren** (ROADMAP.md updaten)

---

## **Vorteile Native Implementation** 💡

### Technisch
- ✅ **Keine Fortran-Dependencies** (gfortran, LAPACK nicht benötigt)
- ✅ **Einfaches Debugging** (C++ Debugger, keine Fortran-C Bridge)
- ✅ **Volle Integration** mit Curcuma-Ecosystem
- ✅ **Thread-Safe** by Design (Eigen thread-safe)

### Pädagogisch
- ✅ **Code-Klarheit** > komplexe Optimierungen
- ✅ **Literatur-Referenzen** bei jeder Formel
- ✅ **Lernbar** für Studenten/Entwickler
- ✅ **Erweiterbar** ohne Fortran-Kenntnisse

### Wissenschaftlich
- ✅ **Vollständige Kontrolle** über Implementierung
- ✅ **Einfache Modifikationen** für Forschung
- ✅ **Transparenz** - jede Formel nachvollziehbar
- ✅ **Validierbar** gegen Fortran-Referenz

---

## **Validierungs-Ziele**

### Energie
- **Ziel**: ±0.5 kcal/mol vs. Fortran
- **Aktuell**: Bonds/Angles ~0.3 kcal/mol (gut!)
- **Problem**: Torsions/Inversions fehlen → große Abweichungen

### Gradienten
- **Ziel**: ±1% vs. Fortran (relative Abweichung)
- **Aktuell**: Bonds/Angles ~2% (akzeptabel)
- **Problem**: CN-Derivatives fehlen → Gradient unvollständig

### Geometrie-Optimierung
- **Ziel**: Identische Minima wie Fortran
- **Aktuell**: Nicht getestet (Gradienten unvollständig)
- **Blocker**: EEQ Charges, Torsions, Inversions

---

## **Alte TODOs - Jetzt beantwortet** ✅

Die ursprünglichen kritischen Rückfragen sind durch externe Library gelöst:

1. ✅ **Parameter-Quelle**: `external/gfnff/src/gfnff_param.f90`
2. ✅ **Formeln**: `external/gfnff/src/gfnff_engrad.F90`
3. ✅ **Topologie**: `external/gfnff/src/gfnff_ini.f90`
4. ✅ **D3/D4 Integration**: In Fortran-Code sichtbar
5. ✅ **Referenz-Tests**: Fortran als Ground Truth

**Alte TODOs 1-9 sind obsolet** - siehe stattdessen `docs/GFNFF_NATIVE_ROADMAP.md`

---

## **Zusammenfassung**

### Stand der Dinge
- **Architektur**: ✅ Komplett (QMInterface, ForceField, MethodFactory)
- **Basis-Terme**: ✅ Bonds/Angles funktional
- **Kritische Terme**: ❌ Torsions/Inversions/EEQ fehlen
- **Validierung**: ⚠️ Nur teilweise (50% Feature-Set)

### Weg nach vorne
1. **Roadmap folgen** (`docs/GFNFF_NATIVE_ROADMAP.md`)
2. **Fortran-Code portieren** (systematisch, Phase für Phase)
3. **Gegen Fortran validieren** (`scripts/validate_gfnff_native.py`)
4. **Dokumentieren & Testen** (CI/CD Integration)
5. **Als Default setzen** (MethodFactory Priority)

### Zeitrahmen
- **Minimal-Funktionalität** (Phases 1-4): 6-8 Wochen
- **Production-Ready** (Phases 1-8): 12-14 Wochen
- **Optimiert & Dokumentiert**: 14-16 Wochen

**Viel Erfolg beim Ersetzen der Fortran-Library!** 🚀

---

*Letzte Aktualisierung: 2025-11-10*
*Nächste Review: Nach Phase 1 Completion*
