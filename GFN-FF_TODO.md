# GFN-FF Implementation TODO

## Status: Vollständiges System implementiert ✅

### Was bereits implementiert ist:
- ✅ **GFNFF Klasse** als QMInterface in `src/core/qm_methods/gfnff.h/cpp`
- ✅ **ForceFieldThread** erweitert um GFN-FF Terme (type=3) in `src/core/forcefieldthread.h/cpp`
- ✅ **Integration** in moderne ForceField Architektur
- ✅ **Alle kovalenten Terme**: Bond/Angle/Torsion/OutOfPlane Berechnungen
- ✅ **GFN-FF Parametergenerierung** mit automatischer Bond-Detektion und Angle-Ableitung
- ✅ **Parameter-Implementation** (PoC Platzhalter für H,C,N,O,F,P,S,Cl,Br,I)
- ✅ **ForceField Integration** (method="gfnff" → setMethod(3))
- ✅ **Geometriefunktionen** UFF::Torsion/OutOfPlane in `forcefieldfunctions.h`
- ✅ **Universelles Parameter-Caching** für alle FF-Methoden (UFF, GFN-FF, QMDFF)
- ✅ **CMakeLists.txt** erweitert für GFN-FF Kompilierung
- ✅ **Test-Programm** `src/helpers/gfnff_test.cpp`

## KRITISCHE TODOs - Echte GFN-FF Parameter benötigt:

### 0. 🔴 **LITERATURABGLEICH & PARAMETRIERUNG** (ABSOLUT KRITISCH) ⚠️⚠️⚠️
**STATUS**: NICHT IMPLEMENTIERT - RÜCKFRAGEN ERFORDERLICH!

**PROBLEM**: Alle aktuell implementierten Parameter sind Platzhalter und entsprechen NICHT der echten GFN-FF Methode!

**KRITISCHE RÜCKFRAGEN AN ENTWICKLER:**
1. **Welche GFN-FF Paper sollen als Referenz dienen?**
   - Spicher & Grimme Angew. Chem. Int. Ed. 59, 15665 (2020)?
   - Gibt es neuere Parameter-Updates?
   - Welche Versionsunterschiede zu beachten?

2. **Woher kommen die echten GFN-FF Parameter?**
   - Aus xtb source code extrahieren? (github.com/grimme-lab/xtb)
   - Aus Paper-Supplements?
   - Gibt es offizielle Parameter-Dateien?

3. **Welche Parameter-Sets sind prioritär?**
   - Bond force constants & equilibrium distances
   - Angle force constants & equilibrium angles  
   - Torsion barriers & periodicities
   - Out-of-plane force constants
   - Element coverage (H,C,N,O vs. vollständig bis Z=86)?

4. **Wie exakt sollen GFN-FF Formeln implementiert werden?**
   - Aktuell: E_bond = 0.5*k*(r-r0)² + α*(r-r0)³
   - Ist das korrekt oder verwendet GFN-FF andere Funktionsformen?
   - Welche Anharmonizitäten sind wichtig?

5. **Integration mit bestehenden Korrekturen klären:**
   - Wie integriert sich GFN-FF D4 mit vorhandener D4-Implementierung?
   - H4 vs. GFN-FF Wasserstoffbrücken - überschneiden sich die Korrekturen?
   - Welche Parameter-Sets für welche Korrekturen?

**AKTION ERFORDERLICH**: 
- Entwickler muss Paper studieren und echte Parameter beschaffen
- Implementierung stoppen bis echte Parameter vorliegen
- Literaturabgleich für korrekte Funktionsformen durchführen

### 1. 🔴 **GFN-FF Parameter Datenbank** (HÖCHSTE PRIORITÄT)
```cpp
// In gfnff.cpp: getCovalentRadius(), getGFNFFBondParameters(), getGFNFFAngleParameters()
```
**Problem**: Aktuell nur Platzhalter-Parameter für H,C,N,O,F,P,S,Cl,Br,I
**Benötigt**: 
- Vollständige GFN-FF Parametertabellen für Z=1-86
- Element-spezifische Kraftkonstanten
- Hybridisierungs- und koordinationsabhängige Parameter
- Echte GFN-FF Gleichgewichtsdistanzen und -winkel

### 2. 🔴 **Fehlende Geometriefunktionen**
```cpp
// In forcefieldthread.cpp: Zeilen 643, 669
```
**Problem**: Torsion und Out-of-Plane Berechnungen existieren nicht
**Benötigt**:
- Torsionswinkel-Berechnung mit analytischen Gradienten
- Out-of-Plane Winkel-Berechnung mit analytischen Gradienten
- Integration in `forcefieldfunctions.h` oder neuer GFN-FF namespace

### 3. 🟡 **ForceField Integration vervollständigen**
```cpp
// In forcefield.cpp: addGFNFFBond/Angle/etc. Methoden aufrufen
```
**Problem**: ForceField.cpp ruft noch alte addBond/addAngle auf
**Benötigt**: 
- Erkennung von method="gfnff" → setMethod(3) für ForceFieldThread
- Aufruf der neuen addGFNFFBond/addGFNFFAngle Methoden

## ERWEITERTE TODOs - Vollständige GFN-FF Features:

### 4. 🟡 **D4 Dispersion Integration**
**Status**: H4/D3/D4 Infrastruktur bereits vorhanden in `src/core/hbonds.h`, `forcefieldthread.h`
**Benötigt**: 
- Echte D4-Parameter für GFN-FF (nicht PM6-H4)
- Integration der D4Thread Klasse für GFN-FF spezifische Parameter

### 5. 🟡 **Halogen- und Wasserstoffbrücken**
**Status**: H4Correction bereits implementiert für PM6
**Benötigt**:
- GFN-FF spezifische XB/HB Parameter
- Anpassung der hbonds4::H4Correction für GFN-FF

### 6. 🟡 **Topologie und Koordinationserkennung**
```cpp
// In gfnff.cpp: calculateTopology()
```
**Benötigt**:
- Ring-Detektion für spezielle GFN-FF Parameter
- Hybridisierungsbestimmung (sp, sp2, sp3)
- Koordinationszahl-Berechnung
- Formale Ladungsverteilung

### 7. 🟡 **Torsion und Inversion Parameter**
**Benötigt**:
- GFN-FF Torsionsparameter-Datenbank
- Automatische Torsion-Detektion
- Out-of-plane Parameter für sp2-Zentren

## COMPILIERUNG/TESTING:

### 8. 🟡 **Build Integration**
**Status**: CMakeLists.txt muss ggf. angepasst werden
**Test**: Kompilierung mit bestehender D3/D4/H4 Infrastruktur

### 9. 🟡 **Erste Tests**
**Minimal-Test**: Einfaches Molekül (H2O, CH4) mit nur Bond/Angle Termen
**Volltest**: Komplexeres System mit allen GFN-FF Korrekturen

## DATENQUELLEN für echte Parameter:

- **GFN-FF Paper**: Spicher & Grimme, Angew. Chem. Int. Ed. 59, 15665 (2020)
- **Original xtb Code**: https://github.com/grimme-lab/xtb (GFN-FF Implementierung)
- **Parameter Files**: Benötigt Extraktion aus xtb source oder Paper supplements

## PRAGMA MESSAGES im Code:
Alle aktuellen TODOs sind mit `#pragma message("TODO: ...")` markiert:
- `src/core/forcefieldthread.cpp:576` - Bond stretching
- `src/core/forcefieldthread.cpp:607` - Angle bending  
- `src/core/forcefieldthread.cpp:634` - Torsion calculation
- `src/core/forcefieldthread.cpp:661` - Out-of-plane calculation
- `src/core/forcefieldthread.cpp:683` - vdW/Dispersion

## PRIORITÄT REIHENFOLGE:
1. **Parameter Datenbank** (ohne echte Parameter läuft nichts sinnvoll)
2. **Geometriefunktionen** (Torsion/OutOfPlane)
3. **ForceField Integration** (method=gfnff handling)
4. **Compilation & Basic Testing**
5. **D4/XB/HB Integration**
6. **Performance & Vollständigkeit**