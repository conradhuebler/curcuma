# EnergyCalculator Big-Bang Refactoring Plan

## Ziel
Vereinheitlichung der QM/MM-Methoden durch polymorphe Architektur mit einheitlichem ComputationalMethod Interface, bei Erhaltung des bestehenden Hierarchiesystems und Threading-Kompatibilität.

## Architektur-Prinzipien
- **Polymorphie** (virtuelle Funktionen) statt Templates für bessere Verständlichkeit
- **Hierarchiesystem beibehalten**: gfn2 kann von TBLite, Ulysses, XTB oder zukünftig native kommen
- **Prioritäten-System**: TBLite > Ulysses > XTB für shared methods (wie aktuell)
- **Big-Bang Refactoring** der Implementierung, aber API/Interface-Kompatibilität
- **Threading**: ForceField threading erhalten, alle anderen Methoden threading-fähig machen

## Aktuelles Hierarchiesystem (beibehalten)
```cpp
// Method Resolution Priority (aus aktuellem SwitchMethod):
"gfn2": TBLite(1) > Ulysses(3) > XTB(2) > Native(zukünftig)
"gfn1": TBLite(1) > XTB(2) > Native(zukünftig)  
"cgfnff": ForceField(0) - explizit
"ugfn2": Ulysses(3) - explizit
"eht": EHT(6) - explizit
"d3": DFT-D3(4) - explizit
"d4": DFT-D4(5) - explizit
```

## Implementierungsplan

### Phase 1: ComputationalMethod Base Interface
**Datei**: `src/core/energy_calculators/computational_method.h`

```cpp
class ComputationalMethod {
public:
    virtual ~ComputationalMethod() = default;
    
    // Core functionality
    virtual bool setMolecule(const Mol& mol) = 0;
    virtual bool updateGeometry(const Matrix& geometry) = 0;
    virtual double calculateEnergy(bool gradient = false, bool verbose = false) = 0;
    
    // Property access
    virtual Matrix getGradient() const = 0;
    virtual Vector getCharges() const = 0;
    virtual Vector getBondOrders() const = 0;
    virtual Position getDipole() const = 0;
    virtual bool hasGradient() const = 0;
    
    // Method information
    virtual std::string getMethodName() const = 0;
    virtual bool isThreadSafe() const = 0;
    virtual void setThreadCount(int threads) = 0;
    
    // Configuration
    virtual void setParameters(const json& params) = 0;
    virtual json getParameters() const = 0;
    virtual bool hasError() const = 0;
};
```

### Phase 2: Method Factory mit Prioritäten-System
**Datei**: `src/core/energy_calculators/method_factory.h/.cpp`

```cpp
class MethodFactory {
public:
    static std::unique_ptr<ComputationalMethod> create(
        const std::string& method_name, 
        const json& config
    );
    
private:
    // Priority resolution for shared methods
    static std::unique_ptr<ComputationalMethod> createGFN2(const json& config);
    static std::unique_ptr<ComputationalMethod> createGFN1(const json& config);
    
    // Explicit method creation
    static std::unique_ptr<ComputationalMethod> createEHT(const json& config);
    static std::unique_ptr<ComputationalMethod> createForceField(const json& config);
    static std::unique_ptr<ComputationalMethod> createUlysses(const json& config);
    
    // Compilation checks
    static bool hasTBLite();
    static bool hasXTB(); 
    static bool hasUlysses();
    static bool hasD3();
    static bool hasD4();
};
```

### Phase 3: Directory Structure
```
src/core/energy_calculators/
├── computational_method.h          # Base interface
├── method_factory.h/.cpp           # Priority-based factory
├── qm_methods/                     # QM method wrappers
│   ├── eht_method.h/.cpp           # Native EHT
│   ├── gfnff_method.h/.cpp         # Native GFN-FF (cgfnff)
│   ├── xtb_method.h/.cpp           # XTB wrapper
│   ├── tblite_method.h/.cpp        # TBLite wrapper
│   ├── ulysses_method.h/.cpp       # Ulysses wrapper
│   └── dispersion_method.h/.cpp    # D3/D4 wrappers
└── ff_methods/
    └── forcefield_method.h/.cpp    # ForceField wrapper (threading beibehalten)
```

### Phase 4: Method Wrappers Implementation

#### QM Method Wrappers
Alle bestehenden QMInterface Implementierungen werden durch ComputationalMethod Wrapper ersetzt:

- **EHTMethod**: Wraps bestehende EHT Klasse
- **XTBMethod**: Wraps XTBInterface, fügt Threading hinzu
- **TBLiteMethod**: Wraps TBLiteInterface, fügt Threading hinzu
- **UlyssesMethod**: Wraps UlyssesInterface, fügt Threading hinzu
- **DispersionMethod**: Wraps DFTD3/DFTD4Interface

#### ForceField Wrapper
**ForceFieldMethod**: Wraps bestehende ForceField Klasse, erhält CxxThreadPool System vollständig.

### Phase 5: EnergyCalculator Refactoring
**Vereinfachung der EnergyCalculator Klasse**:

```cpp
class EnergyCalculator {
private:
    std::unique_ptr<ComputationalMethod> m_method;  // Single interface pointer
    
    // Remove: m_qminterface, m_forcefield, m_eht, etc.
    // Remove: Complex SwitchMethod() logic
    
public:
    EnergyCalculator(const std::string& method, const json& controller);
    
    void setMolecule(const Mol& mol);
    double calculateEnergy(bool gradient = false, bool verbose = false);
    
    // All property access delegated to m_method
    Matrix getGradient() const { return m_method->getGradient(); }
    Vector getCharges() const { return m_method->getCharges(); }
    // etc.
};
```

### Phase 6: Threading Strategy
- **ForceField**: Bestehendes CxxThreadPool System komplett erhalten
- **QM Methods**: Threading nach ForceField-Vorbild implementieren
- **EnergyCalculator**: Thread-safe durch einheitliches Interface
- **Configuration**: Threading per Method konfigurierbar

### Phase 7: API-Erhaltung
**Bestehende APIs bleiben verfügbar**:
- EnergyCalculator Constructor signatures unverändert
- Method strings ("uff", "eht", "gfnff", "gfn2", etc.) kompatibel
- Parameter handling (JSON) konsistent
- Property access methods unverändert

### Phase 8: Migration Steps
1. **Base Interface** und Factory erstellen
2. **Method Wrappers** parallel zu bestehenden Klassen implementieren
3. **EnergyCalculator** komplett ersetzen (Big-Bang)
4. **Tests** für alle Methoden durchlaufen
5. **Alte Interfaces** schrittweise entfernen

## Performance-Überlegungen
- **Virtual Function Overhead**: ~1-3ns pro Aufruf - vernachlässigbar bei chemischen Berechnungen (μs-ms Range)
- **Memory Layout**: Polymorphie hat bessere Cache-Lokalität als Template-Instanziierung  
- **Threading**: Method-spezifische Threading-Kontrolle für optimale Performance

## Testing Strategy
- **Unit Tests**: Für jeden Method Wrapper
- **Integration Tests**: EnergyCalculator mit allen Methoden
- **Performance Tests**: Threading performance validation
- **Regression Tests**: Alle bestehenden Testfälle müssen passieren

## Benefits
- **Unified Architecture**: Single interface für alle Berechnungsmethoden
- **Simplified EnergyCalculator**: Keine komplexe Switch-Logik mehr
- **Better Threading**: Konsistente Threading-Unterstützung
- **Easier Extension**: Neue Methoden einfach hinzufügbar
- **Maintained Compatibility**: Alle bestehenden APIs erhalten

---

**Status**: Ready for implementation
**Estimated Time**: 2-3 development cycles
**Risk Level**: Medium (Big-Bang approach, but with comprehensive testing)