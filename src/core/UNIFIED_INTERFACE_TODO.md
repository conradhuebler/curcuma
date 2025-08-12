# TODO: Unified QM-ForceField Interface

## Vision: Merge QMInterface and ForceField systems

Currently we have two parallel interface systems:
- **QMInterface**: For quantum mechanical methods
- **ForceField**: For classical force field methods

## Proposed Unification

### 1. Create Common Base Interface
```cpp
class ComputationalMethod {
public:
    virtual ~ComputationalMethod() = default;
    virtual void setMolecule(const Mol& mol) = 0;
    virtual double calculateEnergy(bool gradient = false) = 0;
    virtual Matrix getGradient() const = 0;
    virtual Vector getCharges() const = 0;
    virtual bool hasGradient() const = 0;
    // ... other common methods
};
```

### 2. Derive Both Systems from Common Base
```cpp
class QMMethod : public ComputationalMethod {
    // QM-specific implementations
};

class ForceFieldMethod : public ComputationalMethod {
    // FF-specific implementations  
};
```

### 3. Simplify EnergyCalculator
- Single interface pointer instead of m_qminterface + m_forcefield
- Unified dispatch without special cases
- Consistent method handling

## Benefits
- **Simplified EnergyCalculator**: No more dual interface system
- **Consistent API**: Same methods for all computational methods
- **Easier Extension**: New methods inherit common interface
- **Better Polymorphism**: True OOP design

## Implementation Plan
1. Define ComputationalMethod base interface
2. Refactor existing QMInterface to inherit from base
3. Refactor ForceField to inherit from base  
4. Update EnergyCalculator to use unified interface
5. Test all methods work correctly

## Migration Strategy
- Phase 1: Add temporary setMolecule to ForceField (✅ Done)
- Phase 2: Define common interface
- Phase 3: Gradual migration of existing code
- Phase 4: Remove old dual-interface system

---
*Claude Generated: Long-term architectural improvement plan*