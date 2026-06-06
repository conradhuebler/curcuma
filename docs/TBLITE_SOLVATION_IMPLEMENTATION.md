# TBLite Solvation Implementation Analysis
## Comprehensive Guide to Extracting Solvation into Native Curcuma Code

**Date**: November 2025  
**Analysis Scope**: Complete solvation model implementation in TBLite library integration  
**Target**: Extracting and implementing solvation models natively in Curcuma

---

## 1. SOLVATION MODELS IMPLEMENTED

### 1.1 Overview

TBLite provides **three distinct solvation models** for implicit solvation calculations:

| Model | Code | Algorithm Type | Default? | Status |
|-------|------|----------------|----------|--------|
| **CPCM** | 1 | Conductor-like PCM | No | ✅ Functional |
| **GB** | 2 | Generalized Born | Yes | ❌ Disabled (bug) |
| **ALPB** | 3 | Analytical Linearized Poisson-Boltzmann | No | ✅ Functional |

### 1.2 Model Selection in Code

**File**: `/home/user/curcuma/src/core/energy_calculators/qm_methods/tbliteinterface.cpp`  
**Function**: `TBLiteInterface::ApplySolvation()` (lines 258-323)

```cpp
// Solvation model selection (line 268-322)
if (m_solvent_model == 1)      // CPCM
    m_tb_cont = tblite_new_cpcm_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps);
else if (m_solvent_model == 2)  // GB (currently broken)
    m_tb_cont = tblite_new_gb_solvation_epsilon(m_error, m_tblite_mol, m_solvent_eps, 
                                                 m_solvent_gb_version, m_solvent_gb_kernel);
else if (m_solvent_model == 3)  // ALPB
    m_tb_cont = tblite_new_alpb_solvation_solvent(m_error, m_tblite_mol, m_solvent, 
                                                   m_solvent_alpb_version, m_solvent_alpb_reference);
```

---

## 2. CORE SOLVATION CODE STRUCTURE

### 2.1 Main Interface File

**Location**: `/home/user/curcuma/src/core/energy_calculators/qm_methods/tbliteinterface.h`

**Key Components**:
- Parameter definitions (lines 31-50)
- Solvation state variables (lines 102-108)
- `ApplySolvation()` method declaration (line 88)

**Parameter Registry** (Claude Generated 2025):
```cpp
BEGIN_PARAMETER_DEFINITION(tblite)
    PARAM(solvent_model, Int, 0, "Solvent model type (0=none, 1=CPCM, 2=GB, 3=ALPB)", "Solvation")
    PARAM(solvent, String, "none", "Solvent name (e.g., 'water', 'acetone', 'none')", "Solvation")
    PARAM(solvent_epsilon, Double, -1.0, "Solvent dielectric constant (epsilon). -1 = auto-detect", "Solvation")
    PARAM(solvent_gb_version, Int, 0, "Generalized Born model version (0=GBSA, 1=ALPB)", "Solvation")
    PARAM(solvent_gb_kernel, Int, 1, "GB kernel function type", "Solvation")
    PARAM(solvent_alpb_version, Int, 12, "ALPB solvation model version", "Solvation")
    PARAM(solvent_alpb_reference, Int, 1, "ALPB reference state", "Solvation")
END_PARAMETER_DEFINITION
```

### 2.2 Solvation State Variables

**Header File** (tbliteinterface.h, lines 102-108):
```cpp
double m_solvent_eps = -1;           // Dielectric constant (-1 = auto-detect)
int m_solvent_model = 0;             // Model ID (1=CPCM, 2=GB, 3=ALPB)
char* m_solvent = NULL;              // Solvent name string
int m_solvent_gb_version = 0;        // GB version selector
int m_solvent_gb_kernel = 1;         // GB kernel function type
int m_solvent_alpb_version = 12;     // ALPB version
int m_solvent_alpb_reference = 1;    // ALPB reference state
```

**Implementation** (tbliteinterface.cpp, lines 55-70):
```cpp
// Constructor parameter loading via ConfigManager
m_solvent = new char[tmp.length() + 1];
strcpy(m_solvent, tmp.c_str());
m_solvent_eps = m_config.get<double>("solvent_epsilon", -1.0);

// Auto-enable solvation if solvent is specified
if (tmp != "none" && tmp != "" && tmp != "vacuum") {
    m_solvent_model = m_config.get<int>("solvent_model", 2);  // Default: GB
} else {
    m_solvent_model = m_config.get<int>("solvent_model", 0);  // No solvation
}
```

### 2.3 Solvation Application Flow

**Call Stack**:
1. **Calculation()** (line 325) → Initializes calculator
2. **ApplySolvation()** (line 408) → Applied before each SCF step
3. **tblite_calculator_push_back()** → Adds solvation to calculator

**Integration Point** (line 408 in Calculation):
```cpp
ApplySolvation();
tblite_set_calculator_save_integrals(m_ctx, m_tblite_calc, 0);
tblite_get_singlepoint(m_ctx, m_tblite_mol, m_tblite_calc, m_tblite_res);
```

---

## 3. SOLVATION MODELS DETAILED

### 3.1 CPCM (Conductor-like Polarizable Continuum Model)

**Code ID**: 1  
**API Function**: `tblite_new_cpcm_solvation_epsilon()`  
**Status**: ✅ Functional

#### Algorithm Overview:
- Treats solvent as a **dielectric continuum** with conductor-like boundary
- Molecular surface divides space: cavity (molecule) + continuum (solvent)
- Surface charge distribution induced by molecular density
- Electrostatic energy = induced surface charge interaction with solute

#### Parameters Required:
```
- epsilon (ε):  Solvent dielectric constant
- Solute density from SCF calculation
```

#### Implementation Details:

**Dielectric Input**:
- Can specify epsilon directly: `m_solvent_eps > -1`
- Auto-detected from solvent name if not specified
- Example: water (ε=78.4), DMSO (ε=46.7), acetone (ε=20.7)

**Solute Density**:
- Obtained from prior SCF calculation via TBLite calculator
- Surface constructed from molecular geometry
- Continuous self-consistent field (CSCF) iteration

#### Energy Term:
```
E_CPCM = (1/2) * ∫ σ(r) * φ(r) dS

where:
σ(r) = surface charge density
φ(r) = electrostatic potential at surface
```

#### Verbosity Output (Level 2+):
```cpp
CurcumaLogger::param("solvation_method", "CPCM");
CurcumaLogger::param("dielectric_constant", fmt::format("{:.3f}", m_solvent_eps));
```

#### Code Location:
- **Setup**: tbliteinterface.cpp, lines 268-282
- **Error handling**: Automatic retry with increased verbosity (lines 274-281)
- **Integration**: `tblite_calculator_push_back()` adds CPCM to SCF calculator

---

### 3.2 GB (Generalized Born) - CURRENTLY BROKEN

**Code ID**: 2  
**API Function**: `tblite_new_gb_solvation_epsilon()`  
**Status**: ❌ **DISABLED** (exits with error at line 291-292)

#### Important Note:
**THE GB MODEL IS INTENTIONALLY DISABLED IN CURCUMA!**
```cpp
CurcumaLogger::error("GB solvation model not functional");
exit(1);
```
This happens at line 291, preventing any GB calculation.

#### Algorithm Overview (Theoretical):
- **Approximate method** for implicit solvation
- Based on: Still et al., J. Am. Chem. Soc. 112, 6127 (1990)
- Uses **effective Born radii** for each atom
- Pairwise atom-based interaction model

#### Parameters Required:
```
- epsilon (ε):       Solvent dielectric constant
- gb_version:        Model variant (0=standard, 1=modified)
- gb_kernel:         Kernel function type (1=default)
- Born radii (r_i):  Effective radius for each atom (computed internally)
```

#### Algorithm (Theoretical):
```cpp
// Born radius calculation
r_i = r_vdw,i - offset * r_vdw,i

// Generalized Born energy
E_GB = -0.5 * (1/ε_solute - 1/ε_solvent) * sum_ij [q_i * q_j / f_GB(r_ij)]

where:
f_GB(r_ij) = sqrt(r_ij² + r_i * r_j * exp(-r_ij²/(4*r_i*r_j)))
```

#### Implementation Status:
- Function exists in C API but produces errors in TBLite
- Likely issue: Born radius calculation or kernel evaluation
- Known workaround: Use ALPB instead (which works)

#### Code Location:
- **Setup**: tbliteinterface.cpp, lines 284-302
- **Error handling**: Same retry mechanism as CPCM (lines 294-300)
- **Critical**: Exit(1) prevents execution (line 292)

#### Why It's Broken:
- TBLite's GB implementation may require additional parameters not exposed via C API
- The C wrapper functions may not properly initialize GB parameters
- Possible incomplete kernel implementation in TBLite Fortran core

---

### 3.3 ALPB (Analytical Linearized Poisson-Boltzmann)

**Code ID**: 3  
**API Function**: `tblite_new_alpb_solvation_solvent()`  
**Status**: ✅ **Fully Functional**

**This is the RECOMMENDED solvation model for TBLite in Curcuma**

#### Algorithm Overview:
- **Advanced variant** of Generalized Born
- Linearized Poisson-Boltzmann equation solution
- Uses **analytical approximation** instead of iterative solving
- Fitted parameters for accuracy improvements
- Reference: Ehlert et al., J. Chem. Theory Comput. 17, 4250 (2021)

#### Parameters Required:
```
- solvent_name:         String identifier (e.g., "water", "acetone")
- alpb_version:         Integer version (typically 12)
- alpb_reference:       Reference state (typically 1)
- Born radii (computed):  Automatically from atomic coordinates
```

#### Energy Model:
```
E_ALPB = E_electrostatic + E_cavity + E_dispersion_repulsion

where:
- E_electrostatic: Interaction with induced surface charge
- E_cavity:        Non-polar cavity formation energy
- E_dispersion_repulsion: Van der Waals corrections
```

#### Parameter Details:

**ALPB Version** (alpb_version):
- `12`: Default version, most commonly used
- Other versions available but `12` is recommended

**Reference State** (alpb_reference):
- `1`: Reference solvent state (default)
- `0`: Alternative parameterization

**Solvent Names**:
Built-in supported solvents (from SOLVATION.md):
```
water, methanol, ethanol, octanol, acetone, dmso, dmf,
acetonitrile, nitromethane, benzene, toluene, aniline,
chloroform, dichloromethane, thf, dioxane, furane,
diethyl ether, ethyl acetate, hexane, hexadecane,
carbon disulfide, phenol, benzaldehyde, octanol wet
```

#### Key Advantages:
1. **No direct epsilon requirement** - uses built-in solvent parameters
2. **More accurate** than standard GB due to fitting
3. **Smooth gradients** - good for geometry optimization
4. **Fast** - analytical formula, no iterative solving

#### Verbosity Output (Level 2+):
```cpp
CurcumaLogger::param("solvation_method", "ALPB");
CurcumaLogger::param("solvent", std::string(m_solvent));
CurcumaLogger::param("alpb_version", fmt::format("{}", m_solvent_alpb_version));
CurcumaLogger::param("alpb_reference", fmt::format("{}", m_solvent_alpb_reference));
```

#### Code Location:
- **Setup**: tbliteinterface.cpp, lines 304-321
- **Error handling**: Automatic retry with verbosity increase (lines 312-319)
- **Integration**: `tblite_calculator_push_back()` adds ALPB to SCF calculator

---

## 4. SOLVENT PARAMETER SYSTEM

### 4.1 Parameter Sources

**TBLite Integration**:
- Solvent database is **compiled into TBLite library**
- Parameters are **NOT extracted from configuration files**
- Accessed via solvent name string passed to C API

**Parameter Content**:
For each solvent, TBLite stores:
```
- Dielectric constant (ε)
- Born radius offset (for GB)
- Surface tension (for cavity term)
- van der Waals parameters
- Probe radius (for surface calculation)
```

### 4.2 Parameter Files (NOT in Curcuma)

TBLite's Fortran source code contains solvent parameters. Since TBLite is external:
- Parameters are **compiled into the library**
- Individual parameter values are **not directly accessible** to Curcuma
- Only via C API functions that reference solvent by name

**To extract parameters**:
1. Inspect TBLite source code (https://github.com/tblite/tblite)
2. Look for: `solvent_data.f90` or similar
3. Extract Born radii, cavity parameters, dispersion coefficients

### 4.3 Automatic Dielectric Constant Detection

**Implementation** (tbliteinterface.cpp, lines 55-71):
```cpp
// If solvent specified but epsilon not given
if (m_solvent_eps < 0) {
    // CPCM needs explicit epsilon
    // ALPB has built-in solvent database
    // GB uses built-in parameters
}
```

**Usage Example**:
```bash
# ALPB - solvent DB used automatically
curcuma -sp mol.xyz -method gfn2:tblite -solvent water

# CPCM - requires epsilon or solvent name resolution
curcuma -sp mol.xyz -method gfn2:tblite -solvent water -solvent_model 1
```

---

## 5. KEY FUNCTIONS AND C API INTERFACE

### 5.1 TBLite C API Functions Used

```cpp
// CPCM Creation
tblite_container tblite_new_cpcm_solvation_epsilon(
    tblite_error error,
    tblite_structure mol,
    double epsilon
);

// GB Creation (currently broken)
tblite_container tblite_new_gb_solvation_epsilon(
    tblite_error error,
    tblite_structure mol,
    double epsilon,
    int version,
    int kernel
);

// ALPB Creation (working)
tblite_container tblite_new_alpb_solvation_solvent(
    tblite_error error,
    tblite_structure mol,
    const char* solvent,
    int version,
    int reference
);

// Integration with calculator
tblite_calculator_push_back(
    tblite_context ctx,
    tblite_calculator calc,
    tblite_container* cont
);
```

### 5.2 Data Structures

```cpp
// TBLite context - manages calculations
tblite_context m_ctx;

// Molecular structure
tblite_structure m_tblite_mol;

// Solvation container (added to calculator)
tblite_container m_tb_cont;

// Main calculator
tblite_calculator m_tblite_calc;

// Results holder
tblite_result m_tblite_res;
```

### 5.3 Error Handling Pattern

```cpp
// Create solvation model
m_tb_cont = tblite_new_cpcm_solvation_epsilon(...);

// Check for errors
if (tblite_check_context(m_ctx)) {
    tbliteError();           // Print C-level error
    tbliteContextError();    // Print context error
    // Retry with increased verbosity
    tblite_set_context_verbosity(m_ctx, 3);
    tblite_delete_container(&m_tb_cont);
    m_tb_cont = tblite_new_cpcm_solvation_epsilon(...);
}

// Add to calculator
tblite_calculator_push_back(m_ctx, m_tblite_calc, &m_tb_cont);
```

---

## 6. SOLVATION INTEGRATION WITH ULYSSES

### 6.1 Ulysses Solvation Support

**File**: `/home/user/curcuma/src/core/energy_calculators/qm_methods/ulyssesinterface.cpp`

**Model**: GBSA (Generalized Born + Surface Area)

**Implementation**:
```cpp
// UlyssesInterface constructor
m_solvent = m_config.get<std::string>("solvent", "none");

// Calculation function
m_ulysses->setSolvent(m_solvent);
m_ulysses->Calculate(gradient, verbose);
```

**Supported Solvents** (line 72 in header):
```cpp
StringList m_solvents = {
    "acetone", "acetonitrile", "aniline", "benzaldehyde", "benzene",
    "dichloromethane", "chloroform", "carbon disulfide", "dioxane", "dmf",
    "dmso", "ethanol", "diethyl ether", "ethyl acetate", "furane",
    "hexadecane", "hexane", "methanol", "nitromethane", "octanol",
    "phenol", "thf", "toluene", "water", "octanol wet"
};
```

**Validation** (ulyssesinterface.cpp, lines 63-66):
```cpp
if (std::find(m_solvents.begin(), m_solvents.end(), m_solvent) == m_solvents.end()) {
    CurcumaLogger::warn("Solvent '" + m_solvent + "' not supported by Ulysses");
    m_solvent = "none";
}
```

**GBSA Features**:
- Generalized Born electrostatic solvation
- Surface area-dependent cavitation energy
- Non-polar cavity formation term

---

## 7. DOCUMENTATION AND USER INTERFACE

### 7.1 Main Documentation

**File**: `/home/user/curcuma/docs/SOLVATION.md`

**Key Sections**:
1. **Overview** - Available models and providers
2. **Supported Solvents** - 25+ solvents listed
3. **Usage** - CLI examples for each model
4. **Method Availability** - Table of method×solvation combinations
5. **Solvation Models Explained** - Detailed algorithm descriptions
6. **Performance Considerations** - Speed and accuracy trade-offs

### 7.2 CLI Usage Examples

```bash
# Single point with water solvation (default GB)
curcuma -sp molecule.xyz -method gfn2:tblite -solvent water

# Geometry optimization in DMSO
curcuma -opt molecule.xyz -method gfn2:tblite -solvent dmso

# CPCM with custom dielectric
curcuma -sp molecule.xyz -method gfn2:tblite \
        -solvent_epsilon 78.4 -solvent_model 1

# ALPB (most accurate for TBLite)
curcuma -sp molecule.xyz -method gfn2:tblite -solvent water -solvent_model 3

# PM6 via Ulysses with GBSA
curcuma -sp molecule.xyz -method pm6:ulysses -solvent acetone
```

### 7.3 Verbosity Output

**Level 1+**: Final energy with solvation
```
[ENERGY] TBLite Energy: -5.123456 Eh
[PARAM] solvation_model: GB (Generalized Born)
[PARAM] solvent: water
```

**Level 2+**: Detailed solvation setup
```
[INFO] Applying solvation model
[PARAM] solvation_method: CPCM
[PARAM] dielectric_constant: 78.400
```

**Level 3**: Full debug output
```
[INFO] TBLite implicit solvation enabled
[PARAM] alpb_version: 12
[PARAM] alpb_reference: 1
```

---

## 8. ALGORITHM DETAILS FOR NATIVE EXTRACTION

### 8.1 CPCM Algorithm Steps

To implement natively, follow these steps:

```
1. Generate molecular surface
   - Use van der Waals radii for each atom
   - Create triangulated surface mesh
   
2. Compute surface area elements
   - Partition into small patches (e.g., GEPOL format)
   - Each patch: position, normal vector, area
   
3. Initialize electrostatic solver
   - Solve Laplace/Poisson equation
   - Boundary condition: conductor surface
   
4. Compute induced surface charge
   - From molecular density at cavity surface
   - Iterative refinement (IEFPCM)
   
5. Calculate solvation energy
   - E_solv = (1/2) * sum(σ_i * φ_i * A_i)
   
6. Compute gradients
   - ∂E_solv/∂R via surface element derivatives
```

### 8.2 Born Radius Parameters

For GB/GBSA native implementation:

```cpp
// Atom-specific Born radius calculation
double computeBornRadius(int Z, double vdw_radius) {
    // Radii from literature (Still et al., Hawkins et al.)
    static const double born_offsets[] = {
        0.0,    // Dummy
        1.2,    // H
        1.7,    // C
        1.55,   // N
        1.52,   // O
        // ... more elements
    };
    
    // Compute from volume integral
    // r_i = (3*V_i/(4π))^(1/3) - offset
    return vdw_radius - born_offsets[Z];
}
```

### 8.3 Dielectric Constants per Solvent

Common values for implementation:

```cpp
std::map<std::string, double> dielectric_constants = {
    {"water", 78.4},
    {"methanol", 32.7},
    {"ethanol", 24.5},
    {"acetone", 20.7},
    {"dmso", 46.7},
    {"dmf", 37.2},
    {"acetonitrile", 37.5},
    {"benzene", 2.3},
    {"hexane", 1.9},
    {"chloroform", 4.8},
    // ... more solvents
};
```

### 8.4 Surface Tension Parameters

For cavity energy term (GBSA):

```cpp
std::map<std::string, double> surface_tensions = {
    // γ in cal/(mol·Ų) - used in E_cavity = γ * SASA
    {"water", 0.0072},
    {"acetone", 0.0063},
    // ... 
};

// Cavity energy: E_cavity = γ * sum(area_i)
```

---

## 9. CURRENT LIMITATIONS AND WORKAROUNDS

### 9.1 Known Issues

| Issue | Status | Workaround |
|-------|--------|-----------|
| **GB Model Broken** | ❌ Does not work | Use ALPB (model 3) instead |
| **No custom solvent params** | ⚠️ Built-in only | TBLite limited to 25+ solvents |
| **GB version mismatch** | ❌ Unclear | ALPB versions are documented (12, 13) |

### 9.2 Best Practices

1. **For production calculations**: Use ALPB (`-solvent_model 3`)
2. **For custom dielectrics**: Use CPCM with explicit epsilon (`-solvent_model 1 -solvent_epsilon X`)
3. **For Ulysses**: Use GBSA (`-method pm6:ulysses -solvent water`)

### 9.3 Limitations of Current Implementation

- GB model disabled (C API issue)
- No analytical gradients (numerical gradients used)
- Solvent database fixed (can't add new solvents)
- ALPB parameters not customizable

---

## 10. FILES FOR NATIVE EXTRACTION

### 10.1 Reference Implementation Files

**Primary Source**:
- `/home/user/curcuma/src/core/energy_calculators/qm_methods/tbliteinterface.h` (126 lines)
- `/home/user/curcuma/src/core/energy_calculators/qm_methods/tbliteinterface.cpp` (full implementation)

**Documentation**:
- `/home/user/curcuma/docs/SOLVATION.md` (225 lines)
- `/home/user/curcuma/docs/NATIVE_QM_IMPLEMENTATION_STATUS.md` (572 lines)

**Helper/Testing**:
- `/home/user/curcuma/src/helpers/tblite_helper.cpp` (test program)
- `/home/user/curcuma/test_cases/examples/solvation_example.sh` (usage examples)

**Secondary Source** (Ulysses):
- `/home/user/curcuma/src/core/energy_calculators/qm_methods/ulyssesinterface.h/cpp`

### 10.2 Key Code Locations

```
TBLite Integration:
├── Parameter Definition: tbliteinterface.h:31-50
├── Solvation Setup: tbliteinterface.cpp:55-71
├── Application Logic: tbliteinterface.cpp:258-323
├── Error Handling: tbliteinterface.cpp:274-319
└── Integration: tbliteinterface.cpp:408

Documentation:
├── User Guide: docs/SOLVATION.md
├── Algorithm Details: docs/NATIVE_QM_IMPLEMENTATION_STATUS.md
├── Examples: test_cases/examples/solvation_example.sh
└── CLI Tests: test_cases/cli/*/run_test.sh
```

---

## 11. EXTRACTION ROADMAP FOR NATIVE IMPLEMENTATION

### Phase 1: Core Algorithms (Medium Effort)

**CPCM Implementation**:
```
1. Surface mesh generation (molecular surface)
2. Charge density evaluation
3. Poisson equation solver
4. Energy assembly
5. Gradient computation
```

**GB/GBSA Implementation**:
```
1. Born radius calculation
2. Pairwise interaction kernel
3. Surface area calculation (SASA)
4. Energy summation
5. Gradient via chain rule
```

### Phase 2: Parameter Integration (Low Effort)

- Extract dielectric constants (25+ solvents)
- Extract Born radius offsets
- Extract surface tension values
- Implement solvent parameter database

### Phase 3: Optimization (Medium Effort)

- Analytical gradients (vs numerical)
- Fast surface mesh generation
- Caching of Born radii
- Parallel SASA computation

### Phase 4: Validation (High Effort)

- Compare with TBLite results
- Test across diverse molecular types
- Benchmark performance
- Document accuracy limits

---

## 12. SUMMARY TABLE

| Component | Location | Status | Effort |
|-----------|----------|--------|--------|
| **CPCM** | TBLite C API | ✅ Functional | Medium native |
| **GB** | TBLite C API | ❌ Broken | High native |
| **ALPB** | TBLite C API | ✅ Functional | Medium native |
| **GBSA** | Ulysses | ✅ Functional | High native |
| **Parameters** | Compiled TBLite | ⚠️ Fixed | Low native |
| **Solvents** | 25+ built-in | ⚠️ Limited | Low native |
| **Gradients** | TBLite internal | ✅ Analytical | Medium native |

---

**End of Analysis**
