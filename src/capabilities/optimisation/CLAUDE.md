# CLAUDE.md - Optimization Algorithms Directory

## Overview

This directory contains the **core optimization algorithms** for molecular geometry optimization. Focus is on **educational clarity** and **scientific transparency** rather than software engineering abstractions.

## 🧪 **SCIENTIFIC CONTENT** (Core Algorithms)

### Native LBFGS Implementation (`lbfgs.h/.cpp`)
**Mathematical algorithms - the actual computational chemistry:**

#### L-BFGS Algorithm (`LBFGSStep()` lines 263-320)
- **Two-loop recursion** for limited-memory BFGS
- **Scientific core**: Quasi-Newton method using gradient history
- **Literature**: Nocedal & Wright "Numerical Optimization", Chapter 7
- **Key equations**: H₀⁻¹ approximation, α/β coefficients, search direction computation

#### DIIS Algorithm (`DIISStep()` lines 245-261) 
- **Direct Inversion of Iterative Subspace** - acceleration method
- **Scientific core**: Linear combination of previous solutions to minimize residual
- **Literature**: Pulay, P. (1980) Chem. Phys. Lett. 73, 393
- **Key equations**: B-matrix construction, SVD-based coefficient solution

#### RFO Algorithm (`RFOStep()` lines 347-391)
- **Rational Function Optimization** - eigenvector following
- **Scientific core**: Newton-Raphson in eigenvector space of Hessian
- **Literature**: Banerjee et al. (1985) J. Phys. Chem. 89, 52
- **Key equations**: Eigenvalue shifting, mass-weighted coordinates, τ coefficients

#### SR1 Hessian Update (`sr1Update()` lines 390-400)
- **Symmetric Rank-1 update** formula
- **Scientific core**: (s·(y-Hs))/(s·(y-Hs)) * outer_product(y-Hs, y-Hs)
- **Literature**: Broyden (1970), Dennis & Moré (1977)

### Line Search Methods
- **Backtracking line search** with Armijo condition
- **RFO-specific line search** with curvature conditions
- **Step size control** and convergence monitoring

## 🏗️ **ORGANIZATIONAL CODE** (Software Engineering)

### Modern Optimizer System (`modern_optimizer_simple.h/.cpp`)
**Pure software engineering - no computational chemistry:**

#### Strategy Pattern Dispatcher
- **OptimizerType enum**: Type-safe method selection (no magic numbers)
- **parseOptimizerType()**: String-to-enum conversion with error handling
- **Auto-selection logic**: Size-based algorithm choice (>200 atoms → LBFGSpp, etc.)

#### Integration Layer  
- **optimizeWithLBFGSpp()**: Wrapper for external LBFGSpp library
- **optimizeWithInternal()**: Wrapper for internal GPTLBFGS
- **Legacy compatibility**: Bridge to existing CurcumaOpt interface

#### Logging & Results
- **CurcumaLogger integration**: Verbosity-controlled output with colors
- **SimpleOptimizationResult**: Standardized result structure
- **Unit-aware output**: Energy in kJ/mol, distances in Å, time in seconds

## File Organization

```
optimisation/
├── lbfgs.h/.cpp           # 🧪 SCIENCE: Core optimization algorithms (L-BFGS, DIIS, RFO)
├── modern_optimizer_simple.h/.cpp  # 🏗️ ORGANIZATION: Strategy pattern dispatcher
└── CLAUDE.md              # This documentation
```

## Key Principles

### Educational-First Design
- **Algorithm visibility**: Core mathematical methods directly accessible
- **Minimal abstraction**: No unnecessary templates or inheritance hiding the science
- **Clear separation**: Scientific algorithms vs. software organization clearly marked
- **Literature references**: Papers and equations documented in code comments

### Scientific Transparency
- **Direct implementation**: L-BFGS two-loop recursion visible and understandable
- **Mathematical naming**: Variables reflect scientific meaning (gradient, hessian, eigenvalues)
- **Parameter meaning**: Physical/mathematical significance clear from context
- **Debugging access**: Easy inspection of intermediate algorithmic steps

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Future tasks and visions defined by operator/programmer*

## 🚀 **EINHEITLICHES OPTIMIZER-SYSTEM** (Claude Generated)

### **Kommandozeilen-Interface:**
```bash
# 🧪 CURCUMA'S EIGENE OPTIMIERER:
./curcuma -opt input.xyz -optimizer native_lbfgs   # Native L-BFGS two-loop recursion
./curcuma -opt input.xyz -optimizer diis           # Native DIIS acceleration  
./curcuma -opt input.xyz -optimizer rfo            # Native RFO eigenvector following
./curcuma -opt input.xyz -optimizer auto           # Auto-selection (>50 Atome → L-BFGS)

# 🏗️ EXTERNE BIBLIOTHEKEN:
./curcuma -opt input.xyz -optimizer lbfgspp        # External LBFGSpp library
./curcuma -opt input.xyz -optimizer internal       # Legacy internal LBFGS

# Standard (ohne -optimizer): Fallback zu legacy system
./curcuma -opt input.xyz                           # → Legacy CurcumaOpt
```

### **Einheitliche Systematik - FESTGEHALTEN:**
- **`-method`**: **WAS** berechnet wird (Potential: uff, gfnff, eht, etc.) 
- **`-optimizer`**: **WIE** optimiert wird (native_lbfgs, diis, rfo, etc.)
- **Einheitlich für alle Module**: ConfScan, RMSD, etc. können das gleiche System nutzen

### **Educational Output:**
```
🧪 Using native Curcuma optimizer: diis
🧪 SCIENCE: Running native DIIS acceleration method  
Algorithm: Direct Inversion of Iterative Subspace
Citation: Pulay, P. (1980) Chem. Phys. Lett. 73, 393
DIIS history: 10
DIIS start iteration: 5

Iteration 5 (L-BFGS): E = -154.123456 Eh, |grad| = 1.23e-04
Iteration 10 (DIIS): E = -154.124567 Eh, |grad| = 5.67e-05
🧪 Native DIIS converged!

🎉 Optimization successful!
Method: native_diis
Iterations: 12  
Time: 2.345 seconds
Optimized structure saved to: input.opt.xyz
```

## Variable Section

### ✅ **COMPLETED MILESTONES:**
- **Native optimizer cleanup**: All variable naming inconsistencies fixed
- **Einheitliches -optimizer System**: Implementiert und funktionsfähig
- **Educational-First Integration**: Wissenschaftliche Algorithmen direkt aufrufbar
- **Curcuma-eigene Optimierer**: L-BFGS, DIIS, RFO nativ verfügbar
- **Build erfolgreich**: Kompiliert ohne Fehler
- **🧪 ALGORITHMUS KORREKTUREN** (August 2025):
  - **L-BFGS Two-Loop Recursion**: Indexierung korrigiert (Nocedal & Wright konform)
  - **Curvature Condition**: Implementiert - verhindert numerische Instabilität
  - **DIIS B-Matrix**: Korrigiert nach Pulay (1980) - einfache Skalarprodukte
  - **RFO Algorithmus**: Komplett neu nach Banerjee et al. (1985) implementiert
  - **SR1 Hessian Update**: Korrigiert - m_lambda Faktor entfernt

### 🎯 **EINHEITLICHES DESIGN FESTGEHALTEN:**
1. **`-optimizer` Parameter**: Bestimmt Optimierungsalgorithmus (nicht `-method`)
2. **Native First**: Curcuma's eigene Implementierungen bevorzugt  
3. **Educational Output**: Literaturzitate und wissenschaftliche Parameter sichtbar
4. **Automatic Fallback**: Legacy system für unbekannte Optimizer
5. **Erweiterbar**: System kann einfach um weitere Optimizer erweitert werden

---

**Educational Focus**: This directory demonstrates computational chemistry optimization algorithms with minimal software engineering overhead, prioritizing learning and understanding over complex abstractions.