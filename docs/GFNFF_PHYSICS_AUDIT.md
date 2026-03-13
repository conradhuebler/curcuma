# GFN-FF Physics and Parameter Audit

## Overview
This document tracks in-depth investigations into the physical parameters and energy formulations of the GFN-FF implementation in Curcuma compared to the original XTB reference.

## Investigated Issues

### 1. Alpha Parameter in Coulomb Interaction (Jan 2026) ✅
- **Issue**: Discrepancy in whether to use base alpha ($\alpha_{base}$) or charge-corrected alpha ($\alpha_{eequiv}$) for the Coulomb damping.
- **Finding**: Using the charge-corrected $\alpha_{eequiv}$ in the Coulomb energy term actually *increased* the error relative to XTB for some systems.
- **Resolution**: Curcuma currently maintains the base alpha for energy evaluation as it provides better compensation for our two-phase charge calculation.

### 2. C-H Bond Distance (R0) Deviation ✅
- **Issue**: Systematic overestimation of C-H bond stretching energy.
- **Analysis**: GFN-FF uses coordinate-dependent equilibrium distances $r_0(CN)$. A mismatch in coordination number (CN) scaling led to shifted $r_0$ values.
- **Fix**: Re-synchronized the `fat` scaling factors and CN exponent with the Fortran reference.

### 3. EEQ Self-Energy Term ✅
- **Issue**: The atomic self-energy contribution $E_{self} = 0.5 \cdot q_i^2 \cdot (\eta_i + J_{ii})$ was inconsistently calculated.
- **Verification**: Confirmed the factor of $0.5$ and the inclusion of the Gaussian width term matching Fortran `gfnff_engrad.F90`.

## Technical Findings
- **Energy Unit**: Always ensure parameters are converted to Hartree when used in the `ForceFieldThread` kernels.
- **Damping Consistency**: Most discrepancies in 1,3 and 1,4 interactions arise from the exact formulation of the damping prefactors ($f_{qq}$, $f_{at}$).
