# Archived notes: src/helpers/CLAUDE.md
Origin: `src/helpers/CLAUDE.md`, previous version moved here verbatim; short version: [src/helpers/CLAUDE.md](../../src/helpers/CLAUDE.md)
Archived: 2026-10-01. Statements are as old as their date and were not re-validated unless listed below.

## Claims found stale or wrong on 2026-10-01

- File tree: `gfnff_helper.cpp`, `gfnff_term_validator.cpp`, `gfnff_term_validator.py`, `gfnff_reference_validator.py`, `validate_ch3oh.py` were missing.
- "All interface helpers operational": `cli_test.cpp` is not built (target commented out in CMakeLists.txt, includes the non-existent `src/tools/cli_parser.h`); `gfnff_test.cpp` has no CMake target and includes the non-existent `src/core/energy_calculators/qm_methods/gfnff.h`; `imagewrite.cpp`, `storage_bench.cpp`, `polymer_topo.cpp` have no CMake target either.
- "main.cpp: Central dispatcher for all helper programs": `curcuma_helper` only knows `-statistic` and `-allxyz`; every other helper is its own executable.
- "storage_bench: Parameter caching and storage performance benchmarking" / "Storage benchmarking shows 96% speedup with caching": `storage_bench.cpp` compares switch/case and container layouts for element-parameter lookup; it does not exercise the force-field parameter cache, and no source for the 96% figure was found.
- "gfnff_test: Native GFN-FF implementation testing": see above, not buildable.

## Previous file (verbatim)

# CLAUDE.md - Helpers Directory

## Overview

The helpers directory contains standalone utilities, test programs, and development tools. These programs are used for testing, benchmarking, and standalone functionality that doesn't fit into the main application structure.

## Structure

```
helpers/
├── main.cpp              # Helper programs main entry point
├── cli_test.cpp          # Command-line interface testing
├── dftd3_helper.cpp      # DFT-D3 dispersion testing utility
├── dftd4_helper.cpp      # DFT-D4 dispersion testing utility
├── gfnff_test.cpp        # GFN-FF method testing
├── tblite_helper.cpp     # TBLite interface testing
├── ulysses_helper.cpp    # Ulysses interface testing  
├── xtb_helper.cpp        # XTB interface testing
├── parallel_scf.cpp      # Parallel SCF testing utility
├── imagewrite.cpp        # Image generation utilities
├── storage_bench.cpp     # Storage/caching benchmarking
└── polymer_topo.cpp      # Polymer topology utilities
```

## Key Programs

### Interface Testing
- **dftd3_helper**: Standalone DFT-D3 dispersion correction testing
- **dftd4_helper**: Standalone DFT-D4 dispersion correction testing
- **tblite_helper**: TBLite quantum chemistry interface testing
- **ulysses_helper**: Ulysses semi-empirical method testing
- **xtb_helper**: XTB tight-binding method testing
- **gfnff_test**: Native GFN-FF implementation testing

### Development Tools
- **cli_test**: Command-line interface and argument parsing testing
- **parallel_scf**: Multi-threading and parallel computation testing
- **storage_bench**: Parameter caching and storage performance benchmarking
- **imagewrite**: Visualization and image generation utilities

### Specialized Utilities
- **polymer_topo**: Polymer topology generation and analysis
- **main.cpp**: Central dispatcher for all helper programs

### Development Support
- Standalone testing of individual components
- Performance benchmarking and profiling
- Interface validation and debugging
- Prototype development and experimentation

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Testing priorities, benchmarking requirements, and development tool specifications to be defined by operator/programmer*

## Variable Section

### Current Testing Focus
- Native GFN-FF (gfnff) method validation
- Parameter caching system performance verification
- Multi-threading correctness testing

### Development Tools Status
- All interface helpers operational
- Storage benchmarking shows 96% speedup with caching
- Image generation utilities working for analysis visualization

### Testing Requirements
- Systematic validation of all quantum chemistry interfaces
- Performance regression testing for optimization changes
- Memory usage profiling for large molecular systems

---

*This documentation covers all development tools and standalone testing utilities*