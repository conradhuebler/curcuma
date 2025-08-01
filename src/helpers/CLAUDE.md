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
- Native GFN-FF (cgfnff) method validation
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