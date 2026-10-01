# Archived notes: src/CLAUDE.md
Origin: `src/CLAUDE.md`, previous version moved here verbatim (relative links adjusted); short version: [src/CLAUDE.md](../../src/CLAUDE.md)
Archived: 2026-10-01. Statements are as old as their date and were not re-validated unless listed below.

## Claims found stale or wrong on 2026-10-01

- `#ifdef DEBUG_ON` debugging rule and "Check CMakeLists.txt and includes for DEBUG_ON configuration": `DEBUG_ON` occurs in no source file under `src/` and in no `CMakeLists.txt` (grep).
- "Replace deprecated suprafit functions": this is curcuma, not SupraFit; the kept rule follows the root CLAUDE.md wording (deprecated function calls in general).
- Structure tree: `pch.h`, `pch/`, `version.h.in` and `molecules/` were missing.
- "Universal parameter caching system working (96% speedup ...)": no measurement or test backing the figure was found; not carried over.

## Previous file (verbatim)

# CLAUDE.md - Curcuma Source Directory

## Overview

This directory contains the complete Curcuma source code organized into four main subsystems:
- **capabilities/** - High-level molecular modeling tasks and applications
- **core/** - Core computational engines and data structures
- **tools/** - Utilities and file I/O operations
- **helpers/** - Development tools and standalone utilities

## Structure

```
src/
├── capabilities/          # Application-level functionality
├── core/                 # Core computational engines  
├── tools/                # Utilities and I/O
├── helpers/              # Development tools
├── main.cpp              # Main application entry point
└── global_config.h.in    # Global configuration template
```

## Development Standards

### Code Quality
- All new functions marked as "Claude Generated" for traceability
- Doxygen-ready documentation for new and frequently used functions
- Remove completed TODO hashtags when approved
- Use std::cout for debugging within `#ifdef DEBUG_ON #endif`
- Port std::cout to fmt for non-debugging console output
- Replace deprecated suprafit functions when compiler warnings appear

### Error Handling
- Implement comprehensive error handling and logging
- Maintain backward compatibility where possible
- Check CMakeLists.txt and includes for DEBUG_ON configuration

### Architecture Principles
- Each subsystem maintains its own CLAUDE.md with specific information
- Variable sections contain short-term information and current bugs
- Preserved sections contain permanent knowledge and patterns
- Instructions blocks contain future tasks and development visions

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*This section is reserved for operator/programmer instructions and approved future development tasks*

## Variable Section

### Current Development Focus
- GFN-FF implementation complete and operational (`energy_calculators/ff_methods/`) — See [docs/GFNFF_STATUS.md](../GFNFF_STATUS.md)
- Universal parameter caching system working (96% speedup for iterative calculations)
- Enhanced force field multi-threading support

---

*This file provides overview information relevant to ALL subdirectories in src/*