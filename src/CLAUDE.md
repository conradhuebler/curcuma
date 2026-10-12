# CLAUDE.md - src/

Source tree of curcuma. Project-wide rules are in the root [CLAUDE.md](../CLAUDE.md); every
subdirectory listed below has its own CLAUDE.md.

## Layout

- `capabilities/` - user-facing commands (`-opt`, `-md`, `-confsearch`, ...); many derive from `CurcumaMethod`
- `core/` - Molecule, EnergyCalculator, ParameterRegistry/ConfigManager, logger, units; methods in `core/energy_calculators/`
- `tools/` - file-format readers, geometry helpers, TrajectoryWriter, BMT output directories
- `helpers/` - standalone helper programs, most of them optional CMake targets
- `main.cpp` - CLI entry point (`CLI2Json`, command dispatch)
- `global_config.h.in`, `version.h.in` - `configure_file` templates (CMakeLists.txt)
- `pch.h`, `pch/` - precompiled headers (`USE_PCH`)
- `molecules/` - monomer XYZ files (PDMAEMA, PEO, PPO); no code in `src/` references this directory

## Conventions for all of src/

- Mark new functions "Claude Generated"; doxygen-ready comments for new and frequently used functions
- Console output goes through `CurcumaLogger` or fmt, not `std::cout`; plain ASCII only (root CLAUDE.md)
- Diagnostic dumps: verbosity level 3 or an env switch `CURCUMA_*`; never select a numerical path by the print level
- Replace deprecated function calls when the compiler warns about them
- Error handling and logging for every new code path; keep backward compatibility where possible
- Remove a TODO marker once the work is done and approved

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*This section is reserved for operator/programmer instructions and approved future development tasks*

---

Previous version (status notes, removed 2026-10-01): [docs/archive/SRC_NOTES_2026-10.md](../docs/archive/SRC_NOTES_2026-10.md)
