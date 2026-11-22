# CLAUDE.md - src/tools/

Utilities used across curcuma. Mostly header-only `inline` functions; `bmt_utils.cpp`,
`trajectory_writer.cpp` and `pbc_utils.cpp` are compiled. Element data and unit constants are NOT here (`src/core/elements.h`,
`periodic_table.h`, `units.h`).

## GeometryTools (`geometry.h`)

Header-only helpers over `Position`/`Geometry`. Besides the translation/rotation
utilities: `Distance`, `Centroid`, and since Sep 2026 `Angle` (degrees, clamped acos),
`Dihedral` (signed degrees, IUPAC) and `GyrationRadius(const Geometry&)` (unweighted —
the mass-weighted, PBC-aware and per-fragment variants live on `Molecule::GyrationRadius`,
which has the masses). Shared with the qurcuma GUI so the same numbers appear in its
measurement HUD and its live plots. Pinned by the `geometry_tools` ctest.


## Files

| File | Content |
|------|---------|
| `formats.h` (`Files::`) | Readers XYZ/TRJ, SDF, MOL2, VTF, Turbomole `coord`; VTF writers; `LoadFile()` / `LoadMol()` dispatch |
| `geometry.h` (`GeometryTools::`) | `Distance`, `Centroid`, rotation matrices, translate/rotate geometry |
| `general.h` (`Tools::`) | String parsing/conversion, `mean`/`median`/`stdev`/`Histogram`, `CreateList()` range parser; global `RunTimer` |
| `info.h` (`General::`) | `StartUp()` banner (version, git hash) |
| `pbc_utils.h/.cpp` (`PBCUtils::`) | Lattice vectors from a,b,c,alpha,beta,gamma and back, minimum image, PBC distance, cell check; orthorhombic `OrthoCell` MIC (GPU `rint` arithmetic), cutoff-vs-L/2 check, fragment-wise wrap (Nov 2026) |
| `pbc_gather.h` (`PBCGather::`) | Gather form of symmetric pair loops with serial-order replay of energies/virials: thread-count independent, serial bits (Nov 2026) |
| `spatial_cell_list.h` | Cell list for O(N) neighbour queries; `build(coords, cutoff, OrthoCell)` wraps around an orthorhombic cell (minimum-image distances, cutoff < L/2, Nov 2026); build cost does not pay off below roughly 800 atoms |
| `string_similarity.h` (`StringUtils::`) | Levenshtein distance for "Did you mean" method suggestions |
| `trajectory_writer.{h,cpp}` | TrajectoryWriter |
| `bmt_utils.{h,cpp}` (`BMTUtils::`) | BMT output directories |

## ✅ TrajectoryWriter

- Formats `HumanTable`, `CSV`, `JSON`, `DAT`, `VTF`; single or multi-frame, plus statistics summaries (`TrajectoryStatistics`)

## 🤖 BMT output directories (AI-generated, machine-tested; human production testing pending)

- Default: every command writes into `Basename.Keyword.YYYYMMDD_HHMMSS/`; `-no_bmt` writes to CWD instead
- `createBMTDir`, `writeMetadata` (`metadata.json`), `outputPath`, `processOutputDir()`/`processOutputPath()`, `processBakFiles`, `collectBakFiles`, `stripExtension`
- `-bak f1 f2` copies listed files back to CWD; `bak` and `no_bmt` are global flags (`main.cpp`)
- Output files must go through `BMTUtils::outputPath()` or, in `CurcumaMethod` subclasses, `CurcumaMethod::outputPath()` (root CLAUDE.md, mandatory)

## Traps

- `LoadFile()` picks the reader by substring (`.xyz`, `.mol2`, ..., `coord`), first match wins; only `LoadFile()` reads `.json`
- `createBMTDir()` creates the directory only in a `C17` non-Windows build; a second run in the same second gets suffix `_2`, `_3`, ...
- Code without a capability object (GFN-FF/D4/EEQ diagnostics) uses `BMTUtils::processOutputPath()`: the first BMT directory of the process, plain CWD with `-no_bmt` or when the command makes none (single-structure `-sp`)
- `-bmt false` is not read anywhere; only `-no_bmt` disables BMT
- `Tools::CreateList("a:b,c")` silently drops malformed tokens and returns nothing for a descending range (`1:-1`); callers resolve such forms first (see `analysis.cpp`)

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Tool development priorities and utility function requirements to be defined by operator/programmer*

---

Previous version (removed sections and status notes, 2026-10-01): [docs/archive/TOOLS_NOTES_2026-10.md](../../docs/archive/TOOLS_NOTES_2026-10.md)
