# CLAUDE.md - src/tools/

Utilities used across curcuma. Mostly header-only `inline` functions; `bmt_utils.cpp` and
`trajectory_writer.cpp` are compiled. Element data and unit constants are NOT here (`src/core/elements.h`,
`periodic_table.h`, `units.h`).

## Files

| File | Content |
|------|---------|
| `formats.h` (`Files::`) | Readers XYZ/TRJ, SDF, MOL2, VTF, Turbomole `coord`; VTF writers; `LoadFile()` / `LoadMol()` dispatch |
| `geometry.h` (`GeometryTools::`) | `Distance`, `Centroid`, rotation matrices, translate/rotate geometry |
| `general.h` (`Tools::`) | String parsing/conversion, `mean`/`median`/`stdev`/`Histogram`, `CreateList()` range parser; global `RunTimer` |
| `info.h` (`General::`) | `StartUp()` banner (version, git hash) |
| `pbc_utils.h` (`PBCUtils::`) | Lattice vectors from a,b,c,alpha,beta,gamma and back, minimum image, PBC distance, cell check |
| `spatial_cell_list.h` | Cell list for O(N) neighbour queries, no PBC; build cost does not pay off below roughly 800 atoms |
| `string_similarity.h` (`StringUtils::`) | Levenshtein distance for "Did you mean" method suggestions |
| `trajectory_writer.{h,cpp}`, `trajectory_helpers.h` | TrajectoryWriter and its JSON schema helpers |
| `bmt_utils.{h,cpp}` (`BMTUtils::`) | BMT output directories |

## ✅ TrajectoryWriter

- Formats `HumanTable`, `CSV`, `JSON`, `DAT`, `VTF`; single or multi-frame, plus statistics summaries (`TrajectoryStatistics`)
- `trajectory_helpers.h` converts geometry-command results (bond, angle, torsion) into its JSON schema

## 🤖 BMT output directories (AI-generated, machine-tested; human production testing pending)

- Default: every command writes into `Basename.Keyword.YYYYMMDD_HHMMSS/`; `-no_bmt` writes to CWD instead
- `createBMTDir`, `writeMetadata` (`metadata.json`), `outputPath`, `processBakFiles`, `collectBakFiles`, `stripExtension`
- `-bak f1 f2` copies listed files back to CWD; `bak` and `no_bmt` are global flags (`main.cpp`)
- Output files must go through `BMTUtils::outputPath()` or, in `CurcumaMethod` subclasses, `CurcumaMethod::outputPath()` (root CLAUDE.md, mandatory)

## Traps

- `LoadFile()` picks the reader by substring (`.xyz`, `.mol2`, ..., `coord`), first match wins; only `LoadFile()` reads `.json`
- `createBMTDir()` creates the directory only in a `C17` non-Windows build; a second run in the same second gets suffix `_2`, `_3`, ...
- `-bmt false` is not read anywhere; only `-no_bmt` disables BMT
- `Tools::CreateList("a:b,c")` silently drops malformed tokens and returns nothing for a descending range (`1:-1`); callers resolve such forms first (see `analysis.cpp`)

## Instructions Block

**PRESERVED - DO NOT EDIT BY CLAUDE**

*Tool development priorities and utility function requirements to be defined by operator/programmer*

---

Previous version (removed sections and status notes, 2026-10-01): [docs/archive/TOOLS_NOTES_2026-10.md](../../docs/archive/TOOLS_NOTES_2026-10.md)
