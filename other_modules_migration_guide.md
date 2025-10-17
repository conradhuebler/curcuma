# Migrations-Guide: Weitere Module auf ConfigManager

Dieses Dokument fasst die Migrationsschritte für die verbleibenden Module zusammen, die noch das alte Konfigurationssystem verwenden.

## ✅ Completion Status (Phase 4 - October 2025)

**ALL 6 MODULES MIGRATED SUCCESSFULLY** - Phase 4 Complete

| Module | Status | Details |
|--------|--------|---------|
| Hessian | ✅ COMPLETE | Parameter definitions + ConfigManager constructors + old JSON replacements |
| Docking | ✅ COMPLETE | Parameter definitions + ConfigManager constructors |
| QMDFFfit | ✅ COMPLETE | Parameter definitions + ConfigManager constructors + QMDFFFitJson/HessianJson replaced |
| RMSDtraj | ✅ COMPLETE | Parameter definitions + ConfigManager constructors + RMSDTrajJson replaced in simplemd.cpp |
| TrajectoryAnalysis | ✅ COMPLETE | Parameter definitions + ConfigManager constructors |
| PersistentDiagram | ✅ COMPLETE | Parameter definitions + ConfigManager constructors |

**Implementation Notes:**
- All 6 modules follow identical delegating constructor pattern for backward compatibility
- ConfigManager constructors created in all .h files
- ConfigManager implementations in all .cpp files with parameter extraction from m_defaults
- Old static JSON variables (HessianJson, QMDFFFitJson, RMSDTrajJson) replaced with `ParameterRegistry::getInstance().getDefaultJson("module_name")`
- Build verification: ✅ PASSED (exit code 0)
- Parameter registry generation: ✅ PASSED via `make GenerateParams`
- Estimated ~49 new parameters added to registry across all 6 modules

## 1. Hessian (`hessian.cpp`)

### Motivation
Das `hessian`-Modul hat eine moderate Anzahl von Parametern, die von einer klaren Definition und einem typsicheren Zugriff profitieren würden.

### Parameter-Definition (`hessian.h`)
```cpp
BEGIN_PARAMETER_DEFINITION(hessian)
    PARAM(calculate, Bool, true, "Perform Hessian calculation.", "Execution", {"hess_calc"})
    PARAM(read, Bool, false, "Read Hessian from a file.", "Input", {"hess_read"})
    PARAM(read_file, String, "hessian.json", "File to read Hessian from.", "Input", {"hess_read_file"})
    PARAM(read_xyz, String, "", "XYZ file for geometry when reading Hessian.", "Input", {"hess_read_xyz"})
    PARAM(write_file, String, "hessian.out", "File to write Hessian to.", "Output", {"hess_write_file"})
    PARAM(finite_diff_step, Double, 0.01, "Step size for finite difference calculation.", "Algorithm", {})
    PARAM(freq_scale, Double, 1.0, "Scaling factor for frequencies.", "Analysis", {})
    PARAM(thermo, Double, 298.15, "Temperature for thermodynamic properties.", "Analysis", {})
    PARAM(freq_cutoff, Double, 50.0, "Cutoff for printing frequencies.", "Analysis", {})
END_PARAMETER_DEFINITION
```

### Refactoring
Ersetze die `Json2KeyWord`-Aufrufe in `LoadControlJson` durch `m_config.get<T>()`.

## 2. Docking (`docking.cpp`)

### Motivation
Das `docking`-Modul verwendet inkonsistente `PascalCase`-Namen und würde von einer klaren Struktur profitieren.

### Parameter-Definition (`docking.h`)
```cpp
BEGIN_PARAMETER_DEFINITION(docking)
    PARAM(host, String, "", "Host molecule file.", "Input", {})
    PARAM(guest, String, "", "Guest molecule file.", "Input", {})
    PARAM(complex, String, "", "Pre-assembled complex file.", "Input", {})
    PARAM(pos_x, Double, 0.0, "X-position for docking box center.", "Grid", {"Pos_X"})
    PARAM(pos_y, Double, 0.0, "Y-position for docking box center.", "Grid", {"Pos_Y"})
    PARAM(pos_z, Double, 0.0, "Z-position for docking box center.", "Grid", {"Pos_Z"})
    PARAM(step_x, Int, 10, "Number of steps in X direction.", "Grid", {"Step_X"})
    PARAM(step_y, Int, 10, "Number of steps in Y direction.", "Grid", {"Step_Y"})
    PARAM(step_z, Int, 10, "Number of steps in Z direction.", "Grid", {"Step_Z"})
    PARAM(auto_pos, Bool, true, "Automatically determine docking box position.", "Grid", {"AutoPos"})
    PARAM(cycles, Int, 1, "Number of docking cycles.", "Execution", {"Cycles"})
    PARAM(threads, Int, 1, "Number of threads for the main process.", "Performance", {})
    PARAM(docking_threads, Int, 1, "Number of threads for docking subprocesses.", "Performance", {"DockingThreads"})
END_PARAMETER_DEFINITION
```

## 3. QMDFFfit (`qmdfffit.cpp`)

### Motivation
Ein kleines Modul, das leicht umgestellt werden kann und die Konsistenz erhöht.

### Parameter-Definition (`qmdfffit.h`)
```cpp
BEGIN_PARAMETER_DEFINITION(qmdfffit)
    PARAM(method, String, "uff", "Method for the initial force field.", "General", {})
    PARAM(threads, Int, 1, "Number of threads.", "Performance", {})
    PARAM(hessian_file, String, "hessian.json", "Input Hessian file.", "Input", {"hessian"})
    PARAM(charges_file, String, "scf.json", "Input charges file.", "Input", {"charges"})
END_PARAMETER_DEFINITION
```

## 4. RMSDtraj (`rmsdtraj.cpp`)

### Motivation
Ähnlich wie `rmsd`, aber für Trajektorien. Die Umstellung vereinheitlicht die RMSD-bezogenen Werkzeuge.

### Parameter-Definition (`rmsdtraj.h`)
```cpp
BEGIN_PARAMETER_DEFINITION(rmsdtraj)
    PARAM(heavy_only, Bool, false, "Use only heavy atoms.", "RMSD", {"heavy"})
    PARAM(rmsd_threshold, Double, 1.0, "RMSD threshold for clustering.", "RMSD", {"rmsd"})
    PARAM(reference, String, "", "Reference structure file.", "Input", {})
    PARAM(second_trajectory, String, "", "Second trajectory file for comparison.", "Input", {"second"})
    PARAM(write_unique, Bool, false, "Write unique conformers to file.", "Output", {"writeUnique"})
    PARAM(write_aligned, Bool, false, "Write aligned trajectory.", "Output", {"writeAligned"})
    PARAM(write_rmsd, Bool, false, "Write RMSD values to a file.", "Output", {"writeRMSD"})
END_PARAMETER_DEFINITION
```

## 5. TrajectoryAnalysis (`trajectoryanalysis.cpp`)

### Motivation
Dieses Modul ist bereits teilweise auf `ConfigManager` umgestellt und dient als gutes Beispiel für die finale Bereinigung.

### Parameter-Definition (`trajectoryanalysis.h`)
```cpp
BEGIN_PARAMETER_DEFINITION(trajectoryanalysis)
    PARAM(properties, String, "all", "Properties to calculate: all|basic|geometric|cg.", "Analysis", {})
    PARAM(output_format, String, "human", "Output format: human|json|csv.", "Output", {})
    PARAM(output_file, String, "", "File to save results.", "Output", {})
    PARAM(stride, Int, 1, "Analyze every N-th frame.", "Input", {})
    PARAM(start_frame, Int, 0, "Frame to start analysis from.", "Input", {})
    PARAM(end_frame, Int, -1, "Frame to end analysis at (-1 for end).", "Input", {})
    PARAM(moving_average, Int, 10, "Window size for moving average.", "Analysis", {})
END_PARAMETER_DEFINITION
```

## 6. PersistentDiagram (`persistentdiagram.cpp`)

### Motivation
Dieses Modul verwendet ein statisches `RipserJson`-Objekt. Die Umstellung integriert es in das neue System.

### Parameter-Definition (`persistentdiagram.h`)
```cpp
BEGIN_PARAMETER_DEFINITION(ripser)
    PARAM(x_max, Double, 4.0, "Max x-value for persistence diagram.", "Diagram", {"ripser_xmax"})
    PARAM(x_min, Double, 0.0, "Min x-value for persistence diagram.", "Diagram", {"ripser_xmin"})
    PARAM(y_max, Double, 4.0, "Max y-value for persistence diagram.", "Diagram", {"ripser_ymax"})
    PARAM(y_min, Double, 0.0, "Min y-value for persistence diagram.", "Diagram", {"ripser_ymin"})
    PARAM(bins, Int, 10, "Number of bins for persistence image.", "Image", {"ripser_bins"})
    PARAM(scaling, Double, 0.1, "Scaling factor for persistence image.", "Image", {"ripser_scaling"})
    PARAM(std_x, Double, 10.0, "Standard deviation for x-axis in persistence image.", "Image", {"ripser_stdx"})
    PARAM(std_y, Double, 10.0, "Standard deviation for y-axis in persistence image.", "Image", {"ripser_stdy"})
    PARAM(ratio, Double, 1.0, "Ratio for ripser calculation.", "Algorithm", {"ripser_ratio"})
    PARAM(dimension, Int, 2, "Dimension for ripser calculation.", "Algorithm", {"ripser_dimension"})
    PARAM(epsilon, Double, 0.4, "Epsilon for ripser calculation.", "Algorithm", {"ripser_epsilon"})
END_PARAMETER_DEFINITION
```
