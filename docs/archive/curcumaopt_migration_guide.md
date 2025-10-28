# Migrations-Guide: Umstellung von CurcumaOpt auf ConfigManager

Dieses Dokument beschreibt die schrittweise Migration des `curcumaopt`-Moduls auf das neue, zentralisierte `ConfigManager`-System.

## Motivation

Das `curcumaopt`-Modul ist für alle Optimierungs- und Single-Point-Aufgaben zentral. Es leidet unter den typischen Problemen des alten Konfigurationssystems:
- **Inkonsistente Namensgebung**: Ein Mix aus `PascalCase` (`SinglePoint`), `camelCase` (`printOutput`) und `UPPER_SNAKE_CASE` (`LBFGS_m`).
- **Flache Struktur für hierarchische Daten**: Parameter für den LBFGS-Algorithmus sind mit einem Präfix versehen (`LBFGS_m`, `LBFGS_eps_abs`), anstatt in einer natürlichen Hierarchie (`lbfgs.m`) zu leben.
- **Hohe Parameteranzahl**: Über 25 `Json2KeyWord`-Aufrufe machen die `LoadControlJson`-Methode unübersichtlich.

Die Umstellung auf den `ConfigManager` wird diese Probleme lösen und eine saubere, hierarchische Konfiguration ermöglichen.

## Phase 1: Parameter-Definition in `curcumaopt.h`

Alle für `curcumaopt` relevanten Parameter werden in einem `PARAMETER_DEFINITION_BLOCK` in der Header-Datei `src/capabilities/curcumaopt.h` deklariert.

**Anleitung:** Füge den folgenden Block in die `CurcumaOpt`-Klasse in `curcumaopt.h` ein.

```cpp
// In src/capabilities/curcumaopt.h, innerhalb der class CurcumaOpt

private:
    // vvvvvvvvvvvv PARAMETER DEFINITION BLOCK vvvvvvvvvvvv
    BEGIN_PARAMETER_DEFINITION(opt)

    // --- General Settings ---
    PARAM(method, String, "uff", "Energy calculation method (e.g., uff, gfn2).", "General", {})
    PARAM(threads, Int, 1, "Number of parallel threads.", "Performance", {})
    PARAM(charge, Double, 0.0, "Total charge of the system.", "General", {"Charge"})
    PARAM(spin, Double, 1.0, "Total spin multiplicity.", "General", {"Spin"})
    PARAM(single_point, Bool, false, "Perform a single point calculation instead of optimization.", "General", {"SinglePoint"})

    // --- Convergence Criteria ---
    PARAM(max_iterations, Int, 500, "Maximum number of optimization iterations.", "Convergence", {"maxiter"})
    PARAM(energy_change, Double, 1e-6, "Convergence criterion for energy change (Eh).", "Convergence", {"dE"})
    PARAM(rmsd_change, Double, 1e-3, "Convergence criterion for RMSD change (Angstrom).", "Convergence", {"dRMSD"})
    PARAM(gradient_norm, Double, 1e-4, "Convergence criterion for the gradient norm.", "Convergence", {"GradNorm"})
    PARAM(convergence_count, Int, 3, "Number of consecutive steps criteria must be met.", "Convergence", {"ConvCount"})

    // --- Optimizer Algorithm ---
    PARAM(optimizer, Int, 0, "Optimizer algorithm (0: LBFGS++, 1: Internal LBFGS).", "Algorithm", {"optimethod"})
    PARAM(optimize_h, Bool, false, "Constrain heavy atoms and optimize only hydrogens.", "Algorithm", {"optH"})

    // --- LBFGS++ Parameters (Nested) ---
    PARAM(lbfgs.m, Int, 6, "Number of corrections to store for LBFGS.", "LBFGS", {"LBFGS_m"})
    PARAM(lbfgs.epsilon_abs, Double, 1e-5, "Absolute convergence tolerance for LBFGS.", "LBFGS", {"LBFGS_eps_abs"})
    PARAM(lbfgs.epsilon_rel, Double, 1e-5, "Relative convergence tolerance for LBFGS.", "LBFGS", {"LBFGS_eps_rel"})
    PARAM(lbfgs.past, Int, 0, "Distance for delta-based convergence test.", "LBFGS", {})
    PARAM(lbfgs.delta, Double, 1e-5, "Delta for convergence test.", "LBFGS", {})
    PARAM(lbfgs.line_search, Int, 0, "Line search algorithm.", "LBFGS", {"LBFGS_LST"})
    PARAM(lbfgs.max_linesearch, Int, 20, "Max line search iterations.", "LBFGS", {"LBFGS_ls_iter"})

    // --- Internal Optimizer Parameters ---
    PARAM(diis_history, Int, 5, "History size for DIIS.", "Internal-Opt", {"diis_hist"})
    PARAM(diis_start, Int, 2, "Iteration to start DIIS.", "Internal-Opt", {"diis_start"})
    PARAM(rfo_lambda, Double, 0.1, "Lambda for RFO algorithm.", "Internal-Opt", {"lambda"})

    END_PARAMETER_DEFINITION
    // ^^^^^^^^^^^^ PARAMETER DEFINITION BLOCK ^^^^^^^^^^^^
```

## Phase 2: Umstellung der Implementierung in `curcumaopt.cpp`

### 2.1. `ConfigManager` als Member-Variable hinzufügen

Füge in `curcumaopt.h` eine `ConfigManager`-Instanz als private Member-Variable hinzu:

```cpp
#include "src/core/config_manager.h"

class CurcumaOpt : public CurcumaMethod {
    // ...
private:
    ConfigManager m_config;
    // ... andere Member
};
```

### 2.2. Konstruktor und `LoadControlJson` refaktorisieren

**Konstruktor (Nachher):**
```cpp
CurcumaOpt::CurcumaOpt(const json& controller, bool silent)
    : CurcumaMethod(ParameterRegistry::getInstance().getDefaultJson("opt"), controller, silent),
      m_config("opt", controller) // ConfigManager initialisieren
{
    UpdateController(controller); 
}
```

**`LoadControlJson` (Nachher):**

Die Methode wird durch saubere `m_config.get<T>()`-Aufrufe ersetzt. Die hierarchischen LBFGS-Parameter werden elegant über die Punktnotation angesprochen.

```cpp
void CurcumaOpt::LoadControlJson()
{
    // Einfacher, typsicherer Zugriff
    m_threads = m_config.get<int>("threads");
    m_method = m_config.get<std::string>("method");
    m_singlepoint = m_config.get<bool>("single_point");
    m_maxiter = m_config.get<int>("max_iterations");
    m_optimethod = m_config.get<int>("optimizer");

    // ... alle anderen Parameter entsprechend ersetzen ...
}
```

**Verwendung der hierarchischen Parameter:**

In der `LBFGSOptimise`-Methode wird der Zugriff auf die LBFGS-Parameter wie folgt umgestellt:

**Vorher:**
```cpp
LBFGSParam<double> param;
param.m = Json2KeyWord<int>(m_defaults, "LBFGS_m");
param.epsilon = Json2KeyWord<double>(m_defaults, "LBFGS_eps_abs");
// ...
```

**Nachher:**
```cpp
LBFGSParam<double> param;
// Sauberer Zugriff über die Punktnotation
param.m = m_config.get<int>("lbfgs.m");
param.epsilon = m_config.get<double>("lbfgs.epsilon_abs");
// ...
```

### 2.3. Bereinigung

-   Nach der vollständigen Umstellung können das globale `CurcumaOptJson`-Objekt und das `m_defaults`-Member entfernt werden.
-   Die `printHelp()`-Funktion kann durch die automatisch generierte Version ersetzt werden.

## Zusammenfassung der Vorteile für `curcumaopt`

-   **Strukturierte Parameter**: Die LBFGS-Parameter sind nun logisch unter `lbfgs` gruppiert, was die Lesbarkeit der Konfiguration und des Codes verbessert.
-   **Konsistente Namensgebung**: Alle Parameter folgen nun dem `snake_case`-Standard.
-   **Reduzierter Boilerplate**: Die über 25 `Json2KeyWord`-Aufrufe werden durch eine saubere und sichere API ersetzt.
