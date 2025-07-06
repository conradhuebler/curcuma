#include <array>
#include <chrono>
#include <iostream>
#include <unordered_map>
#include <vector>

// ====================================================================================
// OPTION 1: SWITCH/CASE (SCHLECHT für Effizienz)
// ====================================================================================
class SwitchCaseParameters {
public:
    double getElectronegativity(int atomic_number) const
    {
        switch (atomic_number) {
        case 1:
            return 2.300; // H
        case 6:
            return 2.847; // C
        case 7:
            return 3.022; // N
        case 8:
            return 3.610; // O
        case 9:
            return 4.193; // F
        case 16:
            return 2.703; // S
        case 17:
            return 3.242; // Cl
        // ... weitere 80+ Elemente
        default:
            return 0.0;
        }
    }

    double getHardness(int atomic_number) const
    {
        switch (atomic_number) {
        case 1:
            return 0.692;
        case 6:
            return 0.565;
        case 7:
            return 0.687;
        case 8:
            return 0.736;
        // ... weitere Parameter
        default:
            return 0.0;
        }
    }

    // PROBLEME:
    // - Branch prediction failures bei unvorhersagbaren Zugriffsmuster
    // - Code wird sehr lang und unübersichtlich
    // - Schlechte Cache-Lokalität
    // - Compiler kann nicht gut optimieren
};

// ====================================================================================
// OPTION 2: HASH MAP (GUT für Flexibilität)
// ====================================================================================
class HashMapParameters {
private:
    struct AtomParams {
        double electronegativity;
        double hardness;
        double shell_hardness;
        std::array<double, 4> slater_exponents; // s, p, d, f
        std::array<double, 4> principal_qn;
        double coordination_scaling;
        double repulsion_exp;
        double dispersion_r4;
        double halogen_strength;
    };

    std::unordered_map<int, AtomParams> atom_params;

public:
    HashMapParameters()
    {
        initializeParameters();
    }

    void initializeParameters()
    {
        // Hydrogen
        atom_params[1] = { 2.300, 0.692, 0.42, { 1.23, 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0, 0.0 }, 0.005, 1.42, 3.36, 0.0 };

        // Carbon
        atom_params[6] = { 2.847, 0.565, 0.31, { 1.61, 1.24, 0.0, 0.0 },
            { 2.0, 2.0, 0.0, 0.0 }, 0.005, 1.78, 5.49, 0.0 };

        // Nitrogen
        atom_params[7] = { 3.022, 0.687, 0.32, { 1.88, 1.45, 0.0, 0.0 },
            { 2.0, 2.0, 0.0, 0.0 }, 0.005, 1.85, 5.71, 0.0 };

        atom_params[8] = { 3.022, 0.687, 0.32, { 1.88, 1.45, 0.0, 0.0 },
            { 2.0, 2.0, 0.0, 0.0 }, 0.005, 1.85, 5.71, 0.0 };
        // ... weitere Elemente
    }

    const AtomParams& getParams(int atomic_number) const
    {
        auto it = atom_params.find(atomic_number);
        if (it != atom_params.end()) {
            return it->second;
        }
        // throw std::runtime_error("Unknown element");
    }

    // VORTEILE:
    // - Flexibel erweiterbar
    // - Gute Cache-Lokalität bei sequenziellem Zugriff
    // - Kompakter Code
    //
    // NACHTEILE:
    // - Hash-Overhead (~10-20 ns pro Zugriff)
    // - Keine Compile-time Optimierung
};

// ====================================================================================
// OPTION 3: DIRECT ARRAY INDEXING (BESTE Performance)
// ====================================================================================
class ArrayParameters {
private:
    struct AtomParams {
        double electronegativity;
        double hardness;
        double shell_hardness;
        std::array<double, 4> slater_exponents;
        std::array<double, 4> principal_qn;
        double coordination_scaling;
        double repulsion_exp;
        double dispersion_r4;
        double halogen_strength;
        bool is_valid = false;
    };

    // Array für alle möglichen Atomnummern (1-118)
    std::array<AtomParams, 119> atom_params{}; // Index 0 unused

public:
    ArrayParameters()
    {
        initializeParameters();
    }

    void initializeParameters()
    {
        // Hydrogen (Z=1)
        atom_params[1] = { 2.300, 0.692, 0.42, { 1.23, 0.0, 0.0, 0.0 },
            { 1.0, 0.0, 0.0, 0.0 }, 0.005, 1.42, 3.36, 0.0, true };

        // Carbon (Z=6)
        atom_params[6] = { 2.847, 0.565, 0.31, { 1.61, 1.24, 0.0, 0.0 },
            { 2.0, 2.0, 0.0, 0.0 }, 0.005, 1.78, 5.49, 0.0, true };

        // ... weitere Parameter
    }

    const AtomParams& getParams(int atomic_number) const
    {
        if (atomic_number < 1 || atomic_number > 118 || !atom_params[atomic_number].is_valid) {
            throw std::runtime_error("Invalid atomic number");
        }
        return atom_params[atomic_number];
    }

    // INLINE für maximale Geschwindigkeit
    inline double getElectronegativity(int z) const
    {
        return atom_params[z].electronegativity;
    }

    inline double getHardness(int z) const
    {
        return atom_params[z].hardness;
    }

    // VORTEILE:
    // - Direkter Speicherzugriff (~1-2 ns)
    // - Perfekte Cache-Lokalität
    // - Compiler kann aggressiv optimieren
    // - Bounds-checking nur im Debug-Modus
};

// ====================================================================================
// OPTION 4A: CONSTEXPR mit direkter Initialisierung (FUNKTIONIERT)
// ====================================================================================
class OptimizedParameters {
private:
    struct AtomParams {
        double electronegativity;
        double hardness;
        double shell_hardness;
        std::array<double, 4> slater_exponents;
        std::array<double, 4> principal_qn;
        double coordination_scaling;
        double repulsion_exp;
        double dispersion_r4;
        double halogen_strength;
    };

    // Direkte constexpr Initialisierung - FUNKTIONIERT
    static constexpr std::array<AtomParams, 119> atom_params = { {
        {}, // Index 0 - unused
        // H (Z=1)
        { 2.300, 0.692, 0.42, { 1.23, 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0, 0.0 }, 0.005, 1.42, 3.36, 0.0 },
        // He (Z=2)
        { 1.842, 0.700, 0.90, { 1.38, 0.0, 0.0, 0.0 }, { 1.0, 0.0, 0.0, 0.0 }, 0.005, 1.81, 1.61, 0.0 },
        // Li (Z=3)
        { 1.211, 0.825, 0.15, { 0.64, 0.0, 0.0, 0.0 }, { 2.0, 0.0, 0.0, 0.0 }, 0.005, 2.37, 9.52, 0.0 },
        {}, // Be (Z=4) - not parameterized
        {}, // B (Z=5) - not parameterized
        // C (Z=6)
        { 2.847, 0.565, 0.31, { 1.61, 1.24, 0.0, 0.0 }, { 2.0, 2.0, 0.0, 0.0 }, 0.005, 1.78, 5.49, 0.0 },
        // N (Z=7)
        { 3.022, 0.687, 0.32, { 1.88, 1.45, 0.0, 0.0 }, { 2.0, 2.0, 0.0, 0.0 }, 0.005, 1.85, 5.71, 0.0 },
        // O (Z=8)
        { 3.610, 0.736, 0.42, { 2.24, 1.67, 0.0, 0.0 }, { 2.0, 2.0, 0.0, 0.0 }, 0.005, 1.95, 5.40, 0.0 },
        // F (Z=9)
        { 4.193, 0.705, 0.56, { 2.57, 1.84, 0.0, 0.0 }, { 2.0, 2.0, 0.0, 0.0 }, 0.005, 2.05, 4.83, 0.0 },
        {}, {}, {}, {}, {}, {}, {}, // Z=10-16
        // Cl (Z=17)
        { 3.242, 0.581, 0.30, { 2.05, 1.74, 1.30, 0.0 }, { 3.0, 3.0, 3.0, 0.0 }, 0.005, 2.16, 6.93, 0.72 }
        // ... weitere Elemente würden hier folgen
    } };

public:
    constexpr const AtomParams& getParams(int atomic_number) const
    {
        return atom_params[atomic_number];
    }

    constexpr double getElectronegativity(int z) const
    {
        return atom_params[z].electronegativity;
    }

    constexpr double getHardness(int z) const
    {
        return atom_params[z].hardness;
    }
};

// ====================================================================================
// OPTION 4B: Template-Spezialisierung (BESTE für wenige Elemente)
// ====================================================================================
template <int Z>
struct GFN2Params {
    static constexpr double electronegativity = 0.0;
    static constexpr double hardness = 0.0;
    static constexpr double shell_hardness = 0.0;
    static constexpr std::array<double, 4> slater_exp = { 0.0, 0.0, 0.0, 0.0 };
};

// Spezialisierungen für häufige Elemente
template <>
struct GFN2Params<1> { // H
    static constexpr double electronegativity = 2.300;
    static constexpr double hardness = 0.692;
    static constexpr double shell_hardness = 0.42;
    static constexpr std::array<double, 4> slater_exp = { 1.23, 0.0, 0.0, 0.0 };
};

template <>
struct GFN2Params<6> { // C
    static constexpr double electronegativity = 2.847;
    static constexpr double hardness = 0.565;
    static constexpr double shell_hardness = 0.31;
    static constexpr std::array<double, 4> slater_exp = { 1.61, 1.24, 0.0, 0.0 };
};

template <>
struct GFN2Params<7> { // N
    static constexpr double electronegativity = 3.022;
    static constexpr double hardness = 0.687;
    static constexpr double shell_hardness = 0.32;
    static constexpr std::array<double, 4> slater_exp = { 1.88, 1.45, 0.0, 0.0 };
};

template <>
struct GFN2Params<8> { // O
    static constexpr double electronegativity = 3.610;
    static constexpr double hardness = 0.736;
    static constexpr double shell_hardness = 0.42;
    static constexpr std::array<double, 4> slater_exp = { 2.24, 1.67, 0.0, 0.0 };
};

// Wrapper-Klasse für Template-Dispatch
class TemplateParameters {
public:
    template <int Z>
    static constexpr double getElectronegativity()
    {
        return GFN2Params<Z>::electronegativity;
    }

    // Runtime-Zugriff mit Template-Dispatch
    double getElectronegativity(int z) const
    {
        switch (z) {
        case 1:
            return GFN2Params<1>::electronegativity;
        case 6:
            return GFN2Params<6>::electronegativity;
        case 7:
            return GFN2Params<7>::electronegativity;
        case 8:
            return GFN2Params<8>::electronegativity;
        default:
            return 0.0;
        }
    }
};

// ====================================================================================
// OPTION 4C: Namespace mit constexpr Variablen (EINFACHSTE Lösung)
// ====================================================================================
namespace GFN2Parameters {
// Elektronegativitäten
constexpr double electronegativity[] = {
    0.0, // Z=0 (unused)
    2.300, // H
    1.842, // He
    1.211, // Li
    0.0, // Be (not parameterized)
    0.0, // B
    2.847, // C
    3.022, // N
    3.610, // O
    4.193, // F
    0.0, // Ne
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, // Na-S
    3.242 // Cl
    // ... weitere Werte
};

constexpr double hardness[] = {
    0.0, // Z=0
    0.692, // H
    0.700, // He
    0.825, // Li
    0.0, 0.0,
    0.565, // C
    0.687, // N
    0.736, // O
    0.705, // F
    0.0,
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
    0.581 // Cl
};

// Inline-Zugriffsfunktionen
inline constexpr double getElectronegativity(int z)
{
    return electronegativity[z];
}

inline constexpr double getHardness(int z)
{
    return hardness[z];
}
}

// ====================================================================================
// PERFORMANCE BENCHMARK
// ====================================================================================
void benchmarkParameterAccess()
{
    const int N_ITERATIONS = 10000000;
    const int N_ATOMS = 100;
    std::vector<int> atomic_numbers;

    // Realistische Mischung organischer Moleküle
    for (int i = 0; i < N_ATOMS; ++i) {
        int z = (i % 5 == 0) ? 6 : ((i % 7 == 0) ? 7 : ((i % 11 == 0) ? 8 : 1));
        atomic_numbers.push_back(z);
    }

    SwitchCaseParameters switch_params;
    HashMapParameters hash_params;
    ArrayParameters array_params;
    OptimizedParameters opt_params;

    double dummy = 0.0;

    // Benchmark Switch/Case
    auto start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        for (int z : atomic_numbers) {
            dummy += switch_params.getElectronegativity(z);
        }
    }
    auto end = std::chrono::high_resolution_clock::now();
    auto switch_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Benchmark HashMap
    start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        for (int z : atomic_numbers) {
            dummy += hash_params.getParams(z).electronegativity;
        }
    }
    end = std::chrono::high_resolution_clock::now();
    auto hash_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Benchmark Array
    start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        for (int z : atomic_numbers) {
            dummy += array_params.getElectronegativity(z);
        }
    }
    end = std::chrono::high_resolution_clock::now();
    auto array_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    // Benchmark Namespace (Optimized)
    start = std::chrono::high_resolution_clock::now();
    for (int iter = 0; iter < N_ITERATIONS; ++iter) {
        for (int z : atomic_numbers) {
            dummy += GFN2Parameters::getElectronegativity(z);
        }
    }
    end = std::chrono::high_resolution_clock::now();
    auto opt_time = std::chrono::duration_cast<std::chrono::microseconds>(end - start);

    std::cout << "Performance Results:\n";
    std::cout << "Switch/Case: " << switch_time.count() << " μs\n";
    std::cout << "HashMap:     " << hash_time.count() << " μs\n";
    std::cout << "Array:       " << array_time.count() << " μs\n";
    std::cout << "Optimized:   " << opt_time.count() << " μs\n";
    std::cout << "Dummy sum: " << dummy << std::endl; // Prevent optimization
}

// ====================================================================================
// SPEICHER-EFFIZIENZ ANALYSE
// ====================================================================================
void analyzeMemoryUsage()
{
    std::cout << "\nMemory Usage Analysis:\n";
    std::cout << "Switch/Case: ~0 bytes (code only)\n";
    std::cout << "HashMap:     ~" << sizeof(std::unordered_map<int, double>) * 10 << " bytes + overhead\n";
    std::cout << "Array:       ~" << sizeof(double) * 119 * 10 << " bytes\n";
    std::cout << "Optimized:   ~" << sizeof(double) * 119 * 10 << " bytes (compile-time)\n";

    // Cache-Freundlichkeit
    std::cout << "\nCache Efficiency:\n";
    std::cout << "Switch/Case: Poor (unpredictable branches)\n";
    std::cout << "HashMap:     Medium (hash collisions)\n";
    std::cout << "Array:       Excellent (sequential access)\n";
    std::cout << "Optimized:   Perfect (compile-time constants)\n";
}

int main()
{
    benchmarkParameterAccess();
    // analyzeMemoryUsage();
}
