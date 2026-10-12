// tests/test_confstat.cpp
#include "src/capabilities/confstat.h"
#include <cassert>
#include <cmath>
#include <iostream>

// Helper function for floating point comparison
bool isClose(double a, double b, double rtol = 1e-10)
{
    return std::abs(a - b) <= rtol * std::max(std::abs(a), std::abs(b));
}

// Test runner class
class TestRunner {
public:
    void run_test(const std::string& name, std::function<void()> test)
    {
        try {
            test();
            std::cout << "✓ " << name << " passed\n";
            m_passed++;
        } catch (const std::exception& e) {
            std::cerr << "✗ " << name << " failed: " << e.what() << "\n";
            m_failed++;
        }
    }

    void summary() const
    {
        std::cout << "\nTest Summary:\n";
        std::cout << "Passed: " << m_passed << "\n";
        std::cout << "Failed: " << m_failed << "\n";
        std::cout << "Total:  " << (m_passed + m_failed) << "\n";
    }

    int get_exit_code() const
    {
        return m_failed > 0 ? 1 : 0;
    }

private:
    int m_passed = 0;
    int m_failed = 0;
};

// Test cases
void test_two_state_system()
{
    json config{
        { "Temp", 298.15 },
        { "Cutoff", 10.0 },
        { "Method", "none" }
    };

    ConfStat stats(config);
    std::vector<double> energies = { -100.0, -99.9 };
    std::vector<int> degeneracies = { 1, 1 };

    stats.setEnergiesWithDegeneracy(energies, degeneracies);
    stats.start();

    // Calculate expected values
    double deltaE = 0.1 * 2625.5; // kJ/mol
    double RT = 8.314 * 298.15 / 1000.0; // kJ/mol
    double expected_pop = 100.0 / (1.0 + exp(-deltaE / RT));

    assert(isClose(stats.getFirstPopulation(), expected_pop));
}

void test_degenerate_system()
{
    json config{
        { "Temp", 298.15 },
        { "Cutoff", 10.0 }
    };

    ConfStat stats(config);
    std::vector<double> energies = { -100.0, -99.9 };
    std::vector<int> degeneracies = { 1, 2 };

    stats.setEnergiesWithDegeneracy(energies, degeneracies);
    stats.start();

    double total_pop = stats.getTotalPopulation();
    assert(isClose(total_pop, 100.0));
}

void test_temperature_dependence()
{
    std::vector<double> energies = { -100.0, -99.9 };
    std::vector<int> degeneracies = { 1, 1 };

    json config_low{ { "Temp", 100.0 } };
    json config_high{ { "Temp", 1000.0 } };

    ConfStat stats_low(config_low);
    ConfStat stats_high(config_high);

    stats_low.setEnergiesWithDegeneracy(energies, degeneracies);
    stats_high.setEnergiesWithDegeneracy(energies, degeneracies);

    stats_low.start();
    stats_high.start();

    double pop_diff_low = std::abs(stats_low.getFirstPopulation() - stats_low.getSecondPopulation());
    double pop_diff_high = std::abs(stats_high.getFirstPopulation() - stats_high.getSecondPopulation());

    assert(pop_diff_low > pop_diff_high);
}

int main()
{
    TestRunner runner;

    runner.run_test("Two State System", test_two_state_system);
    runner.run_test("Degenerate System", test_degenerate_system);
    runner.run_test("Temperature Dependence", test_temperature_dependence);

    runner.summary();
    return runner.get_exit_code();
}
