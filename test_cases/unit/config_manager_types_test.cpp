/*
 * Test of the type handling in ConfigManager::get (docs/debt_2026-11/PARAMS.md M8).
 * Claude Generated (Nov 2026).
 *
 *  - key missing                    -> get<T>(key) throws std::runtime_error (not ConfigTypeError),
 *                                      get<T>(key, default) returns the default
 *  - key present, wrong value type  -> get<T>(key) throws ConfigTypeError naming key and value,
 *                                      get<T>(key, default) returns the default (and logs a warning)
 *  - key present, convertible       -> the converted value (an integral double read as int)
 */
#include "src/core/config_manager.h"
#include "src/core/parameter_registry.h"
#include "generated/parameter_registry.h"

#include <iostream>

static int g_failed = 0;
#define CHECK(cond, text)                                                          \
    do {                                                                           \
        if (!(cond)) {                                                             \
            std::cerr << "FAIL: " << text << " (" #cond ")\n";                     \
            ++g_failed;                                                            \
        }                                                                          \
    } while (0)

int main()
{
    initialize_generated_registry();

    nlohmann::json input = { { "opt", { { "threads", "abc" }, { "write_trajectory", 0.0 }, { "max_iterations", 12.0 } } } };
    ConfigManager config("opt", input);

    // wrong type: string where an int is needed
    bool typed = false;
    try {
        config.get<int>("threads");
    } catch (const ConfigTypeError& e) {
        typed = std::string(e.what()).find("\"abc\"") != std::string::npos
            && std::string(e.what()).find("'threads'") != std::string::npos;
    } catch (...) {
    }
    CHECK(typed, "get<int>(threads) with value \"abc\" throws ConfigTypeError naming key and value");
    CHECK(config.get<int>("threads", 7) == 7, "get<int>(threads, 7) returns the default for \"abc\"");

    // wrong type: number where a bool is needed
    bool bool_typed = false;
    try {
        config.get<bool>("write_trajectory");
    } catch (const ConfigTypeError&) {
        bool_typed = true;
    } catch (...) {
    }
    CHECK(bool_typed, "get<bool>(write_trajectory) with value 0.0 throws ConfigTypeError");

    // convertible: an integral double read as int
    CHECK(config.get<int>("max_iterations") == 12, "get<int>(max_iterations) reads 12.0 as 12");

    // missing key: runtime_error, not ConfigTypeError; default without warning
    bool missing_plain = false;
    try {
        config.get<int>("no_such_parameter_at_all");
    } catch (const ConfigTypeError&) {
    } catch (const std::runtime_error&) {
        missing_plain = true;
    }
    CHECK(missing_plain, "a missing key throws a plain runtime_error");
    CHECK(config.get<int>("no_such_parameter_at_all", 5) == 5, "a missing key returns the default");

    std::cout << (g_failed ? "FAILED" : "all checks passed") << "\n";
    return g_failed ? 1 : 0;
}
