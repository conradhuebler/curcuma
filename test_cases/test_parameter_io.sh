#!/bin/bash
# Comprehensive test suite for parameter export/import system
# Claude Generated - October 2025

set -e  # Exit on error

TEST_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CURCUMA="$TEST_DIR/../release/curcuma"
TEMP_DIR="$TEST_DIR/temp_param_test"

# Colors
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
NC='\033[0m' # No Color

# Counters
TESTS_RUN=0
TESTS_PASSED=0
TESTS_FAILED=0

# Setup
setup() {
    echo "Setting up test environment..."
    mkdir -p "$TEMP_DIR"
    cd "$TEMP_DIR"

    # Create minimal test molecule
    cat > test_molecule.xyz <<EOF
3
Test molecule
C 0.0 0.0 0.0
H 1.0 0.0 0.0
H 0.0 1.0 0.0
EOF
}

# Cleanup
cleanup() {
    echo "Cleaning up..."
    rm -rf "$TEMP_DIR"
}

# Test functions
test_export_defaults() {
    echo -n "Test 1: Export default config... "
    TESTS_RUN=$((TESTS_RUN + 1))

    $CURCUMA -export-config analysis 2>&1 | sed -n '/^{/,/^}/p' > default_config.json

    if [ -f default_config.json ] && grep -q "output_format" default_config.json; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

test_export_run_basic() {
    echo -n "Test 2: Export run config (basic)... "
    TESTS_RUN=$((TESTS_RUN + 1))

    $CURCUMA -analysis test_molecule.xyz -output_format json -export_run run1.json 2>/dev/null

    if [ -f run1.json ] && grep -q '"output_format": "json"' run1.json; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC} - run1.json not created or missing output_format"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        [ -f run1.json ] && echo "Content:" && cat run1.json
    fi
}

test_export_run_complex() {
    echo -n "Test 3: Export run config (complex parameters)... "
    TESTS_RUN=$((TESTS_RUN + 1))

    $CURCUMA -analysis test_molecule.xyz \
        -properties basic \
        -window 20 \
        -topological_colormap viridis \
        -export_run run2.json 2>/dev/null

    if [ -f run2.json ] && \
       grep -q '"properties": "basic"' run2.json && \
       grep -q '"window": 20' run2.json && \
       grep -q '"topological_colormap": "viridis"' run2.json; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC} - Parameters not correctly exported"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        [ -f run2.json ] && echo "Content:" && cat run2.json
    fi
}

test_import_config() {
    echo -n "Test 4: Import custom config... "
    TESTS_RUN=$((TESTS_RUN + 1))

    # Create custom config
    cat > custom.json <<EOF
{
  "analysis": {
    "output_format": "csv",
    "properties": "geometric",
    "window": 50
  }
}
EOF

    # Run with import - export to verify
    $CURCUMA -analysis test_molecule.xyz \
        -import_config custom.json \
        -export_run import_test.json 2>/dev/null

    if [ -f import_test.json ] && grep -q '"output_format": "csv"' import_test.json; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC} - Config not imported correctly"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        [ -f import_test.json ] && echo "Content:" && cat import_test.json
    fi
}

test_cli_priority() {
    echo -n "Test 5: CLI arguments override imported config... "
    TESTS_RUN=$((TESTS_RUN + 1))

    # custom.json has output_format: csv (from previous test)
    # CLI specifies output_format: json (should win)
    $CURCUMA -analysis test_molecule.xyz \
        -import_config custom.json \
        -output_format json \
        -export_run priority_test.json 2>/dev/null

    if [ -f priority_test.json ] && grep -q '"output_format": "json"' priority_test.json; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC} - CLI priority not working"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        [ -f priority_test.json ] && echo "Content:" && cat priority_test.json
    fi
}

test_roundtrip() {
    echo -n "Test 6: Full roundtrip (export → modify → import → export)... "
    TESTS_RUN=$((TESTS_RUN + 1))

    # Step 1: Run with custom params, export
    $CURCUMA -analysis test_molecule.xyz \
        -properties topology \
        -window 100 \
        -export_run roundtrip1.json 2>/dev/null

    # Step 2: Import and run again, export
    $CURCUMA -analysis test_molecule.xyz \
        -import_config roundtrip1.json \
        -export_run roundtrip2.json 2>/dev/null

    # Step 3: Compare (should be identical)
    if [ -f roundtrip1.json ] && [ -f roundtrip2.json ]; then
        val1=$(grep '"properties"' roundtrip1.json | cut -d':' -f2)
        val2=$(grep '"properties"' roundtrip2.json | cut -d':' -f2)

        if [ "$val1" = "$val2" ]; then
            echo -e "${GREEN}PASS${NC}"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC} - Roundtrip altered values"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            echo "Round 1: $val1"
            echo "Round 2: $val2"
        fi
    else
        echo -e "${RED}FAIL${NC} - Roundtrip files not created"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

test_list_modules() {
    echo -n "Test 7: List available modules... "
    TESTS_RUN=$((TESTS_RUN + 1))

    $CURCUMA -list-modules 2>/dev/null | grep -q "analysis (25 parameters)"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

test_help_module() {
    echo -n "Test 8: Module-specific help... "
    TESTS_RUN=$((TESTS_RUN + 1))

    $CURCUMA -help-module analysis 2>/dev/null | grep -q "output_format"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC}"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

test_invalid_config() {
    echo -n "Test 9: Handle invalid JSON gracefully... "
    TESTS_RUN=$((TESTS_RUN + 1))

    echo "{ invalid json" > invalid.json

    $CURCUMA -analysis test_molecule.xyz -import_config invalid.json 2>&1 \
        | grep -q "Failed to parse JSON"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC} - Should report parse error"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

test_missing_config() {
    echo -n "Test 10: Handle missing config file... "
    TESTS_RUN=$((TESTS_RUN + 1))

    $CURCUMA -analysis test_molecule.xyz -import_config nonexistent.json 2>&1 \
        | grep -q "Could not open config file"

    if [ $? -eq 0 ]; then
        echo -e "${GREEN}PASS${NC}"
        TESTS_PASSED=$((TESTS_PASSED + 1))
    else
        echo -e "${RED}FAIL${NC} - Should report missing file"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

test_meta_param_removal() {
    echo -n "Test 11: Meta-parameters removed from exported config... "
    TESTS_RUN=$((TESTS_RUN + 1))

    # Create config with import path
    cat > import_source.json <<EOF
{
  "analysis": {
    "output_format": "csv",
    "properties": "basic"
  }
}
EOF

    # Import and export - meta params should NOT appear in export
    $CURCUMA -analysis test_molecule.xyz \
        -import_config import_source.json \
        -export_run meta_test.json 2>/dev/null

    if [ -f meta_test.json ]; then
        if ! grep -q "import_config" meta_test.json && \
           ! grep -q "export_run" meta_test.json; then
            echo -e "${GREEN}PASS${NC}"
            TESTS_PASSED=$((TESTS_PASSED + 1))
        else
            echo -e "${RED}FAIL${NC} - Meta-parameters should be removed"
            TESTS_FAILED=$((TESTS_FAILED + 1))
            cat meta_test.json
        fi
    else
        echo -e "${RED}FAIL${NC} - Export file not created"
        TESTS_FAILED=$((TESTS_FAILED + 1))
    fi
}

# Main test runner
main() {
    echo "================================"
    echo "Parameter I/O System Test Suite"
    echo "================================"
    echo ""

    if [ ! -f "$CURCUMA" ]; then
        echo -e "${RED}Error: curcuma executable not found at $CURCUMA${NC}"
        echo "Please build curcuma first: cd release && make curcuma"
        exit 1
    fi

    setup

    echo "Running tests..."
    echo ""

    test_export_defaults
    test_export_run_basic
    test_export_run_complex
    test_import_config
    test_cli_priority
    test_roundtrip
    test_list_modules
    test_help_module
    test_invalid_config
    test_missing_config
    test_meta_param_removal

    echo ""
    echo "================================"
    echo "Test Results"
    echo "================================"
    echo "Total:  $TESTS_RUN"
    echo -e "Passed: ${GREEN}$TESTS_PASSED${NC}"
    echo -e "Failed: ${RED}$TESTS_FAILED${NC}"
    echo ""

    cleanup

    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "${GREEN}✅ All tests passed!${NC}"
        exit 0
    else
        echo -e "${RED}❌ Some tests failed!${NC}"
        exit 1
    fi
}

# Run if executed directly
if [ "${BASH_SOURCE[0]}" = "${0}" ]; then
    main "$@"
fi
