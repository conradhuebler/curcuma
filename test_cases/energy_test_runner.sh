#!/bin/bash

# Energy Test Runner Script
# Comprehensive test suite for Curcuma energy calculation methods
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>

set -e  # Exit on any error

# Colors for output
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Test configuration
TEST_MOLECULE="AAA-bGal/A.xyz"
CURCUMA_BINARY="../curcuma"
TEST_OUTPUT_DIR="energy_test_results"
TIMESTAMP=$(date +"%Y%m%d_%H%M%S")

echo -e "${BLUE}=== Curcuma Energy Method Test Suite ===${NC}"
echo "Timestamp: $(date)"
echo "Test molecule: ${TEST_MOLECULE}"
echo "Curcuma binary: ${CURCUMA_BINARY}"
echo

# Create output directory
mkdir -p "${TEST_OUTPUT_DIR}"

# Function to test a single method
test_method() {
    local method=$1
    local expected_energy=$2
    local tolerance=$3
    
    echo -n "Testing ${method}... "
    
    # Run curcuma with the method
    local start_time=$(date +%s.%N)
    local output_file="${TEST_OUTPUT_DIR}/${method}_${TIMESTAMP}.out"
    local error_file="${TEST_OUTPUT_DIR}/${method}_${TIMESTAMP}.err"
    
    if timeout 60s ${CURCUMA_BINARY} -sp ${TEST_MOLECULE} -method ${method} > "${output_file}" 2> "${error_file}"; then
        local end_time=$(date +%s.%N)
        local runtime=$(echo "${end_time} - ${start_time}" | bc -l)
        
        # Extract energy from output
        local calculated_energy=$(grep "Single Point Energy" "${output_file}" | tail -1 | grep -o "\-[0-9]*\.[0-9]*" | head -1)
        
        if [[ -n "${calculated_energy}" ]]; then
            # Calculate absolute difference
            local energy_diff=$(echo "${calculated_energy} - ${expected_energy}" | bc -l | sed 's/^-//')
            
            # Check if within tolerance
            local within_tolerance=$(echo "${energy_diff} <= ${tolerance}" | bc -l)
            
            if [[ "${within_tolerance}" == "1" ]]; then
                echo -e "${GREEN}PASS${NC} Energy: ${calculated_energy} Eh (Expected: ${expected_energy}) Time: ${runtime}s"
                echo "${method},PASS,${calculated_energy},${expected_energy},${energy_diff},${runtime}" >> "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv"
            else
                echo -e "${RED}FAIL${NC} Energy: ${calculated_energy} Eh (Expected: ${expected_energy}, Diff: ${energy_diff}, Tol: ${tolerance})"
                echo "${method},FAIL,${calculated_energy},${expected_energy},${energy_diff},${runtime}" >> "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv"
            fi
        else
            echo -e "${RED}FAIL${NC} Could not extract energy from output"
            echo "${method},FAIL,N/A,${expected_energy},N/A,${runtime}" >> "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv"
        fi
    else
        echo -e "${RED}FAIL${NC} Command failed or timed out"
        echo "${method},FAIL,ERROR,${expected_energy},N/A,N/A" >> "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv"
    fi
}

# Check if curcuma binary exists
if [[ ! -f "${CURCUMA_BINARY}" ]]; then
    echo -e "${RED}Error: Curcuma binary not found at ${CURCUMA_BINARY}${NC}"
    echo "Please build curcuma first or adjust the CURCUMA_BINARY path"
    exit 1
fi

# Check if test molecule exists
if [[ ! -f "${TEST_MOLECULE}" ]]; then
    echo -e "${RED}Error: Test molecule ${TEST_MOLECULE} not found${NC}"
    echo "Please ensure A.xyz exists in the current directory"
    exit 1
fi

# Initialize summary CSV
echo "method,result,calculated_energy,expected_energy,energy_diff,time_seconds" > "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv"

echo -e "${YELLOW}Running comprehensive energy method tests...${NC}"
echo

# Test quantum methods (high precision expected)
echo -e "${BLUE}=== Quantum Methods ===${NC}"
test_method "gfn2"   "-165.75590286" "0.000001"    # TBLite GFN2-xTB (reference)
test_method "ugfn2"  "-165.75593207" "0.000001"    # Ulysses ugfn2
test_method "gfn1"   "-165.7"        "0.1"         # GFN1-xTB (approximate)
test_method "ipea1"  "-165.7"        "0.1"         # iPEA1-xTB (approximate)

echo
echo -e "${BLUE}=== Semi-empirical Methods ===${NC}"
test_method "pm6"    "-328.41882169" "0.000001"    # Ulysses PM6 (reference)
test_method "pm3"    "-320.0"        "5.0"         # PM3 (approximate)
test_method "am1"    "-315.0"        "5.0"         # AM1 (approximate)

echo
echo -e "${BLUE}=== Force Field Methods ===${NC}"
test_method "gfnff"  "-20.29485162"  "0.00001"     # External GFN-FF (reference)
test_method "uff"    "-20.0"         "2.0"         # UFF (approximate)

echo
echo -e "${BLUE}=== Native Methods ===${NC}"
test_method "eht"    "-150.0"        "10.0"        # Extended Hückel Theory (approximate)

echo
echo -e "${YELLOW}=== Test Summary ===${NC}"

# Generate summary report
total_tests=$(tail -n +2 "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv" | wc -l)
passed_tests=$(tail -n +2 "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv" | grep ",PASS," | wc -l)
failed_tests=$(tail -n +2 "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv" | grep ",FAIL," | wc -l)

echo "Total tests: ${total_tests}"
echo -e "Passed: ${GREEN}${passed_tests}${NC}"
echo -e "Failed: ${RED}${failed_tests}${NC}"

if [[ ${failed_tests} -gt 0 ]]; then
    echo
    echo -e "${RED}Failed tests:${NC}"
    tail -n +2 "${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv" | grep ",FAIL," | cut -d',' -f1 | while read method; do
        echo "  - ${method}"
    done
fi

echo
echo "Detailed results saved in: ${TEST_OUTPUT_DIR}/"
echo "Summary CSV: ${TEST_OUTPUT_DIR}/summary_${TIMESTAMP}.csv"

# Create permanent symlink to latest results
ln -sf "summary_${TIMESTAMP}.csv" "${TEST_OUTPUT_DIR}/latest_summary.csv"

echo -e "${BLUE}Test suite completed!${NC}"

# Exit with non-zero code if any tests failed
exit ${failed_tests}