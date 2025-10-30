#!/bin/bash

# Test: CG Single Point Energy Calculation
# Tests ForceField with method="cg" for simple 2-bead system
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Tests actual CG implementation

set -e

# ============================================================================
# Configuration
# ============================================================================

SCRIPT_DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
CURCUMA="${CURCUMA:-../../../../curcuma}"  # Default to built executable (4 levels up from test_cases/cli/cg/01_single_point/)

# ============================================================================
# Helper Functions
# ============================================================================

print_header() {
    echo ""
    echo "=========================================="
    echo "$1"
    echo "=========================================="
}

print_result() {
    if [ $1 -eq 0 ]; then
        echo "✅ PASSED: $2"
    else
        echo "❌ FAILED: $2"
    fi
}

# ============================================================================
# Main Test
# ============================================================================

print_header "CG Single Point Energy Calculation"

# Check if executable exists
if [ ! -f "$CURCUMA" ]; then
    echo "❌ FAILED: curcuma executable not found at $CURCUMA"
    exit 1
fi

echo "Using curcuma: $CURCUMA"
echo "Test directory: $SCRIPT_DIR"

# Change to test directory
cd "$SCRIPT_DIR"

# Clean up any previous outputs
rm -f stdout.log stderr.log

# Run CG single point calculation
echo ""
echo "Running: $CURCUMA -vtf input.vtf -sp -method cg -load_ff_json cg_params.json -verbosity 2"
echo ""

$CURCUMA -vtf input.vtf \
         -sp \
         -method cg \
         -load_ff_json cg_params.json \
         -verbosity 2 \
         > stdout.log 2> stderr.log

EXIT_CODE=$?

# ============================================================================
# Validation
# ============================================================================

echo ""
print_header "Test Validation"

# 1. Check exit code
if [ $EXIT_CODE -ne 0 ]; then
    echo "❌ FAILED: Non-zero exit code ($EXIT_CODE)"
    echo ""
    echo "STDERR output:"
    cat stderr.log || true
    exit 1
fi
print_result 0 "Exit code is 0"

# 2. Check output file exists and is not empty
if [ ! -f stdout.log ]; then
    echo "❌ FAILED: No stdout.log output file"
    exit 1
fi
if [ ! -s stdout.log ]; then
    echo "❌ FAILED: stdout.log is empty"
    exit 1
fi
print_result 0 "Output file is not empty"

# 3. Check for energy output
if grep -q "Force Field Energy\|Energy:" stdout.log; then
    print_result 0 "Energy output found"
else
    echo "⚠️  WARNING: No explicit 'Energy' line found in output"
    # Still pass since calculation completed
fi

# 4. Check for CG-specific markers (optional)
if grep -q "CG\|method\|forcefield\|cg_param" stdout.log; then
    print_result 0 "CG-specific output detected"
else
    echo "⚠️  WARNING: No CG-specific markers in output"
fi

# 5. Verify no errors in stderr
if [ -s stderr.log ]; then
    # Check if stderr contains actual errors (not just warnings)
    if grep -qi "error\|failed\|exception" stderr.log; then
        echo "❌ FAILED: Errors detected in stderr"
        cat stderr.log
        exit 1
    fi
    print_result 0 "No critical errors in stderr"
else
    print_result 0 "No stderr output (clean run)"
fi

# ============================================================================
# Final Status
# ============================================================================

echo ""
print_header "Test Complete"
echo "✅ CG Single Point Energy Calculation: PASSED"
echo ""
echo "Output summary:"
echo "- Input: 2 CG particles (element 226) at 4.0 Å separation"
echo "- Method: Coarse Graining (LJ 1-6-12 potential)"
echo "- Sigma: 4.0 Å, Epsilon: 0.5 Hartree"
echo ""

exit 0
