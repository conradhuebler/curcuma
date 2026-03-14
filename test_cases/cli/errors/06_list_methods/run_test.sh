#!/bin/bash
# Phase 2 Test: List available methods with --methods flag
# Expected: Help message showing all available computational methods

CURCUMA="${PROJECT_ROOT}/release/curcuma"

# Test --methods flag
output=$("$CURCUMA" -methods 2>&1)

# Check for quantum methods section
if ! echo "$output" | grep -q "Quantum Methods"; then
    echo "FAILED: Missing 'Quantum Methods' section"
    exit 1
fi

# Check for force fields section
if ! echo "$output" | grep -q "Force Fields"; then
    echo "FAILED: Missing 'Force Fields' section"
    exit 1
fi

# Check for at least one method (should list available methods based on compilation)
if ! echo "$output" | grep -q "^\s*-"; then
    echo "FAILED: Expected method listing"
    exit 1
fi

exit 0
