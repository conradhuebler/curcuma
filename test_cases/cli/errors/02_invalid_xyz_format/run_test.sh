#!/bin/bash
# Phase 1 Test: Invalid XYZ format error handling
# Expected: Clear error message with filename and line number

CURCUMA="${PROJECT_ROOT}/release/curcuma"

# Test with invalid XYZ format (first line not a number)
if ! "$CURCUMA" -sp invalid.xyz -method uff 2>&1 | grep -q "Invalid XYZ format.*line 1.*Expected atom count"; then
    echo "FAILED: Expected detailed parsing error with line number"
    exit 1
fi

exit 0
