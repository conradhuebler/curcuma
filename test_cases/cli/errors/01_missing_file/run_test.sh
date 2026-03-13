#!/bin/bash
# Phase 1 Test: Missing file error handling
# Expected: Clear error message "File not found", NO HANG

CURCUMA="${PROJECT_ROOT}/release/curcuma"

# Test with non-existent file - should fail fast with clear error
if ! "$CURCUMA" -sp nonexistent.xyz -method uff 2>&1 | grep -q "File not found: nonexistent.xyz"; then
    echo "FAILED: Expected 'File not found' error message"
    exit 1
fi

exit 0
