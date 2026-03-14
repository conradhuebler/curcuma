#!/bin/bash
# Phase 1 Test: Permission denied error handling
# Expected: Clear error message mentioning permissions

CURCUMA="${PROJECT_ROOT}/release/curcuma"

# Create a file with no read permissions
echo "2" > no_read.xyz
echo "Test" >> no_read.xyz
echo "H 0.0 0.0 0.0" >> no_read.xyz
echo "H 1.0 0.0 0.0" >> no_read.xyz
chmod 000 no_read.xyz

# Cleanup function
cleanup() {
    chmod 644 no_read.xyz 2>/dev/null || true
    rm -f no_read.xyz 2>/dev/null || true
}
trap cleanup EXIT

# Test with no-read-permission file
if ! "$CURCUMA" -sp no_read.xyz -method uff 2>&1 | grep -qE "(Failed to open|permission|Permission)"; then
    echo "FAILED: Expected permission-related error message"
    exit 1
fi

exit 0
