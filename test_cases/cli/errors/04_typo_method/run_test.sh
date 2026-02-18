#!/bin/bash
# Phase 2 Test: Method name typo handling with fuzzy matching
# Expected: Error message suggests closest matching methods

CURCUMA="${PROJECT_ROOT}/release/curcuma"

# Create a simple test molecule in current directory (tests run from build tree)
cat > typo_test.xyz << 'MOLEOF'
2
Test molecule
H 0.0 0.0 0.0
H 1.0 0.0 0.0
MOLEOF

# Test with typo method name (gfn3 instead of gfn2)
output=$("$CURCUMA" -sp typo_test.xyz -method gfn3 2>&1)

# Check for the error message
if ! echo "$output" | grep -q "Unknown computational method"; then
    echo "FAILED: Expected unknown method error"
    exit 1
fi

# Check for suggestions (fuzzy matching)
if ! echo "$output" | grep -q "Did you mean"; then
    echo "FAILED: Expected fuzzy matching suggestions"
    exit 1
fi

# Just verify we don't silently fall back to UFF
if echo "$output" | grep -q "Using UFF"; then
    echo "FAILED: Should not silently fall back to UFF"
    exit 1
fi

exit 0
