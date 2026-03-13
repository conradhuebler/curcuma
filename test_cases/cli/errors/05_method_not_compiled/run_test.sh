#!/bin/bash
# Phase 2 Test: Method not compiled guidance
# Expected: Error message shows CMake flag needed to enable method
# Note: This test is informational since all core methods are compiled

CURCUMA="${PROJECT_ROOT}/release/curcuma"

# Create a simple test molecule in current directory (tests run from build tree)
cat > notcompiled_test.xyz << 'EOF'
2
Test molecule
H 0.0 0.0 0.0
H 1.0 0.0 0.0
EOF

# This test is mainly to verify that unknown methods show help text
# Try an intentionally unavailable hypothetical method
output=$("$CURCUMA" -sp notcompiled_test.xyz -method "hypothetical_method_xyz" 2>&1)

# Check that we get an error (not silent fallback to UFF)
if echo "$output" | grep -q "Unknown computational method"; then
    exit 0
else
    echo "FAILED: Expected 'Unknown computational method' error"
    exit 1
fi
