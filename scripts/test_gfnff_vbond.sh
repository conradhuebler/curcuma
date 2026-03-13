#!/bin/bash
# Test script for verifying GFN-FF vbond parameters

set -ex

echo "Building Curcuma with GFN-FF vbond parameter verification test..."

# Build the project
cd release
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo -DUSE_XTB=true -DUSE_TBLITE=true -DUSE_D3=true -DUSE_D4=true ..
make -j4 test_gfnff_vbond

# Run the test
echo "Running GFN-FF vbond parameter verification test..."
./test_cases/test_gfnff_vbond

echo "Test completed successfully!"