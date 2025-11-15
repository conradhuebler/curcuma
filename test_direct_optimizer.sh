#!/bin/bash
# Direct test of new optimizer infrastructure

cd /home/user/curcuma/build

echo "=========================================="
echo "Testing New Optimizer Infrastructure"
echo "=========================================="
echo ""

# Test with water molecule
echo "Test 1: LBFGSPP optimizer with UFF on water (H2O)"
echo "--------------------------------------------------"
./curcuma -opt /home/user/curcuma/test_cases/cli/curcumaopt/01_default_uff_opt/input.xyz \
    -method uff -MaxIter 20 -verbosity 2

echo ""
echo ""
echo "Test 2: Testing with larger molecule (A.xyz, 117 atoms)"
echo "--------------------------------------------------------"  
./curcuma -opt /home/user/curcuma/test_cases/AAA-bGal/A.xyz \
    -method uff -MaxIter 10 -verbosity 2

echo ""
echo "=========================================="
echo "Tests completed!"
echo "=========================================="
