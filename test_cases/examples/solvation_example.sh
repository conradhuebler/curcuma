#!/bin/bash
# Solvation Examples for Curcuma
# Demonstrates implicit solvation with GFN and MNDO methods

echo "=== Curcuma Solvation Examples ==="
echo ""

# Create test molecule (water)
cat > water.xyz << EOF
3
Water molecule
O      0.000000    0.000000    0.119262
H      0.000000    0.763239   -0.477047
H      0.000000   -0.763239   -0.477047
EOF

echo "1. Gas Phase (No Solvation)"
curcuma -sp water.xyz -method gfn2:tblite -verbosity 1 | grep "Energy:"
echo ""

echo "2. Water Solvation (Default: GB Model)"
curcuma -sp water.xyz -method gfn2:tblite -solvent water -verbosity 1 | grep -E "solvation|Energy:"
echo ""

echo "3. DMSO Solvation"
curcuma -sp water.xyz -method gfn2:tblite -solvent dmso -verbosity 1 | grep -E "solvent|Energy:"
echo ""

echo "4. Custom Dielectric Constant (epsilon=10.0)"
curcuma -sp water.xyz -method gfn2:tblite -solvent_epsilon 10.0 -solvent_model 1 -verbosity 1 | grep -E "dielectric|Energy:"
echo ""

echo "5. CPCM Model (Explicit)"
curcuma -sp water.xyz -method gfn2:tblite -solvent water -solvent_model 1 -verbosity 1 | grep -E "CPCM|Energy:"
echo ""

echo "6. ALPB Model (Most Accurate)"
curcuma -sp water.xyz -method gfn2:tblite -solvent water -solvent_model 3 -verbosity 1 | grep -E "ALPB|Energy:"
echo ""

# Cleanup
rm -f water.xyz

echo "=== Solvation Examples Complete ==="
echo ""
echo "Note: These examples require TBLite compilation (USE_TBLITE=ON)"
echo "For Ulysses solvation, use methods like: pm6:ulysses, gfn2:ulysses"
