#!/bin/bash
# H-H PES Scan: Bond Energy Analysis
# Scan from 0.5 Å to 3.0 Å in 0.1 Å steps

XTB_BIN="/home/conrad/Downloads/xtb-6.6.1/bin/xtb"
CURCUMA_BIN="/home/conrad/src/claude_curcuma/curcuma/release/curcuma"
SCAN_DIR="/tmp/HH_PES_scan"

mkdir -p $SCAN_DIR
cd $SCAN_DIR

# Create header for results
cat > results.txt << 'EOF'
r_ang	r_bohr	E_total_xtb	E_bond_xtb	E_rep_xtb	E_disp_xtb	E_total_cur	E_bond_cur	E_rep_cur	E_disp_cur	factor
EOF

echo "Starting H-H PES scan..."
echo ""

# Loop from 0.5 Å to 3.0 Å in 0.1 Å steps
for r_ang in $(seq 0.5 0.1 3.0); do
    r_bohr=$(echo "scale=8; $r_ang / 0.529177" | bc)

    # Generate H-H XYZ file
    cat > HH_${r_ang}.xyz << XYZEOF
2

H  0.0 0.0 0.0
H  ${r_ang} 0.0 0.0
XYZEOF

    echo -n "Distance: $r_ang Å ($r_bohr Bohr) ... "

    # Run XTB
    $XTB_BIN HH_${r_ang}.xyz --gfnff --sp > HH_${r_ang}.xtb.out 2>&1
    E_total_xtb=$(grep "TOTAL ENERGY" HH_${r_ang}.xtb.out | awk '{print $3}')
    E_bond_xtb=$(grep "bond energy" HH_${r_ang}.xtb.out | awk '{print $3}' | head -1)
    E_rep_xtb=$(grep "repulsion energy" HH_${r_ang}.xtb.out | awk '{print $3}')
    E_disp_xtb=$(grep "dispersion energy" HH_${r_ang}.xtb.out | awk '{print $3}')

    # Run Curcuma
    $CURCUMA_BIN -sp HH_${r_ang}.xyz -method cgfnff -verbosity 3 > HH_${r_ang}.cur.out 2>&1
    E_total_cur=$(grep "cgfnff Final Energy" HH_${r_ang}.cur.out | awk '{print $4}')
    E_bond_cur=$(grep "Bond 0-1:" HH_${r_ang}.cur.out | grep "E=" | awk -F'E=' '{print $2}' | awk '{print $1}')

    # Extract Repulsion and Dispersion from Force Field Energy
    # Need to parse the verbose output more carefully
    E_rep_cur=$(grep "thread_repulsion_energy" HH_${r_ang}.cur.out | awk '{print $2}' | sed 's/Eh//')
    E_disp_cur=$(grep "thread_dispersion_energy" HH_${r_ang}.cur.out | awk '{print $2}' | sed 's/Eh//')

    # Default to 0 if not found
    E_rep_cur=${E_rep_cur:-0}
    E_disp_cur=${E_disp_cur:-0}

    # Calculate scaling factor
    if [ ! -z "$E_bond_xtb" ] && [ ! -z "$E_bond_cur" ]; then
        factor=$(echo "scale=2; $E_bond_xtb / $E_bond_cur" | bc)
    else
        factor="N/A"
    fi

    # Write results
    printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "$r_ang" "$r_bohr" "$E_total_xtb" "$E_bond_xtb" "$E_rep_xtb" "$E_disp_xtb" \
        "$E_total_cur" "$E_bond_cur" "$E_rep_cur" "$E_disp_cur" "$factor" >> results.txt

    echo "XTB: $E_total_xtb Eh (bond: $E_bond_xtb), Cur: $E_total_cur Eh (bond: $E_bond_cur), factor: $factor"
done

echo ""
echo "Results saved to: $SCAN_DIR/results.txt"
echo ""
echo "=== PES SCAN RESULTS ==="
cat results.txt
