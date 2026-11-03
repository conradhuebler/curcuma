#!/bin/bash
# GPU plugin ABI version check through the curcuma CLI (S0) - Claude Generated (Nov 2026).
# usage: check_gpu_plugin_abi_cli.sh <curcuma> <dir with stub libcurcuma_cuda.so without version> <dir with stub of the
#        right version> <structure.xyz> <work dir>
# Stand-ins have no GPU code and no create_* entry points (test_cases/unit/gpu_plugin_abi_stub/stub_plugin.cpp).
#  (a) old plugin (no version symbol): refusal message, counted fallback, CPU energy equal to the run without -gpu
#  (b) old plugin with -gpu_strict true: exit code 3
#  (c) plugin of the right version: no refusal message (loaded; then declines for lack of an entry point), CPU energy
BIN=$(realpath "$1"); OLD=$(realpath "$2"); OK=$(realpath "$3"); XYZ=$(realpath "$4"); WORK=$5
fail=0
rm -rf "$WORK"; mkdir -p "$WORK"; cd "$WORK" || exit 2
cp "$XYZ" in.xyz
energy() { sed -n 2p "$1" | awk '{print $3}'; }

env -u CURCUMA_PLUGIN_DIR "$BIN" -sp in.xyz -method gfnff -no_bmt -dump_gradient g_cpu.txt -threads 1 > cpu.log 2>&1
e_cpu=$(energy g_cpu.txt)

CURCUMA_PLUGIN_DIR="$OLD" "$BIN" -sp in.xyz -method gfnff -gpu cuda -no_bmt -dump_gradient g_old.txt -threads 1 > old.log 2>&1
rc=$?
if [ $rc -ne 0 ]; then echo "FAIL (a): exit code $rc"; fail=1; fi
grep -q "has plugin ABI version none" old.log || { echo "FAIL (a): no 'has plugin ABI version none' in the output"; fail=1; }
grep -q "rebuild the GPU plugin from the same source tree" old.log || { echo "FAIL (a): no rebuild hint"; fail=1; }
grep -q "ABI/variant check" old.log || { echo "FAIL (a): fallback not counted in the summary"; fail=1; }
e_old=$(energy g_old.txt)
[ -n "$e_cpu" ] && [ "$e_old" == "$e_cpu" ] || { echo "FAIL (a): energy '$e_old' != CPU '$e_cpu'"; fail=1; }

CURCUMA_PLUGIN_DIR="$OLD" "$BIN" -sp in.xyz -method gfnff -gpu cuda -gpu_strict true -no_bmt -threads 1 > strict.log 2>&1
rc=$?
[ $rc -eq 3 ] || { echo "FAIL (b): -gpu_strict exit code $rc, expected 3"; fail=1; }

CURCUMA_PLUGIN_DIR="$OK" "$BIN" -sp in.xyz -method gfnff -gpu cuda -no_bmt -dump_gradient g_ok.txt -threads 1 > ok.log 2>&1
rc=$?
if [ $rc -ne 0 ]; then echo "FAIL (c): exit code $rc"; fail=1; fi
grep -q "plugin ABI version" ok.log && { echo "FAIL (c): plugin of the right version was refused"; fail=1; }
grep -q "missing entry point 'curcuma_cuda_create_gfnff'" ok.log || { echo "FAIL (c): accepted plugin was not asked for its entry point"; fail=1; }
e_ok=$(energy g_ok.txt)
[ "$e_ok" == "$e_cpu" ] || { echo "FAIL (c): energy '$e_ok' != CPU '$e_cpu'"; fail=1; }

echo "E(cpu) = $e_cpu, E(old plugin) = $e_old, E(right version) = $e_ok"
[ $fail -eq 0 ] && echo "gpu_plugin_abi_cli: all checks passed"
exit $fail
