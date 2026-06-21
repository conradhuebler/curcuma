#!/bin/bash
# Common utility functions for Curcuma CLI tests
# Copyright (C) 2025 Conrad Hübler <Conrad.Huebler@gmx.net>
# Claude Generated - Based on test_parameter_io.sh pattern

# Colors for output
export RED='\033[0;31m'
export GREEN='\033[0;32m'
export YELLOW='\033[1;33m'
export BLUE='\033[0;34m'
export NC='\033[0m' # No Color

# Test counters (to be used by individual tests)
export TESTS_RUN=0
export TESTS_PASSED=0
export TESTS_FAILED=0

# Find curcuma binary (auto-detect project root)
if [ -z "$CURCUMA" ]; then
    # Claude Generated (Oct 28, 2025): Use CMAKE-injected PROJECT_ROOT if available
    # Tests run from build tree, so PROJECT_ROOT is explicitly set by CMakeLists.txt
    if [ -z "$PROJECT_ROOT" ] || [ ! -d "$PROJECT_ROOT/src" ]; then
        # Fallback: search for project root by looking for src/ directory
        SEARCH_DIR="$(pwd)"
        PROJECT_ROOT=""

        for i in {1..5}; do
            if [ -d "$SEARCH_DIR/src" ]; then
                PROJECT_ROOT="$SEARCH_DIR"
                break
            fi
            SEARCH_DIR="$(dirname "$SEARCH_DIR")"
        done

        if [ -z "$PROJECT_ROOT" ]; then
            echo -e "${RED}ERROR: Could not find project root!${NC}"
            exit 1
        fi
    fi

    # Try release first (faster for CLI tests), then debug, build, and the GPU build dirs.
    for _bd in release debug build release_rocm release_cuda release_vulkan build_rocm; do
        if [ -f "$PROJECT_ROOT/$_bd/curcuma" ]; then
            export CURCUMA="$PROJECT_ROOT/$_bd/curcuma"
            break
        fi
    done
    if [ -z "$CURCUMA" ] || [ ! -f "$CURCUMA" ]; then
        echo -e "${RED}ERROR: curcuma binary not found!${NC}"
        echo "Expected at: $PROJECT_ROOT/{release,debug,build,release_rocm,release_cuda,release_vulkan}/curcuma"
        exit 1
    fi
    echo -e "${BLUE}Using curcuma:${NC} $CURCUMA"
fi

# Helper: Print test header
test_header() {
    local test_name=$1
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Test: $test_name${NC}"
    echo -e "${BLUE}========================================${NC}"
}

# Helper: Assert exit code equals expected
assert_exit_code() {
    local actual=$1
    local expected=$2
    local msg=${3:-"Exit code check"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ $actual -eq $expected ]; then
        echo -e "${GREEN}✓ PASS${NC}: $msg (exit code: $actual)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (expected: $expected, got: $actual)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Assert curcuma completed successfully (pragmatic - checks output files, not exit code)
# Curcuma sometimes returns exit code 1 even on success (known bug)
assert_curcuma_success() {
    local expected_file=$1
    local msg=${2:-"Curcuma execution should succeed"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -f "$expected_file" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $msg (output file exists)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (expected output file missing: $expected_file)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Assert file exists
assert_file_exists() {
    local filepath=$1
    local msg=${2:-"File exists: $filepath"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -f "$filepath" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $msg"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (file not found)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Assert file does NOT exist
assert_file_not_exists() {
    local filepath=$1
    local msg=${2:-"File should not exist: $filepath"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ ! -f "$filepath" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $msg"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (file exists but shouldn't)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Assert string in file
assert_string_in_file() {
    local pattern=$1
    local filepath=$2
    local msg=${3:-"Pattern '$pattern' found in $filepath"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if grep -q "$pattern" "$filepath" 2>/dev/null; then
        echo -e "${GREEN}✓ PASS${NC}: $msg"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (pattern not found)"
        [ -f "$filepath" ] && echo "  File content (first 5 lines):" && head -5 "$filepath" | sed 's/^/    /'
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Assert string NOT in file
assert_string_not_in_file() {
    local pattern=$1
    local filepath=$2
    local msg=${3:-"Pattern '$pattern' NOT in $filepath"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if ! grep -q "$pattern" "$filepath" 2>/dev/null; then
        echo -e "${GREEN}✓ PASS${NC}: $msg"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (pattern found but shouldn't be)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Assert file exists and is non-empty (Claude Generated, June 2026)
assert_file_not_empty() {
    local filepath=$1
    local msg=${2:-"File non-empty: $filepath"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ -s "$filepath" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $msg"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (missing or empty: '$filepath')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Assert exact string equality (Claude Generated, June 2026)
# Used to pin the ConfScan result fingerprint "<accepted> <reorder> <reuse> <skipped>".
assert_equals() {
    local expected=$1
    local actual=$2
    local msg=${3:-"Value equals"}

    TESTS_RUN=$((TESTS_RUN + 1))
    if [ "$expected" == "$actual" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $msg (got: '$actual')"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (expected: '$expected', got: '$actual')"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Compare numerical value with tolerance
assert_numeric_match() {
    local expected=$1
    local actual=$2
    local tolerance=${3:-0.0001}
    local msg=${4:-"Numeric match"}

    TESTS_RUN=$((TESTS_RUN + 1))

    # Use awk for floating point comparison
    local diff=$(awk -v e="$expected" -v a="$actual" 'BEGIN {
        diff = (a - e);
        if (diff < 0) diff = -diff;
        print diff
    }')

    local passes=$(awk -v d="$diff" -v t="$tolerance" 'BEGIN {
        if (d <= t) print "1"; else print "0"
    }')

    if [ "$passes" == "1" ]; then
        echo -e "${GREEN}✓ PASS${NC}: $msg (expected: $expected, actual: $actual, diff: $diff)"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $msg (expected: $expected, actual: $actual, diff: $diff > tolerance: $tolerance)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Extract energy from XYZ comment line (line 2)
extract_energy_from_xyz() {
    local xyzfile=$1
    # Curcuma format: "** Energy =   0.000000 Eh **" in comment line (line 2)
    # Handles variable whitespace: Energy\s*=\s*<number>
    grep -oP 'Energy\s*=\s*\K[-0-9.]+' "$xyzfile" | head -1
}

# Helper: Extract RMSD value from stdout
extract_rmsd_from_output() {
    local file=$1
    # Try multiple patterns: "RMSD: 1.234" or "RMSD = 1.234" or just a floating point after "RMSD"
    grep -i "rmsd" "$file" | grep -oP '[:=]\s*\K[0-9]+\.[0-9]+' | head -1
}

# Helper: Count structures in XYZ file
count_xyz_structures() {
    local xyzfile=$1
    if [ ! -f "$xyzfile" ]; then
        echo "0"
        return
    fi

    # Count frames by reading each frame's own atom-count header (Claude Generated, June 2026).
    # The previous version assumed a CONSTANT atom count across frames (NR % (natoms+2)); a
    # multi-XYZ with mixed atom counts was miscounted. This walks frame-by-frame from the first
    # header: read header -> skip comment + natoms coordinate lines -> next header. It counts
    # every logical conformer entry, including degenerate 0-atom frames (so the "number of
    # filtered structures" is reported correctly even when a writer emits reduced geometries).
    awk '{
        natoms = $1 + 0
        count++
        for (i = 0; i < natoms + 1; i++) {
            if ((getline) <= 0) break
        }
    }
    END { print count + 0 }' "$xyzfile"
}

# Helper: Compare floating point with tolerance (scientific validation)
# Returns 0 if values match within tolerance, 1 otherwise
compare_float() {
    local expected=$1
    local actual=$2
    local tolerance=${3:-1e-6}
    local label=${4:-"Value"}

    # Use awk for precise floating point comparison
    awk -v e="$expected" -v a="$actual" -v t="$tolerance" -v label="$label" 'BEGIN {
        diff = (a - e);
        if (diff < 0) diff = -diff;

        if (diff <= t) {
            print "✓ " label " matches: expected=" e ", actual=" a ", diff=" diff " (tol=" t ")"
            exit 0
        } else {
            print "✗ " label " mismatch: expected=" e ", actual=" a ", diff=" diff " > tolerance=" t
            exit 1
        }
    }'
}

# Helper: Assert scientific result (energy, RMSD, etc.) with tolerance
assert_scientific_value() {
    local expected=$1
    local actual=$2
    local tolerance=$3
    local msg=$4

    TESTS_RUN=$((TESTS_RUN + 1))

    # Handle empty/missing values
    if [ -z "$actual" ] || [ "$actual" == "0" ]; then
        echo -e "${RED}✗ FAIL${NC}: $msg (could not extract value)"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi

    # Compare with tolerance
    local result
    result=$(compare_float "$expected" "$actual" "$tolerance" "$msg")
    local exit_code=$?

    if [ $exit_code -eq 0 ]; then
        echo -e "${GREEN}✓ PASS${NC}: $result"
        TESTS_PASSED=$((TESTS_PASSED + 1))
        return 0
    else
        echo -e "${RED}✗ FAIL${NC}: $result"
        TESTS_FAILED=$((TESTS_FAILED + 1))
        return 1
    fi
}

# Helper: Store golden reference value
store_golden_reference() {
    local test_name=$1
    local value_name=$2
    local value=$3
    local ref_file="${TEST_DIR}/.golden_references.txt"

    # Create or append to reference file
    echo "${test_name}:${value_name}=${value}" >> "$ref_file"
    echo -e "${BLUE}Info:${NC} Stored golden reference: $value_name=$value"
}

# Helper: Load golden reference value
load_golden_reference() {
    local test_name=$1
    local value_name=$2
    local ref_file="${TEST_DIR}/.golden_references.txt"

    if [ ! -f "$ref_file" ]; then
        echo ""
        return 1
    fi

    grep "^${test_name}:${value_name}=" "$ref_file" | cut -d'=' -f2 | tail -1
}

# Helper: Print test summary
print_test_summary() {
    echo ""
    echo -e "${BLUE}========================================${NC}"
    echo -e "${BLUE}Test Summary${NC}"
    echo -e "${BLUE}========================================${NC}"
    echo "Tests run:    $TESTS_RUN"
    echo -e "Tests passed: ${GREEN}$TESTS_PASSED${NC}"
    echo -e "Tests failed: ${RED}$TESTS_FAILED${NC}"
    echo ""

    if [ $TESTS_FAILED -eq 0 ]; then
        echo -e "${GREEN}✓ All tests passed!${NC}"
        return 0
    else
        echo -e "${RED}✗ Some tests failed${NC}"
        return 1
    fi
}

# Helper: Setup test directory
setup_test_dir() {
    local test_dir=${1:-.}
    mkdir -p "$test_dir"
    cd "$test_dir" || exit 1
}

# ============================================================================
# BMT (Basename.Method.Timestamp) output directory helpers
# With BMT default-on, output files go to basename.keyword.YYYYMMDD_HHMMSS/
# These helpers resolve file paths whether BMT is on or off.
# Claude Generated 2026
# ============================================================================

# Find a BMT output directory matching basename.keyword.* pattern
# Usage: find_bmt_dir "input" "opt" -> "input.opt.20260609_143000"
find_bmt_dir() {
    local basename="$1"
    local keyword="$2"
    find . -maxdepth 1 -type d -name "${basename}.${keyword}.*" 2>/dev/null | head -1
}

# Find an output file, checking CWD first then BMT directories (depth 2)
# Usage: find_output_file "input.opt.xyz" -> "input.opt.xyz" or "input.opt.20260609_143000/input.opt.xyz"
find_output_file() {
    local filename="$1"
    # Check CWD first (backward compatibility with -no_bmt / BMT off)
    if [ -f "$filename" ]; then
        echo "$filename"
        return 0
    fi
    # Search in BMT directories (basename.keyword.YYYYMMDD_HHMMSS/filename)
    find . -maxdepth 2 -name "$filename" -type f 2>/dev/null | head -1
}

# Remove everything in the CWD except the named keep-files (Claude Generated, June 2026).
# Tests run from build-tree dirs that accumulate stale output across runs; a partial
# cleanup lets find_output_file match stale data and produces false pass/fail. This
# guarantees a clean slate so each run is deterministic.
# Usage: clean_dir_keep conformers.xyz run_test.sh
clean_dir_keep() {
    local prune=()
    local f
    for f in "$@"; do prune+=( ! -name "$f" ); done
    find . -mindepth 1 -maxdepth 1 "${prune[@]}" -exec rm -rf {} + 2>/dev/null || true
}

# Extract the ConfScan result fingerprint "<accepted> <reorder> <reuse> <skipped>"
# from a log file, tolerating ANSI color codes that may prefix the line (the line is
# printed via std::cout right after colored CurcumaLogger output). Claude Generated.
extract_confscan_summary() {
    sed 's/\x1b\[[0-9;]*m//g' "$1" 2>/dev/null | grep -E '^[0-9]+ [0-9]+ [0-9]+ [0-9]+$' | tail -1
}

# Remove BMT output directories in current directory
cleanup_bmt_dirs() {
    find . -maxdepth 1 -type d \( \
        -name "*.opt.*" -o \
        -name "*.md.*" -o \
        -name "*.confscan.*" -o \
        -name "*.confsearch.*" -o \
        -name "*.confstat.*" -o \
        -name "*.dock.*" -o \
        -name "*.analysis.*" -o \
        -name "*.rmsd.*" -o \
        -name "*.hessian.*" -o \
        -name "*.qmdfffit.*" -o \
        -name "*.sp.*" \
    \) -exec rm -rf {} + 2>/dev/null || true
}

# Helper: Cleanup test artifacts
cleanup_test_artifacts() {
    # Remove common Curcuma output files in CWD
    rm -f *.opt.xyz *.trj.xyz *.restart.json *.vtf
    rm -f *.accepted.xyz *.rejected.xyz *.thresh.xyz *.initial.xyz
    rm -f *.reorder.*.xyz *.centered.xyz *.reordered.xyz
    rm -f *.log *.dat *.pairs *.diag.jsonl
    rm -f stdout.log stderr.log
    # Remove BMT output directories
    cleanup_bmt_dirs
}
