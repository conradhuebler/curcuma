# Curcuma Error Messages & Solutions

A comprehensive guide to common Curcuma errors and how to fix them.

## File Errors

### "File not found: <filename>"

**Cause**: The input file does not exist at the specified path

**Solution**:
1. Check the file path spelling carefully
2. Verify the file exists: `ls <filename>`
3. Use absolute paths if working from different directories
4. Check file extension matches actual format (`.xyz`, `.mol2`, etc.)

**Example**:
```bash
# Wrong
./curcuma -sp my_molecule.xy -method uff  # Wrong extension!

# Correct
./curcuma -sp my_molecule.xyz -method uff
```

### "Failed to open file: <filename> (check permissions)"

**Cause**: File exists but cannot be read

**Solution**:
1. Check file permissions: `ls -l <filename>`
2. Make file readable: `chmod 644 <filename>`
3. Check if file is locked by another process
4. Ensure you have read access to the directory

**Example**:
```bash
# Check permissions
ls -l my_molecule.xyz
# -rw-r--r-- (good)
# ---------- (bad - no read permissions)

# Fix permissions
chmod 644 my_molecule.xyz
```

### "Invalid XYZ format in <filename> (line 1): Expected atom count (integer)"

**Cause**: First line of XYZ file doesn't contain a valid number of atoms

**Solution**:
1. Check the first line of your file: `head -1 <filename>`
2. First line must be a single integer (number of atoms)
3. Make sure there's no extra text or whitespace
4. Verify file is actual XYZ format, not another format with wrong extension

**Valid XYZ Format**:
```
2
Comment line (can be empty)
H 0.0 0.0 0.0
H 1.0 0.0 0.0
```

**Invalid Examples**:
```
H2              # ✗ Not a number
two atoms       # ✗ Text instead of number
2 atoms         # ✗ Extra text
```

### "Invalid coordinate in <filename> (line <N>): Expected 'Element X Y Z'"

**Cause**: Coordinate line has wrong format or invalid numbers

**Solution**:
1. Check line N of your file: `sed -n '<N>p' <filename>`
2. Each coordinate line must have exactly 4 columns: `Element X Y Z`
3. Element must be valid chemical symbol (H, C, N, O, etc.)
4. X, Y, Z must be valid numbers (not text like "inf", "nan", etc.)
5. Watch for missing columns or extra whitespace

**Valid Coordinate Lines**:
```
H 0.0 0.0 0.0
C 1.5 2.3 -0.5
N -1.2 0.8 1.1
```

**Invalid Coordinate Lines**:
```
H 0.0 0.0         # ✗ Missing Z coordinate
0.0 0.0 0.0 H     # ✗ Wrong order
H X Y Z           # ✗ Non-numeric values
H 0.0 0.0 inf     # ✗ Invalid number
```

## Method Errors

### "Unknown computational method: '<method>'"

**Cause**: Method name is not recognized or not available

**Solution**:
1. Check spelling carefully (case-insensitive)
2. Run `curcuma --methods` to see all available methods
3. Check method suggestions if typo is detected
4. Verify method is compiled into your build

**Example - Common Typos**:
```bash
# Typo detected - will show suggestions
./curcuma -sp input.xyz -method gfn3  # Did you mean: gfn2, gfn1?
./curcuma -sp input.xyz -method uf    # Did you mean: uff, uff-d3?
```

### "Did you mean one of these? - <suggestions>"

**What it means**: You entered an unknown method, but here are similar method names

**How to fix**:
1. Choose one of the suggested methods
2. Run `curcuma --methods` for full list
3. Read method documentation for correct spelling

**Example**:
```bash
$ ./curcuma -sp input.xyz -method gfn3
[ERROR] Unknown computational method: 'gfn3'
[ERROR] Did you mean one of these?
[ERROR]   - gfn2
[ERROR]   - gfn1
[ERROR]   - gfnff
```

### "Method '<method>' requires <provider> which was not compiled in this build"

**Cause**: Method exists but requires a library not included in compilation

**Solution**:
1. Find the CMake flag to enable (shown in error message)
2. Reconfigure CMake with the flag enabled
3. Rebuild Curcuma

**Example**:
```bash
# Error says: Method 'ipea1' requires TBLite which was not compiled
# Fix:
cd build
cmake .. -DUSE_TBLITE=ON
make -j4
```

**Common Compilation Flags**:
- `-DUSE_TBLITE=ON` - Tight-binding methods (gfn2, gfn1, ipea1)
- `-DUSE_XTB=ON` - Extended tight-binding (gfn-ff, gfn1, gfn2)
- `-DUSE_ULYSSES=ON` - Semi-empirical methods (pm6, am1, pm3, etc.)
- `-DUSE_D3=ON` - DFT-D3 dispersion corrections
- `-DUSE_D4=ON` - DFT-D4 dispersion corrections

## Parameter Errors

### "Parameter '<param>' must be positive, got <value>. Using default: <default>"

**Cause**: Parameter needs a positive value but got zero or negative

**Solution**:
1. Use a positive value instead
2. Common positive parameters: `stride`, `max_iterations`, `temperature`
3. See parameter documentation for allowed ranges

**Example**:
```bash
# Wrong - stride must be positive
./curcuma -rmsd ref.xyz target.xyz -stride -1

# Correct
./curcuma -rmsd ref.xyz target.xyz -stride 1
```

### "Parameter '<param>' must be in range [<min>, <max>], got <value>"

**Cause**: Value outside allowed range

**Solution**:
1. Use a value within the specified range
2. Common ranges:
   - Verbosity: 0-3
   - Threads: 1 to number of CPU cores
   - Temperature: positive values (K)
   - Convergence thresholds: typically 1e-6 to 1e-4

**Example**:
```bash
# Wrong - verbosity only 0-3
./curcuma -sp input.xyz -method uff -verbosity 5

# Correct
./curcuma -sp input.xyz -method uff -verbosity 2
```

### "Invalid value for '<param>': '<value>'. Valid choices: [<options>]"

**Cause**: Choice parameter got invalid option

**Solution**:
1. Use one of the valid options shown
2. Check documentation for parameter meanings
3. Values are case-insensitive usually

**Example**:
```bash
# Wrong - invalid optimizer
./curcuma -opt input.xyz -method uff -optimizer steepest_descent

# Correct
./curcuma -opt input.xyz -method uff -optimizer lbfgs
```

## Structure Errors

### "Failed to load molecular structures"

**This error typically includes more detail**. Look for additional error messages that specify which file(s) failed.

**Common causes**:
1. **Wrong file format** - Check file extension matches content
2. **Invalid coordinates** - Non-numeric coordinate values
3. **Missing atoms** - XYZ first line doesn't match actual atom count
4. **Encoding issues** - File saved in wrong text encoding

**Diagnostic steps**:
```bash
# 1. Check file format
file molecule.xyz

# 2. Check first few lines
head -5 molecule.xyz

# 3. Check for non-ASCII characters
file -i molecule.xyz  # Should show "text/plain; charset=us-ascii"

# 4. Count atoms vs declared count
wc -l molecule.xyz
head -1 molecule.xyz
# atom_count = (total_lines - 2) / num_structures
```

## Troubleshooting Checklist

### Before Running Calculation
- [ ] File exists and is readable: `ls -l <file>`
- [ ] File format matches extension
- [ ] First line of XYZ is a number
- [ ] Coordinate lines have 4 columns (Element X Y Z)
- [ ] Method is spelled correctly or listed in `curcuma --methods`
- [ ] All parameters are valid values

### If Calculation Hangs
- [ ] Press Ctrl+C to interrupt
- [ ] Check that input file actually exists
- [ ] Try with verbosity: `curcuma -sp input.xyz -method uff -verbosity 2`
- [ ] Check disk space with `df -h`
- [ ] Check system resources with `top`

### If Results Look Wrong
- [ ] Verify correct method was used: check output "Using method: ..."
- [ ] Check atom count: `head -1 input.xyz`
- [ ] Try with a smaller molecule first
- [ ] Verify input file is not corrupted

## Getting Help

### Built-in Help
```bash
# List available methods
curcuma --methods

# Get general help
curcuma -help

# Get help for specific capability
curcuma -help analysis

# View this error guide
# See docs/ERROR_MESSAGES.md
```

### Common Commands for Verification
```bash
# Check Curcuma version
curcuma --version

# List available methods by category
curcuma --methods

# Test with simple example (2 hydrogens)
cat > test.xyz << 'EOF'
2
Test
H 0.0 0.0 0.0
H 1.0 0.0 0.0
EOF
curcuma -sp test.xyz -method uff -verbosity 2
```

---

**Note**: All error messages include the parameter name and context to help you fix the problem quickly. Read the full error message carefully - it usually tells you exactly what needs to be fixed!
