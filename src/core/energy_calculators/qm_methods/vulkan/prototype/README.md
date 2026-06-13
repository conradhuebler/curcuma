# Vulkan Jacobi eigensolver — standalone validation harness

`jacobi_validation.cpp` is the standalone prototype used to validate the FP64 two-sided
cyclic Jacobi symmetric eigensolver (the shaders in `../shaders/`) against Eigen, on a real
Vulkan device, before it was integrated into the native-xTB Vulkan path. Kept as a
re-runnable correctness reference; not built as part of curcuma.

## Build & run

```bash
# 1) compile the shaders to SPIR-V next to the sources (uses glslc from the Vulkan SDK)
cd ../shaders
for s in angles col row vec; do glslc --target-env=vulkan1.1 -fshader-stage=compute $s.comp -o /tmp/vkjac_$s.spv; done
# (the harness loads /tmp/vkjac/shaders/<name>.spv — adjust the paths in the .cpp or copy there)

# 2) build & run the harness (Eigen from external/)
g++ -O2 -std=c++17 -I<curcuma>/external/eigen-3.4.0 jacobi_validation.cpp -o jac -lvulkan
./jac            # auto sweep count per size
./jac 20         # force 20 sweeps
```

## Reference result (AMD Radeon 890M, RADV, shaderFloat64)

```
n=   4   |eps_err|=4.4e-16  recon=4.4e-16  ortho=4.4e-16  OK
n=  24   |eps_err|=2.0e-14  recon=1.4e-14  ortho=7.3e-15  OK
n= 100   |eps_err|=2.7e-13  recon=5.0e-14  ortho=2.4e-14  OK
n= 128   |eps_err|=3.0e-13  recon=6.3e-14  ortho=3.1e-14  OK
```

Eigenvalues vs Eigen `SelfAdjointEigenSolver` to ~1e-13, reconstruction / orthogonality to
~1e-14 — well inside the 1e-8 Eh SCF gate.
