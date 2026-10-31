#!/bin/bash
# Build the Linux GPU variant (ISA=avx512 (default) or avx2; CUDA + ROCm plugins) directly on an ubuntu-22.04
# runner, without a container, and pack it with package_linux.sh. Claude Generated.
# Usage: linux_host_gpu_build.sh configure | make <target> | package  (one workflow step each,
# so a lost runner shows which part was running).
# CUDA from NVIDIA's apt repository, ROCm from repo.radeon.com; run on GitHub Actions only.
set -euo pipefail
SRC="$PWD"
BUILD="$SRC/build_gpu"
OUT="$SRC/out"
ROCM_ARCHS="gfx1100;gfx1150;gfx1200"
ISA="${ISA:-avx512}"
# GPU=0: CPU-only build, no CUDA/ROCm toolchains installed, no plugins.
GPU="${GPU:-1}"
case "$ISA" in avx512) ISA_FLAG=-DUSE_AVX512=ON;; avx2) ISA_FLAG=-DUSE_AVX2=ON;; *) echo "ISA must be avx512 or avx2"; exit 1;; esac
ROCM_REPO_VERSION="${ROCM_REPO_VERSION:-7.2.1}"

report() { echo "[$1] $(df -h / | awk 'NR==2{print "disk free: "$4}') $(free -m | awk 'NR==2{print "mem used MB: "$3}')"; }

export PATH=/usr/local/cuda-12.8/bin:/opt/rocm/bin:$PATH
case "${1:?configure|make|package}" in
configure)
report start
sudo apt-get update
sudo apt-get install -y g++-12 gfortran-12 cmake patchelf wget gnupg ca-certificates

if [ "$GPU" = 1 ]; then
# CUDA toolkit (nvcc + cudart + cublas + cusolver)
wget -q https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get install -y --no-install-recommends cuda-nvcc-12-8 cuda-cudart-dev-12-8 \
    libcublas-dev-12-8 libcusolver-dev-12-8 cuda-nvrtc-dev-12-8 libcusparse-dev-12-8 libnvjitlink-dev-12-8
report cuda

# ROCm (hipcc, runtime, rocBLAS, rocSOLVER)
wget -qO- https://repo.radeon.com/rocm/rocm.gpg.key | sudo gpg --dearmor -o /etc/apt/keyrings/rocm.gpg
for v in "$ROCM_REPO_VERSION" latest; do
    echo "deb [arch=amd64 signed-by=/etc/apt/keyrings/rocm.gpg] https://repo.radeon.com/rocm/apt/$v jammy main" \
        | sudo tee /etc/apt/sources.list.d/rocm.list
    sudo apt-get update && break
done
# Ubuntu ships its own rocminfo, rocm-cmake and rocm-device-libs (version 5.0.0), which
# outrank the vendor's 7.x packages by version number. Priority above 1000 makes the vendor
# repository win.
printf 'Package: *\nPin: release o=repo.radeon.com\nPin-Priority: 1001\n' | sudo tee /etc/apt/preferences.d/rocm-pin
sudo apt-get update
# rocm-hip-sdk is the vendor's meta package (compiler, device libraries, runtime, rocBLAS,
# rocSOLVER). Ubuntu's own rocm-device-libs (5.0) must not be mixed in.
sudo apt-get install -y --no-install-recommends rocm-hip-sdk
BC=$(find /opt/rocm* -name ockl.bc 2>/dev/null | head -n 3)
echo "ROCm device library: ${BC:-NOT FOUND}"
[ -n "$BC" ] || { ls /opt; ls /opt/rocm* | head -n 30; exit 1; }
report rocm
fi

mkdir -p "$BUILD"; cd "$BUILD"
if [ "$GPU" = 1 ]; then
    GPU_FLAGS=(-DUSE_CUDA=ON -DUSE_ROCM=ON "-DROCM_GPU_ARCH=${ROCM_ARCHS}" -DCMAKE_PREFIX_PATH=/opt/rocm -DCUDAToolkit_ROOT=/usr/local/cuda-12.8)
else
    GPU_FLAGS=(-DUSE_CUDA=OFF -DUSE_ROCM=OFF)
fi
CXX=g++-12 CC=gcc-12 FC=gfortran-12 cmake "$SRC" -DCMAKE_BUILD_TYPE=Release \
    $ISA_FLAG "${GPU_FLAGS[@]}" \
    -DUSE_XTB=OFF -DUSE_TBLITE=OFF -DUSE_D3=OFF -DUSE_D4=OFF -DUSE_Plumed=OFF > cmake.log 2>&1 \
    || { tail -n 60 cmake.log; exit 1; }

;;
make)
    cd "$BUILD"
    # Two jobs at most, at low priority: nvcc, hipcc and gfnff_method.cpp need several GB each
    # and a full-speed -j4 starved the runner (lost communication).
    JOBS=$(awk -v n="$(nproc)" '/MemAvailable/{m=int($2/4000000); m=(m<n?m:n); print (m<2?(m<1?1:m):2)}' /proc/meminfo)
    ( while sleep 120; do report "make ${2:-all}"; done ) &
    nice -n 19 make -j"$JOBS" "${2:-all}" > "make_${2:-all}.log" 2>&1 \
        || { grep -n -B3 -A6 -m6 -E " error|Error [0-9]|undefined reference|cannot find|fatal" "make_${2:-all}.log" | cut -c1-300; exit 1; }
    tail -n 3 "make_${2:-all}.log";;
package)
    # ldd resolves the dependencies of the CUDA/ROCm libraries (e.g. libcusparse -> libnvJitLink)
    # only through the library path.
    export LD_LIBRARY_PATH=/usr/local/cuda-12.8/lib64:/opt/rocm/lib:${LD_LIBRARY_PATH:-}
mkdir -p "$OUT"
bash "$SRC/scripts/ci/package_linux.sh" "$BUILD" "$OUT" "$ISA"
;;
esac
