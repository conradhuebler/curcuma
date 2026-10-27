#!/bin/bash
# Build the Linux GPU variant (AVX512, CUDA + ROCm plugins) directly on an ubuntu-22.04
# runner, without a container, and pack it with package_linux.sh. Claude Generated.
# CUDA from NVIDIA's apt repository, ROCm from repo.radeon.com; run on GitHub Actions only.
set -euo pipefail
SRC="$PWD"
BUILD="$SRC/build_gpu"
OUT="$SRC/out"
ROCM_ARCHS="gfx1100;gfx1150;gfx1200"
ROCM_REPO_VERSION="${ROCM_REPO_VERSION:-7.2.1}"

report() { echo "[$1] $(df -h / | awk 'NR==2{print "disk free: "$4}') $(free -m | awk 'NR==2{print "mem used MB: "$3}')"; }

report start
sudo apt-get update
sudo apt-get install -y g++-12 gfortran-12 cmake patchelf wget gnupg ca-certificates

# CUDA toolkit (nvcc + cudart + cublas + cusolver)
wget -q https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-keyring_1.1-1_all.deb
sudo dpkg -i cuda-keyring_1.1-1_all.deb
sudo apt-get update
sudo apt-get install -y --no-install-recommends cuda-nvcc-12-8 cuda-cudart-dev-12-8 \
    libcublas-dev-12-8 libcusolver-dev-12-8 cuda-nvrtc-dev-12-8 libcusparse-dev-12-8
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

export PATH=/usr/local/cuda-12.8/bin:/opt/rocm/bin:$PATH
mkdir -p "$BUILD"; cd "$BUILD"
CXX=g++-12 CC=gcc-12 FC=gfortran-12 cmake "$SRC" -DCMAKE_BUILD_TYPE=Release \
    -DUSE_AVX512=ON -DUSE_CUDA=ON -DUSE_ROCM=ON "-DROCM_GPU_ARCH=${ROCM_ARCHS}" \
    -DCMAKE_PREFIX_PATH=/opt/rocm -DCUDAToolkit_ROOT=/usr/local/cuda-12.8 \
    -DUSE_XTB=OFF -DUSE_TBLITE=OFF -DUSE_D3=OFF -DUSE_D4=OFF -DUSE_Plumed=OFF > cmake.log 2>&1 \
    || { tail -n 60 cmake.log; exit 1; }

# One make job per 4 GB of RAM (nvcc, hipcc and gfnff_method.cpp are memory hungry).
JOBS=$(awk -v n="$(nproc)" '/MemAvailable/{m=int($2/4000000); print (m<n?(m<1?1:m):n)}' /proc/meminfo)
( while sleep 120; do report build; done ) &
make -j"$JOBS" > make.log 2>&1 || { grep -n -B3 -A6 -m6 -E " error|Error [0-9]|undefined reference|cannot find|fatal" make.log | cut -c1-300; exit 1; }
tail -n 3 make.log

mkdir -p "$OUT"
bash "$SRC/scripts/ci/package_linux.sh" "$BUILD" "$OUT/curcuma-gpu"
tar -C "$OUT" -czf "$OUT/curcuma-linux-gpu.tar.gz" curcuma-gpu
ls -l "$OUT"
