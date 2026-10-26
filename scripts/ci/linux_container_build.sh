#!/bin/bash
# Build one Linux variant inside an AlmaLinux 8 container (glibc 2.28, so the result runs
# on current distributions) and pack it into a tar.gz that carries its own libraries.
#   linux_container_build.sh gpu   AVX512 + CUDA + ROCm plugins   (image rocm/dev-almalinux-8)
#   linux_container_build.sh cpu   AVX2, no AVX512, no GPU plugins (image almalinux:8)
# Run with the repository mounted at /src and an output directory mounted at /out, see
# .github/workflows/ci-variants.yml. Claude Generated; only run on GitHub Actions, never locally.
set -euo pipefail
VARIANT="${1:?gpu|cpu}"
SRC=/src
BUILD=/tmp/build
ROCM_ARCHS="gfx1100;gfx1101;gfx1102;gfx1150;gfx1151;gfx1200;gfx1201"

dnf -y install dnf-plugins-core epel-release
dnf config-manager --set-enabled powertools || dnf config-manager --set-enabled crb
dnf -y install gcc-toolset-12-gcc-c++ gcc-toolset-12-gcc-gfortran make git wget patchelf \
    python3-pip openblas-devel lapack-devel
pip3 install cmake
# shellcheck disable=SC1091
source /opt/rh/gcc-toolset-12/enable

if [ "$VARIANT" = gpu ]; then
    dnf config-manager --add-repo https://developer.download.nvidia.com/compute/cuda/repos/rhel8/x86_64/cuda-rhel8.repo
    dnf -y install cuda-toolkit-12-8
    export PATH=/usr/local/cuda/bin:/opt/rocm/bin:$PATH
    CMAKE_VARIANT=(-DUSE_AVX512=ON -DUSE_CUDA=ON -DUSE_ROCM=ON
                   "-DROCM_GPU_ARCH=${ROCM_ARCHS}" -DCMAKE_PREFIX_PATH=/opt/rocm)
else
    CMAKE_VARIANT=(-DUSE_AVX2=ON)
fi

mkdir -p "$BUILD"
cd "$BUILD"
cmake "$SRC" -DCMAKE_BUILD_TYPE=Release "${CMAKE_VARIANT[@]}" \
    -DUSE_XTB=OFF -DUSE_TBLITE=OFF -DUSE_D3=OFF -DUSE_D4=ON -DUSE_Plumed=OFF > cmake.log 2>&1 \
    || { tail -n 60 cmake.log; exit 1; }
make -j"$(nproc)" > make.log 2>&1 || { tail -n 80 make.log; exit 1; }
tail -n 3 make.log

bash "$SRC/scripts/ci/package_linux.sh" "$BUILD" "/out/curcuma-$VARIANT"
tar -C /out -czf "/out/curcuma-linux-$VARIANT.tar.gz" "curcuma-$VARIANT"
ls -l /out
