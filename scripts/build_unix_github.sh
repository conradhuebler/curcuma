#!/bin/bash
set -ex

export CXX="g++-12"
export CC="gcc-12"
export FC="gfortran-12"

mkdir -p release_xtb
cd release_xtb

cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. \
    -DUSE_XTB=false \
    -DUSE_TBLITE=false \
    -DUSE_D3=false \
    -DUSE_D4=false \
    -DCMAKE_INSTALL_PREFIX=~/curcuma_xtb

make -j$(nproc)
make install
cd ..
