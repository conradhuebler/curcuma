#!/bin/bash
set -ex

export CXX="g++-9"
export CC="gcc-9"
export FC="gfortran-9"

git submodule init
git submodule update --recursive
# check submodules, seems not to work automatically

cd external

for i in $(ls -d */|grep -v 'tb' |grep -v 'cpp' |grep -v 'd3' ); do cd $i; git checkout master || true; git submodule init; git submodule update --recursive; git pull; cd ..; done

cd xtb
git checkout main
git pull
cd ..

cd tblite
git checkout main
git pull
cd ..

cd cpp-d4
git checkout main
git pull
cd ..

cd simple-dftd3
git checkout main
git pull
cd ..

cd ..


mkdir -p release_xtb
cd release_xtb

cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. -DUSE_XTB=true -DUSE_TBLITE=true -DUSE_D3=true -DUSE_D4=true -DCMAKE_INSTALL_PREFIX=~/curcuma_xtb -DLIBS_DIR=$(find /usr -name '*cblas.so' |head -n1 |sed 's/libcblas.so//g') -DINCLUDE_DIR=$(find /usr -name 'cblas.h' |head -n1 |sed 's/cblas.h//g')
make 
make install
cd ..



