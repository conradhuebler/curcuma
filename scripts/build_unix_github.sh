#!/bin/bash
set -ex

export CXX="g++-9"
export CC="gcc-9"
export FC="gfortran-9"
git submodule init
git submodule update --recursive
# check submodules, seems not to work automatically

cd external

for i in $(ls -d */|grep -v 'xtb'); do cd $i; git checkout master || true ; git submodule init; git submodule update --recursive; git pull; cd ..; done

cd xtb
git checkout main
git pull
cd ..
cd ..

mkdir -p release_no_xtb
cd release_no_xtb
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. -DCMAKE_INSTALL_PREFIX=~/curcuma_no_xtb
make 
make install
cd ..

mkdir -p release_xtb
cd release_xtb
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo .. -DCOMPILE_XTB=true -DCMAKE_INSTALL_PREFIX=~/curcuma_xtb
make 
make install
cd ..



