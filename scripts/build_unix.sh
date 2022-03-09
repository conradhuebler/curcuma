#!/bin/bash
set -ex

#export CXX="g++-9"
#export CC="gcc-9"
git submodule init
git submodule update --recursive
# check submodules, seems not to work automatically

cd external
for i in $(ls -d */|grep -v 'xtb'); do cd $i; git checkout master || true; git submodule init; git submodule update --recursive; git pull; cd ..; done
cd xtb
git checkout v6.4.1 
git pull
cd ..
cd ..

mkdir -p release
cd release
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make 
