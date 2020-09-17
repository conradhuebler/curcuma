#!/bin/bash
set -ex

export CXX="g++-9"
export CC="gcc-9"
git submodule init
git submodule update --recursive
# check submodules, seems not to work automatically

cd external
for i in $(ls -d */); do cd $i; git checkout master; git submodule init; git submodule update --recursive; cd ..; done
cd ..

mkdir -p release
cd release
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..
make
