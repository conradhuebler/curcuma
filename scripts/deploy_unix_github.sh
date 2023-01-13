#!/bin/bash
set -ex

mkdir -p package
cp misc/* package
cd package
mkdir -p lib
cp ../release_xtb/curcuma .
cp ../release_xtb/curcuma_helper .
cp ../release_xtb/external/xtb/libxtb.so* lib

cd ..
wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage
cd package


../linuxdeployqt-continuous-x86_64.AppImage curcuma.desktop -appimage
cp curcuma*.AppImage curcuma-nightly-x86_64.AppImage
