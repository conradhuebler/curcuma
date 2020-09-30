#!/bin/bash
set -ex

mkdir -p package
cp misc/* package
cd package
mkdir lib
cp ../release_xtb/curcuma .
cp ../release_xtb/curcuma_helper .
cp ../release_xtb/libxtb.so.6 lib


wget -c -nv "https://github.com/probonopd/linuxdeployqt/releases/download/continuous/linuxdeployqt-continuous-x86_64.AppImage"
chmod a+x linuxdeployqt-continuous-x86_64.AppImage

./linuxdeployqt-continuous-x86_64.AppImage curcuma.desktop -appimage
cp curcuma*.AppImage curcuma-nightly-x86_64.AppImage
