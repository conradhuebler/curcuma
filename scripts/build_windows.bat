echo on

SET project_dir="%cd%"
echo Building curcuma using MinGW ...
git submodule update --init --recursive
git pull --recurse-submodules
mkdir build_windows
cd build_windows

cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release .. -DUSE_TBLITE=OFF -DUSE_XTB=OFF -DUSE_Plumed=OFF
cmake --build . --config Release

mkdir -p ../bin
copy curcuma.exe ../bin
REM Kopiere die benötigten DLL-Dateien in den bin-Ordner
copy %MINGW_DIR%\bin\libgcc_s_seh-1.dll ../bin
copy %MINGW_DIR%\bin\libgfortran-5.dll ../bin
copy %MINGW_DIR%\bin\libgomp-1.dll ../bin
copy %MINGW_DIR%\bin\libquadmath-0.dll ../bin
copy %MINGW_DIR%\bin\libstdc++-6.dll ../bin
copy %MINGW_DIR%\bin\libwinpthread-1.dll ../bin
cd ..

cd bin

REM Packe den Inhalt des bin-Verzeichnisses in eine ZIP-Datei
REM powershell Compress-Archive -Path * -DestinationPath ../curcuma-windows.zip

echo Done! The executable and required DLLs are in the bin folder, and the ZIP file is created.