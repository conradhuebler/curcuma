@echo off
setlocal enabledelayedexpansion

echo ==============================================================================
echo Building curcuma for Windows using MinGW
echo ==============================================================================

set project_dir=%cd%

REM 1. Update dependencies
echo Updating dependencies...
git submodule update --init --recursive

REM 2. Create and enter build directory
if not exist build_windows mkdir build_windows
cd build_windows

REM 3. Configure CMake
echo Configuring project...
cmake -G "MinGW Makefiles" -DCMAKE_BUILD_TYPE=Release .. -DUSE_TBLITE=OFF -DUSE_XTB=OFF -DUSE_Plumed=OFF
if %errorlevel% neq 0 (
    echo [ERROR] CMake configuration failed!
    exit /b %errorlevel%
)

REM 4. Build project
echo Building project...
cmake --build . --config Release
if %errorlevel% neq 0 (
    echo [ERROR] Build failed!
    exit /b %errorlevel%
)

REM 5. Prepare bin directory
cd ..
if not exist bin mkdir bin

REM 6. Copy executable
echo Copying executable...
copy /Y build_windows\curcuma.exe bin\curcuma.exe

REM 7. Dynamic DLL Collection (PowerShell)
echo Collecting required MinGW DLLs...
powershell -Command ^
    "$mingwBin = if ('$env:MINGW_DIR') { '$env:MINGW_DIR\bin' } else { 'C:\msys64\mingw64\bin' }; ^
     $dlls = Get-ChildItem -Path $mingwBin -Filter '*.dll' | Where-Object { ^
        $_.Name -match 'libgcc_s_seh' -or ^
        $_.Name -match 'libgfortran' -or ^
        $_.Name -match 'libgomp' -or ^
        $_.Name -match 'libquadmath' -or ^
        $_.Name -match 'libstdc\+\+' -or ^
        $_.Name -match 'libwinpthread' ^
     }; ^
     foreach ($dll in $dlls) { Copy-Item $dll.FullName -Destination 'bin\' -Force }"

if %errorlevel% neq 0 (
    echo [WARN] Some DLLs could not be copied, but continuing...
)

REM 8. Create ZIP archive
echo Creating distribution ZIP...
powershell -Command "Compress-Archive -Path bin\* -DestinationPath curcuma-windows.zip -Force"

echo.
echo ==============================================================================
echo Done!
echo Executable and DLLs are in: %project_dir%\bin
echo ZIP archive created: %project_dir%\curcuma-windows.zip
echo ==============================================================================
endlocal
