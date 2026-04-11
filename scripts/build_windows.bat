@echo off
setlocal enabledelayedexpansion

echo ==============================================================================
echo Building curcuma for Windows using MinGW
echo ==============================================================================

set project_dir=%cd%

REM 1. FetchContent handles dependencies (no submodules needed)
echo Configuring dependencies via CMake FetchContent...

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

REM 7. Dynamic DLL Collection
echo Collecting required MinGW DLLs...
powershell -File scripts\collect_dlls.ps1

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