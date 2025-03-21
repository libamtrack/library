@echo off
setlocal enabledelayedexpansion

:: Define directories
set BUILD_DIR=build
set INSTALL_DIR=usr

:: Remove the previous installation directory
if exist %INSTALL_DIR% (
    rmdir /s /q %INSTALL_DIR%
    if exist %INSTALL_DIR% (
        echo [ERROR] Failed to remove %INSTALL_DIR%.
        exit /b 1
    )
)

:: Clean previous build artifacts
@REM cmake --build %BUILD_DIR% --target clean
@REM if %errorlevel% neq 0 (
@REM     echo [ERROR] Cleaning build directory failed.
@REM     exit /b %errorlevel%
@REM )

:: Configure CMake with the specified build and install directories
cmake -S . -B %BUILD_DIR% ^
    -DGSL_INCLUDE_DIR="vcpkg_installed/x64-windows/include" ^
    -DGSL_LIBRARY="vcpkg_installed/x64-windows/lib/gsl.lib" ^
    -DGSL_CBLAS_LIBRARY="vcpkg_installed/x64-windows/lib/gslcblas.lib" ^
    -DBUILD_EXAMPLES=OFF ^
    -G "Ninja" ^
    -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR%
if %errorlevel% neq 0 (
    echo [ERROR] CMake configuration failed.
    exit /b %errorlevel%
)

:: Build the project
cmake --build %BUILD_DIR% --parallel
if %errorlevel% neq 0 (
    echo [ERROR] Build process failed.
    exit /b %errorlevel%
)

:: Install the compiled files to the usr directory
cmake --install %BUILD_DIR%
if %errorlevel% neq 0 (
    echo [ERROR] Installation failed.
    exit /b %errorlevel%
)

echo [SUCCESS] Build and installation completed successfully.
exit /b 0
