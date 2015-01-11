@echo off
setlocal

set "SEAMASS_TOOLSET=intel32"
mkdir build
mkdir build\%SEAMASS_TOOLSET%

set "SEAMASS_BUILD=debug"
mkdir build\%SEAMASS_TOOLSET%\%SEAMASS_BUILD%
pushd build\%SEAMASS_TOOLSET%\%SEAMASS_BUILD%
cmake -G"Visual Studio 11 2012" -T"Intel C++ Compiler XE 15.0" -C"%~dp0..\seamass-windeps\build.cmake" %* "%~dp0"
if %errorlevel% neq 0 goto eof
popd

set "SEAMASS_BUILD=release"
mkdir build\%SEAMASS_TOOLSET%\%SEAMASS_BUILD%
pushd build\%SEAMASS_TOOLSET%\%SEAMASS_BUILD%
cmake -G"Visual Studio 11 2012" -T"Intel C++ Compiler XE 15.0" -C"%~dp0..\seamass-windeps\build.cmake" %* "%~dp0"
if %errorlevel% neq 0 goto eof

:eof
popd