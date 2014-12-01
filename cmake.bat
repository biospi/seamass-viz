@pushd "%~dp0"
@mkdir build
@cd build
@cmake -C"..\..\seamass-windeps\build.cmake"  %* ..
@popd