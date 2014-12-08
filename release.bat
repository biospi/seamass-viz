@pushd "%~dp0"
@mkdir release
@pushd release
@cmake -C"..\..\seamass-windeps\release.cmake" %* ..
@popd
@popd