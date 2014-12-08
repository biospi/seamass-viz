@pushd "%~dp0"
@mkdir debug
@pushd debug
@cmake -C"..\..\seamass-windeps\debug.cmake" %* ..
@popd
@popd