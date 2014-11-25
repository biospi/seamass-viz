@pushd "%~dp0"
@mkdir build
@cd build
@cmake -G"Visual Studio 10 2010 Win64" -T"Intel C++ Compiler XE 14.0" -DBoost_USE_STATIC_LIBS=true -DBOOST_ROOT="..\..\seamass-windeps\boost_1_57_0" -DOSGEO4W_ROOT="..\..\seamass-windeps\spatialindex-src-1.8.5\stage" %* ..
@popd.