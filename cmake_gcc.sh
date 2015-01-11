#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                   

mkdir $DIR/build
mkdir $DIR/build/gcc

mkdir $DIR/build/gcc/debug
pushd $DIR/build/gcc/debug
cmake -DCMAKE_BUILD_TYPE=Debug $@ ../../..
popd

mkdir $DIR/build/gcc/release
pushd $DIR/build/gcc/release
cmake -DCMAKE_BUILD_TYPE=Release $@ ../../..
popd
