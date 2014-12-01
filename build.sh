#!/bin/bash                                                     
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"                   
mkdir $DIR/build
pushd $DIR/build
cmake $@ ..
make
popd
