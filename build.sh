#!/usr/bin/env bash

clitool="cmake"
version=3.23.1
cmd="$clitool .."
dockerimage="rikorose/gcc-cmake:gcc-9"

if [ ! -r ./build ]; then
    mkdir build
fi

if ! $clitool --version | grep $version > /dev/null; then
    echo "$clitool $version not installed"
    echo "running build via docker..."
#    docker run --rm -v "$PWD:/src" -w "/src" $dockerimage bash -c 'cd build && $cmd && make && ctest -V'
    docker run --rm -v "$PWD:/src" -w "/src" $dockerimage bash -c 'ls && cd build && ls'
else
    echo "$clitool $version installed!"
    echo "running build..."
    cd build
    $cmd
    make
    ctest -V
fi
