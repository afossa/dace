#!/bin/bash

# create a clean build directory
rm -rf build
mkdir build && cd build

# get the Julia root path
JULIA_ROOT_PATH=$1

# configure the project
if [ -z "$1" ]; then
    echo "Building dace without Julia interface"

    cmake .. \
    -DCMAKE_INSTALL_PREFIX="${HOME}/.local" \
    -DCMAKE_BUILD_TYPE=Release \
    -DWITH_PTHREAD=ON \
    -DWITH_ALGEBRAICMATRIX=ON \
    -DCMAKE_CXX_STANDARD=20 \
    -DWITH_JULIA=OFF \
    -LA
else
    echo "Building dace with Julia interface"

    CXX_WRAP_PREFIX_PATH="/artifacts/6d40e7ea7d9da1197bf07d3970257e8a9260fe32"

    cmake .. \
        -DCMAKE_INSTALL_PREFIX="${JULIA_ROOT_PATH}/dev/DACE_jll/override" \
        -DCMAKE_PREFIX_PATH="${JULIA_ROOT_PATH}${CXX_WRAP_PREFIX_PATH}" \
        -DCMAKE_INSTALL_RPATH="${JULIA_ROOT_PATH}${CXX_WRAP_PREFIX_PATH}/lib:${JULIA_ROOT_PATH}/juliaup/julia-1.10.4+0.x64.linux.gnu/lib" \
        -DCMAKE_BUILD_TYPE=Release \
        -DWITH_PTHREAD=ON \
        -DWITH_ALGEBRAICMATRIX=ON \
        -DCMAKE_CXX_STANDARD=20 \
        -DWITH_JULIA=ON \
        -LA
fi
