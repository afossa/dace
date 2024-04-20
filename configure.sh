#!/bin/bash

# create a clean build directory
rm -rf build
mkdir build && cd build

# get the installation and prefix paths
JULIA_ROOT_PATH=$1
CXX_WRAP_PREFIX_PATH="/artifacts/5209ca23f516fb3391b885eef717e49b4ee0a268"

# configure the project
cmake .. \
    -DCMAKE_INSTALL_PREFIX="${JULIA_ROOT_PATH}/dev/DACE_jll/override" \
    -DCMAKE_PREFIX_PATH="${JULIA_ROOT_PATH}${CXX_WRAP_PREFIX_PATH}" \
    -DCMAKE_INSTALL_RPATH="${JULIA_ROOT_PATH}${CXX_WRAP_PREFIX_PATH}/lib:${JULIA_ROOT_PATH}/juliaup/julia-1.10.2+0.x64.linux.gnu/lib" \
    -DCMAKE_BUILD_TYPE=Release \
    -DWITH_PTHREAD=ON \
    -DWITH_ALGEBRAICMATRIX=ON \
    -DCMAKE_CXX_STANDARD=17 \
    -DWITH_JULIA=ON \
    -LA
