#!/bin/bash

# create a clean build directory
rm -rf build
mkdir build && cd build

# configure the project
cmake .. \
    -DCMAKE_INSTALL_PREFIX=/workspace/programs/julia/dev/DACE_jll/override \
    -DCMAKE_PREFIX_PATH=/workspace/programs/julia/artifacts/5209ca23f516fb3391b885eef717e49b4ee0a268 \
    -DCMAKE_INSTALL_RPATH=/workspace/programs/julia/artifacts/5209ca23f516fb3391b885eef717e49b4ee0a268/lib:/workspace/programs/julia/juliaup/julia-1.10.2+0.x64.linux.gnu/lib \
    -DCMAKE_BUILD_TYPE=Release \
    -DWITH_PTHREAD=ON \
    -DWITH_ALGEBRAICMATRIX=ON \
    -DCMAKE_CXX_STANDARD=17 \
    -DWITH_JULIA=ON \
    -LA
