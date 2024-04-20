#!/bin/bash

# move to the build directory
cd build

# build and install the library
VERBOSE=ON cmake --build . --config Release --target install -- -j8
