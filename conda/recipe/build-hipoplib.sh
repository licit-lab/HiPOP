#!/usr/bin/env bash
set -ex

cd cpp

cmake \
    -DCMAKE_PREFIX_PATH=$PREFIX\
    -DCMAKE_INSTALL_PREFIX=$PREFIX\
    -DBUILD_TESTS=ON\
    -DCMAKE_BUILD_TYPE=Release\
    -B build

cmake --build build --parallel "${CPU_COUNT}"
cmake --install build --config Release

cd build && ctest