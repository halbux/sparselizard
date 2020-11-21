#!/usr/bin/env bash

cmake -S . -B build -DCMAKE_BUILD_TYPE=Debug -DBUILD_EXAMPLES=NO \
    && cmake --build build -j $(nproc) \
    && cd build/simulations/default \
    && ./default

