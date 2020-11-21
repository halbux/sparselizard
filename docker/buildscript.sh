#!/usr/bin/env bash

COLOR_NC='\e[0m' # No Color
COLOR_RED='\e[0;31m'
BUILD_DIR=docker_build

error_exit(){
    echo -e ${COLOR_RED}Error: $1${COLOR_NC}
    exit 1
}

cmake -S . -B ${BUILD_DIR} -DCMAKE_BUILD_TYPE=Debug -DBUILD_EXAMPLES=YES \
    && cmake --build ${BUILD_DIR} -j $(nproc) \
    && for f in $(find ${BUILD_DIR}/examples -type f -executable); do
        echo -e ${COLOR_RED}Running $f${COLOR_NC}
        pushd $(dirname $f) \
            && ./$(basename $f) || error_exit $f
        popd
    done \
        && echo -e ${COLOR_RED}Ok${COLOR_NC}
