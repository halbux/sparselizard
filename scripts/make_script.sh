#!/bin/sh

GCC=${GCC_BINARY:-gcc-11}
GXX=${GXX_BINARY:-g++-11}
SLlibs=${PETSC_DIR:-/root/SLlibs}
PETSC_DIR=$SLlibs/petsc
echo "GCC_BINARY is set to ${GCC}"
echo "GXX_BINARY is set to ${GXX}"

pwd=$(pwd)
echo "$pwd"
echo $SLlibs
echo $PETSC_DIR

mkdir build
cd build
cmake -DCMAKE_C_COMPILER=$GCC -DCMAKE_CXX_COMPILER=$GXX -DPETSC_PATH=${PETSC_DIR} ..
cmake --build . -j$(nproc)

