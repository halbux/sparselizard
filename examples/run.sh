#!/bin/bash
#
# Simple wrapper to set library paths without admin
#
# ~~~
# How to use this file
#
# > cd [example]
# > make -f ../Makefile.example
# > bash ../run.sh
# ~~~


# LOAD REQUIRED LIBRARIES:

if [ "$(uname)" == "Linux" ]; then
LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux2-c-opt/lib:$HOME/Desktop/sparselizard":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
elif [ "$(uname)" == "Darwin"  ]; then
DYLD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-darwin-c-opt/lib":$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH
fi


# Run executable
./main
