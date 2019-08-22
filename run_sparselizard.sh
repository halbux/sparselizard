#!/bin/bash

# THIS SCRIPT RUNS THE SPARSELIZARD EXECUTABLE

SCRIPT=$(readlink -f "$0")
SCRIPTDIR=$(dirname "$SCRIPT")
EXTERNALLIBSDIR=$SCRIPTDIR/external_libs/libs

if [ "$(uname)" == "Linux" ]; then
LD_LIBRARY_PATH="$EXTERNALLIBSDIR/petsc/arch-linux2-c-opt/lib":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH
elif [ "$(uname)" == "Darwin"  ]; then
DYLD_LIBRARY_PATH="$EXTERNALLIBSDIR/petsc/arch-darwin-c-opt/lib":$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH
fi


./sparselizard;
