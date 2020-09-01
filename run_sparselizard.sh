#!/bin/bash

# THIS SCRIPT RUNS THE SPARSELIZARD EXECUTABLE


if [ "$(uname)" == "Linux" ]; then
# With or without the GMSH API:
if [ -d ~/SLlibs/gmsh ]; then
    LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux-c-opt/lib":"$HOME/SLlibs/gmsh/lib":$LD_LIBRARY_PATH
else
    LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux-c-opt/lib":$LD_LIBRARY_PATH
fi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH    
elif [ "$(uname)" == "Darwin"  ]; then
DYLD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-darwin-c-opt/lib":$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH
fi


./sparselizard $@;
