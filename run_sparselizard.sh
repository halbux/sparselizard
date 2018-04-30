# THIS SCRIPT RUNS THE SPARSELIZARD EXECUTABLE


if [ "$(uname)" == "Linux" ]; then
LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux2-c-opt/lib":"$HOME/SLlibs/slepc/arch-linux2-c-opt/lib":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH    
elif [ "$(uname)" == "Darwin"  ]; then
DYLD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-darwin-c-opt/lib":"$HOME/SLlibs/slepc/arch-darwin-c-opt/lib":$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH
fi


./sparselizard;
