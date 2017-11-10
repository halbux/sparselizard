# THIS SCRIPT RUNS THE SPARSELIZARD EXECUTABLE


LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux2-c-opt/lib":"$HOME/SLlibs/openblas/install/lib":"$HOME/SLlibs/fftw/install/lib":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH


./sparselizard;
