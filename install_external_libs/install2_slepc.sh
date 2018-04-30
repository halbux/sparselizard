# !!! THE GFORTRAN, GCC, AND G++ COMPILERS MUST BE AVAILABLE !!!
#
# THIS SCRIPT INSTALLS IN ~/SLlibs THE SLEPC LIBRARY.


#!/bin/bash


########## ALL LIBRARIES REQUIRED BY SPARSELIZARD ARE PUT IN ~/SLlibs.

rm -rf ~/SLlibs/slepc;
mkdir ~/SLlibs;
cd ~/SLlibs;


########## COMPILING SLEPC :

echo '__________________________________________';
echo 'FETCHING SLEPC';
if [ "$(uname)" == "Linux" ]; then
wget http://slepc.upv.es/download/distrib/slepc-3.8.1.tar.gz;
elif [ "$(uname)" == "Darwin"  ]; then
curl http://slepc.upv.es/download/distrib/slepc-3.8.1.tar.gz -o slepc.tar.gz;
fi
tar -xf slepc*.tar.gz;
rm slepc*.tar.gz;
mv slepc* slepc;
cd slepc;

echo '__________________________________________';
echo 'CONFIGURING SLEPC';

export SLEPC_DIR=$(pwd);
if [ "$(uname)" == "Linux" ]; then
export PETSC_DIR=~/SLlibs/petsc;
export PETSC_ARCH=arch-linux2-c-opt;
elif [ "$(uname)" == "Darwin"  ]; then
export PETSC_DIR=~/SLlibs/petsc;
export PETSC_ARCH=arch-darwin-c-opt;
fi

./configure;
echo '__________________________________________';
echo 'COMPILING SLEPC';
make SLEPC_DIR=$PWD;
make test; 

cd ..;



