# !!! THE GFORTRAN, GCC, AND G++ COMPILERS MUST BE AVAILABLE !!!
#
# THIS SCRIPT INSTALLS IN ~/SLlibs THE OPENBLAS LIBRARY.


#!/bin/bash


########## ALL LIBRARIES REQUIRED BY SPARSELIZARD ARE PUT IN ~/SLlibs.

rm -rf ~/SLlibs/openblas;
mkdir ~/SLlibs;
cd ~/SLlibs;


########## COMPILING OPENBLAS :

echo '__________________________________________';
echo 'FETCHING THE LATEST OPENBLAS VERSION FROM GIT';
git clone https://github.com/xianyi/OpenBLAS.git openblas;
echo '__________________________________________';
echo 'COMPILING OPENBLAS';
cd openblas;
make -j4;
mkdir install;
make PREFIX=$(pwd)/install install;
cd ..;

