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
# Create symlink in case libblas is called instead of libopenblas:
cd install/lib;
ln -s libopenblas.so libblas.so;
ln -s libopenblas.so libblas.so.0;
ln -s libopenblas.so libblas.so.1;
ln -s libopenblas.so libblas.so.2;
ln -s libopenblas.so libblas.so.3;
ln -s libopenblas.so libblas.so.4;
ln -s libopenblas.so libblas.so.5;
ln -s libopenblas.a libblas.a;
cd ../../..;

