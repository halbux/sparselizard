# !!! THE GCC, AND G++ COMPILERS MUST BE AVAILABLE !!!
#
# THIS SCRIPT INSTALLS IN ~/SLlibs ALL EXTERNAL LIBRARIES REQUIRED.
# AFTER A SUCCESSFULL RUN ALL WHAT NEEDS TO BE DONE IS WRITE YOUR 
# main.cpp USING THE SPARSELIZARD LIBRARY THEN RUN make -j4.



#!/bin/bash


########## ALL LIBRARIES REQUIRED BY SPARSELIZARD ARE PUT IN ~/SLlibs.

rm -rf ~/SLlibs;
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


########## COMPILING FFTW :

echo '__________________________________________';
echo 'FETCHING FFTW';
mkdir fftw;
cd fftw;
fftwdir=$(pwd);
mkdir install;

if [ "$(uname)" == "Linux" ]; then
wget http://www.fftw.org/fftw-3.3.7.tar.gz;
elif [ "$(uname)" == "Darwin"  ]; then
curl http://www.fftw.org/fftw-3.3.7.tar.gz -o fftw-3.3.7.tar.gz;
fi

tar -xf *.tar.gz;
cd fftw*;
echo '__________________________________________';
echo 'CONFIGURING FFTW';
./configure --prefix $fftwdir/install;
echo '__________________________________________';
echo 'COMPILING FFTW';
make -j4;
make install;
cd ../..;


########## COMPILING PETSC WITH OPTIONS --download-mumps --download-scalapack --download-mpich --with-debugging=0 :

echo '__________________________________________';
echo 'FETCHING THE LATEST PETSC VERSION FROM GIT';
git clone -b maint https://bitbucket.org/petsc/petsc petsc;
echo '__________________________________________';
echo 'CONFIGURING PETSC';
cd petsc;

if [ "$(uname)" == "Linux" ]; then
PETSC_DIR=$(pwd);
PETSC_ARCH=arch-linux2-c-opt;
elif [ "$(uname)" == "Darwin"  ]; then
PETSC_DIR=$(pwd);
PETSC_ARCH=arch-darwin-c-opt;
fi

./configure --download-mumps --download-scalapack --download-mpich --with-blaslapack-dir=~/SLlibs/openblas/install --with-debugging=0;
echo '__________________________________________';
echo 'COMPILING PETSC';
make $PETSC_DIR $PETSC_ARCH all;
make $PETSC_DIR $PETSC_ARCH test;    

cd ..;

