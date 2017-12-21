# !!! THE GFORTRAN, GCC, AND G++ COMPILERS MUST BE AVAILABLE !!!
#
# THIS SCRIPT INSTALLS IN ~/SLlibs THE FFTW LIBRARY.


#!/bin/bash


########## ALL LIBRARIES REQUIRED BY SPARSELIZARD ARE PUT IN ~/SLlibs.

rm -rf ~/SLlibs/fftw;
mkdir ~/SLlibs;
cd ~/SLlibs;


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


