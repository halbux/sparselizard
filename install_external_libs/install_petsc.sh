#!/bin/bash

# !!! THE GFORTRAN, GCC, AND G++ COMPILERS MUST BE AVAILABLE !!!
#
# THIS SCRIPT INSTALLS IN ~/SLlibs THE PETSC LIBRARY.



########## ALL LIBRARIES REQUIRED BY SPARSELIZARD ARE PUT IN ~/SLlibs.

rm -rf ~/SLlibs/petsc;
mkdir ~/SLlibs;
cd ~/SLlibs;


########## DOWNLOAD PETSC :

echo '__________________________________________';
echo 'FETCHING THE LATEST PETSC VERSION FROM GIT';
git clone -b main https://gitlab.com/petsc/petsc.git petsc;


########## CONFIGURE PETSC (SELECT THE APPROPRIATE CONFIGURATION OPTIONS BELOW) :

echo '__________________________________________';
echo 'CONFIGURING PETSC';
cd petsc;

if [ "$(uname)" == "Linux" ]; then
PETSC_DIR=$(pwd);
PETSC_ARCH=arch-linux-c-opt;
elif [ "$(uname)" == "Darwin"  ]; then
PETSC_DIR=$(pwd);
PETSC_ARCH=arch-darwin-c-opt;
fi

# Metis is recommended but not mandatory. It can provide a major speedup for MUMPS during resolution.
./configure --with-openmp --with-mpi=0 --with-shared-libraries=1 --with-mumps-serial=1 --download-mumps --download-openblas --download-metis --download-slepc --with-debugging=0 --with-scalar-type=real --with-x=0 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3';


########## COMPILE PETSC :

# In case petsc appends a 2 make a link:
if [ ! -d arch-linux-c-opt ]; then
    ln -s arch-linux2-c-opt arch-linux-c-opt;
fi

echo '__________________________________________';
echo 'COMPILING PETSC';
make $PETSC_DIR $PETSC_ARCH all;

cd ..;
