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
git clone -b maint https://gitlab.com/petsc/petsc.git petsc;


########## CONFIGURE PETSC (SELECT THE APPROPRIATE CONFIGURATION OPTIONS BELOW) :

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

# The configuration below does not add support for additional mesh formats but does not require mpi.
./configure --with-mpi=0 --with-mumps-serial=1 --download-mumps --download-openblas --download-slepc --with-debugging=0 --with-scalar-type=real COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3';

# The configuration below adds support for .exo and .med mesh formats (mpi is needed and it is therefore added to the configuration options).
# Support for cgns can be added by manually installing cgns then providing the cgns folder to petsc.
#
# ADDITIONAL STEPS TO PERFORM FOR THIS CONFIGURATION OF PETSC:
# --> INSTALL CMAKE, AUTOTOOLS AND AUTOCONF BEFORE RUNNING THE CONFIGURE COMMAND BELOW  (on Ubuntu type: sudo apt-get install cmake autotools-dev autoconf)
# --> IN THE MAKEFILE REMOVE '-I ~/SLlibs/petsc/include/petsc/mpiuni' (otherwise the wrong mpi.h header is selected during make) 
#
#./configure --download-mpich --download-mumps --download-scalapack --download-openblas --download-slepc --download-med --download-hdf5 --download-zlib --download-netcdf --download-pnetcdf --download-exodusii --with-scalar-type=real --with-debugging=0 COPTFLAGS='-O3' CXXOPTFLAGS='-O3' FOPTFLAGS='-O3';


########## COMPILE PETSC :

echo '__________________________________________';
echo 'COMPILING PETSC';
make $PETSC_DIR $PETSC_ARCH all;
make $PETSC_DIR $PETSC_ARCH test;    

cd ..;


