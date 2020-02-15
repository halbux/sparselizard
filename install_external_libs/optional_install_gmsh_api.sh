#!/bin/bash



########## ALL LIBRARIES REQUIRED BY SPARSELIZARD ARE PUT IN ~/SLlibs.

rm -rf ~/SLlibs/gmsh;
mkdir ~/SLlibs;
cd ~/SLlibs;


########## DOWNLOAD THE GMSH API :

echo '__________________________________________';
echo 'FETCHING THE GMSH API';
rm gmsh*.tgz;
wget http://gmsh.info/bin/Linux/gmsh-4.5.2-Linux64-sdk.tgz;
tar -xf *.tgz;
rm gmsh*.tgz;
mv gmsh* gmsh;

