#!/bin/bash

cd ..;


# LOAD REQUIRED LIBRARIES:

if [ "$(uname)" == "Linux" ]; then
# With or without the GMSH API:
if [ -d ~/SLlibs/gmsh ]; then
    LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux-c-opt/lib":"$HOME/SLlibs/gmsh/lib":$LD_LIBRARY_PATH
else
    LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux-c-opt/lib":$LD_LIBRARY_PATH
fi
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH    
elif [ "$(uname)" == "Darwin"  ]; then
DYLD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-darwin-c-opt/lib":$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH
fi


# COMPILE SPARSELIZARD:

make -j4;


# LOOP ON ALL EXAMPLES:

echo ''; 
for i in $(ls examples); 
do 

# Skip this script as well as all too heavy examples:
if [ $i == "validate.sh" ] || [ $i == "electromagnetic-waveguide-time-resolution-3d" ] || [ $i == "superconductor-3d" ] || [ $i == "nonlinear-vonkarman-vortex-2d" ] || [ $i == "thermoacoustic-elasticity-axisymmetry-2d" ] || [ $i == "nonlinear-natural-convection-hpfem-2d" ]
then
continue
fi

# Copy the example to the main directory: 
cp examples/$i/main.cpp ./;
cp examples/$i/*.msh ./ 2>/dev/null;
cp examples/$i/*.nas ./ 2>/dev/null;
cp examples/$i/*.txt ./ 2>/dev/null;

# Compile the current example
makeout=$(make -j4);

# Run the current example:
out=$(./sparselizard);

# Get the last character in the output:
out="${out:$((${#out}-1)):1}"

# If the last character is 1 the current example was run successfully:
if [ $out == "1" ]
then
echo 'SUCCESS AT' $i;
sleep 2;
else
echo 'FAILED AT' $i;
echo '';
exit 1;
fi

done

echo ''; 
echo '';
echo 'ALL OK!';


# CLEAN:

rm *.msh;
rm *.nas;
rm *.pos;
rm *.vtk;
rm *.vtu;
rm *.pvd;
rm *.csv;
rm *.slz;
