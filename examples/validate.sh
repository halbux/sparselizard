# THIS SCRIPTS ENSURES SPARSELIZARD COMPILES AND ALL EXAMPLES RUN SUCCESSFULLY.

#!/bin/bash

cd ..;


# LOAD REQUIRED LIBRARIES:

if [ "$(uname)" == "Linux" ]; then
LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux2-c-opt/lib":$LD_LIBRARY_PATH
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
if [ $i == "validate.sh" ] || [ $i == "electromagnetic-waveguide-time-resolution-3d" ] || [ $i == "superconductor-3d" ] || [ $i == "nonlinear-vonkarman-vortex-2d" ] || [ $i == "thermoacoustic-elasticity-axisymmetry-2d" ]
then
continue
fi

# Copy the example to the main directory: 
cp examples/$i/main.cpp ./;
cp examples/$i/*.msh ./ 2>/dev/null;

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

rm *.pos;
rm *.vtk;
rm *.vtu;
rm *.pvd;
rm *.csv;








