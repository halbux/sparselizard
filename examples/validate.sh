# THIS SCRIPTS ENSURES SPARSELIZARD COMPILES AND ALL EXAMPLES RUN SUCCESSFULLY.

#!/bin/bash

cd ..;


# LOAD REQUIRED LIBRARIES:

if [ "$(uname)" == "Linux" ]; then
LD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-linux2-c-opt/lib":"$HOME/SLlibs/slepc/arch-linux2-c-opt/lib":"$HOME/SLlibs/openblas/install/lib":"$HOME/SLlibs/fftw/install/lib":$LD_LIBRARY_PATH
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH    
elif [ "$(uname)" == "Darwin"  ]; then
DYLD_LIBRARY_PATH="$HOME/SLlibs/petsc/arch-darwin-c-opt/lib":"$HOME/SLlibs/slepc/arch-darwin-c-opt/lib":"$HOME/SLlibs/openblas/install/lib":"$HOME/SLlibs/fftw/install/lib":$DYLD_LIBRARY_PATH
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH
fi


# LOOP ON ALL EXAMPLES:

for i in $(ls examples); 
do 

if [ $i == "validate.sh" ]
then
continue
fi

# Copy the example to the main directory: 
cp examples/$i/main.cpp ./;
gmsh examples/$i/*.geo -3;
mv examples/$i/*.msh ./;

# Compile the current example
make -j4;

# Run the current example:
out=$(./sparselizard);

# Get the last character in the output:
out="${out:$((${#out}-1)):1}"

# If the last character is 1 the current example was run successfully:
if [ $out == "1" ]
then
echo '______________________________________________________';
echo 'SUCCESS AT' $i;
sleep 2;
else
echo '______________________________________________________';
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








