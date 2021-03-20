#!/bin/sh

dir=$1
echo "Copying simulations/default to simulations/$dir"
cp -r simulations/default simulations/$dir

echo "Fixing simulations/$dir/CMakeLists.txt"
sed -i.bak "s/default /$dir /g" ./simulations/$dir/CMakeLists.txt
rm simulations/$dir/CMakeLists.txt.bak

echo "Adding subdirectory to simulations/CMakeLists.txt"
sed -i.bak -e "s/add_subdirectory[(]$dir[)][\r\n]*//g" ./simulations/CMakeLists.txt
sed -i.bak2 -e "/^[[:space:]]*$/d" ./simulations/CMakeLists.txt
rm simulations/CMakeLists.txt.bak
rm simulations/CMakeLists.txt.bak2

echo "add_subdirectory($dir)" >> simulations/CMakeLists.txt
