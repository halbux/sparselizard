#!/bin/sh

GMSH_VERSION=4.8.4

mkdir temp
cd temp
wget https://gmsh.info/bin/Linux/gmsh-$GMSH_VERSION-Linux64.tgz -O gmsh-binary.tgz
tar zxvf gmsh-binary.tgz --strip-components 1
cp bin/gmsh /usr/bin/
cp -r share/doc/gmsh /usr/share/doc/
cp share/man/man1/gmsh.1 /usr/share/man/man1/
cd ..
rm -rf temp