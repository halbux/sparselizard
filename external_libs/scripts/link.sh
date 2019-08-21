#!/bin/bash

if [ ! -d $1 ]; then
  echo "The first argument to this script should be a directory with existing Sparselizard external libs. It should contain subfolders for each library (i.e. 'petsc' should be a subdirectory)";
  exit;
fi;

SCRIPT=$(readlink -f "$0")
SCRIPTDIR=$(dirname "$SCRIPT")
ln -s $1 ${SCRIPTDIR}/../libs