#!/bin/sh

pwd=$(pwd)
echo "Current Working directory: $pwd"

cd $pwd
./scripts/make_script.sh

cd build/simulations/mcvt
./mcvt Parrot "$@"
sleep 1m