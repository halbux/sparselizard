#!/bin/sh

pwd=$(pwd)
echo "Current Working directory: $pwd"

# raise the memory limit for the stack trace
# since calculating the spanning tree is done recursively and
# it ends up to be a deep recursion in dense models
ulimit -s 30000

cd $pwd
./scripts/make_script.sh

cd build/simulations/mcvt
./mcvt Parrot "$@"
sleep 1m