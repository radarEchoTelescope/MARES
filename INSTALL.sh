#!/bin/bash

# Use this script to install the example (i.e. create the executable)
# Use argument 0 to re-build programs and libraries that have been modified sice last build
# Use argument 1 for complete re-build (removing everythinh from build dir and start over)

ERROR_FLAG=0

# Handeling the argument passed to the script
MODE=$1
if [ -z "$1" ]; then
    echo "Useage " $0 " 0 <re-build modified bin / libs> 1 <re-build all>"
    ERROR_FLAG=1
fi
if [ $ERROR_FLAG == 1 ]; then
    exit 1
fi

initial=$PWD

# # We grab Iceraytracing
# if [ -d lib/IceRayTracing ]
# then
#   echo "We have IceRayTracing already, checking if up to date"
#   cd lib/IceRayTracing && git pull
# else
#   echo "We don't have IceRayTracing yet, let's grab it"
#   cd lib && git clone https://github.com/radarEchoTelescope/IceRayTracing
# fi

cd $initial

# We grab libconfig
#if [ -d lib/libconfig ]
#then
 # echo "We have libconfig already, checking if up to date"
 # cd lib/libconfig && git pull
#else
 # echo "We don't have libconfig yet, let's grab it"
 # cd lib && git clone git@github.com:hyperrealm/libconfig.git
#fi

#cd $initial

# This forces a complete re-build
if [ $MODE == 1 ]; then
    rm -rf ./build/*
fi

# Now we run cmake
mkdir -p build
cd build
cmake ../

# Now we make
make 
