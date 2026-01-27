#!/bin/bash

# Use this script to install the example (i.e. create the executable)
# Use argument 0 to re-build programs and libraries that have been modified sice last build
# Use argument 1 for complete re-build (removing everythinh from build dir and start over)

ERROR_FLAG=0

# Handeling the argument passed to the script
MODE=$1
if [ -z "$1" ]; then
    echo "Useage " $0 " 0 <re-build modified lib(s) and bin(s)> 1 <re-build all> / 2 <update>"
    ERROR_FLAG=1
fi
if [ $ERROR_FLAG == 1 ]; then
    exit 1
fi

initial=$PWD

# Add Iceraytracing
if [ ! -d extern/IceRayTracing ]
then
    echo "We don't have IceRayTracing yet, let's grab it"
    cd extern && git clone https://github.com/radarEchoTelescope/IceRayTracing
    cd $initial
fi

# MODE 2 upgrades all dependencies and re-builds
if [ $MODE == 2 ]; then 
    echo "Checking if IceRayTracing is up to date"
    cd extern/IceRayTracing && git pull
    cd $initial
# If we have user-installed libconfig, update that too
    if [ -d extern/libconfig ]; then
        echo "Checking if libconfig is up to date"
        cd extern/libconfig && git pull
        cd $initial
    fi
    MODE=1
fi

# MODE 1 forces a complete re-build
if [ $MODE == 1 ]; then
    rm -rf ./build/*
fi

# Regardless of above changes, now we run cmake
mkdir -p build
cd build
cmake ../
make 
