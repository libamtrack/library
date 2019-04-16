#!/bin/bash
echo "installing emsdk"

if [ -z "$1" ]
then
	EMSDK_INSTALL_PATH=$LIB_AT_PATH/..
else
    EMSDK_INSTALL_PATH=$1
fi

cd $EMSDK_INSTALL_PATH	

git clone https://github.com/emscripten-core/emsdk.git --depth 1
cd emsdk
./emsdk install latest
./emsdk activate latest
source ./emsdk_env.sh

cd $GSL_PATH