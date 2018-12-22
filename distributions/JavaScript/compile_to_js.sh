#!/bin/bash

#!/bin/bash
if [ "$1" != "" ]; then
    echo "WASM parameter set to $1"
    WASM=$1
else
    echo "WASM parameter EMPTY !!! Default set to 1"
    WASM=1
fi

#emcmake="/home/osboxes/emsdk/emscripten/1.37.36/./emcmake"
#emmake="/home/osboxes/emsdk/emscripten/1.37.36/./emmake"
#emcc="/home/osboxes/emsdk/emscripten/1.37.36/./emcc"

cd ../../
mkdir _build
cd _build
cp ../distributions/JavaScript/libgsl.a .
ls .
emcmake cmake ..
emmake make


emcc libat.a libgsl.a -o libat.html -s WASM=$WASM -s EXPORT_ALL=1 \
-s EXTRA_EXPORTED_RUNTIME_METHODS='["ccall", "cwrap"]'

rm ../distributions/JavaScript/output/*
cp libat.a ../distributions/JavaScript/output/
cp libat.html ../distributions/JavaScript/output/
cp libat.wasm ../distributions/JavaScript/output/ 2>/dev/null || : #ignore error when build with -s WASM=0
cp libat.js ../distributions/JavaScript/output/

cd ..
rm -r _build
