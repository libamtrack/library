#!/bin/bash

pwd

cd /io/gsl-2.6
mkdir gsl
./configure --prefix=/io/gsl
make -j2
make install

/opt/python/cp36-cp36m/bin/python3 -m pip install /io/cBinder/
CFLAGS=-std=c99 /opt/python/cp36-cp36m/bin/python3 -m cBinder pyamtrack -f /io/libamtrack/src/ -f /io/libamtrack/include/ -d /io/generated -es /io/symbols_to_export.txt  -mono libAT compile -i /io/gsl/include/ -i /io/libamtrack/include/ -b /io/gsl/lib/ -l gsl -l gslcblas -l m
