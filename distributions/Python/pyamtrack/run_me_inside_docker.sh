#!/bin/bash

pwd
/opt/python/cp36-cp36m/bin/python3 -m pip install /io/cBinder/
CFLAGS=-std=c99 /opt/python/cp36-cp36m/bin/python3 -m cBinder -v pyamtrack -f ../library/src/ -f ../library/include/ -d generated -es ../cBinder/tests/libamtrack/libamtrack_export_symbols_full.txt  -mono libAT compile -i /io/gsl/include/ -i ../library/include/ -b /io/gsl/lib/ -l gsl -l gslcblas -l m -e 'std=c99'