#!/bin/bash

pwd

CFLAGS=-std=c99 /opt/python/cp36-cp36m/bin/python3 /io/cBinder/cBinder/main.py -v pyamtrack -f ../library/src/ -f ../library/include/ -d generated -es ../cBinder/tests/libamtrack/libamtrack_export_symbols_full.txt  -mono libAT compile -i /io/gsl/include/ -i ../library/include/ -b /io/gsl/lib/ -l gsl -l gslcblas -l m -e 'std=c99'