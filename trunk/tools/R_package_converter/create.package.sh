#!/bin/bash
rm -rf package.Rcheck
cp ../../include/*.h package/src/
cp ../../src/*.c package/src/
rm package/src/AT_Wrapper_R.c
rm package/src/AT_Wrapper_R.h
R --no-save < read.header.R
R --no-save < convert.header.to.cwrapper.R
R --no-save < convert.header.to.Rwrapper.R
R --no-save < convert.doxygen.to.Rd.R
R CMD build package
R CMD check libamtrack_0.3.tar.gz
rm package/src/*.h
rm package/src/*.c
cd ../..
