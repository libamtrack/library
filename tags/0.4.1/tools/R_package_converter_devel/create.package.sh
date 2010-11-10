#!/bin/bash
rm -rf package.Rcheck
cp ../../include/*.h package/src/
cp ../../src/*.c package/src/
rm package/src/AT_Wrapper_R.c
rm package/src/AT_Wrapper_R.h
cp hardcoded_documentation/* package/man/
cp hardcoded_wrappers/*.R package/R/
cp hardcoded_wrappers/hardcoded_wrapper.c package/src/
cp hardcoded_wrappers/hardcoded_wrapper.h package/src/
R --no-save < read.header.R
R --no-save < convert.header.to.cwrapper.R
R --no-save < convert.header.to.Rwrapper.R
R --no-save < convert.doxygen.to.Rd.R
R CMD check ./package --no-latex # this command produce directory with R package
rm package/src/*.h
rm package/src/*.c
cd package.Rcheck
tar -zcf libamtrack_0.5.0.tar.gz libamtrack/ # packaging of real package
mv libamtrack_0.5.0.tar.gz ../
rm -rf package.Rcheck
cd ../..
