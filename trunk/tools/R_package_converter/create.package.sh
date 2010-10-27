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
R CMD build package # this command produce .tar.gz file which is not a R package
R CMD check libamtrack_0.3.tar.gz --no-latex # this command produce directory with R package
rm libamtrack_0.3.tar.gz # let us remove this file as this is not really a R package
rm package/src/*.h
rm package/src/*.c
cd libamtrack.Rcheck
tar -zcf libamtrack_0.3.tar.gz libamtrack/ # packaging of real package
mv libamtrack_0.3.tar.gz ../
rm -rf libamtrack.Rcheck
cd ../..
