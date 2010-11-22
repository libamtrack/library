#!/bin/bash
rm -rf package.Rcheck
cp ../../../include/*.h package/src/
cp ../../../src/*.c package/src/
cp hardcoded_documentation/* package/man/
cp hardcoded_wrappers/hardcoded_wrapper.R package/R/
cp hardcoded_wrappers/hardcoded_wrapper.c package/src/
cp hardcoded_wrappers/hardcoded_wrapper.h package/src/
R --no-save < ../../../tools/automatic_wrapper_generator/collect.doxygen.information.R
R --no-save < ../../../tools/automatic_wrapper_generator/R.generate.C.wrapper.R
R --no-save < ../../../tools/automatic_wrapper_generator/R.generate.C.wrapper.R
R --no-save < ../../../tools/automatic_wrapper_generator/R.generate.C.wrapper.R
mv AT_R_Wrapper.* package/src
mv libamtrack.R package/R
R CMD check ./package --no-latex # this command produce directory with R package
rm package/src/*.h
rm package/src/*.c
cd package.Rcheck
tar -zcf libamtrack_0.4.1.tar.gz libamtrack/ # packaging of real package
mv libamtrack_0.4.1.tar.gz ../
rm -rf package.Rcheck
cd ../../..
