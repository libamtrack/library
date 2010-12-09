#!/bin/bash
rm -rf package.Rcheck
cp ../../../include/*.h package/src/
cp ../../../src/*.c package/src/
cp hardcoded_documentation/* package/man/
cp hardcoded_wrappers/*.R package/R/
cp hardcoded_wrappers/hardcoded_wrapper.c package/src/
cp hardcoded_wrappers/hardcoded_wrapper.h package/src/
R --no-save < ../../../tools/automatic_wrapper_generator/collect.doxygen.information.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing collect.doxygen.information.R"
  exit 1
fi
R --no-save < ../../../tools/automatic_wrapper_generator/R.generate.C.wrapper.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.C.wrapper.R"
  exit 1
fi
R --no-save < ../../../tools/automatic_wrapper_generator/R.generate.R.wrapper.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.R.wrapper.R"
  exit 1
fi
R --no-save < ../../../tools/automatic_wrapper_generator/R.generate.Rd.documentation.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.Rd.documentation.R"
  exit 1
fi
mv AT_R_Wrapper.* package/src
mv libamtrack.R package/R
R CMD check ./package --no-latex # this command produce directory with R package
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R CMD check ./package --no-latex"
  exit 1
fi
rm package/src/*.h
rm package/src/*.c
cd package.Rcheck
tar -zcf libamtrack_0.5.0.tar.gz libamtrack/ # packaging of real package 
mv libamtrack_0.5.0.tar.gz ../
rm -rf package.Rcheck
cd ../../..
