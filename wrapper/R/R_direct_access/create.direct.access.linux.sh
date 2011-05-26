#!/bin/bash
echo "*************************************************"
echo "Batch file for R wrappers compiling under Windows"
echo "*************************************************"

echo "Auto-generate wrappers and documentation"
R --no-save < ../../../tools/automatic_wrapper_generator/collect.doxygen.information.R
R --no-save < ../../../tools/automatic_wrapper_generator/R.generate.C.wrapper.R
R --no-save --args nopackage < ../../../tools/automatic_wrapper_generator/R.generate.R.wrapper.R

echo "Copy sources to local directory"
cp -r ../../../include/*.h .
cp -r ../../../src/*.c .

echo "Compile libamtrack including R wrappers"
make all

echo "Delete temporary files"
rm -rf *.o
rm -rf *.c
rm -rf *.h
rm -rf .RData
rm -rf functions.sdd
rm -rf ../../../tools/automatic_wrapper_generator/*.Rout

echo "Done."