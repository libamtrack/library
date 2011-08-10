#!/bin/bash

echo "################################################"
echo "This script creates the R package for libamtrack"
echo "################################################"
echo

# *** Create new temporary folder from template package structure ***
echo "Copy package template..."
cp -r package libamtrack

# *** Copy sources of libamtrack ***
echo "Copying libamtrack sources..."
cp ../../../include/*.h libamtrack/src/
cp ../../../src/*.c libamtrack/src/

# *** Copy hardcoded documentation ***
echo "Copying hardcoded documentation..."
cp hardcoded_documentation/* libamtrack/man/

# *** Copy hardcoded wrappers ***
echo "Copying hardcoded wrappers..."
cp hardcoded_wrappers/*.R libamtrack/R/
cp hardcoded_wrappers/hardcoded_wrapper.c libamtrack/src/
cp hardcoded_wrappers/hardcoded_wrapper.h libamtrack/src/

echo
echo "Running R script to parse doxygen information from sources..."
R --no-save CMD BATCH ../../../tools/automatic_wrapper_generator/collect.doxygen.information.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing collect.doxygen.information.R"
  exit 1
fi

echo "Running R script to create C headers from parsed doxygen information..."
R --no-save CMD BATCH ../../../tools/automatic_wrapper_generator/R.generate.C.wrapper.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.C.wrapper.R"
  exit 1
fi

echo "Running R script to create R headers from parsed doxygen information..."
R --no-save CMD BATCH ../../../tools/automatic_wrapper_generator/R.generate.R.wrapper.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.R.wrapper.R"
  exit 1
fi

echo "Running R script to create Rd documentation from parsed doxygen information..."
R --no-save CMD BATCH ../../../tools/automatic_wrapper_generator/R.generate.Rd.documentation.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.Rd.documentation.R"
  exit 1
fi

# *** Move resulting file (wrapper C and R), dynamically coded documentation to temporary folder
echo "Moving results from R scripts into package structure..."
mv AT_R_Wrapper.* libamtrack/src
mv libamtrack.R libamtrack/R
mv *.Rd libamtrack/man

# *** Run autoconf
echo
echo "Running autoconf to create configure script..."
autoconf libamtrack/configure.ac >libamtrack/configure
chmod 755 libamtrack/configure

# *** Build test package ***
echo
echo "Running package check..."
R CMD check ./libamtrack --no-manual
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R CMD check ./package --no-manual"
  exit 1
fi

# *** Build binary distribution ***
echo
echo "Building binary distribution package..."
sudo R CMD INSTALL --build libamtrack

# *** Build tarball source package ***
echo
echo "Building source package tarball..."
R CMD build libamtrack

# *** Remove transient files ***
echo
echo "Removing transient files and folder..."
rm -f -r libamtrack
rm -f -r libamtrack.Rcheck
rm -f functions.sdd
rm -f *.Rout
rm -f *.RData

echo "Done!"
