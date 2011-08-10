#!/bin/bash

# *** Create new temporary folder from template package structure ***
cp -r package libamtrack

# *** Copy sources of libamtrack ***
cp ../../../include/*.h libamtrack/src/
cp ../../../src/*.c libamtrack/src/

# *** Copy hardcoded documentation ***
cp hardcoded_documentation/* libamtrack/man/

# *** Copy hardcoded wrappers ***
cp hardcoded_wrappers/*.R libamtrack/R/
cp hardcoded_wrappers/hardcoded_wrapper.c libamtrack/src/
cp hardcoded_wrappers/hardcoded_wrapper.h libamtrack/src/

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

# *** Move resulting file (wrapper C and R), dynamically coded documentation to temporary folder
mv AT_R_Wrapper.* libamtrack/src
mv libamtrack.R libamtrack/R
mv *.Rd libamtrack/man

# *** Run autoconf
autoconf libamtrack/configure.ac >libamtrack/configure
chmod 755 libamtrack/configure

# *** Build test package ***
R CMD check ./libamtrack --no-manual
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R CMD check ./package --no-manual"
  exit 1
fi

# *** Build binary distribution ***
sudo R CMD INSTALL --build libamtrack

# *** Build tarball source package ***
R CMD build libamtrack

# *** Remove transient files ***
rm -f -r libamtrack
rm -f -r libamtrack.Rcheck
rm -f functions.sdd

echo Done.