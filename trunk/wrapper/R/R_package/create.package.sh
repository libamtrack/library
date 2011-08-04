#!/bin/bash

# *** Create new temporary folder from template package structure ***
cp package package.tmp -r

# *** Copy sources of libamtrack ***
cp ../../../include/*.h package.tmp/src/
cp ../../../src/*.c package.tmp/src/

# *** Copy hardcoded documentation ***
cp hardcoded_documentation/* package.tmp/man/

# *** Copy hardcoded wrappers ***
cp hardcoded_wrappers/*.R package.tmp/R/
cp hardcoded_wrappers/hardcoded_wrapper.c package.tmp/src/
cp hardcoded_wrappers/hardcoded_wrapper.h package.tmp/src/

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
mv AT_R_Wrapper.* package.tmp/src
mv libamtrack.R package.tmp/R
mv *.Rd package.tmp/man

# *** Build test package ***
R CMD check ./package.tmp --no-manual
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R CMD check ./package --no-manual"
  exit 1
fi

# *** Build binary distribution ***
sudo R CMD INSTALL --build package.tmp

# Delete temporary files and folders

# Uncomment if source package is not required
# rm package.tmp -r -f

# Create tar ball
mv package.tmp libamtrack
tar -zcf libamtrack_0.5-1.tar.gz libamtrack -X exclude.txt

rm libamtrack -r -f
rm functions.sdd -f
rm ../../../tools/automatic_wrapper_generator/*.Rout -f

echo Done.
