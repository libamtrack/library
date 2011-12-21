#!/bin/bash

echo "################################################"
echo "This script creates the R package for libamtrack"
echo "################################################"
echo 

HELP_ARG="FALSE"
FAST_ARG="FALSE"
INSTALL_ARG="FALSE"
NOCLEAN_ARG="FALSE"

for var in "$@"
do
    if [ "$var" == "--help" ] ; then
       HELP_ARG="TRUE"
    fi

    if [ "$var" == "--fast" ] ; then
       FAST_ARG="TRUE"
       echo "Fast execution chosen."
    fi

    if [ "$var" == "--install" ] ; then
       INSTALL_ARG="TRUE"
       echo "Will install package after compilation."
    fi

    if [ "$var" == "--noclean" ] ; then
       NOCLEAN_ARG="TRUE"
       echo "Will not remove transient files and folders."
    fi
done

if [ $HELP_ARG == "TRUE" ] ; then
   echo "This script compiles libamtrack to an R package source tarball"
   echo
   echo "Use --fast to skip svn information update (really slow)."
   echo "Use --install to install package after compilation."
   echo "Use --noclean to leave transient files after compilation for debugging."
   echo
   exit
fi

echo "Use --help for information"
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
cp hardcoded_documentation/*.Rd libamtrack/man/

# *** Copy hardcoded wrappers ***
echo "Copying hardcoded wrappers..."
cp hardcoded_wrappers/*.R libamtrack/R/
cp hardcoded_wrappers/hardcoded_wrapper.c libamtrack/src/
cp hardcoded_wrappers/hardcoded_wrapper.h libamtrack/src/

# *** Run autoconfigure (to update svn version) ***
if [ -d $FAST_ARG ] ; then
	echo "Running autoreconf and configure in main folder (to update svnversion information)..."
	cd ../../..
	autoreconf --force --install
	chmod 755 configure
	./configure
	cd wrapper/R/R_package
fi


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

echo "Running R script to add metainformation (date, version, etc.) to R package description..."
R --no-save CMD BATCH ../../../tools/automatic_wrapper_generator/R.add.metainfo.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.add.metainfo.R"
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

# *** Create namespace for package
echo "Create NAMESPACE file for package..."
R --no-save CMD BATCH ../../../tools/automatic_wrapper_generator/R.create.package.namespace.R
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.create.package.namespace.R"
  exit 1
fi

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
#echo
#echo "Building binary distribution package..."
#sudo R CMD INSTALL --build libamtrack

# *** Build tarball source package ***
echo
echo "Building source package tarball..."
R CMD build libamtrack

# *** Remove transient files ***
if [ ! $NOCLEAN_ARG == "TRUE" ] ; then
	echo
	echo "Removing transient files and folder..."
	rm -f -r libamtrack
	rm -f -r libamtrack.Rcheck
	rm -f *.sdd
	rm -f *.Rout
	rm -f *.RData
	# These files should not even be there, TODO: check RD creation script
	rm -f ./package/man/*.Rd
fi

# *** Install package if chose ***
if [ $INSTALL_ARG == "TRUE" ] ; then
    sudo R CMD INSTALL libamtrack_*.tar.gz
fi
echo "Done!"
