#!/bin/bash

echo "################################################"
echo "This script creates the R package for libamtrack"
echo "################################################"
echo 

ROOT_DIR=../..
THIS_DIR_FROM_ROOT=distributions/R
SCRIPT_DIR=scripts

HELP_ARG="FALSE"
FAST_ARG="FALSE"
INSTALL_ARG="FALSE"
NOCLEAN_ARG="FALSE"
NOCOMPILE_ARG="FALSE"
SYNCWITH_ARG=""

# set bound checking environmental variable
export R_C_BOUNDS_CHECK=yes

for var in "$@"
do
    if [ "$var" == "--help" ] ; then
       HELP_ARG="TRUE"
    fi

    if [ "$var" == "--fast" ] ; then
       FAST_ARG="TRUE"
       echo "Fast execution chosen. Will not update svn information nor run R examples."
    fi

    if [ "$var" == "--install" ] ; then
       INSTALL_ARG="TRUE"
       echo "Will install package after compilation."
    fi

    if [ "$var" == "--noclean" ] ; then
       NOCLEAN_ARG="TRUE"
       echo "Will not remove transient files and folders."
    fi

    if [ "$var" == "--nocompile" ] ; then
       NOCLEAN_ARG="TRUE"
       echo "Will not remove transient files and folders."
    fi

    if [ "$var" == "--sync-with" ] ; then
       echo "Will sync with remote package dir." 
    fi
done

if [ $HELP_ARG == "TRUE" ] ; then
   echo "This script compiles libamtrack to an R package source tarball"
   echo
   echo "Use --fast to skip svn information update and R examples."
   echo "Use --install to install package after compilation."
   echo "Use --noclean to leave transient files after compilation for debugging."
   echo
   exit
fi

echo "Use --help for information"
echo

# Clean up residues from earlier runs if applicable
if [ -d "./libamtrack" ]; then
  rm ./libamtrack -rf
fi
if [ -d "./libamtrack.Rcheck" ]; then
  rm ./libamtrack.Rcheck -rf
fi
if [ ls ./*.Rout &> /dev/null ]; then
    rm *.Rout
fi
if [ ls ./*.sdd &> /dev/null ]; then
    rm *.sdd
fi

# *** Create new temporary folder from template package structure ***
echo "Copy package template..."
cp -r package libamtrack
# Remove .svn subdirs
find libamtrack/ -iname .svn -print | xargs rm -rf

# *** Copy sources of libamtrack ***
echo "Copying libamtrack sources..."
cp $ROOT_DIR/include/*.h libamtrack/src/
cp $ROOT_DIR/src/*.c libamtrack/src/

# *** Clean temporary files (e.g. from gedit that will confuse the collect.doxygen.information.R script)
if [ ls libamtrack/src/*.h~ &> /dev/null ]; then
    rm $ROOT_DIR/include/*.h~
fi
if [ ls libamtrack/src/*.c~ &> /dev/null ]; then
    rm libamtrack/src/*.c~
fi

 
# *** Copy hardcoded documentation ***
echo "Copying hardcoded documentation..."
cp hardcoded_documentation/*.Rd libamtrack/man/

# *** Copy hardcoded wrappers ***
echo "Copying hardcoded wrappers..."
cp hardcoded_wrappers/*.R libamtrack/R/
cp hardcoded_wrappers/hardcoded_wrapper.c libamtrack/src/
cp hardcoded_wrappers/hardcoded_wrapper.h libamtrack/src/

# *** Run autoconfigure (to update svn version) ***
if [ ! -d $FAST_ARG ] ; then
	echo "Running autoreconf and configure in main folder (to update svnversion information)..."
	cd $ROOT_DIR
	autoreconf --force --install
	chmod 755 configure
	./configure
	cd $THIS_DIR_FROM_ROOT
fi

# *** Copy generated file ***
cp $ROOT_DIR/config.h libamtrack/src/

echo
echo "Running R script to parse doxygen information from sources..."
Rscript --no-save $ROOT_DIR/scripts/collect.doxygen.information.R $ROOT_DIR/include libamtrack >collect.doxygen.information.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing collect.doxygen.information.R"
  exit 1
fi

echo "Running R script to create C headers from parsed doxygen information..."
Rscript --no-save scripts/R.generate.C.wrapper.R scripts libamtrack >R.generate.C.wrapper.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.C.wrapper.R"
  exit 1
fi

echo "Running R script to create R headers from parsed doxygen information..."
Rscript --no-save scripts/R.generate.R.wrapper.R scripts/R.type.conversion.R libamtrack >R.generate.R.wrapper.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.R.wrapper.R"
  exit 1
fi

echo "Running R script to add metainformation (date, version, etc.) to R package description..."
Rscript --no-save scripts/R.add.metainfo.R $ROOT_DIR libamtrack >R.add.metainfo.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.add.metainfo.R"
  exit 1
fi

echo "Running R script to create Rd documentation from parsed doxygen information..."
Rscript --no-save scripts/R.generate.Rd.documentation.R ./hardcoded_documentation/ libamtrack >R.generate.Rd.documentation.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.Rd.documentation.R"
  exit 1
fi

# *** Move resulting file (wrapper C and R), dynamically coded documentation to temporary folder
echo "Moving results from R scripts into package structure..."
mv AT_R_Wrapper.* libamtrack/src
mv libamtrack.R libamtrack/R
mv *.Rd libamtrack/man

# *** Truncate Rd files (esp. the automatically generated ones to 80 characters line width)
cd libamtrack/man
echo "Running R script to truncate Rd files to 80 characters line width..."
Rscript --no-save ../../scripts/R.truncate.Rd.files.R >../../R.truncate.Rd.files.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.truncate.Rd.files.R"
  exit 1
fi
cd ../..

# *** Create namespace for package
echo "Create NAMESPACE file for package..."
Rscript --no-save scripts/R.create.package.namespace.R libamtrack AT >R.create.package.namespace.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.create.package.namespace.R"
  exit 1
fi

# *** Run autoconf
echo
echo "Running autoconf to create configure script..."
cd libamtrack
autoconf
chmod 755 configure
cd ..

#autoconf libamtrack/configure.ac >libamtrack/configure
#chmod 755 libamtrack/configure

# *** Build test package ***
echo
echo "Running package check..."
if [ ! -d $FAST_ARG ] ; then
   R CMD check ./libamtrack --no-manual 
   if [ "$?" -ne "0" ]; then
     echo "Problem with executing R CMD check ./package --no-manual"
   exit 1
   fi
else
   R CMD check ./libamtrack --no-manual --no-examples
   if [ "$?" -ne "0" ]; then
     echo "Problem with executing R CMD check ./package --no-manual --no-examples"
   exit 1
   fi
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
