#!/bin/bash

echo "################################################"
echo "This script creates the R package for libamtrack"
echo "################################################"
echo 

HELP_ARG="FALSE"
FAST_ARG="FALSE"
INSTALL_ARG="FALSE"
NOCLEAN_ARG="FALSE"

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
cp ../../../include/*.h libamtrack/src/
cp ../../../src/*.c libamtrack/src/

# *** Clean temporary files (e.g. from gedit that will confuse the collect.doxygen.information.R script)
if [ ls ../../../include/*.h~ &> /dev/null ]; then
    rm ../../../include/*.h~
fi
if [ ls ../../../src/*.c~ &> /dev/null ]; then
    rm ../../../src/*.c~
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
	cd ../../..
	autoreconf --force --install
	chmod 755 configure
	./configure
	cd wrapper/R/R_package
fi

# *** Copy generated file ***
cp ../../../config.h libamtrack/src/

echo
echo "Running R script to parse doxygen information from sources..."
Rscript --no-save ../../../tools/automatic_wrapper_generator/collect.doxygen.information.R ../../../include libamtrack >collect.doxygen.information.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing collect.doxygen.information.R"
  exit 1
fi

echo "Running R script to create C headers from parsed doxygen information..."
Rscript --no-save ../../../tools/automatic_wrapper_generator/R.generate.C.wrapper.R ../../../tools/automatic_wrapper_generator libamtrack >R.generate.C.wrapper.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.C.wrapper.R"
  exit 1
fi

echo "Running R script to create R headers from parsed doxygen information..."
Rscript --no-save ../../../tools/automatic_wrapper_generator/R.generate.R.wrapper.R ../../../tools/automatic_wrapper_generator/R.type.conversion.R libamtrack >R.generate.R.wrapper.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.R.wrapper.R"
  exit 1
fi

echo "Running R script to add metainformation (date, version, etc.) to R package description..."
Rscript --no-save ../../../tools/automatic_wrapper_generator/R.add.metainfo.R ../../../ libamtrack >R.add.metainfo.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.add.metainfo.R"
  exit 1
fi

echo "Running R script to create Rd documentation from parsed doxygen information..."
Rscript --no-save ../../../tools/automatic_wrapper_generator/R.generate.Rd.documentation.R ./hardcoded_documentation/ libamtrack >R.generate.Rd.documentation.Rout 2>&1
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
#cd ./libamtrack/man
#for f in *.Rd
#do
# fold -w 80 -s $f > $f.output
#done
#rm *.Rd
#rename 's/\.output$//' *.output
#rm *.output
#cd ../..

# *** Create namespace for package
echo "Create NAMESPACE file for package..."
Rscript --no-save ../../../tools/automatic_wrapper_generator/R.create.package.namespace.R libamtrack AT >R.create.package.namespace.Rout 2>&1
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
