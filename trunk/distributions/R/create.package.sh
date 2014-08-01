#!/bin/bash

echo "################################################"
echo "This script creates the R package for libamtrack"
echo "################################################"
echo 
echo "Use --help for information"
echo

ROOT_DIR=$(pwd)
WORK_DIR=$ROOT_DIR/distributions/R
SCRIPT_DIR=$WORK_DIR/scripts

HELP_ARG="FALSE"
FAST_ARG="FALSE"
INSTALL_ARG="FALSE"
NOCLEAN_ARG="FALSE"
NOCOMPILE_ARG="FALSE"
SYNC_ARG="FALSE"


# set bound checking environmental variable
export R_C_BOUNDS_CHECK=yes

while test $# -gt 0; do
        case "$1" in
                -h|--help)
			echo "This script compiles libamtrack to an R package source tarball"
			echo
			echo "Use '--fast' to skip svn information update and R examples."
			echo "Use '--noclean' to leave transient files after compilation for debugging."
			echo "Use '--sync syncdir' to sync source package with directory syncdir."
			echo
			exit 0
                        ;;
                --fast)
                   	FAST_ARG="TRUE"
       			echo "Fast execution chosen. Will not update svn information nor run R examples."
			shift
                        ;;
                --sync)
                        shift
                        if test $# -gt 0; then
                                SYNC_ARG="TRUE"
				SYNC_DIR=$1
                                echo "Sync chosen with" $SYNC_DIR
                        else
                                echo "Error: no sync dir specified"
                                exit 1
                        fi
                        shift
			;;
        esac
done
echo


# Clean up residues from earlier runs if applicable
if [ -d "$WORK_DIR/libamtrack" ]; then
  rm $WORK_DIR/libamtrack -rf
fi
if [ -d "$WORK_DIR/libamtrack.Rcheck" ]; then
  rm $WORK_DIR/libamtrack.Rcheck -rf
fi
if [ ls $WORK_DIR/*.Rout &> /dev/null ]; then
    rm $WORK_DIR/*.Rout
fi
if [ ls $WORK_DIR/*.sdd &> /dev/null ]; then
    rm $WORK_DIR/*.sdd
fi


# *** Create new temporary folder from template package structure ***
echo "Copy package template..."
cp -r $WORK_DIR/package $WORK_DIR/libamtrack


# Remove .svn subdirs
find $WORK_DIR/libamtrack/ -iname .svn -print | xargs rm -rf


# *** Copy sources of libamtrack ***
echo "Copying libamtrack sources..."
cp $ROOT_DIR/include/*.h $WORK_DIR/libamtrack/src/
cp $ROOT_DIR/src/*.c $WORK_DIR/libamtrack/src/


# *** Clean temporary files (e.g. from gedit that will confuse the collect.doxygen.information.R script)
if [ ls $WORK_DIR/libamtrack/src/*.h~ &> /dev/null ]; then
    rm $ROOT_DIR/include/*.h~
fi
if [ ls $WORK_DIR/libamtrack/src/*.c~ &> /dev/null ]; then
    rm $WORK_DIR/libamtrack/src/*.c~
fi

 
# *** Copy hardcoded documentation ***
echo "Copying hardcoded documentation..."
cp $WORK_DIR/hardcoded_documentation/*.Rd $WORK_DIR/libamtrack/man/


# *** Copy hardcoded wrappers ***
echo "Copying hardcoded wrappers..."
cp $WORK_DIR/hardcoded_wrappers/*.R $WORK_DIR/libamtrack/R/
cp $WORK_DIR/hardcoded_wrappers/hardcoded_wrapper.c $WORK_DIR/libamtrack/src/
cp $WORK_DIR/hardcoded_wrappers/hardcoded_wrapper.h $WORK_DIR/libamtrack/src/


# *** Run autoconfigure (to update svn version) ***
if [ $FAST_ARG == "FALSE" ] ; then
	echo "Running autoreconf and configure in main folder (to update svnversion information)..."
	autoreconf --force --install
	chmod 755 configure
	./configure
fi


# *** Copy generated file ***
cp $ROOT_DIR/config.h $WORK_DIR/libamtrack/src/

echo
echo "Running R script to parse doxygen information from sources..."
Rscript --no-save $ROOT_DIR/scripts/collect.doxygen.information.R $WORK_DIR $ROOT_DIR/include >$WORK_DIR/collect.doxygen.information.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing collect.doxygen.information.R"
  exit 1
fi

echo "Running R script to create C headers from parsed doxygen information..."
Rscript --no-save $SCRIPT_DIR/R.generate.C.wrapper.R $SCRIPT_DIR libamtrack >$WORK_DIR/R.generate.C.wrapper.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.C.wrapper.R"
  exit 1
fi

echo "Running R script to create R headers from parsed doxygen information..."
Rscript --no-save $SCRIPT_DIR/R.generate.R.wrapper.R $SCRIPT_DIR/R.type.conversion.R libamtrack >$WORK_DIR/R.generate.R.wrapper.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.R.wrapper.R"
  exit 1
fi

echo "Running R script to add metainformation (date, version, etc.) to R package description..."
Rscript --no-save $SCRIPT_DIR/R.add.metainfo.R $ROOT_DIR libamtrack $WORK_DIR >$WORK_DIR/R.add.metainfo.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.add.metainfo.R"
  exit 1
fi

echo "Running R script to create Rd documentation from parsed doxygen information..."
Rscript --no-save $SCRIPT_DIR/R.generate.Rd.documentation.R $WORK_DIR libamtrack hardcoded_documentation >$WORK_DIR/R.generate.Rd.documentation.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.generate.Rd.documentation.R"
  exit 1
fi

# *** Move resulting file (wrapper C and R), dynamically coded documentation to temporary folder
echo "Moving results from R scripts into package structure..."
mv $ROOT_DIR/AT_R_Wrapper.* $WORK_DIR/libamtrack/src
mv $ROOT_DIR/libamtrack.R $WORK_DIR/libamtrack/R
mv $ROOT_DIR/*.Rd $WORK_DIR/libamtrack/man


# *** Truncate Rd files (esp. the automatically generated ones to 80 characters line width)
echo "Running R script to truncate Rd files to 80 characters line width..."
Rscript --no-save $SCRIPT_DIR/R.truncate.Rd.files.R $WORK_DIR/libamtrack/man >$WORK_DIR/R.truncate.Rd.files.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.truncate.Rd.files.R"
  exit 1
fi


# *** Create namespace for package
echo "Create NAMESPACE file for package..."
Rscript --no-save $SCRIPT_DIR/R.create.package.namespace.R $WORK_DIR libamtrack AT >$WORK_DIR/R.create.package.namespace.Rout 2>&1
if [ "$?" -ne "0" ]; then
  echo "Problem with executing R.create.package.namespace.R"
  exit 1
fi


# *** Run autoconf
echo
echo "Running autoconf to create configure script..."
cd $WORK_DIR/libamtrack
autoconf
chmod 755 configure
cd $ROOT_DIR


# *** Copy package directory for sync ***
if [ $SYNC_ARG == "TRUE" ] ; then
  echo "Syncing..."
  cp -rp $WORK_DIR/libamtrack $ROOT_DIR/libamtrack.pkg
  rm -rf $ROOT_DIR/libamtrack.pkg/autom4te.cache
  rsync -avh $ROOT_DIR/libamtrack.pkg/* $SYNC_DIR
  rm -rf $ROOT_DIR/libamtrack.pkg
fi


# *** Build test package ***
echo
echo "Running package check..."
if [ $FAST_ARG == "FALSE" ] ; then
   R CMD check $WORK_DIR/libamtrack --no-manual 
   if [ "$?" -ne "0" ]; then
     echo "Problem with executing R CMD check ./package --no-manual"
   exit 1
   fi
else
   R CMD check $WORK_DIR/libamtrack --no-manual --no-examples
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
R CMD build $WORK_DIR/libamtrack


# *** Remove transient files ***
if [ $NOCLEAN_ARG == "FALSE" ] ; then
	echo
	echo "Removing transient files and folder..."
	rm -f -r $WORK_DIR/libamtrack
	rm -f -r $ROOT_DIR/libamtrack.Rcheck
	rm -f $ROOT_DIR/*.sdd
	rm -f $WORK_DIR/*.Rout
	rm -f $WORK_DIR/*.RData
fi

echo "Done!"
