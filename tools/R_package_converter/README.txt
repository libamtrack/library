################################
R package creator for libamtrack
################################
Created 2010/10/13 by F. Klein

REQUIREMENTS:
- R (>2.10.1) has to be installed
- Under Windows make sure that your R installation is included in the PATH environment variable.
- GSL (+devel) has to be installed

More general information on the creation of R package can be found here: http://cran.r-project.org/doc/manuals/R-exts.pdf


USAGE:
1. Download the latest version of libamtrack from https://libamtrack.sourceforge.net/svnroot/libamtrack/trunk. At least the folders "/include", "/src" 
and "/tools/R_package_converter" and the Makefile in the main directory have to be present. The folder structure necessary to compile the package as well a a
a mandatory description is already found under "/tools/R_package_converter" and does not have to be provided by the user.

2. The names of all functions that should be included in the package have to be included in the "NAMESPACE" text file.

3. Open console / command line and navigate to "/tools/R_package_converter"

4. Extract information on the present functions in libamtrack from the header files by "R CMD BATCH read.header.R". This will create a file "functions.sdd".

5. Create a script with R wrapper routines by "R CMD BATCH convert.header.to.Rwrapper.R". It will use the information stored in "functions.sdd" and create
a file "wrapper.R" under "/tools/R_package_converter/package/R".

6. Create the include and source files containing the C wrapper routines by "R CMD BATCH convert.header.to.cwrapper.R". This will create "Rwrapper.h" and
"Rwrapper.c" under "/tools/R_package_converter/package/C".

7. Create documentation in RD format by "R CMD BATCH convert.doxygen.to.Rd.R". This will create *.Rd documentation files under "/tools/R_package_converter/package/man".

8. A Makevars file for compiling and linking against a shared object is provided in "/tools/R_package_converter/package/src".  
	
9. Run the command "R CMD check package" to check for possible errors

# till now it works (LG)

10. If all errors are fixed, use "R CMD build package" to build package