################################
R package creator for libamtrack
################################
Created 2010/10/13 by F. Klein

REQUIREMENTS:
- R (>2.10.1) has to be installed
- Under Windows make sure that your R installation is included in the PATH environment variable.
- Also under Windows, some additional tools have to be installed. They are collected and provided here: http://www.murdoch-sutherland.com/Rtools/
- Again under Windows, check if the required paths have been added to your PATH variable.
- GSL (+devel) has to be installed

More general information on the creation of R package can be found here: http://cran.r-project.org/doc/manuals/R-exts.pdf
A very helpful description on how to deal with R package creation under Windows has been provided by P. Rossi:  http://gsbwww.uchicago.edu/fac/peter.rossi/research/bayes%20book/bayesm/Making%20R%20Packages%20Under%20Windows.pdf)

USAGE:
1. Download the latest version of libamtrack from https://libamtrack.sourceforge.net/svnroot/libamtrack/trunk. At least the folders "/include", "/src" 
and "/tools/R_package_converter" and the Makefile in the main directory have to be present. The folder structure necessary to compile the package as well a a
a mandatory description is already found under "/tools/R_package_converter" and does not have to be provided by the user.

2. The names of all functions that should be included in the package have to be included in the "NAMESPACE" text file.

3. Open console / command line, navigate to the main directory of libamtrack and build libamtrack using "make static".

4. Navigate to "/tools/R_package_converter"

5. Extract information on the present functions in libamtrack from the header files by "R CMD BATCH read.header.R". This will create a file "functions.sdd".

6. Create a script with R wrapper routines by "R CMD BATCH convert.header.to.Rwrapper.R". It will use the information stored in "functions.sdd" and create
a file "wrapper.R" under "/tools/R_package_converter/package/R".

7. Create the include and source files containing the C wrapper routines by "R CMD BATCH convert.header.to.cwrapper.R". This will create "Rwrapper.h" and
"Rwrapper.c" under "/tools/R_package_converter/package/C".

8. Create documentation in RD format by "R CMD BATCH convert.doxygen.to.Rd.R". This will create *.Rd documentation files under "/tools/R_package_converter/package/man".

9. A Makevars file for compiling and linking against a shared object is provided in "/tools/R_package_converter/package/src". Under Windows, please delete the "Makevars" file and rename "Makevars.windows" to "Makevars". Moreover, adjust the paths to your installation of MinGW and GSL.
	
10. Run the command "R CMD check package" to check for possible errors

11. If all errors are fixed, use "R CMD build package" to build package