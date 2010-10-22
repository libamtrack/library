REM ************************************************
REM Batch file for R package compiling under Windows
REM ************************************************

REM *** Remove old folders ***
del package.Rcheck /Q

REM *** Copy sources of libamtrack ***
copy ..\..\include\*.h package\src\
copy ..\..\src\*.c package\src\
del package\src\AT_Wrapper_R.c
del package\src\AT_Wrapper_R.h

REM *** Auto-generate wrappers and documentation ***
R --no-save < read.header.R
R --no-save < convert.header.to.cwrapper.R
R --no-save < convert.header.to.Rwrapper.R
R --no-save < convert.doxygen.to.Rd.R

REM *** Build test package ***
R CMD check ./package

REM *** Copy GSL to libs ***
copy C:\Programme\GnuWin32\bin\libgsl.dll package.Rcheck\libamtrack\libs
copy C:\Programme\GnuWin32\bin\libgslcblas.dll package.Rcheck\libamtrack\libs

REM Pack into zip file
cd package.Rcheck
zip -r libamtrack.zip libamtrack
move libamtrack.zip ..
cd ..
