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
R CMD BATCH read.header.R
R CMD BATCH convert.header.to.cwrapper.R
R CMD BATCH convert.header.to.Rwrapper.R
R CMD BATCH convert.doxygen.to.Rd.R

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

REM Delete temporary files
del package\R\wrapper.R
del package\src\*.o
del package\src\*.c
del package\src\*.h
del package\src\*.dll