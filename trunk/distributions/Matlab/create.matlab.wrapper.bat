REM ********************************************
REM Batch file for Matlab wrappers under Windows
REM ********************************************

REM *** Auto-generate wrappers ***
R CMD BATCH ..\..\tools\automatic_wrapper_generator\collect.doxygen.information.R
R CMD BATCH ..\..\tools\automatic_wrapper_generator\Matlab.generate.C.header.R
R CMD BATCH ..\..\tools\automatic_wrapper_generator\Matlab.generate.Matlab.wrapper.R

REM *** Copy sources to local directory
copy ..\..\src\*.c .

REM *** Compile libamtrack ***
make all

REM Delete temporary files
rename libamtrack.h xx
del *.o
del *.c
del *.h
del .RData
del functions.sdd
del ..\..\tools\automatic_wrapper_generator\*.Rout

rename xx libamtrack.h 

copy C:\Programme\GnuWin32\bin\libgsl.dll .
copy C:\Programme\GnuWin32\bin\libgslcblas.dll .

echo "Done."