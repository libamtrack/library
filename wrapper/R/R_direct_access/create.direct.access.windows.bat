REM *************************************************
REM Batch file for R wrappers compiling under Windows
REM *************************************************

REM *** Auto-generate wrappers and documentation ***
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\collect.doxygen.information.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.C.wrapper.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.R.wrapper.R

REM *** Copy sources to local directory
copy ..\..\..\include\*.h .
copy ..\..\..\src\*.c .

REM *** Compile libamtrack including R wrappers ***
make all

REM Delete temporary files
del *.o
del *.c
del *.h
del .RData
del functions.sdd
del ..\..\..\tools\automatic_wrapper_generator\*.Rout

copy C:\Programme\GnuWin32\bin\libgsl.dll .
copy C:\Programme\GnuWin32\bin\libgslcblas.dll .

echo "Done."