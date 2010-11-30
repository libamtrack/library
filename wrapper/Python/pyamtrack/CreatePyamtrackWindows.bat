REM ********************************************
REM Batch file for Python wrappers under Windows
REM ********************************************

REM *** Auto-generate wrappers and documentation ***
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\collect.doxygen.information.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\Python.generate.wrapper.R

REM *** Copy sources to local directory
copy ..\..\..\include\*.h .
copy ..\..\..\src\*.c .

REM *** Compile libamtrack including R wrappers ***
make all

REM Delete temporary files
del *.o
del *.c
del *.h
del functions.sdd
del ..\..\..\tools\automatic_wrapper_generator\*.Rout

copy C:\Programme\GnuWin32\bin\libgsl.dll package.Rcheck\libamtrack\libs
copy C:\Programme\GnuWin32\bin\libgslcblas.dll package.Rcheck\libamtrack\libs

echo "Done."