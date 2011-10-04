REM *************************************************
REM Batch file for R wrappers compiling under Windows
REM *************************************************

REM *** Auto-generate wrappers and documentation ***
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\collect.doxygen.information.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.C.wrapper.R
R --slave --no-save --args nopackage <..\..\..\tools\automatic_wrapper_generator\R.generate.R.wrapper.R

REM *** Copy sources to local directory
copy ..\..\..\include\*.h .
copy ..\..\..\src\*.c .
copy ..\R_package\hardcoded_wrappers\hardcoded_wrapper.* .

REM *** Merge R wrapper
copy libamtrack.R+hardcoded_wrapper.R libamtrack.R

REM *** Compile libamtrack including R wrappers ***
make all

REM Uncomment to leave temporary file for debugging
REM goto EOF

REM TODO REMOVE 'PACKAGE="libamtrack"' from R wrapper!
grep

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