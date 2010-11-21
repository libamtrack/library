REM ************************************************
REM Batch file for R package compiling under Windows
REM ************************************************

REM *** Remove old folders ***
del package.Rcheck /Q

REM *** Copy sources of libamtrack ***
copy ..\..\..\include\*.h package\src\
copy ..\..\..\src\*.c package\src\

REM *** Copy hardcoded documentation ***
copy .\hardcoded_documentation .\package\man\

REM *** Copy hardcoded wrappers ***
copy .\hardcoded_wrappers\hardcoded_wrapper.R .\package\R\
copy .\hardcoded_wrappers\hardcoded_wrapper.c .\package\src\
copy .\hardcoded_wrappers\hardcoded_wrapper.h .\package\src\

REM *** Auto-generate wrappers and documentation ***
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\collect.doxygen.information.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.C.wrapper.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.R.wrapper.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.Rd.documentation.R

move AT_R_Wrapper.* .\package\src
move libamtrack.R .\package\R

REM *** Build test package ***
R CMD check ./package

REM *** Enable for debugging ***
REM goto :eof

REM *** Copy GSL to libs ***
copy C:\Programme\GnuWin32\bin\libgsl.dll package.Rcheck\libamtrack\libs
copy C:\Programme\GnuWin32\bin\libgslcblas.dll package.Rcheck\libamtrack\libs

REM Pack into zip file
cd package.Rcheck
zip -r libamtrack.zip libamtrack
move libamtrack.zip ..
cd ..

REM Delete temporary files
del package\R\libamtrack.R
del package\R\hardcoded_wrapper.R
del package\src\*.o
del package\src\*.c
del package\src\*.h
del package\src\*.dll
del package\man\*.* /Q
del .RData
del functions.sdd
del package.Rcheck /Q
del ..\..\..\tools\automatic_wrapper_generator\*.Rout

echo "Done."