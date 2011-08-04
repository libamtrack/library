REM ************************************************
REM Batch file for R package compiling under Windows
REM ************************************************

REM *** Create new temporary folder from template package structure ***
xcopy package package.tmp /E /I /EXCLUDE:exclude.txt

REM *** Copy sources of libamtrack ***
copy ..\..\..\include\*.h package.tmp\src\
copy ..\..\..\src\*.c package.tmp\src\

REM *** Copy hardcoded documentation ***
copy .\hardcoded_documentation .\package.tmp\man\

REM *** Copy hardcoded wrappers ***
copy .\hardcoded_wrappers\*.R .\package.tmp\R\
copy .\hardcoded_wrappers\hardcoded_wrapper.c .\package.tmp\src\
copy .\hardcoded_wrappers\hardcoded_wrapper.h .\package.tmp\src\

REM *** Auto-generate wrappers and documentation ***
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\collect.doxygen.information.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.C.wrapper.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.R.wrapper.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.Rd.documentation.R

REM *** Move resulting file (wrapper C and R), dynamically coded documentation to temporary folder
move AT_R_Wrapper.* .\package\src
move libamtrack.R .\package\R
move *.Rd package.tmp/man

REM *** Copy GSL to src directory ***
copy C:\Programme\GnuWin32\bin\libgsl.dll package.tmp\src
copy C:\Programme\GnuWin32\bin\libgslcblas.dll package.tmp\src

REM *** Build test package ***
R CMD check package.tmp

REM *** Enable for debugging ***
REM goto :eof

REM *** Build actual package ***
R CMD INSTALL --build package.tmp

REM Delete temporary files and folders
del package.tmp /Q
del package.tmp.Rcheck /Q
del .RData
del functions.sdd
del ..\..\..\tools\automatic_wrapper_generator\*.Rout

echo "Done."