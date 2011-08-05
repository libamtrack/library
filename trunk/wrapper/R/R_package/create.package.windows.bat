REM ************************************************
REM Batch file for R package compiling under Windows
REM ************************************************

REM *** Create new temporary folder from template package structure ***
xcopy package libamtrack /E /I /EXCLUDE:exclude.txt

REM *** Copy sources of libamtrack ***
copy ..\..\..\include\*.h libamtrack\src\
copy ..\..\..\src\*.c libamtrack\src\

REM *** Copy hardcoded documentation ***
copy .\hardcoded_documentation .\libamtrack\man\

REM *** Copy hardcoded wrappers ***
copy .\hardcoded_wrappers\*.R .\libamtrack\R\
copy .\hardcoded_wrappers\hardcoded_wrapper.c .\libamtrack\src\
copy .\hardcoded_wrappers\hardcoded_wrapper.h .\libamtrack\src\

REM *** Auto-generate wrappers and documentation ***
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\collect.doxygen.information.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.C.wrapper.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.R.wrapper.R
R CMD BATCH ..\..\..\tools\automatic_wrapper_generator\R.generate.Rd.documentation.R

REM *** Move resulting file (wrapper C and R), dynamically coded documentation to temporary folder
move AT_R_Wrapper.* .\libamtrack\src
move libamtrack.R .\libamtrack\R
move *.Rd libamtrack/man

REM *** Copy GSL to src directory ***
copy C:\Programme\GnuWin32\bin\libgsl.dll libamtrack\src
copy C:\Programme\GnuWin32\bin\libgslcblas.dll libamtrack\src

REM *** Build test package ***
R CMD check libamtrack

REM *** Enable for debugging ***
REM goto :eof

REM *** Build actual package ***
R CMD INSTALL --build libamtrack

REM Delete temporary files and folders
del libamtrack /Q
del libamtrack.Rcheck /Q
del .RData
del functions.sdd
del ..\..\..\tools\automatic_wrapper_generator\*.Rout

echo "Done."