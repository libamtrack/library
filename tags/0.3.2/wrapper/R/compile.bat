@echo off
cls

set swigexe="Z:\swig\swig.exe"
set gslinclude="C:\Program Files\GnuWin32\include"
set rinclude="C:\Program Files\R\R-2.9.0\include"
set rlib="C:\Program Files\R\R-2.9.0\bin"
set gsllib="C:\Program Files\GnuWin32\lib"
set gsldllA="C:\Program Files\GnuWin32\bin\libgsl.dll"
set gsldllB="C:\Program Files\GnuWin32\bin\libgslcblas.dll"
set gccexe="C:\Program Files\MinGW\bin\gcc.exe"

:: Generate interface elements using SWIG 
mkdir c-swig-src
set swigwrapper=example_wrap
%swigexe% -dll example -r -o c-swig-src\%swigwrapper%.c -outdir . ..\example.i

::  Compile libamtrack C library + SWIG C wrapper

mkdir obj
del obj\%swigwrapper%.o
:: Compilation
%gccexe% -DDUPA -I.. -I..\..\include -I%rinclude% -I%gslinclude% -O3 -fno-strict-aliasing -c c-swig-src\%swigwrapper%.c -oobj\%swigwrapper%.o
%gccexe% -I.. -I..\..\include -I%gslinclude% -O3 -fno-strict-aliasing -c ../example.c -oobj\example.o
for /f %%a IN ('dir /b ..\..\src\*.c') do %gccexe% -DDUPA -I.. -I..\..\include -I%gslinclude% -O3 -fno-strict-aliasing -c ..\..\src\%%a -oobj\%%~na.o
:: Linking
%gccexe% -shared -L%gsllib% -L%rlib% -mno-cygwin -Wl,--add-stdcall-alias -oexample.dll obj\*.o -lgsl -lgslcblas -lm -lR 

del obj\*.o

copy %gsldllA% .
copy %gsldllB% .