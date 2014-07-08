@echo off
cls

set swigexe="Z:\swig\swig.exe"
set jarsignerexe="C:\Program Files\Java\jdk1.6.0_18\bin\jarsigner.exe"
set javacexe="C:\Program Files\Java\jdk1.6.0_18\bin\javac.exe"
set jarexe="C:\Program Files\Java\jdk1.6.0_18\bin\jar.exe"
set javaincludeA="C:\Program Files\Java\jdk1.6.0_18\include"
set javaincludeB="C:\Program Files\Java\jdk1.6.0_18\include\win32"
set gslinclude="C:\Program Files\GnuWin32\include"
set gsllib="C:\Program Files\GnuWin32\lib"
set gsldllA="C:\Program Files\GnuWin32\bin\libgsl.dll"
set gsldllB="C:\Program Files\GnuWin32\bin\libgslcblas.dll"
set gccexe="C:\Program Files\MinGW\bin\gcc.exe"

:: Generate interface elements using SWIG 
mkdir java-swig-src
mkdir c-swig-src
set swigwrapper=example_wrap
%swigexe% -java -o c-swig-src\%swigwrapper%.c -outdir java-swig-src example.i

::  Compile libamtrack C library + SWIG C wrapper

mkdir obj
del obj\%swigwrapper%.o
:: Compilation
%gccexe% -DSWIG -I. -I..\..\include -I%javaincludeB% -I%javaincludeA% -I%gslinclude% -O3 -fno-strict-aliasing -c c-swig-src\%swigwrapper%.c -oobj\%swigwrapper%.o
%gccexe% -I. -I..\..\include -I%javaincludeB% -I%javaincludeA% -I%gslinclude% -O3 -fno-strict-aliasing -c example.c -oobj\example.o
for /f %%a IN ('dir /b ..\..\src\*.c') do %gccexe% -DDUPA -I.. -I..\..\include -I%javaincludeB% -I%javaincludeA% -I%gslinclude% -O3 -fno-strict-aliasing -c ..\..\src\%%a -oobj\%%~na.o
:: Linking
%gccexe% -shared -L%gsllib% -mno-cygwin -Wl,--add-stdcall-alias -oexample.dll obj\*.o -lgsl -lgslcblas -lm 

del obj\*.o

::  Compile Java GUI

%javacexe% java-swig-src\*.java src\*.java -d bin\

::  Create JAR file

copy bin\*.class .
copy %gsldllA% .
copy %gsldllB% .
%jarexe% cvfm example.jar MANIFEST.MF *.class example.dll libgsl.dll libgslcblas.dll
del *.class
del example.dll libgsl.dll libgslcblas.dll
del bin\*.class

::  Sign JAR file

:: keystore created using command "keytool -genkey -keystore myk -alias jdc", password "libamtrack"

%jarsignerexe% -keystore myk -storepass libamtrack -signedjar examplesigned-Windows.jar example.jar jdc
move examplesigned-Windows.jar webstart\
