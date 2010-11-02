==================================================================================================================================
= Prerequisites

In order to use R interface for libamtrack following things needs to be installed:

1. R environment for statistical computing and graphics (http://www.r-project.org/index.html). On most linux platforms it is 
available in the repositories (please install also -dev version of packages). Windows users can download it from project webpage.

2. GCC compiler (preferably v.4.0+) . Linux user usually have it by default. Windows users can download
 port of GCC for windows - MinGW compiler (Minimalist GNU for Windows), from the webpage http://www.mingw.org/.
 
3. Swig interface compiler, preferably v1.3. On most linux platform it available in the repositories as swig package. 
Windows users can download it from project webpage (visit http://www.swig.org/download.html).

4. GNU Scientific Library (http://www.gnu.org/software/gsl/). On most linux platforms it is avalaible in the repositories (please install
also -dev version of packages). Windows users can download it from here: http://gnuwin32.sourceforge.net/packages/gsl.htm




==================================================================================================================================
= Compilation


==========================================
== Windows systems

Locate file wrapper\R\compile.bat and setup following variables:

set swigexe="Z:\swig\swig.exe"
set gslinclude="C:\Program Files\GnuWin32\include"
set rinclude="C:\Program Files\R\R-2.9.0\include"
set rlib="C:\Program Files\R\R-2.9.0\bin"
set gsllib="C:\Program Files\GnuWin32\lib"
set gsldllA="C:\Program Files\GnuWin32\bin\libgsl.dll"
set gsldllB="C:\Program Files\GnuWin32\bin\libgslcblas.dll"
set gccexe="C:\Program Files\MinGW\bin\gcc.exe"

Start windows command line (start -> Run -> cmd.exe), go to the directory
containing compile.bat file (using cd command) and run it (typing compile.bat).


==========================================
== Linux systems

Locate file wrapper\R\compile.sh and setup following variables:

SWIGEXE=swig
GSLINCLUDE="/usr/include"
RINCLUDE="/usr/share/R/include/"
RLIB="/usr/lib/R/lib"
GSLLIB="/usr/lib"
GSLDLLA="/usr/lib/libgsl.so"
GSLDLLB="/usr/lib/libgslcblas.so"
GCCEXE=gcc

Here if you have everything properly installed you will need only to adjust
R include and library paths (you can locate them by trying to find R.h and 
libR.so / R.dll file using comand locate) and gsl libraries.

Start your favorite console,  go to the directory
containing compile.sh file (using cd command) and run it (typing ./compile.sh).



==================================================================================================================================
= Compilation details


Both compile scripts do the same job, just with slightly different parameters of compiler.

==========================================
== Step 1. Generate interface elements using SWIG


First - please read SWIG tutorial : http://www.swig.org/tutorial.html then please 
inspect file wrapper/example.i
This file includes header files from libamtrack C library and contains list
of function that we will use in interface.

Command 
$SWIGEXE -dll example -r -o c-swig-src/$SWIGWRAPPER.c -outdir . ../example.i
will create interface files which consist of two sets:

A) C interface (*.c) which will be saved into c-swig-src directory and will be later compiled
together with the C libramtrack library
B) R wrapper that contains caller method to C functions and does necessary casting. 
It will be later loaded together with library during R session

Please do not commit generated files into repository.


==========================================
== Step 2. Compile libamtrack C library + SWIG C wrapper


In this step all the *.c files will be compiled to the *.o files (these
will go to obj directory). Finally all *.o files and necessary libraries (GSL) 
will be linked together to create C library file (example.dll or example.so respectively).

Note that library name cannot follow usual naming schema under Linux
(libXXX.so), but have to be named example.so. 

Please also note that there is difference in compilation flags:
Windows: -fno-strict-aliasing 
Linux: -fPIC 
and in linking flags:
Windows: -mno-cygwin -Wl,--add-stdcall-alias 
Linux:  




==================================================================================================================================
= Usage details



In order to use library functions, simply run R session and type:

> dyn.load("example.so")
> source("example.R")

Now library functions can be called with the same functions name
as in C:

> AT_GetNumber()
[1] 137
 

 