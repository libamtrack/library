demo README

0. Rationale

This example was written to provide easy way to check that library
is working.
Source file AT_demo.c can be compiled and linked to libamtrack library
to produce executable file AT_demo.exe. 

Lines starting with "$" character are commands to be executed by reader of this tutorial.

1. Compilation

Our example assume that compiled library is in lib folder. Compilation
is done using Makefile provided in current directory.
First step is to adjust operating system (Linux or Windows) in the Makefile in current directory.

To adjust PATH variable in Windows (which keeps track of the directories
containing executable files that can be run without providing full path), 
run windows command line shell (cmd.exe) and use following command:

$ path %path%;"C:\Program Files\MinGW\bin"
$ path %path%;"C:\Program Files\GnuWin32\bin"

These two command should append directories containing make.exe and gcc.exe to the PATH variable. 
This trick will work only for one cmd session.

Compile:

$ make all

After this step two files should be generated: AT_demo.o (not needed later) and AT_demo.exe.
Compilation schema is following:

AT_demo.c + libamtrack header files (*.h) + gsl header files (*.h) ---(compiling)---->  AT_demo.o 
AT_demo.o + libamtrack library (*.dll or *.so) + gsl libraries (*.dll or *.so) ---(linking)---->  AT_demo.exe

2. Run

When you will have AT_demo.exe file and want to run it under Windows, copy to the same
directory as .exe file all libraries (3 dll files for GSL and libamtrack). 
Under Linux, adjust LD_LIBRARY_PATH variable:

$ export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:../../lib

First let us try to run:

$ ./AT_demo.exe
