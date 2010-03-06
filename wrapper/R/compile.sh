#!/bin/bash

SWIGEXE=swig
GSLINCLUDE="/usr/include"
RINCLUDE="/usr/share/R/include/"
RLIB="/usr/lib/R/lib"
GSLLIB="/usr/lib"
GSLDLLA="/usr/lib/libgsl.so"
GSLDLLB="/usr/lib/libgslcblas.so"
GCCEXE=gcc

# Generate interface elements using SWIG 
mkdir c-swig-src
SWIGWRAPPER=example_wrap
$SWIGEXE -dll example -r -o c-swig-src/$SWIGWRAPPER.c -outdir . ../example.i

# Compile libamtrack C library + SWIG C wrapper

mkdir obj
rm obj/$SWIGWRAPPER.o
# Compilation
for a in ../../src/*.c
do
 $GCCEXE -I.. -I../../include -I$GSLINCLUDE -O3 -fmessage-length=0 -fPIC -c $a -oobj/`basename $a .c`.o
done
$GCCEXE -I.. -I../../include -I$GSLINCLUDE -fPIC -c ../example.c -oobj/example.o
$GCCEXE -DDUPA -I.. -I../../include -I$RINCLUDE -fPIC -c c-swig-src/$SWIGWRAPPER.c -oobj/$SWIGWRAPPER.o

# Linking 
# be carefull here - library should have name example.so (not libexample.so !!!) or example.dll
$GCCEXE -shared -L$GSLLIB -L$RLIB -oexample.so obj/*.o -lgsl -lgslcblas -lm -lR 

rm obj/*.o
rm c-swig-src/$SWIGWRAPPER.c
