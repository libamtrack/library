#!/bin/bash

SWIGEXE=swig
JARSIGNEREXE=jarsigner
JAVACEXE=javac
JAREXE=jar
JAVAINCLUDEA="/usr/lib/jvm/java-6-sun/include/"
JAVAINCLUDEB="/usr/lib/jvm/java-6-sun/include/linux"
GSLINCLUDE="/usr/include"
GSLLIB="/usr/lib"
GSLDLLA="/usr/lib/libgsl.so"
GSLDLLB="/usr/lib/libgslcblas.so"
GCCEXE=gcc

# Generate interface elements using SWIG 
mkdir -p java-swig-src
mkdir -p c-swig-src
SWIGWRAPPER=example_wrap
$SWIGEXE -java -o c-swig-src/$SWIGWRAPPER.c -outdir java-swig-src example.i

# Compile libamtrack C library + SWIG C wrapper

mkdir -p obj
rm -f obj/$SWIGWRAPPER.o
# Compilation
for a in ../../src/*.c
do
 $GCCEXE -I. -I../../include -I$JAVAINCLUDEB -I$JAVAINCLUDEA  -I$GSLINCLUDE -O3 -fmessage-length=0 -fPIC -c $a -oobj/`basename $a .c`.o
done
$GCCEXE -I. -I../../include -I$JAVAINCLUDEB -I$JAVAINCLUDEA -I$GSLINCLUDE -fPIC -c example.c -oobj/example.o
$GCCEXE -DSWIG -I. -I../../include -I$JAVAINCLUDEB -I$JAVAINCLUDEA -fPIC -c c-swig-src/$SWIGWRAPPER.c -oobj/$SWIGWRAPPER.o

# Linking
$GCCEXE -shared -L$GSLLIB -olibexample.so obj/*.o -lgsl -lgslcblas -lm 

rm -f obj/*.o

# Compile Java GUI

$JAVACEXE java-swig-src/*.java src/*.java -d bin/

# Create JAR file

cp bin/*.class .
cp $GSLDLLA .
cp $GSLDLLB .
$JAREXE cvfm example.jar MANIFEST.MF *.class libexample.so libgsl.so libgslcblas.so
rm -f *.class
rm -f libexample.so libgsl.so libgslcblas.so
rm -f bin/*.class

# Sign JAR file
# keystore created using command "keytool -genkey -keystore myk -alias jdc", password "libamtrack"

$JARSIGNEREXE -keystore myk -storepass libamtrack -signedjar examplesigned-Linux.jar example.jar jdc
mv examplesigned-Linux.jar webstart/
