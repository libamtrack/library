#!/bin/bash

SWIGEXE=swig
JARSIGNEREXE=jarsigner
JAVACEXE=javac
JAREXE=jar
JAVAINCLUDE="/System/Library/Frameworks/JavaVM.framework/Headers"
GSLINCLUDE="/usr/local/include/gsl"
SYSINCLUDE="/usr/include/sys"
GSLLIB="/usr/local/lib"
GSLDLLA="/usr/local/lib/libgsl.dylib"
GSLDLLB="/usr/local/lib/libgslcblas.dylib"
GCCEXE=gcc

# Generate interface elements using SWIG 
mkdir java-swig-src
mkdir c-swig-src
SWIGWRAPPER=example_wrap
$SWIGEXE -java -o c-swig-src/$SWIGWRAPPER.c -outdir java-swig-src example.i

# Compile libamtrack C library + SWIG C wrapper

mkdir obj
rm obj/$SWIGWRAPPER.o
# Compilation
for a in ../../src/*.c
do
 $GCCEXE -I. -I../../include -I../../include/sys -I$JAVAINCLUDE -I$SYSINCLUDE -I$GSLINCLUDE -O3 -fmessage-length=0 -fPIC -arch i386 -c $a -oobj/`basename $a .c`.o
done
$GCCEXE -I. -I../../include -I$JAVAINCLUDE -I$SYSINCLUDE -I$GSLINCLUDE -fPIC -arch i386 -c example.c -oobj/example.o
$GCCEXE -DSWIG -I. -I../../include -I$JAVAINCLUDE -I$SYSINCLUDE -fPIC -arch i386 -c c-swig-src/$SWIGWRAPPER.c -oobj/$SWIGWRAPPER.o

# Linking
$GCCEXE -shared -arch i386 -L$GSLLIB -o libexample.dylib obj/*.o -lgsl -lgslcblas -lm 

rm obj/*.o

# Compile Java GUI

$JAVACEXE java-swig-src/*.java src/*.java -d bin/

# Create JAR file

cp bin/*.class .
cp $GSLDLLA .
cp $GSLDLLB .
$JAREXE cvfm example.jar MANIFEST.MF *.class libexample.dylib libgsl.dylib libgslcblas.dylib
rm *.class
rm libexample.dylib libgsl.dylib libgslcblas.dylib
rm bin/*.class

# Sign JAR file
# keystore created using command "keytool -genkey -keystore myk -alias jdc", password "libamtrack"

$JARSIGNEREXE -keystore myk -storepass libamtrack -signedjar examplesigned-Mac.jar example.jar jdc
mv examplesigned-Mac.jar webstart/
