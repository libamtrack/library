==================================================================================================================================
= Prerequisites

In order to compile Java GUI for libamtrack following things needs to be installed:

1. Web browser (preferably Firefox) with Java support installed (Java Java Runtime Environment, see http://java.sun.com/javase/downloads/index.jsp)  

2. GCC compiler (preferably v.4.0+) . Linux user usually have it by default. Windows users can download
 port of GCC for windows - MinGW compiler (Minimalist GNU for Windows), from the webpage http://www.mingw.org/.
 
3. Java JDK, for example Sun Java SE v1.6 (visit http://java.sun.com/javase/downloads/widget/jdk6.jsp)

4. Swig interface compiler, preferably v1.3. On most linux platform it available in the repositories as swig package. 
Windows users can download it from project webpage (visit http://www.swig.org/download.html).

5. GNU Scientific Library (http://www.gnu.org/software/gsl/). On most linux platforms it is avalaible in the repositories (please install
also -dev version of packages). Windows users can download it from here: http://gnuwin32.sourceforge.net/packages/gsl.htm


==================================================================================================================================
= Compilation


==========================================
== Windows systems

Locate file wrapper\Java\make_jar.bat and setup following variables:

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

Start windows command line (start -> Run -> cmd.exe), go to the directory
containing make_jar.bat file (using cd command) and run it (typing make_jar.bat).


==========================================
== Linux systems

Locate file distributions\Java\make_jar.sh and setup following variables:

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

Here if you have everything properly installed you will need only to adjust
java include path (you can locate them by trying to find jni.h file using
comand locate jni.h) and gsl libraries (locate libgsl.so).

Start your favorite console,  go to the directory
containing make_jar.sh file (using cd command) and run it (typing ./make_jar.sh).



==================================================================================================================================
= Compilation details


Both make_jar scripts do the same job, just with slightly different parameters of compiler.

==========================================
== Step 1. Generate interface elements using SWIG


First - please read SWIG tutorial : http://www.swig.org/tutorial.html then please 
inspect file example.i
This file includes header files from libamtrack C library and contains list
of function that we will use in GUI.

Command 
swig -java -o c-swig-src/$SWIGWRAPPER.c -outdir java-swig-src example.i
will create interface files which consist of two sets:

A) C interface (*.c) which will be saved into c-swig-src directory and will be later compiled
together with the C libramtrack library
B) Java interface files (*.java) which will be saved into java-swig-src directory and will be
later compiled with Java libamtrack GUI

Please do not commit generated files into repository.


==========================================
== Step 2. Compile libamtrack C library + SWIG C wrapper


In this step all the *.c files will be compiled to the *.o files (these
will go to obj directory). SWIG C wrappers needs to be compiled
with Java JDK include directory as they contains JNI stubs.

Finally all *.o files and necessary libraries (GSL) will be linked together
to create C library file (example.dll or libamtrack.so respectively).

Please note that there is difference in compilation flags:
Windows: -fno-strict-aliasing 
Linux: -fPIC 
and in linking flags:
Windows: -mno-cygwin -Wl,--add-stdcall-alias 
Linux:  


==========================================
== Step 3. Compile Java GUI


All java source files (*.java) will be compiled in this
step to the java bytecode (*.class):

$JAVACEXE java-swig-src/*.java src/*.java -d bin/

Compilation is done on SWIG generated stubs and on JAVA GUI source files.
Generated *.class files will be save in bin directory. 


==========================================
== Step 3. Create JAR file


Java bytecode files (*.class), manifest file together with libamtrack C library
and with GSL libraries will packed into JAR file. 

$JAREXE cvfm example.jar MANIFEST.MF *.class libexample.so libgsl.so libgslcblas.so
or
%jarexe% cvfm example.jar MANIFEST.MF *.class example.dll libgsl.dll libgslcblas.dll

Previously generated *.class files are deleted in this step.

Generated example.jar file can be now used as local GUI.
Can be started by typing command:
java -jar example.jar


==========================================
== Step 4. Sign JAR file


All JAR files that needs to be used in Java Web Start technology
needs to be signed. I've prepared self-signed certificate, stored
in myk file. During signing process it will ask for password, which
is "libamtrack". If you want to generate your own certificate,
please use command:
keytool -genkey -keystore myk -alias jdc

Signed JAR file will be save as webstart/examplesigned-Linux.jar or
webstart/examplesigned-Windows.jar 


==========================================
== Step 5. Deploying web start application


In order to run generated JAR archive from within web browser
you will need to adjust *.jnpl file. Please inspect
webstart/libamtrack-linux.jnlp or webstart/libamtrack-windows.jnlp
respectively.

In the second line of that file you will find:
<jnlp spec="1.0+" codebase="XXX">

Instead of XX you need to put URL. If you test how the applet works
locally, please put there local URL, like:
file:///home/libamtrack/workspace/AmTrack_reorgan/wrapper/Java/webstart
if you test it on the webserver, put there usuall URL:
http://libamtrack.dkfz.org/libamtrack

This URL should point to location where JAR file is stored.

There is also webstart/HelloWorld.html file with links to jnlp files.
 
