#****************************************************************************
# Makefile for AT_CPPSC
#****************************************************************************

############################################################
# USER SETTINGS
############################################################
## Please specify your OS here
#  OS        = Linux (or)
#  OS        = Windows
OS        = Linux
## Please set the correct path to your GSL installation here:
#  most likely:
#  GSLPATH   = /usr (Linux)
#  GSLPATH   = C:\Program Files\GnuWin32 (Windows)
GSLPATH   = /usr
## Under Windows please specify the path to your MinGW installation
#  leave empty under Linux: MINGWPATH = 
#  for Windows most likely:
#  MINGWPATH = C:\Program Files\MinGW (Windows)
MINGWPATH = 
############################################################

CFLAGS    = -Wall -c -O3
LFLAGS    = -lm -lgsl -lgslcblas

ifeq ($(OS),Linux)
STATICLIB   = libamtrack.a
LIBDIR    = ../../../../lib/
DSEP      = /
INCLDIR   = ../../../../include/
RMCMD     = rm
GCC       = gcc 
GSLINCLUDE = "$(GSLPATH)/include"
GSLLIB     = "$(GSLPATH)/lib"
else
STATICLIB    = libamtrack.lib
LIBDIR    = ..\..\..\..\lib\
DSEP       = \\
INCLDIR    = ..\..\..\..\include\
RMCMD      = del
GCC        = "$(MINGWPATH)\bin\gcc.exe"
GSLINCLUDE = "$(GSLPATH)\include"
GSLLIB     = "$(GSLPATH)\lib"
endif


AT_CPPSC.o:AT_CPPSC.c 
		$(GCC) -I$(GSLINCLUDE) -I$(INCLDIR) $(CFLAGS) AT_CPPSC.c


all:AT_CPPSC.o
		$(GCC) -L$(GSLLIB) -L$(LIBDIR) AT_CPPSC.o $(LFLAGS) -lamtrack -o AT_CPPSC.exe

static:AT_CPPSC.c
		$(GCC) -I$(INCLDIR) -L$(GSLLIB) AT_CPPSC.c $(LIBDIR)$(STATICLIB) $(LFLAGS) -o AT_CPPSC.exe

clean:
		- $(RMCMD) *.o
		- $(RMCMD) AT_CPPSC.exe
