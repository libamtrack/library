#****************************************************************************
# Makefile for AT_demo
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
STATICLIB = libamtrack.a
LIBDIR    = ../../lib/
INCLDIR   = ../../include/
RMCMD     = rm
GCC       = gcc 
GSLINCLUDE = "$(GSLPATH)/include/"
GSLLIB     = "$(GSLPATH)/lib/"
else
STATICLIB = libamtrack.lib
LIBDIR    = ..\..\lib\
INCLDIR    = ..\..\include\
RMCMD      = del
GCC        = "$(MINGWPATH)\bin\gcc.exe"
GSLINCLUDE = "$(GSLPATH)\include"
GSLLIB     = "$(GSLPATH)\lib"
endif

AT_demo.o:AT_demo.c 
		$(GCC) -I$(GSLINCLUDE) -I$(INCLDIR) $(CFLAGS) -fPIC AT_demo.c

all:AT_demo.o
		$(GCC) -L$(GSLLIB) -L$(LIBDIR) AT_demo.o $(LFLAGS) -lamtrack -o AT_demo.exe

static:AT_demo.c
		$(GCC) -I$(INCLDIR) -L$(GSLLIB) AT_demo.c $(LIBDIR)$(STATICLIB) $(LFLAGS)  -o AT_demo.exe


clean:
		- $(RMCMD) *.o
		- $(RMCMD) AT_demo.exe
