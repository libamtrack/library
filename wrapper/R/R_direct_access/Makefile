#****************************************
# Makefile for libamtrack direct R access
#****************************************

############################################################
# USER SETTINGS
############################################################
## Please specify your OS here
#  OS        = Linux (or)
#  OS        = Windows
OS        = Windows
## Please set the correct path to your GSL installation here:
#  most likely:
#  GSLPATH   = /usr (Linux)
#  GSLPATH   = C:\Program Files\GnuWin32 (Windows)
GSLPATH   = C:\Programme\GnuWin32
## Under Windows please specify the path to your MinGW installation
#  leave empty under Linux: MINGWPATH = 
#  for Windows most likely:
#  MINGWPATH = C:\Program Files\MinGW (Windows)
MINGWPATH = "C:\Programme\MinGW"
############################################################

CFLAGS    = -Wall -c -O3 -fPIC
LFLAGS    = -lm -lgsl -lgslcblas

ifeq ($(OS),Linux)
NAMELIB   = libamtrack.so
LIBDIR    = ../../../lib
DSEP      = /
INCLDIR   = ../../../include
RMCMD     = rm
GCC       = gcc 
GSLINCLUDE = "$(GSLPATH)/include"
GSLLIB     = "$(GSLPATH)/lib"
else
NAMELIB    = libamtrack.dll
LIBDIR    = "..\..\..\lib"
DSEP       = \\
INCLDIR    = "..\..\..\include"
RMCMD      = del
GCC        = "gcc.exe"
GSLINCLUDE = "$(GSLPATH)\include"
GSLLIB     = "$(GSLPATH)\bin"
endif

all:
		$(GCC) -I$(GSLINCLUDE) -I$(INCLDIR) -I$(MINGWPATH)/include -I. $(CFLAGS) *.c
		$(GCC) -L$(LIBDIR) -L$(GSLLIB) $(LFLAGS) -shared -o $(NAMELIB) *.o

clean:
		- $(RMCMD) *.o
		- $(RMCMD) $(NAMELIB)
