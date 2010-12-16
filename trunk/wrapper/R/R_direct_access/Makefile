#****************************************
# Makefile for libamtrack direct R access
#****************************************

############################################################
# USER SETTINGS
############################################################
## OS setting detection
TRG       = $(shell uname)
ifeq ($(TRG),Linux)
OS = Linux
else
OS = Windows
endif
## Please set the correct path to your GSL installation here:
#  most likely:
ifeq ($(OS),Linux)
GSLPATH   = /usr
else
#  GSLPATH   = C:\Program Files\GnuWin32 (Windows)
GSLPATH   = C:\Programme\GnuWin32
endif
## Under Windows please specify the path to your MinGW installation
#  leave empty under Linux: MINGWPATH = 
#  for Windows most likely:
#  MINGWPATH = C:\Program Files\MinGW (Windows)
ifeq ($(OS),Linux)
MINGWPATH = 
else
MINGWPATH = "C:\Programme\MinGW"
endif
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
